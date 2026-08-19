"""Supervise the private C++ core and expose its columnar batch stream.

This component owns only process lifetime and pipe transport. Core stdout is
never logged or decoded here; a bounded transport hands every byte exclusively
to the explicitly selected strict protocol reader. Core stderr is drained on a
dedicated thread so diagnostics cannot deadlock stdout, while only a bounded
tail is retained.

Use the supervisor as a context manager and exhaust its selected stream
iterator before leaving the context. A run is successful only after the
protocol trailer, protocol EOF, an exit status of zero, and a non-overflowing
stderr collector have all been
observed.  Leaving early is a failure and terminates the child.
"""

from __future__ import annotations

from collections import deque
import math
import os
from pathlib import Path
import subprocess
import threading
import time
from types import TracebackType
from typing import (
    BinaryIO,
    Callable,
    Deque,
    Dict,
    Iterator,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
)

from .protocol import (
    DecodedBatchView,
    Header,
    PROTOCOL_MAJOR,
    PROTOCOL_MINOR,
    ProtocolReader,
    Trailer,
)


PathLike = Union[str, os.PathLike]
class CoreProcessError(RuntimeError):
    """The core process did not produce one complete, trustworthy run."""

    def __init__(
        self,
        message: str,
        *,
        argv: Sequence[str],
        returncode: Optional[int],
        stderr_tail: bytes,
        stderr_truncated: bool,
    ) -> None:
        self.argv = tuple(argv)
        self.returncode = returncode
        self.stderr_tail = bytes(stderr_tail)
        self.stderr_truncated = stderr_truncated
        detail = message
        if returncode is not None:
            detail += "; core exit status={}".format(returncode)
        if stderr_tail:
            diagnostic = stderr_tail.decode("utf-8", "replace")
            detail += "; stderr tail={!r}".format(diagnostic)
        if stderr_truncated:
            detail += "; stderr exceeded the configured limit"
        super().__init__(detail)


class _ProtocolIdleTimeout(TimeoutError):
    pass


def _positive_timeout(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("{} must be a positive finite number".format(name))
    converted = float(value)
    if not math.isfinite(converted) or converted <= 0.0:
        raise ValueError("{} must be a positive finite number".format(name))
    return converted


def _optional_positive_timeout(name: str, value: Optional[float]) -> Optional[float]:
    if value is None:
        return None
    return _positive_timeout(name, value)


def _positive_int(name: str, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError("{} must be a positive integer".format(name))
    return value


def _normalize_argv(argv: Sequence[PathLike]) -> Tuple[str, ...]:
    if isinstance(argv, (str, bytes, os.PathLike)):
        raise ValueError("argv must be a non-empty sequence of arguments")
    try:
        raw_arguments = tuple(argv)
    except TypeError as error:
        raise ValueError("argv must be a non-empty sequence of arguments") from error
    if not raw_arguments:
        raise ValueError("argv must not be empty")

    normalized = []
    for index, argument in enumerate(raw_arguments):
        try:
            value = os.fspath(argument)
        except TypeError as error:
            raise ValueError("argv[{}] must be a string or path".format(index)) from error
        if not isinstance(value, str):
            raise ValueError("argv[{}] must resolve to text, not bytes".format(index))
        if not value:
            raise ValueError("argv[{}] must not be empty".format(index))
        if "\x00" in value:
            raise ValueError("argv[{}] must not contain NUL".format(index))
        normalized.append(value)
    return tuple(normalized)


class _BoundedStdoutTransport:
    """Backpressured byte transport used only by ``ProtocolReader``."""

    _READ_SIZE = 64 * 1024

    def __init__(
        self,
        pipe: BinaryIO,
        *,
        maximum_buffer_bytes: int,
        idle_timeout_seconds: Optional[float],
    ) -> None:
        self._pipe = pipe
        self._maximum_buffer_bytes = maximum_buffer_bytes
        self._idle_timeout_seconds = idle_timeout_seconds
        self._condition = threading.Condition()
        self._chunks = deque()  # type: Deque[bytes]
        self._buffered_bytes = 0
        self._eof = False
        self._closing = False
        self._error = None  # type: Optional[BaseException]
        self._thread = threading.Thread(
            target=self._pump,
            name="htsim-core-stdout",
            daemon=True,
        )
        self._thread.start()

    def _pump(self) -> None:
        try:
            while True:
                with self._condition:
                    while (
                        not self._closing
                        and self._buffered_bytes >= self._maximum_buffer_bytes
                    ):
                        self._condition.wait()
                    if self._closing:
                        return
                    capacity = self._maximum_buffer_bytes - self._buffered_bytes
                    read_size = min(self._READ_SIZE, capacity)
                chunk = self._pipe.read(read_size)
                if chunk is None:
                    raise OSError("core stdout read made no progress")
                if not isinstance(chunk, bytes):
                    raise OSError("core stdout returned non-bytes data")
                with self._condition:
                    if self._closing:
                        return
                    if not chunk:
                        self._eof = True
                        self._condition.notify_all()
                        return
                    self._chunks.append(chunk)
                    self._buffered_bytes += len(chunk)
                    self._condition.notify_all()
        except BaseException as error:
            with self._condition:
                if not self._closing:
                    self._error = error
                self._condition.notify_all()
        finally:
            try:
                self._pipe.close()
            except Exception:
                pass

    def read(self, size: int = -1) -> bytes:
        if not isinstance(size, int) or isinstance(size, bool) or size == 0 or size < -1:
            raise ValueError("protocol transport reads require a positive size or -1")
        deadline = None
        if self._idle_timeout_seconds is not None:
            deadline = time.monotonic() + self._idle_timeout_seconds

        with self._condition:
            while not self._chunks:
                if self._error is not None:
                    raise OSError("core stdout transport failed") from self._error
                if self._eof or self._closing:
                    return b""
                if deadline is None:
                    self._condition.wait()
                else:
                    remaining = deadline - time.monotonic()
                    if remaining <= 0.0:
                        raise _ProtocolIdleTimeout(
                            "core stdout was idle for {:.3f} seconds".format(
                                self._idle_timeout_seconds
                            )
                        )
                    self._condition.wait(remaining)

            requested = self._buffered_bytes if size == -1 else size
            output = bytearray()
            while self._chunks and len(output) < requested:
                chunk = self._chunks[0]
                needed = requested - len(output)
                if len(chunk) <= needed:
                    output.extend(self._chunks.popleft())
                    self._buffered_bytes -= len(chunk)
                else:
                    output.extend(chunk[:needed])
                    self._chunks[0] = chunk[needed:]
                    self._buffered_bytes -= needed
            self._condition.notify_all()
            return bytes(output)

    def close(self, join_timeout_seconds: float) -> bool:
        with self._condition:
            self._closing = True
            self._chunks.clear()
            self._buffered_bytes = 0
            self._condition.notify_all()
        try:
            self._pipe.close()
        except Exception:
            pass
        self._thread.join(join_timeout_seconds)
        return not self._thread.is_alive()


class _StderrTailCollector:
    _READ_SIZE = 64 * 1024

    def __init__(self, pipe: BinaryIO, limit_bytes: int) -> None:
        self._pipe = pipe
        self._limit_bytes = limit_bytes
        self._tail = bytearray()
        self._total_bytes = 0
        self._error = None  # type: Optional[BaseException]
        self._closing = False
        self._lock = threading.Lock()
        self._thread = threading.Thread(
            target=self._pump,
            name="htsim-core-stderr",
            daemon=True,
        )
        self._thread.start()

    def _pump(self) -> None:
        try:
            while True:
                chunk = self._pipe.read(self._READ_SIZE)
                if chunk is None:
                    raise OSError("core stderr read made no progress")
                if not isinstance(chunk, bytes):
                    raise OSError("core stderr returned non-bytes data")
                if not chunk:
                    return
                with self._lock:
                    self._total_bytes += len(chunk)
                    self._tail.extend(chunk)
                    excess = len(self._tail) - self._limit_bytes
                    if excess > 0:
                        del self._tail[:excess]
        except BaseException as error:
            with self._lock:
                if not self._closing:
                    self._error = error
        finally:
            try:
                self._pipe.close()
            except Exception:
                pass

    @property
    def tail(self) -> bytes:
        with self._lock:
            return bytes(self._tail)

    @property
    def truncated(self) -> bool:
        with self._lock:
            return self._total_bytes > self._limit_bytes

    @property
    def error(self) -> Optional[BaseException]:
        with self._lock:
            return self._error

    def close(self, join_timeout_seconds: float) -> bool:
        with self._lock:
            self._closing = True
        try:
            self._pipe.close()
        except Exception:
            pass
        self._thread.join(join_timeout_seconds)
        return not self._thread.is_alive()

    def join(self, timeout_seconds: float) -> bool:
        self._thread.join(timeout_seconds)
        return not self._thread.is_alive()


class CoreProcess:
    """One-shot, fail-closed supervisor for a core protocol process.

    ``argv`` is always passed directly to :class:`subprocess.Popen` with
    ``shell=False``.  The class is intentionally one-shot: protocol streams
    cannot be resumed or retried after any failure.
    """

    def __init__(
        self,
        argv: Sequence[PathLike],
        *,
        expected_header: Optional[Header] = None,
        read_length: Optional[int] = None,
        paired_end: Optional[bool] = None,
        expected_skipped_fragment_count: Optional[int] = None,
        cwd: Optional[PathLike] = None,
        env: Optional[Mapping[str, str]] = None,
        stderr_limit_bytes: int = 64 * 1024,
        stdout_buffer_bytes: int = 256 * 1024,
        protocol_idle_timeout_seconds: Optional[float] = None,
        exit_timeout_seconds: float = 10.0,
        terminate_timeout_seconds: float = 1.0,
        kill_timeout_seconds: float = 2.0,
    ) -> None:
        self._argv = _normalize_argv(argv)
        if expected_header is not None and not isinstance(expected_header, Header):
            raise ValueError("expected_header must be a protocol Header")
        if read_length is not None:
            _positive_int("read_length", read_length)
        if paired_end is not None and not isinstance(paired_end, bool):
            raise ValueError("paired_end must be a boolean")
        self._expected_header = expected_header
        self._read_length = read_length
        self._paired_end = paired_end
        self._expected_skipped_count = expected_skipped_fragment_count
        self._cwd = None if cwd is None else Path(cwd)
        self._env = None if env is None else self._normalize_env(env)
        self._stderr_limit_bytes = _positive_int(
            "stderr_limit_bytes", stderr_limit_bytes
        )
        self._stdout_buffer_bytes = _positive_int(
            "stdout_buffer_bytes", stdout_buffer_bytes
        )
        self._protocol_idle_timeout = _optional_positive_timeout(
            "protocol_idle_timeout_seconds", protocol_idle_timeout_seconds
        )
        self._exit_timeout = _positive_timeout(
            "exit_timeout_seconds", exit_timeout_seconds
        )
        self._terminate_timeout = _positive_timeout(
            "terminate_timeout_seconds", terminate_timeout_seconds
        )
        self._kill_timeout = _positive_timeout(
            "kill_timeout_seconds", kill_timeout_seconds
        )

        self._state = "new"
        self._entered = False
        self._iterated = False
        self._process = None  # type: Optional[subprocess.Popen]
        self._stdout_transport = None  # type: Optional[_BoundedStdoutTransport]
        self._stderr_collector = None  # type: Optional[_StderrTailCollector]
        self._reader = None  # type: Optional[ProtocolReader]
        self._returncode = None  # type: Optional[int]
        self._kill_used = False
        self._cleanup_issue = None  # type: Optional[str]

    @staticmethod
    def _normalize_env(env: Mapping[str, str]) -> Dict[str, str]:
        if not isinstance(env, Mapping):
            raise ValueError("env must be a string mapping")
        normalized = {}  # type: Dict[str, str]
        for key, value in env.items():
            if not isinstance(key, str) or not isinstance(value, str):
                raise ValueError("env keys and values must be strings")
            if not key or "=" in key or "\x00" in key or "\x00" in value:
                raise ValueError("env contains an invalid key or value")
            normalized[key] = value
        return normalized

    @property
    def argv(self) -> Tuple[str, ...]:
        return self._argv

    @property
    def pid(self) -> Optional[int]:
        return None if self._process is None else self._process.pid

    @property
    def returncode(self) -> Optional[int]:
        if self._process is not None:
            observed = self._process.poll()
            if observed is not None:
                self._returncode = observed
        return self._returncode

    @property
    def stderr_tail(self) -> bytes:
        if self._stderr_collector is None:
            return b""
        return self._stderr_collector.tail

    @property
    def stderr_text(self) -> str:
        return self.stderr_tail.decode("utf-8", "replace")

    @property
    def stderr_truncated(self) -> bool:
        return (
            self._stderr_collector is not None
            and self._stderr_collector.truncated
        )

    @property
    def kill_used(self) -> bool:
        return self._kill_used

    @property
    def succeeded(self) -> bool:
        return self._state == "succeeded"

    @property
    def protocol_version(self) -> Tuple[int, int]:
        if self._reader is None:
            raise self._error(
                "protocol version is unavailable before process startup"
            )
        return PROTOCOL_MAJOR, PROTOCOL_MINOR

    @property
    def header(self) -> Header:
        if self._reader is None:
            raise CoreProcessError(
                "core header is unavailable before process startup",
                argv=self._argv,
                returncode=self.returncode,
                stderr_tail=self.stderr_tail,
                stderr_truncated=self.stderr_truncated,
            )
        return self._reader.header

    @property
    def trailer(self) -> Trailer:
        if not self.succeeded or self._reader is None or self._reader.trailer is None:
            raise CoreProcessError(
                "core trailer is unavailable before successful completion",
                argv=self._argv,
                returncode=self.returncode,
                stderr_tail=self.stderr_tail,
                stderr_truncated=self.stderr_truncated,
            )
        return self._reader.trailer

    def __enter__(self) -> "CoreProcess":
        if self._state != "new" or self._entered:
            raise CoreProcessError(
                "core process supervisors are one-shot context managers",
                argv=self._argv,
                returncode=self.returncode,
                stderr_tail=self.stderr_tail,
                stderr_truncated=self.stderr_truncated,
            )
        self._entered = True
        self._start()
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> bool:
        del exc_value, traceback
        if exc_type is not None:
            self._abort()
            return False
        if not self.succeeded:
            self._abort()
            raise self._error(
                "fragment stream was not fully consumed and verified"
            )
        self._close_resources()
        return False

    def _start(self) -> None:
        try:
            process = subprocess.Popen(
                self._argv,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=None if self._cwd is None else str(self._cwd),
                env=self._env,
                shell=False,
                bufsize=0,
                close_fds=True,
            )
            self._process = process
            if process.stdout is None or process.stderr is None:
                raise RuntimeError("core process pipes were not created")
            self._stderr_collector = _StderrTailCollector(
                process.stderr, self._stderr_limit_bytes
            )
            self._stdout_transport = _BoundedStdoutTransport(
                process.stdout,
                maximum_buffer_bytes=self._stdout_buffer_bytes,
                idle_timeout_seconds=self._protocol_idle_timeout,
            )
            self._reader = ProtocolReader(
                self._stdout_transport,
                expected_header=self._expected_header,
                expected_skipped_fragment_count=self._expected_skipped_count,
            )
            self._validate_execution_shape(self._reader.header)
            self._state = "running"
        except Exception as error:
            self._abort()
            raise self._error(
                "could not start and decode the core stream: {}".format(error)
            ) from error

    def iter_batches(self) -> Iterator[DecodedBatchView]:
        """Yield authenticated immutable columnar batches."""
        if not self._entered or self._state != "running" or self._reader is None:
            raise self._error(
                "batches may be consumed only inside the active context"
            )
        if self._iterated:
            raise self._error("the core protocol stream may be iterated only once")
        self._iterated = True
        return self._iterate_batches()

    def _iterate_batches(self) -> Iterator[DecodedBatchView]:
        assert isinstance(self._reader, ProtocolReader)
        try:
            for batch in self._reader:
                yield batch
            if self._reader.trailer is None:
                raise RuntimeError("protocol reader reached EOF without a trailer")
            self._complete_process()
            self._close_resources()
            if self._cleanup_issue is not None:
                raise RuntimeError(self._cleanup_issue)
            self._state = "succeeded"
        except Exception as error:
            self._abort()
            if isinstance(error, CoreProcessError):
                raise
            raise self._error(
                "core batch stream validation failed: {}".format(error)
            ) from error

    def consume_batches(
        self, consumer: Callable[[DecodedBatchView], object]
    ) -> Trailer:
        """Run a batch consumer under the strict process lifetime."""
        if not callable(consumer):
            raise ValueError("consumer must be callable")
        with self as running:
            for batch in running.iter_batches():
                consumer(batch)
        return self.trailer

    def _validate_execution_shape(self, header: Header) -> None:
        if self._paired_end is not None:
            expected_mates = 2 if self._paired_end else 1
            if header.mates_per_fragment != expected_mates:
                raise ValueError(
                    "protocol header disagrees with paired_end"
                )
        if self._read_length is not None:
            lengths = (header.read_length_r1,)
            if header.mates_per_fragment == 2:
                lengths += (header.read_length_r2,)
            if any(length != self._read_length for length in lengths):
                raise ValueError(
                    "protocol header disagrees with read_length"
                )

    def _complete_process(self) -> None:
        assert self._process is not None
        try:
            self._returncode = self._process.wait(timeout=self._exit_timeout)
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(
                "core did not exit after closing its protocol stream"
            ) from error

        if self._stderr_collector is None:
            raise RuntimeError("core stderr collector is unavailable")
        if not self._stderr_collector.join(self._kill_timeout):
            raise RuntimeError("core stderr collector did not reach EOF")
        if self._stderr_collector.error is not None:
            raise RuntimeError("core stderr collection failed") from (
                self._stderr_collector.error
            )
        if self.stderr_truncated:
            raise RuntimeError(
                "core stderr exceeded {} bytes".format(
                    self._stderr_limit_bytes
                )
            )
        self._reader.validate_core_exit_status(self._returncode)

    def _abort(self) -> None:
        if self._state == "succeeded":
            self._close_resources()
            return
        process = self._process
        if process is not None and process.poll() is None:
            try:
                process.terminate()
            except (OSError, ProcessLookupError):
                pass
            try:
                self._returncode = process.wait(timeout=self._terminate_timeout)
            except subprocess.TimeoutExpired:
                self._kill_used = True
                try:
                    process.kill()
                except (OSError, ProcessLookupError):
                    pass
                try:
                    self._returncode = process.wait(timeout=self._kill_timeout)
                except subprocess.TimeoutExpired:
                    self._cleanup_issue = "core remained alive after kill timeout"
            else:
                self._returncode = process.returncode
        elif process is not None:
            self._returncode = process.returncode
        self._state = "failed"
        self._close_resources()

    def _close_resources(self) -> None:
        if self._stdout_transport is not None:
            if not self._stdout_transport.close(self._kill_timeout):
                if self._cleanup_issue is None:
                    self._cleanup_issue = "stdout transport thread did not stop"
        if self._stderr_collector is not None:
            if not self._stderr_collector.close(self._kill_timeout):
                if self._cleanup_issue is None:
                    self._cleanup_issue = "stderr collector thread did not stop"

    def _error(self, message: str) -> CoreProcessError:
        if self._cleanup_issue:
            message += "; " + self._cleanup_issue
        return CoreProcessError(
            message,
            argv=self._argv,
            returncode=self.returncode,
            stderr_tail=self.stderr_tail,
            stderr_truncated=self.stderr_truncated,
        )


__all__ = ["CoreProcess", "CoreProcessError"]
