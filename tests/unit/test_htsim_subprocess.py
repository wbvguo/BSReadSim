"""Process-lifetime tests for the private C++ core supervisor boundary."""

import hashlib
from dataclasses import replace
from pathlib import Path
import sys
import tempfile
import time
import unittest


from bsreadsim.htsim.subprocess import CoreProcess, CoreProcessError
from bsreadsim.htsim.protocol import (
    PROTOCOL_VERSION,
    RNG_CONTRACT,
    AmbiguityPolicy,
    BaseEncoding,
    CaptureStrand,
    Contig,
    FragmentBatch,
    Header,
    MethylationAllele,
    MethylationContext,
    MethylationSource,
    Technology,
)
from tests.helpers.stream_support import encode_stream


HELPER = r"""
import os
import signal
import sys
import time

mode = sys.argv[1]
payload = bytes.fromhex(sys.argv[2])

def write_all(fd, data):
    view = memoryview(data)
    while view:
        written = os.write(fd, view)
        if written <= 0:
            raise RuntimeError("helper write made no progress")
        view = view[written:]

if "hang" in mode:
    signal.signal(signal.SIGTERM, signal.SIG_IGN)

if mode == "valid":
    write_all(1, payload)
elif mode == "valid-diagnostic":
    write_all(2, b"core diagnostic\n")
    write_all(1, payload)
elif mode == "short":
    write_all(1, payload[:-7])
elif mode == "nonzero":
    write_all(1, payload)
    raise SystemExit(7)
elif mode == "large-stderr":
    write_all(2, b"x" * 400000 + b"TAIL-MARKER")
    write_all(1, payload)
elif mode == "valid-hang":
    write_all(1, payload)
    while True:
        time.sleep(1)
elif mode == "invalid-hang":
    write_all(1, b"NOT-A-PROTOCOL")
    while True:
        time.sleep(1)
elif mode == "closed-stdout-hang":
    write_all(1, payload)
    os.close(1)
    while True:
        time.sleep(1)
else:
    raise SystemExit("unknown helper mode")
"""


HEADER = Header(
    run_id="00000000-0000-4000-8000-0000000000aa",
    core_version="1.2.3",
    rng_contract=RNG_CONTRACT,
    master_seed=17,
    normalized_config_sha256=hashlib.sha256(b"config").digest(),
    technology=Technology.WGBS,
    has_details=False,
    mates_per_fragment=1,
    base_encoding=BaseEncoding.ACGTN_U8,
    ambiguity_policy=AmbiguityPolicy.PRESERVE_N,
    read_length_r1=4,
    read_length_r2=0,
    contigs=(
        Contig(
            name="chr1",
            length=1000,
            reference_sha256=hashlib.sha256(b"reference").digest(),
        ),
    ),
)
BATCH = FragmentBatch(
    first_fragment_ordinal=0,
    contig_indices=(0,),
    reference_starts=(100,),
    reference_ends=(104,),
    template_offsets=(0, 4),
    mate_offsets=(0, 1),
    site_offsets=(0, 1),
    mate_template_starts=(0,),
    mate_template_ends=(4,),
    site_template_offsets=(1,),
    site_probabilities=(0.5,),
    haplotypes=(0,),
    capture_strands=(CaptureStrand.FORWARD,),
    mate_indices=(0,),
    mate_reverse_complements=(0,),
    site_contexts=(MethylationContext.CG_C,),
    methylation_sources=(MethylationSource.CGMAP,),
    site_alleles=(MethylationAllele.SHARED,),
    template_bases=(0, 1, 2, 3),
    details=None,
)
VALID_STREAM = encode_stream(HEADER, (BATCH,), skipped_fragment_count=3)


def helper_argv(mode: str, *extra: str):
    return (sys.executable, "-c", HELPER, mode, VALID_STREAM.hex()) + extra


def supervisor(mode: str, **overrides) -> CoreProcess:
    options = {
        "expected_header": HEADER,
        "read_length": 4,
        "paired_end": False,
        "expected_skipped_fragment_count": 3,
        "stderr_limit_bytes": 4096,
        "stdout_buffer_bytes": 128,
        "exit_timeout_seconds": 0.3,
        "terminate_timeout_seconds": 0.1,
        "kill_timeout_seconds": 1.0,
    }
    options.update(overrides)
    return CoreProcess(helper_argv(mode), **options)


def assert_reaped(testcase: unittest.TestCase, core: CoreProcess) -> None:
    testcase.assertIsNotNone(core.returncode)
    testcase.assertIsNotNone(core._process)
    testcase.assertIsNotNone(core._stdout_transport)
    testcase.assertIsNotNone(core._stderr_collector)
    testcase.assertIsNotNone(core._process.poll())
    testcase.assertFalse(core._stdout_transport._thread.is_alive())
    testcase.assertFalse(core._stderr_collector._thread.is_alive())
    testcase.assertTrue(core._process.stdout.closed)
    testcase.assertTrue(core._process.stderr.closed)


class CoreProcessSuccessTests(unittest.TestCase):
    def test_batches_require_trailer_eof_exit_and_retain_payload(self) -> None:
        core = supervisor("valid")

        with core as running:
            self.assertEqual(running.header, HEADER)
            self.assertEqual(running.protocol_version, PROTOCOL_VERSION)
            batches = tuple(running.iter_batches())
            self.assertEqual(len(batches), 1)
            self.assertEqual(tuple(batches[0].template_bases), (0, 1, 2, 3))
            payload = batches[0].raw_payload
            self.assertTrue(payload.readonly)

        self.assertTrue(core.succeeded)
        self.assertEqual(core.trailer.fragment_batch_count, 1)
        self.assertEqual(core.trailer.skipped_fragment_count, 3)
        self.assertGreater(len(payload), 0)
        assert_reaped(self, core)

    def test_consume_batches_uses_the_same_supervised_boundary(self) -> None:
        observed = []
        core = supervisor("valid-diagnostic")

        trailer = core.consume_batches(observed.append)

        self.assertEqual(len(observed), 1)
        self.assertEqual(trailer.fragment_count, 1)
        self.assertEqual(core.stderr_text, "core diagnostic\n")

    def test_argv_is_not_interpreted_by_a_shell(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            marker = Path(directory) / "shell-injection-marker"
            malicious_argument = "; touch {}".format(marker)
            core = CoreProcess(
                helper_argv("valid", malicious_argument),
                expected_header=HEADER,
                read_length=4,
                paired_end=False,
                expected_skipped_fragment_count=3,
            )
            with core as running:
                self.assertEqual(len(tuple(running.iter_batches())), 1)
            self.assertFalse(marker.exists())
            self.assertEqual(core.argv[-1], malicious_argument)


class CoreProcessFailureTests(unittest.TestCase):
    def test_rejects_wrong_header_truncation_and_nonzero_exit(self) -> None:
        cases = (
            (
                supervisor("valid", expected_header=replace(HEADER, master_seed=18)),
                "start and decode",
            ),
            (supervisor("short"), "batch stream validation"),
            (supervisor("nonzero"), "non-zero status"),
        )
        for core, pattern in cases:
            with self.subTest(pattern=pattern):
                with self.assertRaisesRegex(CoreProcessError, pattern):
                    with core as running:
                        tuple(running.iter_batches())
                self.assertFalse(core.succeeded)
                assert_reaped(self, core)

    def test_stderr_is_drained_but_overflow_rejects_and_retains_tail(self) -> None:
        core = supervisor("large-stderr", stderr_limit_bytes=1024)
        with self.assertRaisesRegex(CoreProcessError, "stderr exceeded") as raised:
            with core as running:
                tuple(running.iter_batches())
        self.assertTrue(core.stderr_truncated)
        self.assertEqual(len(core.stderr_tail), 1024)
        self.assertTrue(core.stderr_tail.endswith(b"TAIL-MARKER"))
        self.assertEqual(raised.exception.stderr_tail, core.stderr_tail)
        assert_reaped(self, core)

    def test_consumer_exception_terminates_then_kills_ignoring_child(self) -> None:
        class ConsumerFailure(Exception):
            pass

        core = supervisor("valid-hang")
        started = time.monotonic()
        with self.assertRaises(ConsumerFailure):
            with core as running:
                for _batch in running.iter_batches():
                    raise ConsumerFailure("processing failed")
        self.assertLess(time.monotonic() - started, 2.0)
        self.assertTrue(core.kill_used)
        assert_reaped(self, core)

    def test_decoder_idle_timeout_kills_a_hung_invalid_child(self) -> None:
        core = supervisor("invalid-hang", protocol_idle_timeout_seconds=0.1)
        started = time.monotonic()
        with self.assertRaisesRegex(CoreProcessError, "start and decode"):
            with core:
                self.fail("invalid protocol unexpectedly entered the context")
        self.assertLess(time.monotonic() - started, 2.0)
        self.assertTrue(core.kill_used)
        assert_reaped(self, core)

    def test_child_that_closes_stdout_but_hangs_hits_exit_timeout(self) -> None:
        core = supervisor("closed-stdout-hang")
        started = time.monotonic()
        with self.assertRaisesRegex(CoreProcessError, "did not exit"):
            with core as running:
                tuple(running.iter_batches())
        self.assertLess(time.monotonic() - started, 2.0)
        self.assertTrue(core.kill_used)
        assert_reaped(self, core)

    def test_early_context_exit_is_a_failure_and_reaps_child(self) -> None:
        core = supervisor("valid-hang")
        with self.assertRaisesRegex(CoreProcessError, "not fully consumed"):
            with core as running:
                next(running.iter_batches())
        self.assertTrue(core.kill_used)
        assert_reaped(self, core)


class CoreProcessArgumentTests(unittest.TestCase):
    def test_invalid_commands_and_bounds_fail_before_launch(self) -> None:
        cases = (
            (("command\x00bad",), {}, "NUL"),
            ("command --flag", {}, "sequence"),
            ((), {}, "empty"),
            ((sys.executable,), {"stderr_limit_bytes": 0}, "positive integer"),
            (
                (sys.executable,),
                {"protocol_idle_timeout_seconds": float("nan")},
                "positive finite",
            ),
            ((sys.executable,), {"expected_header": object()}, "protocol Header"),
            ((sys.executable,), {"read_length": 0}, "positive integer"),
            ((sys.executable,), {"paired_end": "false"}, "boolean"),
        )
        for argv, options, pattern in cases:
            with self.subTest(argv=argv, options=options):
                with self.assertRaisesRegex(ValueError, pattern):
                    CoreProcess(argv, **options)

    def test_iteration_requires_active_context_and_is_one_shot(self) -> None:
        core = supervisor("valid")
        with self.assertRaisesRegex(CoreProcessError, "active context"):
            core.iter_batches()

        with core as running:
            iterator = running.iter_batches()
            self.assertEqual(next(iterator).fragment_count, 1)
            with self.assertRaisesRegex(CoreProcessError, "only once"):
                running.iter_batches()
            with self.assertRaises(StopIteration):
                next(iterator)

        with self.assertRaisesRegex(CoreProcessError, "one-shot"):
            with core:
                pass


if __name__ == "__main__":
    unittest.main()
