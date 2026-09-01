"""Side-effect-limited preparation of one reproducible simulation run.

This component sits between schema normalization and process launch.  It owns
materializing omitted seeds and hashing every immutable input.  It does not
create output directories, start the C++ core, or interpret any biological
file format.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import secrets
import stat
from collections.abc import Callable, Iterable, Mapping

from .config import LoadedRunConfig, UINT64_MAX


DEFAULT_HASH_CHUNK_SIZE = 1024 * 1024


class PreparationError(ValueError):
    """A run cannot be made reproducible and safe to launch."""


@dataclass(frozen=True)
class FileDigest:
    """Identity of one input role at the instant its exact bytes were hashed."""

    role: str
    path: Path
    size_bytes: int
    sha256: str
    declared_sha256: str | None = None


@dataclass(frozen=True)
class PreparedRun:
    """Effective config plus verified immutable input identities."""

    config: LoadedRunConfig
    files: tuple[FileDigest, ...]

    def file_for_role(self, role: str) -> FileDigest:
        """Return the unique file assigned to *role*."""
        for file_digest in self.files:
            if file_digest.role == role:
                return file_digest
        raise KeyError(role)


EntropySource = Callable[[int], int]


SEED_DERIVATION_CONTRACT = "sha256-domain-separated-v1"
_SEED_DERIVATION_PREFIX = b"BSReadSim/stage-seed/v1\0"
_SEED_DOMAINS = ("mutation", "phasing", "methylation")


def derive_stage_seed(master_seed: int, domain: str) -> int:
    """Derive one stable unsigned 64-bit stage seed from the master seed."""
    if (
        isinstance(master_seed, bool)
        or not isinstance(master_seed, int)
        or master_seed < 0
        or master_seed > UINT64_MAX
    ):
        raise PreparationError("master seed must be an unsigned 64-bit integer")
    if domain not in _SEED_DOMAINS:
        raise PreparationError("stage seed domain is unsupported")
    digest = hashlib.sha256(
        _SEED_DERIVATION_PREFIX
        + master_seed.to_bytes(8, byteorder="little")
        + b"\0"
        + domain.encode("ascii")
    ).digest()
    return int.from_bytes(digest[:8], byteorder="little")


def materialize_master_seed(
    config: LoadedRunConfig,
    *,
    entropy: EntropySource = secrets.randbits,
) -> LoadedRunConfig:
    """Return the effective normalized config with an explicit u64 seed.

    Explicit seeds are preserved without consulting the entropy source.  An
    omitted seed is generated exactly once and inserted as a decimal string,
    after which canonical JSON and its SHA-256 identity are recomputed.
    """
    if not isinstance(config, LoadedRunConfig):
        raise PreparationError("config must be a LoadedRunConfig")
    if config.master_seed is not None:
        return config

    generated_seed = entropy(64)
    if (
        isinstance(generated_seed, bool)
        or not isinstance(generated_seed, int)
        or generated_seed < 0
        or generated_seed > UINT64_MAX
    ):
        raise PreparationError("entropy source must return an unsigned 64-bit integer")

    return config.with_master_seed(generated_seed)


def materialize_run_seeds(
    config: LoadedRunConfig,
    *,
    entropy: EntropySource = secrets.randbits,
) -> LoadedRunConfig:
    """Materialize the master seed and derive every omitted stage seed."""
    effective = materialize_master_seed(config, entropy=entropy)
    if effective.master_seed is None:  # Guaranteed by materialization.
        raise PreparationError("master seed materialization failed")
    normalized_seeds = effective.normalized.get("seeds")
    if not isinstance(normalized_seeds, Mapping):
        raise PreparationError("config stage seeds are missing")
    unresolved = tuple(
        domain for domain in _SEED_DOMAINS if normalized_seeds.get(domain) is None
    )
    if not unresolved:
        return effective
    return effective.with_resolved_seeds(
        {
            domain: derive_stage_seed(effective.master_seed, domain)
            for domain in unresolved
        }
    )


def prepare_run(
    config: LoadedRunConfig,
    *,
    entropy: EntropySource = secrets.randbits,
    hash_chunk_size: int = DEFAULT_HASH_CHUNK_SIZE,
) -> PreparedRun:
    """Materialize all seeds and hash referenced input/model files.

    Model hashes declared by the config are verified before a core process can
    be launched.  Files referenced by multiple roles are read only once, while
    the returned records retain every semantic role for the manifest.
    """
    if isinstance(hash_chunk_size, bool) or not isinstance(hash_chunk_size, int):
        raise PreparationError("hash_chunk_size must be a positive integer")
    if hash_chunk_size <= 0:
        raise PreparationError("hash_chunk_size must be a positive integer")

    effective_config = materialize_run_seeds(config, entropy=entropy)
    descriptors = tuple(_iter_file_descriptors(effective_config.normalized))
    cache: dict[Path, tuple[int, str]] = {}
    identities = []

    for role, path_text, declared_sha256 in descriptors:
        path = Path(path_text)
        if not path.is_absolute():
            raise PreparationError(
                "normalized input path must be absolute for {}: {}".format(role, path)
            )
        cached = cache.get(path)
        if cached is None:
            cached = _hash_regular_file(path, hash_chunk_size)
            cache[path] = cached
        size_bytes, actual_sha256 = cached
        if declared_sha256 is not None and actual_sha256 != declared_sha256:
            raise PreparationError(
                "{} SHA-256 mismatch for {}: expected {}, observed {}".format(
                    role, path, declared_sha256, actual_sha256
                )
            )
        identities.append(
            FileDigest(
                role=role,
                path=path,
                size_bytes=size_bytes,
                sha256=actual_sha256,
                declared_sha256=declared_sha256,
            )
        )

    return PreparedRun(config=effective_config, files=tuple(identities))


def snapshot_prepared_file(
    file_digest: FileDigest,
    *,
    maximum_size: int,
    chunk_size: int = DEFAULT_HASH_CHUNK_SIZE,
) -> bytes:
    """Read and revalidate one prepared file into an immutable byte snapshot.

    This closes the preparation-to-use path replacement window for small model
    artifacts. A replacement with byte-identical content is harmless and is
    accepted; size/digest drift or identity drift during the read fails before
    process launch.
    """
    if not isinstance(file_digest, FileDigest):
        raise PreparationError("file_digest must be a prepared FileDigest")
    for name, value in (("maximum_size", maximum_size), ("chunk_size", chunk_size)):
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise PreparationError("{} must be a positive integer".format(name))
    if (
        isinstance(file_digest.size_bytes, bool)
        or not isinstance(file_digest.size_bytes, int)
        or file_digest.size_bytes < 0
    ):
        raise PreparationError("prepared file size is invalid")
    if file_digest.size_bytes > maximum_size:
        raise PreparationError(
            "prepared file exceeds the snapshot size limit: {}".format(
                file_digest.path
            )
        )
    if (
        not isinstance(file_digest.sha256, str)
        or len(file_digest.sha256) != 64
        or any(character not in "0123456789abcdef" for character in file_digest.sha256)
    ):
        raise PreparationError("prepared file SHA-256 is invalid")
    if (
        file_digest.declared_sha256 is not None
        and file_digest.declared_sha256 != file_digest.sha256
    ):
        raise PreparationError("prepared file digest disagrees with its declaration")

    path = file_digest.path
    if not isinstance(path, Path) or not path.is_absolute():
        raise PreparationError("prepared file path must be absolute")
    digest = hashlib.sha256()
    payload = bytearray()
    try:
        initial = path.stat()
        if not stat.S_ISREG(initial.st_mode):
            raise PreparationError("input is not a regular file: {}".format(path))
        with path.open("rb") as input_file:
            before = os.fstat(input_file.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise PreparationError("input is not a regular file: {}".format(path))
            if (initial.st_dev, initial.st_ino) != (before.st_dev, before.st_ino):
                raise PreparationError(
                    "input changed before it could be snapshotted: {}".format(path)
                )
            while True:
                chunk = input_file.read(chunk_size)
                if not chunk:
                    break
                if len(chunk) > maximum_size - len(payload):
                    raise PreparationError(
                        "prepared file exceeds the snapshot size limit: {}".format(
                            path
                        )
                    )
                payload.extend(chunk)
                digest.update(chunk)
            after = os.fstat(input_file.fileno())
    except PreparationError:
        raise
    except OSError as error:
        raise PreparationError(
            "cannot snapshot prepared file {}: {}".format(path, error)
        ) from error

    before_identity = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    )
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if before_identity != after_identity:
        raise PreparationError(
            "prepared file changed while it was snapshotted: {}".format(path)
        )
    if len(payload) != file_digest.size_bytes or after.st_size != len(payload):
        raise PreparationError(
            "prepared file size changed before use: {}".format(path)
        )
    observed_sha256 = digest.hexdigest()
    if observed_sha256 != file_digest.sha256:
        raise PreparationError(
            "prepared file SHA-256 changed before use: {}".format(path)
        )
    return bytes(payload)


def _iter_file_descriptors(
    config: Mapping[str, object]
) -> Iterable[tuple[str, str, str | None]]:
    yield "reference", str(config["reference"]), None

    inputs = config["inputs"]
    if not isinstance(inputs, Mapping):  # Protected by the schema; fail closed here too.
        raise PreparationError("normalized inputs section is not an object")
    for name in (
        "vcf",
        "cgmap",
        "bed_methyl",
        "methbg",
        "methbed",
        "methdb",
        "asm",
        "asm_bed",
    ):
        if name in inputs:
            yield "input.{}".format(name), str(inputs[name]), None

    if config["technology"] in ("TBS", "WES", "TS"):
        tbs = config["tbs"]
        if not isinstance(tbs, Mapping):
            raise PreparationError("normalized target section is not an object")
        yield "input.tbs-bed", str(tbs["bed"]), None

    model_containers = (
        ("model.coverage", config, "coverage"),
        ("model.quality", config["sequencing"], "quality"),
        ("model.error", config["sequencing"], "error"),
    )
    for role, parent, key in model_containers:
        if not isinstance(parent, Mapping):
            raise PreparationError("normalized model section is not an object")
        container = parent[key]
        if not isinstance(container, Mapping):
            raise PreparationError("normalized model declaration is not an object")
        artifact = container.get("artifact")
        if artifact is None:
            continue
        if not isinstance(artifact, Mapping):
            raise PreparationError("normalized model artifact is not an object")
        yield role, str(artifact["path"]), str(artifact["sha256"])


def _hash_regular_file(path: Path, chunk_size: int) -> tuple[int, str]:
    digest = hashlib.sha256()
    size_bytes = 0
    try:
        initial = path.stat()
        if not stat.S_ISREG(initial.st_mode):
            raise PreparationError("input is not a regular file: {}".format(path))
        with path.open("rb") as input_file:
            before = os.fstat(input_file.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise PreparationError("input is not a regular file: {}".format(path))
            if (initial.st_dev, initial.st_ino) != (before.st_dev, before.st_ino):
                raise PreparationError(
                    "input changed before it could be hashed: {}".format(path)
                )
            while True:
                chunk = input_file.read(chunk_size)
                if not chunk:
                    break
                size_bytes += len(chunk)
                digest.update(chunk)
            after = os.fstat(input_file.fileno())
    except PreparationError:
        raise
    except OSError as error:
        raise PreparationError("cannot hash input {}: {}".format(path, error)) from error

    before_identity = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    )
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if before_identity != after_identity or size_bytes != after.st_size:
        raise PreparationError("input changed while it was being hashed: {}".format(path))
    return size_bytes, digest.hexdigest()
