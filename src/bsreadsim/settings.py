"""The settings of one bsreadsim command.

A ``Settings`` holds every choice of a run, build, or validation once it is
resolved: defaults filled in, paths absolute, and seeds known (``cli`` resolves
it from the command line). The rest derives from it: the htsim command line
(``htsim``), the read simulator and its models (``reads``), and the manifest
(``output``). Constructing one checks the options that Python owns; htsim
checks the rest and how they combine.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re

from .errors import BSReadSimError

BISULFITE = ("WGBS", "RRBS", "TBS")
WHOLE_GENOME = ("WGBS", "WGS")
DEFAULT_CUT_SITES = ("C|CGG",)  # MspI
_PREFIX = re.compile(r"[A-Za-z0-9._-]{1,128}")


@dataclass(frozen=True)
class Settings:
    technology: str = "WGBS"
    reference: Path | None = None
    # variant and methylation inputs
    vcf: Path | None = None
    cgmap: Path | None = None
    bed_methyl: Path | None = None
    methbg: Path | None = None
    methbed: Path | None = None
    methdb: Path | None = None
    asm: Path | None = None
    asm_bed: Path | None = None
    # amount and shape of fragments
    reads: int | None = 1_000_000  # read records; None when depth is set
    depth: float | None = None
    paired_end: bool = True
    read_length: tuple[int, ...] = (100,)  # every read, or read 1 and read 2
    insert_min: int = 100
    insert_mean: int = 400
    insert_max: int = 1000
    insert_sd: float = 25.0
    max_ambiguous_fraction: float = 0.05
    # sampling
    sampling: str = "uniform"
    gc_profile: Path | None = None
    cut_sites: tuple[str, ...] = ()
    rrbs_candidates: Path | None = None
    targets: Path | None = None
    center_sd: float = 50.0
    # de novo variants
    mutation_rate: float = 0.0
    indel_fraction: float = 0.15
    indel_extension_probability: float = 0.15
    homozygous_only: bool = False
    # methylation
    beta_cg: tuple[float, float] = (0.5, 0.5)
    beta_chg: tuple[float, float] = (0.01, 0.05)
    beta_chh: tuple[float, float] = (0.01, 0.05)
    cpg_only: bool = False
    pool_meth: bool = False
    meth_model: str = "bernoulli"
    # chemistry and sequencing
    conversion_rate: float = 0.998
    directional: bool = True
    phred: int | None = 40
    quality_model: Path | None = None
    error_rate: float | None = 0.005
    error_model: Path | None = None
    # seeds
    seed: int = 0
    seed_mut: int = 0
    seed_phase: int = 0
    seed_meth: int = 0
    # execution and output
    threads: int = 1
    output: Path | None = None  # run output directory
    prefix: str = "sim"
    format: str = "fastq.gz"
    gzip_level: int = 4  # as bcl2fastq and fastp: much faster than 6, a few % larger
    fragment_summary: bool = False
    fragment_realization: bool = False
    save_methdb: bool = False
    save_vcf: bool = False

    def __post_init__(self) -> None:
        def require(condition: bool, message: str) -> None:
            if not condition:
                raise BSReadSimError(message)

        require(_PREFIX.fullmatch(self.prefix) is not None,
                "--prefix must be 1-128 characters from A-Z, a-z, 0-9, '.', '_', '-'")
        require(0 <= self.gzip_level <= 9, "--gzip-level must be in [0, 9]")
        require(1 <= self.threads <= 256, "--threads must be in [1, 256]")
        require(0 <= self.conversion_rate <= 1, "--conversion-rate must be in [0, 1]")
        require(self.phred is None or 0 <= self.phred <= 93, "--phred must be in [0, 93]")
        require(self.error_rate is None or 0 <= self.error_rate <= 1,
                "--error-rate must be in [0, 1]")
        require(self.alignments or not self.fragment_summary,
                "--fragment-summary/--fragment-realization requires --format bam or sam")

    @property
    def bisulfite(self) -> bool:
        return self.technology in BISULFITE

    @property
    def alignments(self) -> bool:
        """Reads are written as annotated alignments (BAM or SAM)."""
        return self.format in ("bam", "sam")

    @property
    def fragments(self) -> int | None:
        if self.reads is None:
            return None
        return self.reads // (2 if self.paired_end else 1)

    @property
    def read_lengths(self) -> tuple[int, int]:
        """The lengths of read 1 and read 2."""
        return self.read_length[0], self.read_length[-1]

    def as_dict(self) -> dict:
        """JSON-compatible values of every setting; seeds are decimal strings."""
        return {
            name: (str(value) if isinstance(value, Path) or name.startswith("seed")
                   else list(value) if isinstance(value, tuple) else value)
            for name, value in dataclasses.asdict(self).items()
        }

    @property
    def sha256(self) -> str:
        text = json.dumps(self.as_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode("utf-8")).hexdigest()


DEFAULTS = Settings()  # every default, for the command-line help and its resolution


def derive_stage_seed(master_seed: int, domain: str) -> int:
    """A stable unsigned 64-bit seed for one stage, derived from the master seed."""
    digest = hashlib.sha256(
        b"BSReadSim/stage-seed/v1\0" + master_seed.to_bytes(8, "little") + b"\0"
        + domain.encode("ascii")).digest()
    return int.from_bytes(digest[:8], "little")
