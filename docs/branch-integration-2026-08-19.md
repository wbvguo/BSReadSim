# Branch integration ledger (2026-08-19)

Status: merge record for development branches reconciled into `dev` after the
BSReadSim 0.3.0 consolidation.

## Closure policy

Every branch created from the development line must eventually become an
ancestor of `dev`. A closure merge records that ancestry only after every
branch-specific change has one explicit disposition:

- `ported`: adapted to the current runtime and retained in the active tree;
- `superseded`: the current runtime already provides the intended contract in
  a different implementation;
- `held/rejected`: retained in Git history but not promoted because it conflicts
  with a released contract or did not pass its own correctness/performance gate.

A closure merge does not claim that every experiment became production code.
It guarantees that no topic history is orphaned and that rejected work remains
auditable from `dev`.

## `backup/pre-squash-2026-08-19`

Disposition: `ported` as history.

The branch tip and the squashed 0.3.0 `dev` tip had identical Git trees. It was
merged without content changes to reconnect the complete pre-squash ancestry
before either topic branch was closed.

## `codex/variable-insert-gc-profile`

Disposition: mixed, with the released feature ported.

| Source commit | Disposition | Current integration |
| --- | --- | --- |
| `44fa4d9`, `60a05c2`, `289e0fd`, `4921c47` | superseded | Their fixed-insert GC audit, sampler, and paired benchmark intent was already present in the replayed 0.3.0 history. |
| `825555f` | ported | Variable-insert GC rejection sampler. |
| `e95c2ad` | ported | Approximate mean-insert calibration integrated with the current WGBS/RRBS/TBS capability boundary. |
| `2c6fcce` | ported | Variable-insert target validation, audit behavior, and contract. |
| `5650f56` | ported | Benchmark requires the native FASTQ acceleration path. |
| `05e2957` | held as historical evidence | Results identify their old source tree and are not release evidence for the reconciled 0.3.x tree. The closure merge retains them in history; current evidence must be regenerated. |

The fixed-insert path remains the exact target-distribution contract. Variable
inserts use the documented mean-insert plug-in approximation, report unreachable
target mass, and may reweight the emitted insert-length distribution.

## `codex/paper-aligned-core`

Disposition: mixed, with the scientific correction ported and experimental
representations retained as auditable history.

| Source commits | Disposition | Reason |
| --- | --- | --- |
| `d96e200` | ported | Conversion now occurs once per physical fragment before mate derivation. Overlapping mates share template-offset conversion draws. |
| `a24592e` | superseded | The current sequencing-model implementation already samples exact cumulative integer weights with right-biased cumulative search. The incompatible probability-artifact v2 is not silently substituted for the released v1 artifact. |
| `0aeaa9c`, `a2f9f3f`, `8875b65` | held/rejected | These commits apply to the branch-only packed MethDB record. The current sparse catalog uses a checked 16-byte typed record; narrowing code from the absent representation is not copied. |
| `f35d916`, `c05a713`, `7c357f5`, `7a22267`, `5d113c0`, `5e2a8a3`, `63e5951` | held/rejected | The complete `uint16_t` haplotype arrays and locality path increased memory and reduced throughput against that branch's sparse baseline. Current typed sparse haplotype projection remains active. |
| `0c4fdda` | superseded | Current production FASTQ keeps the third line exactly `+`; lossless Truth is versioned JSONL or truth BAM. Inline FASTQ annotations would change the released output contract. |
| benchmark and conformance commits | held as historical evidence | Their measurements remain tied to exact historical binaries. They are preserved by ancestry but are not current 0.3.x release evidence. |

The paper branch's non-model ownership boundary remains represented in the
current runtime: C++ owns reference, variants, haplotype projection, MethDB, and
fragment generation; Python owns fragment-level conversion, quality, errors,
and output. The closure merge does not promote unmeasured historical bit layouts
over the current typed boundary.

## Required release gates

Before the reconciled tree is described as release-ready:

1. build the bundled-HTSlib C++ target and pass CTest;
2. pass the Python unit and installed-package checks;
3. prove worker/chunk invariance and unchanged fixed-insert target-GC output;
4. validate shared conversion events for overlapping mates;
5. rerun variable-insert GC accuracy, insert-drift, and paired throughput
   measurements on the reconciled source and record exact source/binary hashes.

## Final consolidation and current evidence

The accepted runtime and documented branch dispositions were consolidated
locally into one development commit over the local `origin/dev` snapshot. The
topic worktrees and local topic branches were removed after their accepted
content, rejected decisions, and measurements were represented in the final
tree. This intentionally replaces topic-merge ancestry with one reviewable
development delta; no remote branch was rewritten or pushed.

Current one-million-read evidence is retained under `data/benchmarks/` for
WGBS producer/FASTQ components, production/debug/truth-BAM output policies,
and balanced WGBS/RRBS/TBS synthetic-corpus throughput. Reports bind the
pre-squash runtime commit, core and native-extension digests, input digests,
counts, output hashes, individual timings, and limitations. The pre-squash
runtime tree is behaviorally identical to the consolidated runtime; only the
benchmark harness, evidence, tests, and documentation were added afterward.
