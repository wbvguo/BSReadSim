# CGmap context pool v1

Status: normative optional methylation-level randomization contract.

`cgmap_pool=true` changes the C++ MethDB construction rule; it does not move
any methylation work into Python. It requires one verified CGmap or bedMethyl
input. Parser and reference checks remain those in
[cgmap-profile-v1.md](cgmap-profile-v1.md) or
[bed-methyl-profile-v1.md](bed-methyl-profile-v1.md).

## Pool construction

For each contig, C++ builds three ordered pools from normalized input rows with
a defined value: CG, CHG, and CHH. CGmap `na` rows do not enter a pool;
bedMethyl rows are always defined. C- and
G-oriented contexts share their context-class pool; this mapping is performed
with the six typed protocol context values. Values retain their binary32
representation and their increasing input position order.

When a class pool is non-empty, every MethDB site of that class samples one
value from it with replacement.  This includes reference-equivalent,
variant-changed, and inserted sites.  The resulting protocol source is
`POOLED_CGMAP`.  When the contig has no defined values for a site's class, that
site keeps its independently addressed Beta fallback and source `BETA`.
At least one defined value must exist in the complete input, so enabling the
option with an all-`na` CGmap source fails closed.

Pooling replaces the direct position-specific CGmap overlay.  If ASM is also
present, ASM is applied afterward and retains precedence at its exact typed
targets: `ASM > pooled CGmap > Beta`.

## Addressed sampling

The algorithm identifier is `cgmap-context-pool-v1`.  For each site:

```text
key              = derive_key(seed, Stage::methylation_level, contig_index)
entity_ordinal   = the frozen uint64 SiteEntity
local_index      = UINT64_MAX
upper_exclusive  = the uint32 number of values in the context-class pool
index            = bounded_integer(key, entity_ordinal,
                                   local_index, upper_exclusive)
```

Reference sites use their zero-based `uint32` coordinate as the tag-zero
`SiteEntity`; variant and insertion sites use the existing tagged `uint64`
identity.  `UINT64_MAX` is outside every address consumed by the released Beta
sampler's bounded attempt layout and is reserved for this selection.  No
mutable RNG state, call-order dependence, or wider coordinate/counter type is
introduced.

## Width and ownership boundary

CGmap positions, per-contig row counts, and each pool size are checked
`uint32_t`.  Site identities, global counts, and Philox counter halves remain
`uint64_t`.  The component does not construct `posidx`, pack `rseq`, reinterpret
`HaplotypeMask`, or introduce a u96 representation.  Only one contig's three
value pools are resident while its MethDB is built.
