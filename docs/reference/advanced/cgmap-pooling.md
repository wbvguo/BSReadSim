# CGmap context pooling

By default, a CGmap or bedMethyl value applies to the same genomic position.
`--cgmap-pool` instead treats the supplied values as empirical distributions
for CG, CHG, and CHH contexts.

```bash
bsreadsim run wgbs \
  --reference reference.fa \
  --output runs/pooled \
  --fragments 100000 \
  --cgmap sample.cgmap.gz \
  --cgmap-pool \
  --mutation-rate 0 \
  --seed-meth 7 \
  --seed 42
```

For each contig, every methylation site samples with replacement from the
values observed for its context class. C- and G-oriented forms share the same
class. CGmap rows whose value is `na` do not enter a pool; bedMethyl rows always
provide a value.

If a contig has no observed value for a context class, sites in that class use
the configured Beta fallback. At least one defined value must exist somewhere
in the input. When ASM is also supplied, its site-specific allele probabilities
take precedence over pooled values.

Pooling deliberately breaks the positional link between an input row and the
same genomic coordinate. Use it only when the input should represent a
context-specific distribution rather than a position-specific profile.
