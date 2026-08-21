# Methylation state model interface v1

Status: stable Python process-stage boundary.

The model consumes fixed MethDB site probabilities plus addressed identity and
returns one boolean state per site. The built-in model is
`BernoulliStateModel` with contract `bernoulli-site-v1`; each site is sampled
independently through the frozen Philox methylation-state domain.

`--methylation-model bilstm` is accepted as a forward-compatible request. The
current release emits a `RuntimeWarning` that the BiLSTM plugin is unavailable
and records an effective fallback to Bernoulli. It never silently claims that a
correlated model ran.

A future correlated model must implement the same typed interface but may
sample a fragment jointly or conditionally. It must declare a model contract,
artifact identity, deterministic addressing rule, and worker/chunk invariance;
relabeling independent Bernoulli draws as BiLSTM is forbidden.
