# Predictive BHC

Rich Savage's BHC implementation, but including some routines to compute the
posterior predictive distribution.

## Assumptions made

- `logEvidence` is the log-odds ratio of $r_k$, since the tree is cut in
    `WriteOutClusterLabels` where `logEvidence < 0`, i.e. where $r_k < 0.5$.
