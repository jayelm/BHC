# Predictive BHC

Rich Savage's BHC implementation, but including some routines to compute the
posterior predictive distribution.

TODO: Figure out how the hyperparameters are used.

## Assumptions made

- `logEvidence` is the log-odds ratio of $r_k$, since in `WriteOutClusterLabels`,
    the tree is cut where `logEvidence < 0`, which must be where $r_k < 0.5$.

## Chain of computation in multinomial case

- `bhc` finds the optimal Beta global hyperparameter using `FindOptimalHyperparameter` (later named `cc`), then passes all info into `RunBhcWrapper`.
- `RunBhcWrapper` runs the C++ function `bhcWrapper_multinomial` in `bhc.cpp`
    - This seems to be older code.
    - It seems like newer code to work with the updated `bhc` stuff was
    implemented in `MultinomialDataSet.cpp` but is not used
- `bhcWrapper_multinomial` hardwires the Dirichlet process hyperparameter `alp = 0.001` hyperparameter to 0.001 and passes this, and the Beta global hyperparameter `cc`, to `bayeslink_binf` in `multinomial_bayeslink_binf.cpp`.
- `bayeslink_binf` does a lot of things:
    - Calculates a 2D array of individual hyperparameters using `CalculateHyperparameters` in `multinomial_CalculateHyperparameters.cpp`. The rows are the number of discrete values of the multinomial (in the binary case, 2), and the columns are the number of columns in the original dataset. Some odd procedure is used for this, but has something to do with the global Beta hyperparameter `cc`.
    - Computes log-evidence for single datapoints using `binevidence` in
    `multinomial_binevidence.cpp`. This is used since there's only one
    observation?
        - `binevidence` appears to calculate the probability under the relevant
            hypothesis using Appendix B. Multinomial-Dirichlet
