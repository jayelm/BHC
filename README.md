# Predictive BHC

Rich Savage's BHC implementation, but including some routines to compute the
posterior predictive distribution.

Things to do to compute posterior predictive:

- Function for getting all $r_k$s of nodes (**Done**)
- Function for computing all $\omega_k$s of nodes (**Done**)
    - Can be done recursively, probably; look at $\omega_k$ formula (just sum of weights of nodes up to root)
- Re-calculate Dirichlet hyperparameters (**Done**)
- Function for getting indices of the data represented by each node (**Done**)
- Function for computing the posterior predictive distribution given a node's subset of data and hyperparameters;
    - More specifically, a function for returning the log-likelihood of a data point according to the posterior predictive distribution of a node (I don't need to do sampling)
    - Follow (and verify) math in https://people.eecs.berkeley.edu/~stephentu/writeups/dirichlet-conjugate-prior.pdf
        - Or find an R version on the internet somewhere???
        - Will require clarifying the prior parameters of the multinomial
- A toplevel function to make predictions from the tree by summing up over all
    nodes $\omega_k$

To not mess around with export, it might be worth defining some of these
functions within toplevel, export-ready functions (or at least composed into
very easy-to-understand pieces)

## Catches when transition to new dataset

- MultinomialDataset doesn't have many virtual methods of the abstract class
    Dataset implemented
- MultinomialDataset accepts a raw vector of ints; make it accept a pointer
    (both in the .h and .cpp files). Seems to work even then
- In the multinomial case, FindOptimalHyperparameter is ran, which runs RunBhcWrapper; make sure these functions take the numThreads/randomised arguments

## Assumptions made

- `logEvidence` is the log-odds ratio of $r_k$, since in `WriteOutClusterLabels`,
    the tree is cut where `logEvidence < 0`, which must be where $r_k < 0.5$.
    - Further evidence that `logEvience` is actually a *log evidence ratio* is that `logEvidence`s as constructed in `ConstructDendrogramObject.R` are assigned from `mergeWeight`, where a comment mentions that `mergeWeight`s are the log-evidence ratios (and that the math in bayeslink_bin checks out)

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
