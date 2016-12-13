# vicar 0.1.6

This version mostly changes the documentation. Additions include:

* Many examples in the main functions in vicar. Including examples in `vruv4`, `ruv3`, `ruvb`, `mouthwash`, and `backwash`.
* Three new vignettes. One providing instructions and examples on how to customize your factor analysis, one providing instructions and
examples on how to customize your prior specification in `ruvb`, and one giving a sample analysis using the functions in vicar and other confounder
adjustment packages. Type `vignette(package = "vicar")` for the list of available vignettes.
* New functions `fa_tester` and `fa_tester_ruvb` for testing whether a user-specified factor analysis is appropriate for the functions in vicar.
* An example simulated dataset, "sim_gtex", based on the characteristics of the GTEx data. This is so that you can test our your confounder adjustment methods on some good test data. Type `data(sim_gtex)` to access the data or `?sim_gtex` for seeing details about the data.
* I've removed `vruv2` from being exported as it was supplanted by `ruv3`.

# vicar 0.1.5

This version added the function `backwash`. This is very similar in spirit to
`mouthwash`, except that rather than estimate the confounders by maximum
likelihood, `backwash` does so using a more Bayesian approach. `backwash`
returns a variational approximation to the posterior.

# vicar 0.1.4

This version added the function `mouthwash` to adjust for hidden
confounding when one does not have control genes. It applies the same
prior from [ashr](https://github.com/stephens999/ashr/tree/general) to
a factor-augmented regression framework.

# vicar 0.1.3

A lot of changes have occurred since my last news update. The biggest changes are:

* RUV3 is a method that can be considered both a version of RUV2 and RUV4. I implemented this in the function `ruv3`.
* `ruvimpute` is a generic function for using matrix imputation for confounder adjustment.
* `ruvb` is a special Bayesian version of RUV. It is highly customize-able, as you can tweak the Bayesian factor analysis and the prior specifications of the effects.
* I no longer recommend `vruv2` as this is now subsumed by `ruv3`. I'll probably remove `vruv2` in the future.

I provide reasonable defaults for all new methods.

# vicar 0.1.2

* `vruv2` now works pretty well and is recommended for use. This is a very different way to do variance inflation in RUV2 than what was previously implemented.
* The previous implementation is now in the function `vruv2_old`, but it may be removed at any time.
* I included `ash_ruv2` that is a wrapper for `vruv2` and `ashr::ash.workhorse`.
* Some new factor analyses are available under the hood, but none of them are recommended for general use: `pca_ruv2`, `qmle_ruv2`, and `pca_2step`. In the future, I plan on only saving `pca_2step`.

# vicar 0.1.1

* Added `vruv2`, a variance-inflated version of RUV2, but it doesn't work too well yet.
* The main function for variance-inflated RUV4 is now `vruv4`. I thought that `vicarius_ruv4` was too verbose. In the future, as I create new calibrated versions of confounder adjustment methods, the function name will just have a "v" in front of the name of the confounder adjustment method.
* To get the standard errors of betahat, `vruv4` now multiplies the estimated variances by ([X, Z]'[X, Z])^{-1} rather than (X'X)^{-1}.
*  I plan on adding a separate function for variance inflation without confounder adjustment and removing future capabilities of using `vruv4` when `k = 0`. Keep this in mind.
* `limmashrink = TRUE` is now the default.
