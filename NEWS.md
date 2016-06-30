# vicar 0.1.1

* Added `vruv2`, a variance-inflated version of RUV2, but it doesn't work too well yet.
* The main function for variance-inflated RUV4 is now `vruv4`. I thought that `vicarius_ruv4` was too verbose. In the future, as I create new calibrated versions of confounder adjustment methods, the function name will just have a "v" in front of the name of the confounder adjustment method.
* To get the standard errors of betahat, `vruv4` now multiplies the estimated variances by ([X, Z]'[X, Z])^{-1} rather than (X'X)^{-1}.
*  I plan on adding a separate function for variance inflation without confounder adjustment and removing future capabilities of using `vruv4` when `k = 0`. Keep this in mind.
* `limmashrink = TRUE` is now the default.
