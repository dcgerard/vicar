# vicar 0.1.1

* Added `r vruv2`, a variance-inflated version of RUV2.
* The main function for variance-inflated RUV4 is now `r vruv4`. I though that `r vicarius_ruv4` was too verbose. In the future, as I create new calibrated versions of confounder adjustment methods, the function name will just have a "v" in front of the name of the confounder adjustment method.
* `r vruv4` no longer works when `k = 0`. I plan on adding a separate function for variance inflation without confounder adjustment.
* `r limmashrink = TRUE` is now the default.
