Various Ideas for Confounder Adjustment in Regression
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/dcgerard/vicar/workflows/R-CMD-check/badge.svg)](https://github.com/dcgerard/vicar/actions)
[![Codecov test
coverage](https://codecov.io/gh/dcgerard/vicar/branch/master/graph/badge.svg)](https://codecov.io/gh/dcgerard/vicar?branch=master)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Description

Let

Y = XB + ZA + E,

for

-   Y an n by p matrix of gene expression data with n samples and p
    genes,
-   X an n by q matrix of q covariates,
-   B a q by p matrix of unobserved coefficients for the observed
    covariates,
-   Z an n by k matrix of hidden confounders,
-   A a k by p matrix of hidden coefficients for the hidden confounders,
    and
-   E an n by p matrix of independent normal errors with column
    variances s1,…,sp.

Not accounting for the hidden covariates, Z, can reduce power and result
in poor control of false discovery rate. This package provides a suite
of functions to adjust for hidden confounders, both when one has and
does not have access to control genes.

The functions `mouthwash()` and `backwash()` can adjust for hidden
confounding when one does not have access to control genes. They do so
via non-parametric empirical Bayes methods that use the powerful
methodology of Adaptive Shrinkage (Stephens 2016) within the
factor-augmented regression framework described in Wang et al. (2017).
`backwash()` is a slightly more Bayesian version of `mouthwash()`. These
methods are described in Gerard and Stephens (2020).

When one has control genes, there are many approaches to take. Such
methods include RUV2 (J. A. Gagnon-Bartsch and Speed 2012), RUV4 (J.
Gagnon-Bartsch, Jacob, and Speed 2013), and CATE (Wang et al. 2017).
This package adds to the field of confounder adjustment with control
genes by

1.  Implementing a version of CATE that is calibrated using control
    genes similarly to the method in J. Gagnon-Bartsch, Jacob, and
    Speed (2013). The function is called `vruv4()`.
2.  Introduces RUV3, a version of RUV that can be considered both RUV2
    and RUV4. The function is called `ruv3()`.
3.  Introduces RUV-impute, a more general framework for accounting for
    hidden confounders in regression. The function is called
    `ruvimpute()`
4.  Introduces RUV-Bayes, a Bayesian version of RUV. The function is
    called `ruvb()`.

These additions are described in detail in Gerard and Stephens (2021).

See also the related R packages
[`cate`](https://cran.r-project.org/package=cate) (Wang and Zhao 2015)
and [`ruv`](https://cran.r-project.org/package=ruv) (J. Gagnon-Bartsch
2015).

Check out [NEWS.md](NEWS.md) to see what’s new with each version.

## How to cite

If you use any of the control-gene based methods, please cite:

> Gerard, D., & Stephens, M. 2021. “Unifying and Generalizing Methods
> for Removing Unwanted Variation Based on Negative Controls.”
> *Statistica Sinica*, 31(3), 1-22
> &lt;[doi:10.5705/ss.202018.0345](https://doi.org/10.5705/ss.202018.0345)&gt;.

Or, using BibTex:

``` tex
@article{gerard2021unifying,
  title={Unifying and Generalizing Methods for Removing Unwanted Variation Based on Negative Controls},
  author={Gerard, David and Stephens, Matthew},
  journal={Statistica Sinica},
  doi={10.5705/ss.202018.0345},
  volume={31},
  number={3},
  pages={1--22},
  year={2021}
}
```

If you use either MOUTHWASH or BACKWASH, please cite:

> Gerard, D., & Stephens, M. 2020. “Empirical Bayes shrinkage and false
> discovery rate estimation, allowing for unwanted variation,”
> *Biostatistics*, 21(1), 15-32
> &lt;[doi:10.1093/biostatistics/kxy029](https://doi.org/10.1093/biostatistics/kxy029)&gt;.

Or, using BibTex:

``` tex
@article{gerard2020empirical,
  author = {Gerard, David and Stephens, Matthew},
  title = "Empirical {B}ayes shrinkage and false discovery rate estimation, allowing for unwanted variation",
  journal = {Biostatistics},
  volume = {21},
  number = {1},
  pages = {15--32},
  year = {2020},
  issn = {1465-4644},
  doi = {10.1093/biostatistics/kxy029},
}
```

## Installation

To install, first install `sva` and `limma` from Bioconductor in R:

``` r
install.packages("BiocManager")
BiocManager::install(c("limma", "sva"))
```

Then run in R:

``` r
install.packages("devtools")
devtools::install_github("dcgerard/vicar")
```

If you want some of the tools in `vicar` to be exactly equivalent to
those in `ruv`, you’ll need to install an older version of `ruv` (`ruv`
was updated and now the those equivalencies are not *exactly* the same)

``` r
devtools::install_version("ruv", version = "0.9.6", repos = "http://cran.us.r-project.org")
```

A note about matrix computations in vicar: Some of the methods in the
vicar package such as mouthwash and backwash rely heavily on
matrix-vector operations. The speed of these operations can have a big
impact on vicar’s performance, especially in large-scale data sets. If
you are applying vicar to large data sets, I recommend that you set up R
with optimized BLAS (optionally, LAPACK) libraries, especially if you
have a multicore computer (most modern laptops and desktops are
multicore). See
[here](https://csgillespie.github.io/efficientR/set-up.html#blas-and-alternative-r-interpreters)
and
[here](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Linear-algebra)
for advice and technical details on this. For example, in [our
experiments on a high-performance compute
cluster](https://github.com/pcarbo/mouthwash_sims/blob/master/mouthwash.sbatch)
we set up R with multithreaded OpenBLAS.

## Vignettes

I’ve provided three vignettes to help you get started with vicar. By
default, the vignettes are not built when you use `install_github`. To
build the vignettes during installation, run

``` r
install.packages("devtools")
devtools::install_github("dcgerard/vicar", build_vignettes = TRUE)
```

Note that this will result in a somewhat slower install. The first
vignette, *sample\_analysis*, gives a sample analysis using vicar to
account for hidden confounding. The second vignette, *customFA*, gives a
few instructions on how to incorporate user-defined factor analyses with
the confounder adjustment procedures implemented in vicar. The third
vignette, *custom\_prior*, gives instructions and examples on
incorporating a user-specified prior into `ruvb`. To see these vignettes
after install, run

``` r
utils::vignette("sample_analysis", package = "vicar")
utils::vignette("customFA", package = "vicar")
utils::vignette("custom_prior", package = "vicar")
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bartsch2015ruv" class="csl-entry">

Gagnon-Bartsch, Johann. 2015. *<span class="nocase">ruv</span>: Detect
and Remove Unwanted Variation Using Negative Controls*.
<https://CRAN.R-project.org/package=ruv>.

</div>

<div id="ref-gagnon2012using" class="csl-entry">

Gagnon-Bartsch, Johann A, and Terence P Speed. 2012. “Using Control
Genes to Correct for Unwanted Variation in Microarray Data.”
*Biostatistics* 13 (3): 539–52.
<https://doi.org/10.1093/biostatistics/kxr034>.

</div>

<div id="ref-gagnon2013removing" class="csl-entry">

Gagnon-Bartsch, Johann, Laurent Jacob, and Terence Speed. 2013.
“Removing Unwanted Variation from High Dimensional Data with Negative
Controls.” Technical Report 820, Department of Statistics, University of
California, Berkeley. <http://statistics.berkeley.edu/tech-reports/820>.

</div>

<div id="ref-gerard2020empirical" class="csl-entry">

Gerard, David, and Matthew Stephens. 2020. “<span
class="nocase">Empirical Bayes shrinkage and false discovery rate
estimation, allowing for unwanted variation</span>.” *Biostatistics* 21
(1): 15–32. <https://doi.org/10.1093/biostatistics/kxy029>.

</div>

<div id="ref-gerard2021unifying" class="csl-entry">

———. 2021. “Unifying and Generalizing Methods for Removing Unwanted
Variation Based on Negative Controls.” *Statistica Sinica* 31 (3): 1–22.
<https://doi.org/10.5705/ss.202018.0345>.

</div>

<div id="ref-stephens2016false" class="csl-entry">

Stephens, Matthew. 2016. “<span class="nocase">False discovery rates: a
new deal</span>.” *Biostatistics* 18 (2): 275–94.
<https://doi.org/10.1093/biostatistics/kxw041>.

</div>

<div id="ref-wang2015cate" class="csl-entry">

Wang, Jingshu, and Qingyuan Zhao. 2015. *<span
class="nocase">cate</span>: High Dimensional Factor Analysis and
Confounder Adjusted Testing and Estimation*.
<https://CRAN.R-project.org/package=cate>.

</div>

<div id="ref-wang2017confounder" class="csl-entry">

Wang, Jingshu, Qingyuan Zhao, Trevor Hastie, and Art B. Owen. 2017.
“Confounder Adjustment in Multiple Hypothesis Testing.” *Ann. Statist.*
45 (5): 1863–94. <https://doi.org/10.1214/16-AOS1511>.

</div>

</div>
