Various Ideas for Confounder Adjustment in Regression
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Windows
Build](https://ci.appveyor.com/api/projects/status/github/dcgerard/vicar?branch=master&svg=true)](https://ci.appveyor.com/project/dcgerard/vicar)
[![Linux
Build](https://travis-ci.org/dcgerard/vicar.svg?branch=master)](https://travis-ci.org/dcgerard/vicar)
[![Coverage
Status](https://coveralls.io/repos/github/dcgerard/vicar/badge.svg?branch=master)](https://coveralls.io/github/dcgerard/vicar?branch=master)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Description

Let

Y = XB + ZA + E,

for

  - Y an n by p matrix of gene expression data with n samples and p
    genes,
  - X an n by q matrix of q covariates,
  - B a q by p matrix of unobserved coefficients for the observed
    covariates,
  - Z an n by k matrix of hidden confounders,
  - A a k by p matrix of hidden coefficients for the hidden confounders,
    and
  - E an n by p matrix of independent normal errors with column
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
methods are described in Gerard and Stephens (2018).

When one has control genes, there are many approaches to take. Such
methods include RUV2 (Gagnon-Bartsch and Speed 2012), RUV4
(Gagnon-Bartsch, Jacob, and Speed 2013), and CATE (Wang et al. 2017).
This package adds to the field of confounder adjustment with control
genes by

1.  Implementing a version of CATE that is calibrated using control
    genes similarly to the method in Gagnon-Bartsch, Jacob, and Speed
    (2013). The function is called `vruv4()`.
2.  Introduces RUV3, a version of RUV that can be considered both RUV2
    and RUV4. The function is called `ruv3()`.
3.  Introduces RUV-impute, a more general framework for accounting for
    hidden confounders in regression. The function is called
    `ruvimpute()`
4.  Introduces RUV-Bayes, a Bayesian version of RUV. The function is
    called `ruvb()`.

These additions are described in detail in Gerard and Stephens (2019).

See also the related R packages
[`cate`](https://cran.r-project.org/package=cate) (Wang and Zhao 2015)
and [`ruv`](https://cran.r-project.org/package=ruv) (Gagnon-Bartsch
2015).

Check out [NEWS.md](NEWS.md) to see what’s new with each version.

## How to cite

If you use any of the control-gene based methods, please cite:

> Gerard, D., & Stephens, M. 2019. “Unifying and Generalizing Methods
> for Removing Unwanted Variation Based on Negative Controls”.
> *Statistica Sinica*, in press
> \<[doi:10.5705/ss.202018.0345](https://doi.org/10.5705/ss.202018.0345)\>.

Or, using BibTex:

``` tex
@article{gerard2019unifying,
  title={Unifying and Generalizing Methods for Removing Unwanted Variation Based on Negative Controls},
  author={Gerard, David and Stephens, Matthew},
  journal={Statistica Sinica},
  doi={10.5705/ss.202018.0345},
  volume={},
  number={},
  year={2019},
  pages={in press}
}
```

If you use either MOUTHWASH or BACKWASH, please cite:

> Gerard, D., & Stephens, M. 2018. “Empirical Bayes shrinkage and false
> discovery rate estimation, allowing for unwanted variation”,
> *Biostatistics*, , kxy029,
> \<[doi:10.1093/biostatistics/kxy029](https://doi.org/10.1093/biostatistics/kxy029)\>.

Or, using BibTex:

``` tex
@article{gerard2018empirical,
  author = {Gerard, David and Stephens, Matthew},
  title = {Empirical {B}ayes shrinkage and false discovery rate estimation, allowing for unwanted variation},
  journal = {Biostatistics},
  volume = {},
  number = {},
  pages = {kxy029},
  year = {2018},
  doi = {10.1093/biostatistics/kxy029},
}
```

Otherwise, please cite this package as:

> Gerard, David. 2016. *vicar: Various Ideas for Confounder Adjustment
> in Regression*. <https://github.com/dcgerard/vicar>.

Or, using BibTex:

``` tex
@Manual{gerard2016vicar,
    title = {{vicar}: Various Ideas for Confounder Adjustment in Regression},
    author = {David Gerard},
    year = {2016},
    note = {R package version 0.1-9},
    url = {https://github.com/dcgerard/vicar},
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
was updated and now the those equivalencies are not *exactly* the
same)

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

<div id="refs" class="references">

<div id="ref-bartsch2015ruv">

Gagnon-Bartsch, Johann. 2015. *ruv: Detect and Remove Unwanted Variation
Using Negative Controls*. <https://CRAN.R-project.org/package=ruv>.

</div>

<div id="ref-gagnon2012using">

Gagnon-Bartsch, Johann A, and Terence P Speed. 2012. “Using Control
Genes to Correct for Unwanted Variation in Microarray Data.”
*Biostatistics* 13 (3). Biometrika Trust: 539–52.
<https://doi.org/10.1093/biostatistics/kxr034>.

</div>

<div id="ref-gagnon2013removing">

Gagnon-Bartsch, Johann, Laurent Jacob, and Terence Speed. 2013.
“Removing Unwanted Variation from High Dimensional Data with Negative
Controls.” Technical Report 820, Department of Statistics, University of
California, Berkeley. <http://statistics.berkeley.edu/tech-reports/820>.

</div>

<div id="ref-gerard2018empirical">

Gerard, David, and Matthew Stephens. 2018. “Empirical Bayes Shrinkage
and False Discovery Rate Estimation, Allowing for Unwanted Variation.”
*Biostatistics*, kxy029. <https://doi.org/10.1093/biostatistics/kxy029>.

</div>

<div id="ref-gerard2019unifying">

———. 2019. “Unifying and Generalizing Methods for Removing Unwanted
Variation Based on Negative Controls.” *Statistica Sinica*, in press.
<https://doi.org/10.5705/ss.202018.0345>.

</div>

<div id="ref-stephens2016false">

Stephens, Matthew. 2016. “False discovery rates: a new deal.”
*Biostatistics* 18 (2): 275–94.
<https://doi.org/10.1093/biostatistics/kxw041>.

</div>

<div id="ref-wang2015cate">

Wang, Jingshu, and Qingyuan Zhao. 2015. *cate: High Dimensional Factor
Analysis and Confounder Adjusted Testing and Estimation*.
<https://CRAN.R-project.org/package=cate>.

</div>

<div id="ref-wang2017confounder">

Wang, Jingshu, Qingyuan Zhao, Trevor Hastie, and Art B. Owen. 2017.
“Confounder Adjustment in Multiple Hypothesis Testing.” *Ann.
Statist.* 45 (5). The Institute of Mathematical Statistics: 1863–94.
<https://doi.org/10.1214/16-AOS1511>.

</div>

</div>
