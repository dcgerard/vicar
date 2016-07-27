
<!-- README.md is generated from README.Rmd. Please edit that file -->
VICAR: Variance Inflation for Confounder Adjustment in Regression
=================================================================

[![Build Status](https://travis-ci.org/dcgerard/vicar.svg?branch=master)](https://travis-ci.org/dcgerard/vicar) [![Coverage Status](https://coveralls.io/repos/github/dcgerard/vicar/badge.svg?branch=master)](https://coveralls.io/github/dcgerard/vicar?branch=master)

Description
-----------

Let

Y = XB + ZA + E,

for

-   Y an n by p matrix of gene expression data with n samples and p genes,
-   X an n by q matrix of of q covariates,
-   B a q by p matrix of unobserved coefficients for the observed covariates,
-   Z an n by k matrix of hidden confounders,
-   A a k by p matrix of hidden coefficients for the hidden confounders, and
-   E an n by p matrix of independent normal errors with column variances s1,...,sp.

Not accounting for the hidden covariates, Z, can reduce power and result in poor control of false discovery rate.

There are many approaches available for estimating hidden confounders. This package is among those that use control genes to account to hidden confounding. Other such methods include RUV2 (J. A. Gagnon-Bartsch and Speed 2012), RUV4 (J Gagnon-Bartsch, Jacob, and Speed 2013), and CATE (Wang et al. 2015). This package adds to this field by

1.  Implementing a version of CATE that is calibrated using control genes similarly to the method in J Gagnon-Bartsch, Jacob, and Speed (2013). The function is called `vruv4`.
2.  Introduces RUV3, a version of RUV that can be considered both RUV2 and RUV4. The function is called `ruv3`.
3.  Introduces RUV-impute, a more general framework for accounting for hidden confounders in regression. The function is called `ruvimpute`
4.  Introduces RUV-Bayes, a Bayesian version of RUV. The function is called `ruvb`.

Please cite this package as:

> Gerard, David. 2016. *vicar: Variance Inflation for Confounder Adjustment in Regression*. <https://github.com/dcgerard/vicar>.

Or, using BibTex:

``` tex
@Manual{gerard2016vicar,
    title = {{vicar}: Variance Inflation for Confounder Adjustment in Regression},
    author = {David Gerard},
    year = {2016},
    note = {R package version 0.1.2},
    url = {https://github.com/dcgerard/vicar},
}
```

See also the related R packages [`cate`](https://cran.r-project.org/web/packages/cate/index.html) (Wang and Zhao 2015) and [`ruv`](https://cran.r-project.org/web/packages/ruv/index.html) (Johann Gagnon-Bartsch 2015).

Check out [NEWS.md](NEWS.md) to see what's new with each version.

Installation
------------

To install, run in R:

``` r
install.packages("devtools")
devtools::install_github("dcgerard/vicar")
```

References
----------

Gagnon-Bartsch, J, L Jacob, and TP Speed. 2013. “Removing Unwanted Variation from High Dimensional Data with Negative Controls.” Technical Report 820, Department of Statistics, University of California, Berkeley.

Gagnon-Bartsch, Johann. 2015. *ruv: Detect and Remove Unwanted Variation Using Negative Controls*. <https://CRAN.R-project.org/package=ruv>.

Gagnon-Bartsch, Johann A, and Terence P Speed. 2012. “Using Control Genes to Correct for Unwanted Variation in Microarray Data.” *Biostatistics* 13 (3). Biometrika Trust: 539–52.

Wang, Jingshu, and Qingyuan Zhao. 2015. *cate: High Dimensional Factor Analysis and Confounder Adjusted Testing and Estimation*. <https://CRAN.R-project.org/package=cate>.

Wang, Jingshu, Qingyuan Zhao, Trevor Hastie, and Art B Owen. 2015. “Confounder Adjustment in Multiple Hypotheses Testing.” *ArXiv Preprint ArXiv:1508.04178*.
