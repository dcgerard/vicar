
<!-- README.md is generated from README.Rmd. Please edit that file -->
VICAR: Variance Inflation for Confounder Adjustment in Regression
=================================================================

[![Build Status](https://travis-ci.org/dcgerard/vicar.svg?branch=master)](https://travis-ci.org/dcgerard/vicar)

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

There are many approaches available for estimating hidden confounders, but in practice they are not well calibrated. This package calibrates RUV2 (J. A. Gagnon-Bartsch and Speed 2012), RUV4 (J. Gagnon-Bartsch, Jacob, and Speed 2013), and LEAPP (Sun et al. 2012). If I come up with a way to calibrate SVA (Leek and Storey 2007, Leek and Storey (2008)), then I'll include that as well.

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

Gagnon-Bartsch, Johann A, and Terence P Speed. 2012. “Using Control Genes to Correct for Unwanted Variation in Microarray Data.” *Biostatistics* 13 (3). Biometrika Trust: 539–52.

Leek, Jeffrey T, and John D Storey. 2007. “Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis.” *PLoS Genet* 3 (9): 1724–35.

———. 2008. “A General Framework for Multiple Testing Dependence.” *Proceedings of the National Academy of Sciences* 105 (48). National Acad Sciences: 18718–23.

Sun, Yunting, Nancy R Zhang, Art B Owen, and others. 2012. “Multiple Hypothesis Testing Adjusted for Latent Variables, with an Application to the Agemap Gene Expression Data.” *The Annals of Applied Statistics* 6 (4). Institute of Mathematical Statistics: 1664–88.
