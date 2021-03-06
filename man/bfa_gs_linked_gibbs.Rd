% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bfa_gs_linked_gibbs}
\alias{bfa_gs_linked_gibbs}
\title{Fast Gibbs sampler for Bayesian factor analysis.}
\usage{
bfa_gs_linked_gibbs(Linit, Finit, xi_init, phi_init, zeta_init, Y22init,
  Y21, Y31, Y32, nsamp, burnin, thin, rho_0, alpha_0, beta_0, eta_0, tau_0,
  display_progress)
}
\arguments{
\item{Linit}{A numeric matrix. The initial values of the loadings.}

\item{Finit}{A numeric matrix. The initial values for the factors.}

\item{xi_init}{A numeric vector. The initial values of the precisions.}

\item{phi_init}{A numeric scalar. The initial value of the mean of the precisions.}

\item{zeta_init}{A numeric vector. The initial values of the augmented row precisions.}

\item{Y22init}{A matrix of numerics. The initial value of Y22.}

\item{Y21}{A matrix of numerics.}

\item{Y31}{A matrix of numerics.}

\item{Y32}{A matrix of numerics.}

\item{nsamp}{The number of iterations to run in the Gibbs sampler, not including the burnin.}

\item{burnin}{The number of iterations to burnin.}

\item{thin}{We only collect samples every \code{thin} iterations.}

\item{rho_0}{The prior sample size for column-specific the
precisions.}

\item{alpha_0}{The prior sample size for the mean of the
column-specific precisions.}

\item{beta_0}{The prior mean of the mean of the column-specific
precisions.}

\item{eta_0}{The prior sample size of the expanded parameters.}

\item{tau_0}{The prior mean of of the expanded parameters.}

\item{display_progress}{A logical. If \code{TRUE}, then a progress bar will
be displayed and you'll be able to interrupt the C++ code. If \code{FALSE},
then neither of these capabilities will be provided.}
}
\description{
This is very similar to \code{\link{bfa_gd_gibbs}} except that we link the precisions
of the observations with the precisions of the factors. For some reason, this works
very well in practice.
}
\author{
David Gerard
}
