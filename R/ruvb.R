#' Bayesian version of Removing Unwanted Variation.
#'
#' This function will take posterior draws from a Bayesian factor
#' analysis (or, more generally, a Bayesian imputation approach) to
#' propagate uncertainty when adjusting for unwanted variation. Useful
#' summaries are returned, such as local false sign rates and
#' posterior means.
#'
#' The model is \deqn{Y = XB + ZA + E,} where \eqn{Y} is a matrix of
#' responses (e.g. log-transformed gene expression levels), \eqn{X} is
#' a matrix of covariates, \eqn{B} is a matrix of coefficients,
#' \eqn{Z} is a matrix of unobserved confounders, \eqn{A} is a matrix
#' of unobserved coefficients of the unobserved confounders, and
#' \eqn{E} is the noise matrix where the elements are independent
#' Gaussian and each column shares a common variance. The rows of
#' \eqn{Y} are the observations (e.g. individuals) and the columns of
#' \eqn{Y} are the response variables (e.g. genes).
#'
#' I have three versions of Bayesian factor analyses that I
#' recommend. The first is \code{\link{bfa_gs_linked}}. This version
#' links the variances between the factors and observations and is the
#' version used in Gerard and Stephens (2016). This version appears to
#' work the best in practice and is thus the default. The second,
#' \code{\link{bfa_gs}}, is the same as the first except it does not
#' link the variances between the factors and the observations. The
#' last is \code{bfa_wrapper}, which is just a wrapper for the R
#' package bfa. The main thing about this version is that they do not
#' use a hierarchical prior on the variances.
#'
#' The user can specify their own Bayesian factor analysis (or
#' Bayesian model for missing data) using the \code{fa_func} and
#' \code{fa_args} parameters. To see instructions and examples on how
#' to do this, type the following in R:
#' \code{utils::vignette("customFA", package = "vicar")}. If you see
#' an error in the above code, then this probably means that the
#' vignettes were not built during installation. In which case, see
#' \url{https://github.com/dcgerard/vicar#vignettes}.
#'
#' The user can also specify their own priors for the second step of
#' RUVB. To do so, use the parameters \code{prior_fun} and
#' \code{prior_args}. To see instructions and an example on how to do
#' this, run the following code in R:
#' \code{utils::vignette("custom_prior", package = "vicar")}.  Again,
#' if you see an error in the above code then you probably need to
#' build the vignettes. Go to
#' \url{https://github.com/dcgerard/vicar#vignettes} for
#' instructions. If a prior is not specified, then the default is to
#' use a non-informative uniform prior. Though improper, using this
#' prior will result in a proper posterior no matter the model for the
#' unwanted variation.
#'
#' @inheritParams vruv4
#' @param fa_func A function that takes as input matrices named
#'     \code{Y21}, \code{Y31}, \code{Y32}, and \code{k} and returns a
#'     list, one of whose elements is called \code{Y22_array}. See
#'     \code{\link{bfa_gs_linked}} for an example function. Type in R
#'     \code{utils::vignette("customFA", package = "vicar")} for
#'     instructions and examples.
#' @param fa_args A list of additional parameters to pass to
#'     \code{fa_func}.
#' @param return_mcmc A logical. Should we return the MCMC draws?
#' @param prior_fun A function. This should take as input a matrix
#'     called \code{beta_mat} and output a positive numeric (if
#'     \code{return_log = TRUE}) or any numeric (if \code{return_log =
#'     FALSE}). This is the prior density (or log-density) of
#'     \code{beta_mat}.  Additional arguments may be passed to
#'     \code{prior_fun} through the \code{prior_args} argument. The
#'     default is to use an improper non-informative uniform prior,
#'     which is guaranteed to result in a proper posterior no matter
#'     the model for the unwanted variation. Type in R
#'     \code{utils::vignette("custom_prior", package = "vicar")} for
#'     instructions and examples.
#' @param return_log A logical. Does \code{prior_fun} return the log
#'     of the density (\code{"TRUE"}) or not (\code{"FALSE"})? For
#'     numerical stability reasons, you should probably make
#'     \code{prior_fun} return the log of the density and set this to
#'     \code{"TRUE"}.
#' @param prior_args A list of arguments to pass to \code{prior_fun}.
#' @param pad_na A logical. Should the indexing of the posterior summaries
#'     be the same as the data (\code{TRUE}) or should the control genes
#'     be removed (\code{FALSE})?
#'
#' @return A list with with some or all of the following elements.
#'
#'     \code{means} The posterior means of the betas.
#'
#'     \code{sd} The posterior standard deviations of the betas.
#'
#'     \code{medians} The posterior medians of the betas
#'
#'     \code{upper} The posterior 97.5th percentile of the betas.
#'
#'     \code{lower} The posterior 2.5th percentile of the betas.
#'
#'     \code{lfsr1} The empirical local false sign rate. This just
#'     counts the number of betas that are less than 0 before
#'     calculating lfsr. For the most significant genes, you should probably
#'     use \code{lfsr2}.
#'
#'     \code{lfsr2} The normal approximation for local false sign
#'     rate. This approximates the posterior of each beta by a normal,
#'     then uses this approximation to calculate lfsr.
#'
#'     \code{t} The posterior means divided by the posterior standard
#'     deviations.
#'
#'     \code{svalues1} The svalues from lfsr1. For the most significant
#'     genes, you should probably use \code{svalues2}.
#'
#'     \code{svalues2} The svalues from lfsr2.
#'
#'     \code{betahat_post} An array of the posterior samples of the
#'     betas. Only returned if \code{return_mcmc} is \code{TRUE}.
#'
#'     \code{fa} The raw output from whatever factor analysis is
#'     used. Only returned if \code{return_mcmc} is \code{TRUE}.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{bfa_gs}}, \code{\link{bfl}}, and
#'     \code{bfa_wrapper} for implemented Bayesian factor analyses.
#'
#' @export
#'
#' @examples
#' library(vicar)
#'
#' ## Generate data and controls ---------------------------------------------
#' set.seed(345)
#' n <- 13
#' p <- 101
#' k <- 2
#' q <- 3
#' is_null       <- rep(FALSE, length = p)
#' is_null[1:51] <- TRUE
#' ctl           <- rep(FALSE, length = p)
#' ctl[1:37]     <- TRUE
#'
#' X <- matrix(stats::rnorm(n * q), nrow = n)
#' B <- matrix(stats::rnorm(q * p), nrow = q)
#' B[2, is_null] <- 0
#' Z <- X %*% matrix(stats::rnorm(q * k), nrow = q) +
#'      matrix(rnorm(n * k), nrow = n)
#' A <- matrix(stats::rnorm(k * p), nrow = k)
#' E <- matrix(stats::rnorm(n * p, sd = 1 / 2), nrow = n)
#' Y <- X %*% B + Z %*% A + E
#'
#' ## Fit RUVB ---------------------------------------------------------------
#' ## I use a much smaller number of samples than reccommended for time.
#' ruvbout <- ruvb(Y = Y, X = X, k = k, ctl = ctl, cov_of_interest = 2,
#'                 include_intercept = FALSE,
#'                 fa_args = list(nsamp = 1000))
#' ruvblfsr <- ruvbout$lfsr2
#'
#' ## Compare to CATE/RUV4 ---------------------------------------------------
#' ruv4out <- cate::cate.fit(Y = Y, X.primary = X[, 2, drop = FALSE],
#'                           X.nuis = X[, -2, drop = FALSE], r = k, fa.method = "pc",
#'                           adj.method = "nc", nc = ctl)
#' ruv4p <- ruv4out$beta.p.value
#' ruv4p[ctl] <- NA
#'
#' ## Plot ROC curves --------------------------------------------------------
#' order_ruv4 <- order(ruv4p, na.last = NA)
#' order_ruvb <- order(ruvblfsr, na.last = NA)
#'
#' nnull <- sum(is_null[!ctl])
#' nsig  <- sum(!is_null[!ctl])
#' fpr4 <- cumsum(is_null[order_ruv4]) / nnull
#' tpr4 <- cumsum(!is_null[order_ruv4]) / nsig
#' fprb <- cumsum(is_null[order_ruvb]) / nnull
#' tprb <- cumsum(!is_null[order_ruvb]) / nsig
#'
#' graphics::plot(fprb, tprb, type = "l", xlab = "False Positive Rate",
#'                ylab = "True Positive Rate", main = "ROC Curves",
#'                col = 3)
#' graphics::lines(fpr4, tpr4, col = 4, lty = 2)
#' graphics::legend("bottomright", legend = c("RUVB", "RUV4"), col = c(3, 4),
#'                  lty = c(1, 2))
#' graphics::abline(0, 1, lty = 2)
#'
ruvb <- function(Y, X, ctl, k = NULL, fa_func = bfa_gs_linked,
                 fa_args = list(), cov_of_interest = ncol(X),
                 include_intercept = TRUE, return_mcmc = FALSE,
                 prior_fun = NULL, prior_args = list(),
                 return_log = NULL, pad_na = TRUE) {

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.numeric(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::assert_that(is.numeric(X))
    assertthat::assert_that(is.vector(ctl))
    assertthat::assert_that(is.logical(ctl))
    assertthat::are_equal(ncol(Y), length(ctl))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(all(cov_of_interest >= 1 & cov_of_interest <= ncol(X)))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.function(fa_func))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(is.null(fa_args$Y21))
    assertthat::assert_that(is.null(fa_args$Y31))
    assertthat::assert_that(is.null(fa_args$Y32))
    assertthat::assert_that(is.list(prior_args))

    rotate_out <- rotate_model(Y = Y, X = X, k = k, cov_of_interest =
                               cov_of_interest, include_intercept =
                               include_intercept, limmashrink = FALSE,
                               do_factor = FALSE)

    k <- rotate_out$k

    Y21 <- rotate_out$Y2[, ctl, drop = FALSE]
    Y22 <- rotate_out$Y2[, !ctl, drop = FALSE]
    Y31 <- rotate_out$Y3[, ctl, drop = FALSE]
    Y32 <- rotate_out$Y3[, !ctl, drop = FALSE]
    R22 <- rotate_out$R22

    ncontrols <- sum(ctl)
    ncovariates <- length(cov_of_interest)


    fa_args$Y21 <- Y21
    fa_args$Y31 <- Y31
    fa_args$Y32 <- Y32
    fa_args$k   <- k
    faout <- do.call(what = fa_func, args = fa_args)

    R22inv <- backsolve(R22, diag(nrow(R22)))

    betahat_post <- array(NA, dim = dim(faout$Y22_array))
    for (index in 1:dim(betahat_post)[3]) {
        betahat_post[, , index] <- R22inv %*% (Y22 - faout$Y22_array[, , index])
    }


    ## Create posterior summaries -----------------------------------------------
    return_list <- list()
    if (is.null(prior_fun)) {
        return_list$means   <- apply(betahat_post, c(1, 2), mean)
        return_list$sd      <- apply(betahat_post, c(1, 2), stats::sd)
        return_list$medians <- apply(betahat_post, c(1, 2), stats::median)
        return_list$upper   <- apply(betahat_post, c(1, 2), stats::quantile, c(0.975))
        return_list$lower   <- apply(betahat_post, c(1, 2), stats::quantile, c(0.025))
        return_list$lfsr1   <- apply(betahat_post, c(1, 2), calc_lfsr)
        pless               <- stats::pnorm(q = 0, mean = return_list$means, sd = return_list$sd)
        return_list$lfsr2   <- pmin(pless, 1 - pless)
        return_list$t       <- return_list$means / return_list$sd
    } else if (is.null(return_log)) {
        stop("if prior_fun is not NULL, then return_log needs to be a logical")
    } else {
        if (return_log) {
            lg_vec <- apply(X = betahat_post, MARGIN = 3, FUN = prior_fun_wrapper,
                            prior_fun = prior_fun, prior_args = prior_args)
            g_vec  <- exp(lg_vec - max(lg_vec))
            g_vec  <- g_vec / sum(g_vec)
        } else {
            g_vec <- apply(X = betahat_post, MARGIN = 3, FUN = prior_fun_wrapper,
                           prior_fun = prior_fun, prior_args = prior_args)
            g_vec <- g_vec / sum(g_vec)
        }

        return_list$means   <- apply(betahat_post, c(1, 2), calc_mean_g, g = g_vec)
        return_list$sd      <- sqrt(apply(betahat_post, c(1, 2), calc_mean_g, g = g_vec, r = 2) -
                                    return_list$means ^ 2)
        return_list$medians <- apply(betahat_post, c(1, 2), calc_quantiles_g, g = g_vec,
                                     quant = 0.5)
        return_list$lower   <- apply(betahat_post, c(1, 2), calc_quantiles_g, g = g_vec,
                                     quant = 0.025)
        return_list$upper   <- apply(betahat_post, c(1, 2), calc_quantiles_g, g = g_vec,
                                     quant = 0.975)
        return_list$lfsr1   <- apply(betahat_post, c(1, 2), calc_lfsr_g, g = g_vec)
        pless               <- stats::pnorm(q = 0, mean = return_list$means, sd = return_list$sd)
        return_list$lfsr2   <- pmin(pless, 1 - pless)
        return_list$t       <- return_list$means / return_list$sd
    }

    lfsr1_order <- order(c(return_list$lfsr1))
    svalues1 <- matrix((cumsum(c(return_list$lfsr1)[lfsr1_order]) /
                        (1:(prod(dim(return_list$lfsr1)))))[order(lfsr1_order)],
                       nrow = nrow(return_list$lfsr1), ncol = ncol(return_list$lfsr1))
    return_list$svalues1 <- svalues1

    lfsr2_order <- order(c(return_list$lfsr2))
    svalues2 <- matrix((cumsum(c(return_list$lfsr2)[lfsr2_order]) /
                        (1:(prod(dim(return_list$lfsr2)))))[order(lfsr2_order)],
                       nrow = nrow(return_list$lfsr2), ncol = ncol(return_list$lfsr2))
    return_list$svalues2 <- svalues2


    ## pad with NA's ---------------------------------------------------------
    if (pad_na) {
      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$means
      return_list$means <- temp

      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$sd
      return_list$sd <- temp

      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$medians
      return_list$medians <- temp

      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$lower
      return_list$lower <- temp

      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$upper
      return_list$upper <- temp

      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$lfsr1
      return_list$lfsr1 <- temp

      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$lfsr2
      return_list$lfsr2 <- temp

      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$svalues1
      return_list$svalues1 <- temp

      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$svalues2
      return_list$svalues2 <- temp

      temp <- matrix(NA, nrow = length(cov_of_interest), ncol = length(ctl))
      temp[, !ctl] <- return_list$t
      return_list$t <- temp
    }

    ## Return MCMC output? --------------------------------------------------------------
    if (return_mcmc) {
        return_list$betahat_post <- betahat_post
        return_list$fa <- faout
    }

    class(return_list) <- "ruvb"

    return(return_list)
}

#' Wrapper for \code{prior_fun} so that can be called in apply.
#'
#' @param beta_mat A matrix
#' @param prior_fun A function.
#' @param prior_args A list of arguments
#'
#' @author David Gerard
prior_fun_wrapper <- function(beta_mat, prior_fun, prior_args = list()) {
    assertthat::assert_that(is.matrix(beta_mat))
    assertthat::assert_that(is.function(prior_fun))
    assertthat::assert_that(is.list(prior_args))
    prior_args$beta_mat <- beta_mat
    return(do.call(what = prior_fun, args = prior_args))
}


#' Hierarchical prior density function as described in Gerard and Stephens (2016)
#'
#' @param beta_mat A matrix. The rows are the coefficients of the
#'     difference covariates. The columns are the different genes.
#' @param shape_param A positive numeric. The shape parameter for the
#'     gamma prior.
#' @param rate_param A positive numeric. The rate parameter for the
#'     gamma prior.
#' @param return_log A logical. Should we return the log-density
#'     (\code{"TRUE"}) or the the density (\code{"FALSE"})?
#'
#' @author David Gerard
#'
#' @export
hier_fun <- function(beta_mat, shape_param = 1, rate_param = 1, return_log = TRUE) {
    if (!is.matrix(beta_mat)) {
        beta_mat <- matrix(beta_mat, nrow = 1)
    }

    pstar <- ncol(beta_mat)
    sample_means     <- rowMeans(beta_mat)
    sample_variances <- apply(beta_mat, 1, stats::var)
    inner_pow <- pstar / (pstar + 1) * sample_means ^ 2 + (pstar - 1) * sample_variances +
                                                          shape_param * rate_param
    llike <- lgamma((pstar + shape_param) / 2) + shape_param * log(2) / 2 -
        (pstar + shape_param) * log(inner_pow) / 2 - pstar * log(pi) / 2 - log(pstar + 1) / 2
    if(return_log) {
        return(sum(llike))
    } else {
        return(exp(sum(llike)))
    }
}

#' A basic normal prior density function.
#'
#' @param beta_mat A matrix. The rows are the coefficients of the
#'     difference covariates. The columns are the different genes.
#' @param prior_mean A numeric. The prior mean.
#' @param prior_variance A positive numeric. The prior variance.
#' @param return_log A logical. Should we return the log-density
#'     (\code{"TRUE"}) or the the density (\code{"FALSE"})?
#'
#' @author David Gerard
#'
#' @export
normal_prior <- function(beta_mat, prior_mean = 0, prior_variance = 100, return_log = TRUE) {
    if (!is.matrix(beta_mat)) {
        beta_mat <- matrix(beta_mat, nrow = 1)
    }
    pstar <- ncol(beta_mat)

    llike <- sum(stats::dnorm(beta_mat, mean = prior_mean, sd = sqrt(prior_variance), log = TRUE))
    if (return_log) {
        return(llike)
    } else {
        return(exp(llike))
    }
}


#' Empirical estimate of lfsr based on posterior samples.
#'
#' @param y A vector of posterior draws.
#'
#' @author David Gerard
calc_lfsr <- function(y) {
    nless <- sum(y < 0)
    lfsr <- min(nless, length(y) - nless) / length(y)
    return(lfsr)
}

#' Same as \code{\link{calc_lfsr}} except with a prior specification.
#'
#' @param y A vector of posterior draws.
#' @param g A vector of priors.
#'
#' @author David Gerard
calc_lfsr_g <- function(y, g) {
    which_less <- y < 0
    pjk <- sum(g[which_less]) / sum(g)
    lfsr <- min(pjk, 1 - pjk)
    return(lfsr)
}

#' Calculate moments of pointmass rv's.
#'
#' @inheritParams calc_lfsr_g
#' @param r A positive numeric. The moment to calculate.
#'
#' @author David Gerard
calc_mean_g <- function(y, g, r = 1) {
    mean_val <- sum(g * (y ^ r)) / sum(g)
}

#' Calculate quantiles of pointmass rv's.
#'
#' @inheritParams calc_lfsr_g
#' @param quant The quantile to calculate.
#'
#' @author David Gerard
calc_quantiles_g <- function(y, g, quant = 0.5) {
    qval <- y[order(y)][max(which(cumsum(g[order(y)]) / sum(g) < quant))]
    return(qval)
}

#' Gibbs sampler for Bayesian SVD.
#'
#' This is a modification of the Bayesian approach from Hoff (2007) to
#' allow for heteroscedastic columns. We start the missing values from
#' the RUV4 solution.
#'
#' The rejection sampler in \code{\link[rstiefel]{rbmf.matrix.gibbs}}
#' almost always stalls, so I am removing it from exports for now.
#'
#' @author David Gerard
#'
#' @inheritParams em_miss
#' @param nsamp A positive integer. The number of samples to draw.
#' @param burnin A positive integer. The number of early samples to
#'     discard.
#' @param keep A positive integer. We will same the updates of
#'     \code{Y22} every \code{keep} iteration of the Gibbs sampler.
#' @param print_update A logical. Should we print a text progress bar
#'     to keep track of the Gibbs sampler (\code{TRUE}) or not
#'     (\code{FALSE})?
#' @param plot_update A logical. Should we make some plots to keep
#'     track of the Gibbs sampler (\code{TRUE}) or not (\code{FALSE})?
#'
#' @return A list with the following elements:
#'
#'     \code{Y22_array} A three-dimensional array containing draws of
#'     Y22. \code{Y22[, , i]} contains the \eqn{i}th draw of Y22.
#'
#'     \code{mu_psi_phi} A matrix with three columns. The rows are the
#'     draws from the posterior. The first column is for the mean of
#'     the singular values. The second column is for the precision of
#'     the singular values. The last column is for the mean of the
#'     variances.
#'
#'     \code{neff_y22} The effective sample sizes from
#'     \code{Y22_array} as calculated by
#'     \code{\link[coda]{effectiveSize}}
#'
#'     \code{neff_mu_psi_phi} The effective sample sizes from
#'     \code{mu_psi_phi} as calculated by
#'     \code{\link[coda]{effectiveSize}}
#'
#'
#' @references Hoff, P. D. (2012). Model averaging and dimension
#'     selection for the singular value decomposition. Journal of the
#'     American Statistical Association.
bsvd <- function(Y21, Y31, Y32, k, nsamp = 10000,
                 burnin = round(nsamp / 4), keep = 20,
                 print_update = TRUE,
                 plot_update = FALSE) {

    if (!requireNamespace("rstiefel", quietly = TRUE)) {
        stop("rstiefel needs to be installed to run bsvd. To install, run in R:\n    install.packages(\"rstiefel\")")
    }

    assertthat::are_equal(ncol(Y21), ncol(Y31))
    assertthat::are_equal(nrow(Y31), nrow(Y32))
    ncovs <- nrow(Y21)
    ncontrols <- ncol(Y21)
    n <- nrow(Y21) + nrow(Y31)
    p <- ncol(Y31) + ncol(Y32)

    ## Get initial values and set hyper parameters----------------------
    pcout <- pca_naive(cbind(Y31, Y32), r = k)
    sig_diag <- pcout$sig_diag
    alpha <- t(pcout$alpha)
    Z3 <- pcout$Z
    shrunk_var <- limma::squeezeVar(sig_diag, df = n - ncovs - k)$var.post
    shrunk_var_c <- shrunk_var[1:ncontrols]
    alpha_c    <- alpha[, 1:ncontrols, drop = FALSE]
    Z2 <- tcrossprod(sweep(Y21, 2, 1 / shrunk_var_c, `*`), alpha_c) %*%
        solve(tcrossprod(sweep(alpha_c, 2, 1 / sqrt(shrunk_var_c), `*`)))
    Y22init <- Z2 %*% alpha[, (ncontrols + 1):p, drop = FALSE]
    Zinit <- rbind(Z2, Z3)
    svout <- svd(Zinit %*% alpha, nv = k, nu = k)
    u_init <- svout$u
    delta_init <- svout$d[1:k]
    v_init <- svout$v


    rho_0   <- 2
    alpha_0 <- 2
    eta_0   <- 2
    mu_0    <- mean(delta_init)
    nu_0    <- 0.01
    beta_0  <- 1
    tau_0   <- 1

    mu_init  <- mu_0
    xi_init  <- 1 / shrunk_var
    phi_init <- mean(xi_init)
    psi_init <- 1 / stats::var(delta_init)


    xi_current    <- xi_init
    mu_current    <- mu_init
    phi_current   <- phi_init
    psi_current   <- psi_init
    u_current     <- u_init
    delta_current <- delta_init
    v_current     <- v_init
    Y22_current   <- Y22init
    Y_current     <- rbind(cbind(Y21, Y22_current), cbind(Y31, Y32))
    Yxi           <- sweep(Y_current, 2, xi_current, `*`)

    ## Run the Gibbs sampler --------------------------------------------------

    Y22_array <- array(NA, dim = c(ncovs, p - ncontrols, floor((nsamp - burnin) / keep)))
    keep_index <- 1
    mu_psi_phi <- matrix(NA, nrow = floor((nsamp - burnin) / keep), ncol = 3)
    plot_iters <- round(seq(0.01, 1, by = 0.01) * nsamp)
    if (print_update) {
        cat("Progress:\n")
        pb <- utils::txtProgressBar(style = 3)
    }
    for (gindex in 1:nsamp) {

        if (print_update) {
            utils::setTxtProgressBar(pb = pb, value = gindex / nsamp)
        }


        ## skip the first few iterations so that we can get to something more stable.
        ## Update U --------------------------------------
        Cmat <- Yxi %*% sweep(v_current, 2, delta_current, `*`)
        u_current <- rstiefel::rmf.matrix.gibbs(M = Cmat, X = u_current)

        ## Update V --------------------------------------
        Cmat <- crossprod(Yxi, sweep(u_current, 2, delta_current, `*`))
        dmax <- max(delta_current ^ 2) / 2
        ximax <- max(xi_current)
        new_mult <- sqrt(abs(dmax) * ximax)

        v_current <- rstiefel::rbmf.matrix.gibbs(A = diag(xi_current * (new_mult / ximax)),
                                                 B = diag(-(abs(delta_current) ^ 2) *
                                                          (new_mult / (2 * dmax))),
                                                 C = Cmat,
                                                 X = v_current)

        ## Update delta ----------------------------------
        for (rindex in 1:k) {
            uysigv <- c(crossprod(u_current[, rindex, drop = FALSE], Yxi) %*%
                        v_current[, rindex, drop = FALSE])
            vsigv <- c(crossprod(v_current[, rindex, drop = FALSE] * xi_current,
                                 v_current[, rindex, drop = FALSE]))
            dvar <- 1 / (vsigv + psi_current)
            dmean <- (uysigv + psi_current * mu_current) * dvar
            delta_current[rindex] <- stats::rnorm(n = 1, mean = dmean, sd = sqrt(dvar))
        }

        ## Update xi -------------------------------------
        theta_current <- tcrossprod(sweep(u_current, 2, delta_current, `*`),
                                    v_current)
        rvec       <- colSums((Y_current - theta_current) ^ 2)
        xi_shape   <- (n + rho_0) / 2
        xi_rate    <- (rvec + rho_0 * phi_current) / 2
        xi_current <- sapply(xi_rate, FUN = stats::rgamma, n = 1, shape = xi_shape)

        ## Update phi ------------------------------------
        phi_shape   <- (p * rho_0 + alpha_0) / 2
        phi_rate    <- (alpha_0 * beta_0 + rho_0 * sum(xi_current)) / 2
        phi_current <- stats::rgamma(n = 1, shape = phi_shape, rate = phi_rate)

        ## Update psi ------------------------------------
        psi_shape   <- (k + eta_0) / 2
        psi_rate    <- (eta_0 * tau_0 + sum((delta_current - mu_current) ^ 2)) / 2
        psi_current <- stats::rgamma(n = 1, shape = psi_shape, rate = psi_rate)

        ## Update mu -------------------------------------
        mu_var     <- 1 / (nu_0 + k * psi_current)
        mu_mean    <- (psi_current * sum(delta_current) + nu_0 * mu_0) * mu_var
        mu_current <- stats::rnorm(n = 1, mean = mu_mean, sd = sqrt(mu_var))

        ## Update Y22 ------------------------------------
        Y22_mean <- theta_current[1:ncovs, (ncontrols + 1):p, drop = FALSE]
        Y22_error <- sweep(matrix(stats::rnorm(n = ncovs * (p - ncontrols)), nrow = ncovs), 2,
                           1 / sqrt(xi_current[(ncontrols + 1):p]), `*`)
        Y22_current <- Y22_mean + Y22_error

        Y_current[1:ncovs, (ncontrols + 1):p] <- Y22_current

        if ((gindex - burnin) %% keep == 0 & gindex > burnin) {
            Y22_array[, , keep_index] <- Y22_current
            mu_psi_phi[keep_index, 1] <- mu_current
            mu_psi_phi[keep_index, 2] <- psi_current
            mu_psi_phi[keep_index, 3] <- phi_current
            keep_index <- keep_index + 1
        }

        if (gindex %in% plot_iters) {
            if (plot_update & gindex > burnin + keep) {
                graphics::par(mfrow = c(3, 1))
                graphics::plot(mu_psi_phi[, 1], type = "l", ylab = expression(mu))
                graphics::plot(mu_psi_phi[, 2], type = "l", ylab = expression(psi))
                graphics::plot(mu_psi_phi[, 3], type = "l", ylab = expression(phi))
                graphics::par(mfrow = c(1, 1))
            }
        }
    }
    if (print_update) {
        cat("\nComplete!\n")
    }

    nmpp <- apply(mu_psi_phi, 2, coda::effectiveSize)
    ny22 <- apply(Y22_array, c(1, 2), coda::effectiveSize)
    return(list(Y22_array = Y22_array, mu_psi_phi = mu_psi_phi, neff_y22 = ny22,
                neff_mu_psi_phi = nmpp))
}

#' Simple Bayesian low rank matrix decomposition.
#'
#' "bfl" = "Bayesian factor loading"
#'
#' This is as simple as they come. I put normal priors on the loadings
#' and factors and gamma priors on the precisions. The hyperparameters
#' are set to provide weak prior information by default.
#'
#' The main difference between this version and others is that the
#' factors are a prior assumed to have the same variances as the data
#' observations. This might be distasteful to some.
#'
#' This also has parameter expansion impelemented and is written in
#' compiled code. To see a slower version without paramaeter epansion,
#' go to \code{\link{bfl}}.
#'
#'
#' @inheritParams em_miss
#' @param nsamp A positive integer. The number of samples to draw.
#' @param burnin A positive integer. The number of early samples to
#'     discard.
#' @param thin A positive integer. We will same the updates of
#'     \code{Y22} every \code{keep} iteration of the Gibbs sampler.
#' @param display_progress A logical. Should we print a text progress
#'     bar to keep track of the Gibbs sampler (\code{TRUE}) or not
#'     (\code{FALSE})? Also, if \code{TRUE}, then you can interupt the
#'     C++ code every 1\% of runtime.
#' @param rho_0 A scalar. The prior "sample size" for the precisions.
#' @param alpha_0 A scalar. The prior "sample size" for the mean of
#'     the precisions.
#' @param beta_0 A scalar. The prior mean of the precisions.
#' @param eta_0 A scalar. The prior "sample size" for the parameter
#'     expansion.
#' @param tau_0 A scalar. The prior mean for the parameter expansion.
#'     matrix.
#' @param use_code A character. Should we use the C++ code
#'     (\code{"cpp"}) or the R code (\code{"r"})?
#'
#' @export
#'
#' @author David Gerard
#'
#' @seealso \code{\link{bfa_gs_linked}}.
#'
bfa_gs_linked <- function(Y21, Y31, Y32, k, nsamp = 10000,
                          burnin = round(nsamp / 4), thin = 10,
                          display_progress = TRUE,
                          rho_0 = 0.1, alpha_0 = 0.1,
                          beta_0 = 1, eta_0 = 1,
                          tau_0 = 1, use_code = c("r", "cpp")) {

    assertthat::are_equal(ncol(Y21), ncol(Y31))
    assertthat::are_equal(nrow(Y31), nrow(Y32))
    ncovs <- nrow(Y21)
    ncontrols <- ncol(Y21)
    n <- nrow(Y21) + nrow(Y31)
    p <- ncol(Y31) + ncol(Y32)
    use_code = match.arg(use_code)

    ## Get initial values and set hyper parameters----------------------
    pcout <- pca_naive(cbind(Y31, Y32), r = k)
    sig_diag <- pcout$sig_diag
    alpha <- t(pcout$alpha)
    Z3 <- pcout$Z
    shrunk_var <- limma::squeezeVar(sig_diag, df = n - ncovs - k)$var.post
    shrunk_var_c <- shrunk_var[1:ncontrols]
    alpha_c    <- alpha[, 1:ncontrols, drop = FALSE]
    Z2 <- tcrossprod(sweep(Y21, 2, 1 / shrunk_var_c, `*`), alpha_c) %*%
        solve(tcrossprod(sweep(alpha_c, 2, 1 / sqrt(shrunk_var_c), `*`)))
    Y22init <- Z2 %*% alpha[, (ncontrols + 1):p, drop = FALSE]
    Zinit <- rbind(Z2, Z3)
    fnorm_z <- sum(Zinit ^ 2)
    fnorm_alpha <- sum(alpha ^ 2)

    Linit <- Zinit
    Finit <- alpha
    xi_init  <- 1 / shrunk_var
    phi_init <- mean(xi_init)
    zeta_init <- 1 / limma::squeezeVar(colMeans(Linit ^ 2), df = k)$var.post

    if (use_code == "cpp") {
        bfout <- bfa_gs_linked_gibbs(Linit = Linit, Finit = Finit,
                                     xi_init = xi_init, phi_init = phi_init,
                                     zeta_init = zeta_init, Y22init = Y22init,
                                     Y21 = Y21, Y31 = Y31, Y32 = Y32,
                                     nsamp = nsamp, burnin = burnin,
                                     thin = thin, rho_0 = rho_0,
                                     alpha_0 = alpha_0, beta_0 = beta_0,
                                     eta_0 = eta_0, tau_0 = tau_0,
                                     display_progress = display_progress)
    } else if (use_code == "r") {
        bfout <- bfa_gs_linked_gibbs_r(Linit = Linit, Finit = Finit,
                                       xi_init = xi_init, phi_init = phi_init,
                                       zeta_init = zeta_init, Y22init = Y22init,
                                       Y21 = Y21, Y31 = Y31, Y32 = Y32,
                                       nsamp = nsamp, burnin = burnin,
                                       thin = thin, rho_0 = rho_0,
                                       alpha_0 = alpha_0, beta_0 = beta_0,
                                       eta_0 = eta_0, tau_0 = tau_0,
                                       display_progress = display_progress)
    }

    return(bfout)
}


#' R implementation of \code{\link{bfa_gs_linked_gibbs}}.
#'
#' @inheritParams bfa_gs_linked_gibbs
#'
#' @author David Gerard
#'
bfa_gs_linked_gibbs_r <- function(Linit, Finit, xi_init, phi_init,
                                  zeta_init, Y22init, Y21, Y31, Y32,
                                  nsamp, burnin, thin, rho_0, alpha_0,
                                  beta_0, eta_0, tau_0,
                                  display_progress) {

    ## Get dimensions of matrices
    n         <- nrow(Linit)
    p         <- ncol(Finit)
    nfac      <- ncol(Linit)
    ncovs     <- nrow(Y21)
    ncontrols <- ncol(Y21)

    ## Initialize matrices
    L_current   <- Linit
    F_current   <- Finit
    xi_current  <- xi_init
    phi_current <- phi_init
    zeta_current <- zeta_init
    Y22_current <- Y22init
    Y_current   <- rbind(cbind(Y21, Y22_current), cbind(Y31, Y32))


    ## Run the Gibbs sampler --------------------------------------------------
    nkeeps <- floor(nsamp / thin)

    Y22_array <- array(NA, dim = c(ncovs, p - ncontrols, nkeeps))
    thin_index <- 1
    phi_mat <- matrix(NA, nrow = nkeeps, ncol = 1)
    xi_mat <- matrix(NA, nrow = nkeeps, ncol = p)
    if (display_progress) {
        cat("Progress:\n")
        pb <- utils::txtProgressBar(style = 3)
    }
    for (gindex in 1:(nsamp + burnin)) {
        if (display_progress) {
            utils::setTxtProgressBar(pb = pb, value = gindex / (nsamp + burnin))
        }

        ## Update L ----------------------------------------------------
        Fsig <- sweep(F_current, 2, xi_current, `*`)
        eigen_fsf <- eigen(tcrossprod(Fsig, F_current) +
                           diag(zeta_current, nrow = nfac, ncol = nfac), symmetric = TRUE)

        L_meanmat <- tcrossprod(sweep(Y_current, 2, xi_current, `*`), F_current) %*%
            tcrossprod(sweep(eigen_fsf$vectors, 2, 1 / eigen_fsf$values, `*`),
                       eigen_fsf$vectors)

        col_cov_half <- tcrossprod(sweep(eigen_fsf$vectors, 2, 1 / sqrt(eigen_fsf$values), `*`),
                                   eigen_fsf$vectors)
        L_current <- L_meanmat +
            matrix(stats::rnorm(prod(dim(L_current))), nrow = nrow(L_current)) %*% col_cov_half

        ## Should equal
        ## Y_current %*% diag(xi_current) %*% t(F_current) %*%
        ##     solve(tcrossprod(Fsig, F_current) + diag(nfac))
        ## L_meanmat

        ## Update F -----------------------------------------------------------

        eigen_ll <- eigen(crossprod(L_current) + diag(nfac), symmetric = TRUE)


        F_meanmat <- tcrossprod(sweep(eigen_ll$vectors, 2, 1 / eigen_ll$values, `*`),
                                eigen_ll$vectors) %*% crossprod(L_current, Y_current)
        row_cov_half <- tcrossprod(sweep(eigen_ll$vectors, 2, 1 / sqrt(eigen_ll$values), `*`),
                                   eigen_ll$vectors)
        F_error <- row_cov_half %*% sweep(matrix(stats::rnorm(prod(dim(F_current))),
                                                 nrow = nrow(F_current)),
                                          2, 1 / sqrt(xi_current), `*`)
        F_current <- F_meanmat + F_error

        ## Update xi ----------------------------------------------------------
        theta_current <- L_current %*% F_current
        r_vec         <- colSums((Y_current - theta_current) ^ 2)
        u_vec         <- colSums(F_current ^ 2)
        xi_shape      <- (n + nfac + rho_0) / 2
        xi_rate       <- (r_vec + u_vec + rho_0 * phi_current) / 2
        xi_current    <- sapply(xi_rate, FUN = stats::rgamma, n = 1, shape = xi_shape)

        ## Update phi ---------------------------------------------------------
        phi_shape   <- (p * rho_0 + alpha_0) / 2
        phi_rate    <- (alpha_0 * beta_0 + rho_0 * sum(xi_current)) / 2
        phi_current <- stats::rgamma(n = 1, shape = phi_shape, rate = phi_rate)

        ## Update zeta --------------------------------------------------------
        s_vec        <- colSums(L_current ^ 2)
        zeta_shape   <- (n + eta_0) / 2
        zeta_rate    <- (s_vec + eta_0 * tau_0) / 2
        zeta_current <- sapply(zeta_rate, FUN = stats::rgamma, n = 1, shape = zeta_shape)

        ## Update Y22 ---------------------------------------------------------
        Y22_mean <- theta_current[1:ncovs, (ncontrols + 1):p, drop = FALSE]
        Y22_error <- sweep(matrix(stats::rnorm(n = ncovs * (p - ncontrols)), nrow = ncovs), 2,
                           1 / sqrt(xi_current[(ncontrols + 1):p]), `*`)
        Y22_current <- Y22_mean + Y22_error
        Y_current[1:ncovs, (ncontrols + 1):p] <- Y22_current

        if ((gindex - burnin) %% thin == 0 & gindex > burnin) {
            Y22_array[, , thin_index] <- Y22_current
            phi_mat[thin_index, 1] <- phi_current
            xi_mat[thin_index,] <- xi_current
            thin_index <- thin_index + 1
        }
    }
    if (display_progress) {
        cat("\nComplete!\n")
    }

    return(list(Y22_array = Y22_array, xi = xi_mat, phi = phi_mat))
}




#' Simple Bayesian low rank matrix decomposition.
#'
#' "bfl" = "Bayesian factor loading"
#'
#' This is as simple as they come. I put normal priors on the loadings
#' and factors and gamma priors on the precisions. The hyperparameters
#' are set to provide weak prior information by default.
#'
#' The main difference between this version and others is that the
#' factors are a prior assumed to have the same variances as the data
#' observations. This might be distasteful to some.
#'
#' There is no parameter expansion in this one. To see one with
#' parameter expansion, and a much faster version, see
#' \code{\link{bfa_gs_linked}}.
#'
#' @inheritParams em_miss
#' @param nsamp A positive integer. The number of samples to draw.
#' @param burnin A positive integer. The number of early samples to
#'     discard.
#' @param keep A positive integer. We will same the updates of
#'     \code{Y22} every \code{keep} iteration of the Gibbs sampler.
#' @param print_update A logical. Should we print a text progress bar
#'     to keep track of the Gibbs sampler (\code{TRUE}) or not
#'     (\code{FALSE})?
#' @param plot_update A logical. Should we make some plots to keep
#'     track of the Gibbs sampler (\code{TRUE}) or not (\code{FALSE})?
#' @param rho_0 A scalar. The prior "sample size" for the precisions.
#' @param alpha_0 A scalar. The prior "sample size" for the mean of the precisions.
#' @param beta_0 A scalar. The prior mean of the precisions.
#' @param eta_0 A scalar. The prior "sample size" for the scale of the mean matrix.
#' @param tau_0 A scalar. The prior mean of the scale of the mean matrix.
#'
#' @export
#'
#' @author David Gerard
#'
#' @seealso \code{\link{bfa_gs_linked}}.
#'
bfl <- function(Y21, Y31, Y32, k, nsamp = 10000,
                burnin = round(nsamp / 4), keep = 20,
                print_update = TRUE,
                plot_update = FALSE,
                rho_0 = 0.1, alpha_0 = 0.1,
                beta_0 = 1, eta_0 = 0.1,
                tau_0 = 1) {

    assertthat::are_equal(ncol(Y21), ncol(Y31))
    assertthat::are_equal(nrow(Y31), nrow(Y32))
    ncovs <- nrow(Y21)
    ncontrols <- ncol(Y21)
    n <- nrow(Y21) + nrow(Y31)
    p <- ncol(Y31) + ncol(Y32)

    ## Get initial values and set hyper parameters----------------------
    pcout <- pca_naive(cbind(Y31, Y32), r = k)
    sig_diag <- pcout$sig_diag
    alpha <- t(pcout$alpha)
    Z3 <- pcout$Z
    shrunk_var <- limma::squeezeVar(sig_diag, df = n - ncovs - k)$var.post
    shrunk_var_c <- shrunk_var[1:ncontrols]
    alpha_c    <- alpha[, 1:ncontrols, drop = FALSE]
    Z2 <- tcrossprod(sweep(Y21, 2, 1 / shrunk_var_c, `*`), alpha_c) %*%
        solve(tcrossprod(sweep(alpha_c, 2, 1 / sqrt(shrunk_var_c), `*`)))
    Y22init <- Z2 %*% alpha[, (ncontrols + 1):p, drop = FALSE]
    Zinit <- rbind(Z2, Z3)
    fnorm_z <- sum(Zinit ^ 2)
    fnorm_alpha <- sum(alpha ^ 2)
    Linit <- Zinit * (prod(dim(Zinit)) / fnorm_z)
    Finit <- alpha * (prod(dim(alpha)) / fnorm_alpha)
    psi_init <- fnorm_z * fnorm_alpha / (prod(dim(Zinit)) * prod(dim(alpha)))
    xi_init  <- 1 / shrunk_var
    phi_init <- mean(xi_init)

    xi_current    <- xi_init
    phi_current   <- phi_init
    psi_current   <- psi_init
    L_current     <- Linit
    F_current     <- Finit
    Y22_current   <- Y22init
    Y_current     <- rbind(cbind(Y21, Y22_current), cbind(Y31, Y32))


    ## Run the Gibbs sampler --------------------------------------------------
    Y22_array <- array(NA, dim = c(ncovs, p - ncontrols, floor((nsamp - burnin) / keep)))
    keep_index <- 1
    psi_phi <- matrix(NA, nrow = floor((nsamp - burnin) / keep), ncol = 2)
    xi_mat <- matrix(NA, nrow = floor((nsamp - burnin) / keep), ncol = p)
    plot_iters <- round(seq(0.01, 1, by = 0.01) * nsamp)
    if (print_update) {
        cat("Progress:\n")
        pb <- utils::txtProgressBar(style = 3)
    }
    for (gindex in 1:nsamp) {
        if (print_update) {
            utils::setTxtProgressBar(pb = pb, value = gindex / nsamp)
        }

        ## Update L ----------------------------------------------------
        Fsig <- sweep(F_current, 2, xi_current, `*`)
        eigen_fsf <- eigen(tcrossprod(Fsig, F_current) + diag(k), symmetric = TRUE)

        L_meanmat <- tcrossprod(sweep(Y_current, 2, xi_current, `*`), F_current) %*%
            tcrossprod(sweep(eigen_fsf$vectors, 2, 1 / eigen_fsf$values, `*`),
                       eigen_fsf$vectors)

        col_cov_half <- tcrossprod(sweep(eigen_fsf$vectors, 2, 1 / sqrt(eigen_fsf$values), `*`),
                                   eigen_fsf$vectors)
        L_current <- L_meanmat +
            matrix(stats::rnorm(prod(dim(L_current))), nrow = nrow(L_current)) %*% col_cov_half

        ## Should equal
        ## Y_current %*% diag(xi_current) %*% t(F_current) %*%
        ##     solve(tcrossprod(Fsig, F_current) + diag(k))
        ## L_meanmat

        ## Update F ----------------------------------------------------

        eigen_ll <- eigen(crossprod(L_current) + diag(psi_current, k), symmetric = TRUE)
        F_meanmat <- tcrossprod(sweep(eigen_ll$vectors, 2, 1 / eigen_ll$values, `*`),
                                eigen_ll$vectors) %*% crossprod(L_current, Y_current)
        row_cov_half <- tcrossprod(sweep(eigen_ll$vectors, 2, 1 / sqrt(eigen_ll$values), `*`),
                                   eigen_ll$vectors)
        F_error <- row_cov_half %*% sweep(matrix(stats::rnorm(prod(dim(F_current))),
                                                 nrow = nrow(F_current)),
                                          2, 1 / sqrt(xi_current), `*`)
        F_current <- F_meanmat + F_error

        ## Should be equal
        ## solve(t(L_current) %*% L_current + psi_current * diag(k)) %*% t(L_current) %*% Y_current
        ## F_meanmat

        ## Update xi -----------------------------------------------------
        theta_current <- L_current %*% F_current
        r_vec <- colSums((Y_current - theta_current) ^ 2)
        s_vec <- colSums(F_current ^ 2)
        xi_shape <- (n + k + rho_0) / 2
        xi_rate  <- (r_vec + s_vec + rho_0 * phi_current) / 2
        xi_current <- sapply(xi_rate, FUN = stats::rgamma, n = 1, shape = xi_shape)

        ## Update phi ----------------------------------------------------
        phi_shape   <- (p * rho_0 + alpha_0) / 2
        phi_rate    <- (alpha_0 * beta_0 + rho_0 * sum(xi_current)) / 2
        phi_current <- stats::rgamma(n = 1, shape = phi_shape, rate = phi_rate)

        ## Update psi ----------------------------------------------------
        psi_shape <- (p * k + eta_0) / 2
        psi_rate <- (eta_0 * tau_0 + sum(sweep(F_current, 2, sqrt(xi_current), `*`) ^ 2)) / 2
        ## sum(diag(F_current %*% diag(xi_current) %*% t(F_current)))
        psi_current <- stats::rgamma(n = 1, shape = psi_shape, rate = psi_rate)

        ## Update Y22 ------------------------------------
        Y22_mean <- theta_current[1:ncovs, (ncontrols + 1):p, drop = FALSE]
        Y22_error <- sweep(matrix(stats::rnorm(n = ncovs * (p - ncontrols)), nrow = ncovs), 2,
                           1 / sqrt(xi_current[(ncontrols + 1):p]), `*`)
        Y22_current <- Y22_mean + Y22_error
        Y_current[1:ncovs, (ncontrols + 1):p] <- Y22_current

        if ((gindex - burnin) %% keep == 0 & gindex > burnin) {
            Y22_array[, , keep_index] <- Y22_current
            psi_phi[keep_index, 1] <- psi_current
            psi_phi[keep_index, 2] <- phi_current
            xi_mat[keep_index,] <- xi_current
            keep_index <- keep_index + 1
        }

        if (gindex %in% plot_iters) {
            if (plot_update & gindex > burnin + keep) {
                graphics::par(mfrow = c(2, 1))
                graphics::plot(psi_phi[, 1], type = "l", ylab = expression(psi))
                graphics::plot(psi_phi[, 2], type = "l", ylab = expression(phi))
                graphics::par(mfrow = c(1, 1))
            }
        }
    }
    if (print_update) {
        cat("\nComplete!\n")
    }

    return(list(Y22_array = Y22_array, psi_phi = psi_phi, xi_mat = xi_mat))
}


#' Wrapper for the bfa package.
#'
#' This is just a wrapper for using the \code{\link[bfa]{bfa_gauss}}
#' from the bfa package. Default settings are used.
#'
#' @inheritParams em_miss
#' @param nsamp A positive integer. The number of samples to draw.
#' @param burnin A positive integer. The number of early samples to
#'     discard.
#' @param keep A positive integer. We will same the updates of
#'     \code{Y22} every \code{keep} iteration of the Gibbs sampler.
#' @param print_status How often should bfa print updates?
#'
#' @author David Gerard
#'
#' @seealso \code{\link[bfa]{bfa_gauss}} for the workhorse function.
bfa_wrapper <- function(Y21, Y31, Y32, k, nsamp = 10000, burnin = round(nsamp / 4), keep = 20,
                        print_status = 500) {

    n <- nrow(Y21) + nrow(Y31)
    p <- ncol(Y31) + ncol(Y32)
    ncovs <- nrow(Y21)
    ncontrols <- ncol(Y21)

    Y22NA <- matrix(NA, nrow = ncovs, ncol = p - ncontrols)

    Y <- as.data.frame(rbind(cbind(Y21, Y22NA), cbind(Y31, Y32)))
    v_vec <- paste("V", 1:p, sep = "")
    colnames(Y) <- v_vec
    form1 <- stats::as.formula(paste("~", paste(v_vec, collapse = " + ")))
    trash <- utils::capture.output(bfout <- bfa::bfa_gauss(form1, data = Y, num.factor = k,
                                                           keep.scores = TRUE, thin = keep,
                                                           nburn = burnin, nsim = nsamp,
                                                           factor.scales = TRUE,
                                                           loading.prior = "normal",
                                                           print.status = print_status))

    nmcmc_samp <- dim(bfout$post.loadings)[3]
    Y22_array <- array(NA, dim = c(ncovs, p - ncontrols, nmcmc_samp))
    for(index in 1:nmcmc_samp) {
        if (k == 1) {
            Y22_array[,, index] <- matrix(bfout$post.scores[, 1:ncovs, index], ncol = 1) %*%
                matrix(bfout$post.loadings[(ncontrols + 1):p,, index], nrow = 1)
        } else {
            Y22_array[,, index] <- t(bfout$post.loadings[(ncontrols + 1):p,, index] %*%
                                     bfout$post.scores[, 1:ncovs, index])
        }
    }

    return(list(Y22_array = Y22_array, sigma2 = bfout$post.sigma2))
}



#' Bayesian factor analysis used in Gerard and Stephens (2016).
#'
#' Similar to that of Ghosh and Dunson (2009) but with two key
#' differences: (1) the prior is order invariant (though this makes
#' the factors and factor loadings unidentified), and (2) we place
#' hierarchical priors on the uniquenesses (variances).
#'
#'
#' @inheritParams em_miss
#' @param nsamp A positive integer. The number of samples to draw.
#' @param burnin A positive integer. The number of early samples to
#'     discard.
#' @param thin A positive integer. We will same the updates of
#'     \code{Y22} every \code{thin} iteration of the Gibbs sampler.
#' @param display_progress A logical. Should we print a text progress bar
#'     to keep track of the Gibbs sampler (\code{TRUE}) or not
#'     (\code{FALSE})?
#' @param hetero_factors A logical. Should we assign colum-specific
#'     variances for the factors (\code{TRUE}) or not (\code{FALSE})?
#' @param rho_0 The prior sample size for column-specific the
#'     precisions.
#' @param alpha_0 The prior sample size for the mean of the
#'     column-specific precisions.
#' @param beta_0 The prior mean of the mean of the column-specific
#'     precisions.
#' @param eta_0 The prior sample size of the expanded parameters.
#' @param tau_0 The prior mean of of the expanded parameters.
#' @param delta_0 The prior sample size of the column-specific
#'     precisions of the factors.
#' @param lambda_0 The prior sample size of the mean of the
#'     column-specific precisions of the factors.
#' @param nu_0 The prior mean of the mean of the column-specific
#'     precisions of the factors.
#'
#' @export
#'
#' @author David Gerard
#'
#' @seealso \code{\link{gdfa}} for the slower R implementation.
#'
#' @references Ghosh, Joyee, and David
#'     B. Dunson. "Default prior distributions and efficient posterior
#'     computation in Bayesian factor analysis."
#'     Journal of Computational and Graphical Statistics 18.2 (2009):
#'     306-320.
bfa_gs <- function(Y21, Y31, Y32, k, nsamp = 10000,
                   burnin = round(nsamp / 4), thin = 20,
                   display_progress = TRUE, hetero_factors = TRUE,
                   rho_0 = 0.1, alpha_0 = 0.1, delta_0 = 0.1,
                   lambda_0 = 0.1, nu_0 = 1, beta_0 = 1, eta_0 = 1,
                   tau_0 = 1) {

    assertthat::are_equal(ncol(Y21), ncol(Y31))
    assertthat::are_equal(nrow(Y31), nrow(Y32))
    ncovs <- nrow(Y21)
    ncontrols <- ncol(Y21)
    n <- nrow(Y21) + nrow(Y31)
    p <- ncol(Y31) + ncol(Y32)

    ## Get initial values and set hyper parameters----------------------
    pcout <- pca_naive(cbind(Y31, Y32), r = k)
    sig_diag <- pcout$sig_diag
    alpha <- t(pcout$alpha)
    Z3 <- pcout$Z
    shrunk_var <- limma::squeezeVar(sig_diag, df = n - ncovs - k)$var.post
    shrunk_var_c <- shrunk_var[1:ncontrols]
    alpha_c    <- alpha[, 1:ncontrols, drop = FALSE]
    Z2 <- tcrossprod(sweep(Y21, 2, 1 / shrunk_var_c, `*`), alpha_c) %*%
        solve(tcrossprod(sweep(alpha_c, 2, 1 / sqrt(shrunk_var_c), `*`)))
    Y22init <- Z2 %*% alpha[, (ncontrols + 1):p, drop = FALSE]
    Zinit <- rbind(Z2, Z3)
    fnorm_z <- sum(Zinit ^ 2)
    fnorm_alpha <- sum(alpha ^ 2)
    Linit <- Zinit
    Finit <- alpha

    if (hetero_factors) {
        theta_init <- 1 / limma::squeezeVar(colMeans(Finit ^ 2), df = k)$var.post
        kappa_init <- mean(theta_init)
    } else {
        theta_init <- rep(1, length = p)
        kappa_init <- 1
    }

    xi_init <- 1 / shrunk_var
    phi_init <- mean(xi_init)

    zeta_init <- 1 / limma::squeezeVar(colMeans(Linit ^ 2), df = k)$var.post

    if (display_progress) {
        cat("Progress:\n")
    }
    bfout <- bfa_gd_gibbs(Linit = Linit, Finit = Finit,
                          xi_init = xi_init, phi_init = phi_init,
                          zeta_init = zeta_init, theta_init = theta_init,
                          kappa_init = kappa_init, Y22init = Y22init,
                          Y21 = Y21, Y31 = Y31, Y32 = Y32,
                          nsamp = nsamp, burnin = burnin,
                          thin = thin, rho_0 = rho_0,
                          alpha_0 = alpha_0, delta_0 = delta_0,
                          lambda_0 = lambda_0, nu_0 = nu_0,
                          beta_0 = beta_0, eta_0 = eta_0,
                          tau_0 = tau_0, hetero_factors = hetero_factors,
                          display_progress = display_progress)
    if (display_progress) {
        cat("Complete!\n")
    }

    return(bfout)
}




#' Old version of Bayesian factor analysis.
#'
#' See \code{\link{bfa_gs}} for a description. This is the same thing
#' but a slower R implementation.
#'
#' @inheritParams em_miss
#' @param nsamp A positive integer. The number of samples to draw.
#' @param burnin A positive integer. The number of early samples to
#'     discard.
#' @param keep A positive integer. We will same the updates of
#'     \code{Y22} every \code{keep} iteration of the Gibbs sampler.
#' @param print_update A logical. Should we print a text progress bar
#'     to keep track of the Gibbs sampler (\code{TRUE}) or not
#'     (\code{FALSE})?
#' @param plot_update A logical. Should we make some plots to keep
#'     track of the Gibbs sampler (\code{TRUE}) or not (\code{FALSE})?
#' @param hetero_factors A logical. Should we assign colum-specific
#'     variances for the factors (\code{TRUE}) or not (\code{FALSE})?
#' @param rho_0 The prior sample size for column-specific the
#'     precisions.
#' @param alpha_0 The prior sample size for the mean of the
#'     column-specific precisions.
#' @param beta_0 The prior mean of the mean of the column-specific
#'     precisions.
#' @param eta_0 The prior sample size of the expanded parameters.
#' @param tau_0 The prior mean of of the expanded parameters.
#' @param delta_0 The prior sample size of the column-specific
#'     precisions of the factors.
#' @param lambda_0 The prior sample size of the mean of the
#'     column-specific precisions of the factors.
#' @param nu_0 The prior mean of the mean of the column-specific
#'     precisions of the factors.
#'
#'
#' @author David Gerard
#'
#' @seealso \code{\link{bfa_gs}}
#'
#' @references Ghosh, Joyee, and David
#'     B. Dunson. "Default prior distributions and efficient posterior
#'     computation in Bayesian factor analysis."
#'     Journal of Computational and Graphical Statistics 18.2 (2009):
#'     306-320.
gdfa <- function(Y21, Y31, Y32, k, nsamp = 10000,
                 burnin = round(nsamp / 4), keep = 20,
                 print_update = TRUE, plot_update = FALSE,
                 hetero_factors = FALSE, rho_0 = 0.1, alpha_0 = 0.1,
                 delta_0 = 0.1, lambda_0 = 0.1, nu_0 = 1, beta_0 = 1,
                 eta_0 = 1, tau_0 = 1) {

    assertthat::are_equal(ncol(Y21), ncol(Y31))
    assertthat::are_equal(nrow(Y31), nrow(Y32))
    ncovs <- nrow(Y21)
    ncontrols <- ncol(Y21)
    n <- nrow(Y21) + nrow(Y31)
    p <- ncol(Y31) + ncol(Y32)

    ## Get initial values and set hyper parameters----------------------
    pcout <- pca_naive(cbind(Y31, Y32), r = k)
    sig_diag <- pcout$sig_diag
    alpha <- t(pcout$alpha)
    Z3 <- pcout$Z
    shrunk_var <- limma::squeezeVar(sig_diag, df = n - ncovs - k)$var.post
    shrunk_var_c <- shrunk_var[1:ncontrols]
    alpha_c    <- alpha[, 1:ncontrols, drop = FALSE]
    Z2 <- tcrossprod(sweep(Y21, 2, 1 / shrunk_var_c, `*`), alpha_c) %*%
        solve(tcrossprod(sweep(alpha_c, 2, 1 / sqrt(shrunk_var_c), `*`)))
    Y22init <- Z2 %*% alpha[, (ncontrols + 1):p, drop = FALSE]
    Zinit <- rbind(Z2, Z3)
    fnorm_z <- sum(Zinit ^ 2)
    fnorm_alpha <- sum(alpha ^ 2)
    Linit <- Zinit * (prod(dim(Zinit)) / fnorm_z)
    Finit <- alpha * (prod(dim(alpha)) / fnorm_alpha)

    if (hetero_factors) {
        theta_init <- 1 / limma::squeezeVar(colMeans(Finit ^ 2), df = k)$var.post
        kappa_init <- mean(theta_init)
    } else {
        theta_init <- rep(1, length = p)
        kappa_init <- 1
    }

    xi_init <- 1 / shrunk_var
    phi_init <- mean(xi_init)

    zeta_init <- 1 / limma::squeezeVar(colMeans(Linit ^ 2), df = k)$var.post

    L_current     <- Linit
    F_current     <- Finit
    xi_current    <- xi_init
    phi_current   <- phi_init
    zeta_current  <- zeta_init
    theta_current <- theta_init
    kappa_current <- kappa_init
    Y22_current   <- Y22init
    Y_current     <- rbind(cbind(Y21, Y22_current), cbind(Y31, Y32))

    ## Gibbs Sampler ---------------------------------------------------------
    Y22_array <- array(NA, dim = c(ncovs, p - ncontrols, floor((nsamp - burnin) / keep)))
    xi_mat    <- matrix(NA, nrow = floor((nsamp - burnin) / keep), ncol = p)
    phi_vec   <- rep(NA, floor((nsamp - burnin) / keep))
    keep_index <- 1
    plot_iters <- round(seq(0.01, 1, by = 0.01) * nsamp)
    if (print_update) {
        cat("Progress:\n")
        pb <- utils::txtProgressBar(style = 3)
    }
    for (gindex in 1:nsamp) {
        if (print_update) {
            utils::setTxtProgressBar(pb = pb, value = gindex / nsamp)
        }

        ## Update L ----------------------------------------------------
        Fsig <- sweep(F_current, 2, xi_current, `*`)
        eigen_fsf <- eigen(tcrossprod(Fsig, F_current) +
                           diag(x = zeta_current, nrow = length(zeta_current),
                                ncol = length(zeta_current)), symmetric = TRUE)

        L_meanmat <- tcrossprod(sweep(Y_current, 2, xi_current, `*`), F_current) %*%
            tcrossprod(sweep(eigen_fsf$vectors, 2, 1 / eigen_fsf$values, `*`),
                       eigen_fsf$vectors)

        col_cov_half <- tcrossprod(sweep(eigen_fsf$vectors, 2, 1 / sqrt(eigen_fsf$values), `*`),
                                   eigen_fsf$vectors)
        L_current <- L_meanmat +
            matrix(stats::rnorm(prod(dim(L_current))), nrow = nrow(L_current)) %*% col_cov_half

        ## checks
        ## L_mean_simp <- Y_current %*% diag(xi_current) %*% t(F_current) %*%
        ##     solve(F_current %*% diag(xi_current) %*% t(F_current) +
        ##           diag(zeta_current, nrow = length(zeta_current), ncol = length(zeta_current)))
        ## L_mean_simp
        ## L_meanmat
        ## solve(F_current %*% diag(xi_current) %*% t(F_current) +
        ##       diag(zeta_current, nrow = length(zeta_current), ncol = length(zeta_current)))
        ## col_cov_half %*% col_cov_half

        ## Update F ---------------------------------------------------
        ltl <- crossprod(L_current)
        LY <- crossprod(L_current, Y_current)
        if (k == 1) {
            F_current <- matrix(mapply(LY, xi_current, theta_current, FUN = update_f,
                                       MoreArgs = list(ltl = ltl)), nrow = 1)
            ## F_current <- matrix(apply(rbind(LY, xi_current, theta_current),
            ##                           MARGIN = 2, FUN = update_f2, ltl = ltl),
            ##                     nrow = 1)
        } else {
            F_current <- mapply(split(LY, col(LY)), xi_current, theta_current, FUN = update_f,
                                MoreArgs = list(ltl = ltl))
            ## F_current <- apply(rbind(LY, xi_current, theta_current),
            ##                    MARGIN = 2, FUN = update_f2, ltl = ltl)
        }

        ## Update xi ---------------------------------------------------------
        mean_current <- L_current %*% F_current

        r_vec <- colSums((Y_current - mean_current) ^ 2)
        ## should equal this
        ## diag(t(Y_current - mean_current) %*% (Y_current - mean_current))

        xi_shape <- (n + rho_0) / 2
        xi_rate  <- (r_vec + rho_0 * phi_current) / 2

        xi_current <- sapply(xi_rate, stats::rgamma, n = 1, shape = xi_shape)

        ## update phi --------------------------------------------------------
        phi_shape <- (p * rho_0 + alpha_0) / 2
        phi_rate  <- (alpha_0 * beta_0 + rho_0 * sum(xi_current)) / 2
        phi_current <- stats::rgamma(n = 1, shape = phi_shape, rate = phi_rate)

        ## Update zeta -------------------------------------------------------
        s_vec <- colSums(L_current ^ 2)
        zeta_shape <- (n + eta_0) / 2
        zeta_rate  <- (s_vec + eta_0 * tau_0) / 2
        zeta_current <- sapply(zeta_rate, stats::rgamma, n = 1, shape = zeta_shape)

        if (hetero_factors) {
            ## Update theta --------------------------------------------------
            u_vec <- colSums(F_current ^ 2)
            theta_shape <- (k + delta_0) / 2
            theta_rate  <- (u_vec + delta_0 * kappa_current) / 2
            theta_current <- sapply(theta_rate, stats::rgamma, n = 1, shape = theta_shape)

            ## Update kappa --------------------------------------------------
            kappa_shape <- (p * delta_0 + lambda_0) / 2
            kappa_rate  <- (lambda_0 * nu_0 + delta_0 * sum(theta_current)) / 2
            kappa_current <- stats::rgamma(n = 1, shape = kappa_shape, rate = kappa_rate)
        }
        ## Update Y22 ------------------------------------
        Y22_mean <- mean_current[1:ncovs, (ncontrols + 1):p, drop = FALSE]
        Y22_error <- sweep(matrix(stats::rnorm(n = ncovs * (p - ncontrols)), nrow = ncovs), 2,
                           1 / sqrt(xi_current[(ncontrols + 1):p]), `*`)
        Y22_current <- Y22_mean + Y22_error
        Y_current[1:ncovs, (ncontrols + 1):p] <- Y22_current

        if ((gindex - burnin) %% keep == 0 & gindex > burnin) {
            Y22_array[, , keep_index] <- Y22_current
            xi_mat[keep_index, ] <- xi_current
            phi_vec[keep_index] <- phi_current
            if (plot_update & gindex %in% plot_iters) {
                graphics::plot(phi_vec, type = "l", xlab = "Index", ylab = expression(phi))
            }
            keep_index <- keep_index + 1
        }
    }
    if (print_update) {
        cat("\nComplete!\n")
    }

    return(list(Y22_array = Y22_array, xi_mat = xi_mat))
}


#' Utility function for updating F in \code{gdfa}.
#'
#' @param LY A vector. A column of \code{t(L_current) \%*\% Y_current}
#' @param xi A numeric scalar. An element of xi_current.
#' @param theta A numeric scalar. An element of theta_current.
#' @param ltl t(L_current) \%*\% L_current
#'
#' @seealso \code{\link{gdfa}}
update_f <- function(LY, xi, theta, ltl) {
    k <- nrow(ltl)

    ## LY <- matrix(args[1:k], ncol = 1)
    ## xi <- args[k + 1]
    ## theta <- args[k + 2]

    eigen_ltl <- eigen(ltl * xi + diag(theta, k),
                       symmetric = TRUE)
    fmean <- sweep(eigen_ltl$vectors, 2, 1 / eigen_ltl$values, `*`) %*%
        crossprod(eigen_ltl$vectors, LY) * xi

    fcovhalf <- tcrossprod(sweep(eigen_ltl$vectors, 2, 1 / sqrt(eigen_ltl$values), `*`),
                           eigen_ltl$vectors)

    ## Checks
    ## fmean_simp <- xi * solve(ltl * xi + theta * diag(k)) %*% LY
    ## fmean
    ## fmean_simp
    ## solve(ltl * xi + theta * diag(k))
    ## fcovhalf %*% fcovhalf

    return(fmean + fcovhalf %*% matrix(stats::rnorm(k), ncol = 1))
}


update_f2 <- function(args, ltl) {
    k <- nrow(ltl)

    LY <- matrix(args[1:k], ncol = 1)
    xi <- args[k + 1]
    theta <- args[k + 2]

    eigen_ltl <- eigen(ltl * xi + diag(theta, k),
                       symmetric = TRUE)
    fmean <- sweep(eigen_ltl$vectors, 2, 1 / eigen_ltl$values, `*`) %*%
        crossprod(eigen_ltl$vectors, LY) * xi

    fcovhalf <- tcrossprod(sweep(eigen_ltl$vectors, 2, 1 / sqrt(eigen_ltl$values), `*`),
                           eigen_ltl$vectors)

    return(fmean + fcovhalf %*% matrix(stats::rnorm(k), ncol = 1))
}
