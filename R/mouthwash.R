## Functions for performing MOUTHWASH. This is a re-coding from
## succotashr because I'm slightly better at programming now.

#' MOUTHWASH: Maximize Over Unobservables To Help With Adaptive SHrinkage.
#'
#' This function implements the full MOUTHWASH method. First, it
#' rotates the response and explanatory variables into a part that we
#' use to estimate the confounding variables and the variances, and a
#' part that we use to estimate the coefficients of the observed
#' covariates. This function will implement a factor analysis for the
#' first part then run \code{\link{mouthwash_second_step}} for the
#' second part.
#'
#' The assumed model is \deqn{Y = X\beta + Z\alpha + E.} \eqn{Y} is a
#' \eqn{n} by \code{p} matrix of response variables. For example, each
#' row might be an array of log-transformed gene-expression data.
#' \eqn{X} is a \eqn{n} by \eqn{q} matrix of observed covariates. It
#' is assumed that all but one column of which contains nuisance
#' parameters. For example, the first column might be a vector of ones
#' to include an intercept. \eqn{\beta} is a \eqn{q} by \eqn{p} matrix
#' of corresponding coefficients.  \eqn{Z} is a \eqn{n} by \eqn{k}
#' matrix of confounder variables. \eqn{\alpha} is the corresponding
#' \eqn{k} by \eqn{p} matrix of coefficients for the unobserved
#' confounders. \eqn{E} is a \eqn{n} by \eqn{p} matrix of error
#' terms. \eqn{E} is assumed to be matrix normal with identity row
#' covariance and diagonal column covariance \eqn{\Sigma}. That is,
#' the columns are heteroscedastic while the rows are homoscedastic
#' independent.
#'
#' This function will first rotate \eqn{Y} and \eqn{X} using the QR
#' decomposition. This separates the model into three parts. The first
#' part contains nuisance parameters, the second part contains
#' the coefficients of interest, and the third part contains the
#' confounders. \code{mouthwash} applies a user-provided factor
#' analysis to the third part to estimate the confounding factors,
#' then runs an EM (or coordinate-ascent) algorithm on the second part
#' to estimate the coefficients of interest.
#'
#' There are a couple forms of factor analysis available in this
#' package. The default is PCA with the column-wise residual
#' mean-squares as the estimates of the column-wise variances.
#'
#' For instructions and examples on how to specify your own factor analysis, run the following code in R:
#' \code{utils::vignette("customFA", package = "vicar")}. If it doesn't work, then you probably haven't built
#' the vignettes. To do so, see \url{https://github.com/dcgerard/vicar#vignettes}.
#'
#'
#' @seealso Factor analyses available in the \code{vicar} package:
#'     \code{\link{pca_naive}}, \code{\link{fa_ml}}.
#'
#' @inheritParams vruv4
#' @param grid_seq The grid for the mixing distribution. If
#'     \code{mixing_dist = "uniform"} or \code{"+uniform"}, then these
#'     should be the non-zero limits of the uniform distributions. If
#'     \code{mixing_dist = "sym_uniform"}, then these should be the
#'     right limits of the uniform distributions. If \code{mixing_dist
#'     = "normal"}, then these should be the variances of the mixing
#'     normal distributions.
#' @param pi_init A numeric vector. These are the initial values of
#'     the mixing proportions.
#' @param lambda_seq A numeric vector with elements all greater than
#'     or equal to 1. These are the tuning parameters for the mixing
#'     proportions. This can only be specified if \code{grid_seq} is
#'     also specified.
#' @param mixing_dist A character. Should we use a mixture of uniforms
#'     (\code{"uniform"}), a mixture of uniforms with minimum at 0
#'     (\code{"+uniform"}), a mixture of uniforms symmetric at 0
#'     (\code{"sym_uniform"}), or a mixture of normals
#'     (\code{"normal"})?
#' @param lambda_type A character. Should we apply a penalty on zero
#'     (\code{"zero_conc"}) or no penalty (\code{"uniform"}). Not used
#'     if \code{lambda_seq} is not \code{NULL}.
#' @param lambda0 A numeric greater than or equal to 1. The penalty on
#'     zero if \code{lambda_type = "zero_conc"}.
#' @param pi_init_type How should we initialize the mixing
#'     proportions? By concentrating on zero (\code{"zero_conc"}), by
#'     equal weights on all mixing distributions (\code{"uniform"}),
#'     or by sampling uniformly on the simplex (\code{"random"})?
#' @param scale_var A logical. Should we estimate a variance inflation
#'     parameter (\code{TRUE}) or not (\code{FALSE})?
#' @param plot_update A logical. Should I plot the the path of the
#'     log-likelihood (\code{TRUE}) or not (\code{FALSE})? Only
#'     applicable when \code{mixing_dist} is not \code{"normal"}.
#' @param sprop If \eqn{b} is an effect and \eqn{s} is an estimated
#'     standard error, then we model \eqn{b/s^{sprop}} as
#'     exchangeable. The default is 0. When \code{sprop = 1}, for
#'     identifiability reasons it must be the case that
#'     \code{scale_var = FALSE}.
#' @param var_inflate_pen The penalty to apply on the variance inflation parameter.
#'     Defaults to 0, but should be something non-zero when \code{alpha = 1}
#'     and \code{scale_var = TRUE}.
#' @param subsample A logical. Should we only use a subsample of the genes to estimate
#'     the hidden covariates (\code{TRUE}) or use all of the genes (\code{FALSE})? If
#'     \code{TRUE}, then \code{\link[ashr]{ash}} will be re-run on the residuals (after
#'     subtracting out the contribution from the unobserved confounders) to obtain the
#'     estimated prior.
#' @param num_sub The number of genes to subsample if \code{subsample = TRUE}. Not used if
#'     \code{subsample = FALSE}.
#' @param same_grid A logical. If \code{subsample = FALSE}, should we use the same grid as
#'     when we estimated the unobserved confounders (\code{TRUE}) or the default grid from
#'     \code{\link[ashr]{ash.workhorse}} (\code{FALSE})?
#' @param use_t_adjust A logical. Should we adjust the variance estimates so that the p-values
#'     from the z-statistics match the corresponding p-values from the original
#'     t-statistics (\code{TRUE}) or not (\code{FALSE})?
#' @param detailed_output A logical. Should we return a lot of output (\code{TRUE}) or the standard
#'     output (\code{FALSE}). Most users should only need this set to (\code{FALSE}).
#' @param verbose If \code{verbose = TRUE}, print progress of the algorithm
#'   to the console.
#' @param cov_of_interest A positive integer. The column number of the
#'     covariate in X whose coefficients you are interested in.
#'     The rest are considered nuisance parameters and are regressed
#'     out by OLS.
#' @param X A matrix of numerics. The observed covariates.
#'
#' @return A list with some or all of the following elements.
#'
#'     \code{fitted_g}: The estimated unimodal prior. It is of class
#'     \code{\link[ashr]{unimix}} if \code{mixing_dist} is one of
#'     \code{"uniform"}, \code{"+uniform"}, or
#'     \code{"sym_uniform"}. It is of class
#'     \code{\link[ashr]{normalmix}} if \code{mixing_dist} is
#'     \code{"normal"}.
#'
#'     \code{loglik} The final log-likelihood.
#'
#'     \code{logLR} The likelihood ratio compared to the all-null setting (point-mass on zero).
#'
#'     \code{data} Post-confounder adjusted ashr data.
#'
#'     \code{pi0} The estimate of the proportion of null genes.
#'
#'     \code{z2} The estimated confounders (after rotation)
#'     corresponding the covariates of interest (\code{z2}). Mostly output for
#'     debugging reasons.
#'
#'     \code{xi} The estimated variance inflation parameter.
#'
#'     \code{Zhat} The estimate of the confounders.
#'
#'     \code{alphahat} The estimate of the coefficients of the confounders.
#'
#'     \code{sig_diag} The estimate of the column-specific variances.
#'
#'     \code{result} A data frame with the results from MOUTHWASH. The columns of which are
#'     \itemize{
#'       \item{NegativeProb}{The probability that the effect is negative.}
#'       \item{PositiveProb}{The probability that the effect is positive.}
#'       \item{lfsr}{The local false sign rates of each effect.}
#'       \item{svalue}{The s-values, a measure of significance.}
#'       \item{lfdr}{The local false discovery rates.}
#'       \item{qvalue}{The q-values, a measure of significance.}
#'       \item{PosteriorMean}{The posterior means of the effects.}
#'       \item{PosteriorSD}{The posterior standard deviations of the effects.}
#'       \item{extra}{If \code{detailed_output = TRUE}, this list is returned with some extra output. Mostly used for debugging.}
#'     }
#'
#' @export
#'
#' @seealso \code{\link{backwash}} for a similar method that puts a prior on the
#'     unobserved confounders rather than maximizes over them.
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., and Stephens, M. 2020. "Empirical Bayes shrinkage and false discovery rate estimation, allowing for unwanted variation", \emph{Biostatistics}, 21(1), 15-32 \doi{10.1093/biostatistics/kxy029}}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' library(vicar)
#'
#' ## Generate data ----------------------------------------------------------
#' set.seed(116)
#' n <- 13
#' p <- 101
#' k <- 2
#' q <- 3
#' is_null       <- rep(FALSE, length = p)
#' is_null[1:57] <- TRUE
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
#' ## Fit MOUTHWASH ----------------------------------------------------------
#' mout <- mouthwash(Y = Y, X = X, k = k, cov_of_interest = 2,
#'                   include_intercept = FALSE)
#' mout$pi0 ## mouthwash estimate
#' mean(is_null) ## truth
#'
#' ## plot ordering
#' order_lfdr <- order(mout$result$lfdr)
#' graphics::plot(mout$result$lfdr[order_lfdr], col = is_null[order_lfdr] + 3,
#'                ylab = "lfdr")
#' graphics::legend("topleft", col = c(3, 4), legend = c("non-null", "null"),
#'                  pch = 1)
#'
#' ## Compare to ASH on OLS coefficients -------------------------------------
#' lmout <- limma::lmFit(t(Y), X)
#' betahat_ols <- lmout$coefficients[, 2]
#' sebetahat_ols <- lmout$stdev.unscaled[, 2] * lmout$sigma
#' aout <- ashr::ash.workhorse(betahat = betahat_ols, sebetahat = sebetahat_ols,
#'                             optmethod = "mixEM")
#' ashr::get_pi0(aout) ## ash estimate
#' mean(is_null) ## truth
#'
#' ash_lfdr <- ashr::get_lfdr(aout)
#' aorder_lfdr <- order(ash_lfdr)
#' graphics::plot(ash_lfdr[aorder_lfdr], col = is_null[aorder_lfdr] + 3,
#'                ylab = "lfdr")
#' graphics::legend("topleft", col = c(3, 4), legend = c("non-null", "null"),
#'                  pch = 1)
#'
#' ## ROC Curves -------------------------------------------------------------
#' afpr <- cumsum(is_null[aorder_lfdr]) / sum(is_null)
#' atpr <- cumsum(!is_null[aorder_lfdr]) / sum(!is_null)
#' mfpr <- cumsum(is_null[order_lfdr]) / sum(is_null)
#' mtpr <- cumsum(!is_null[order_lfdr]) / sum(!is_null)
#' graphics::plot(afpr, atpr, type = "l", xlab = "False Positive Rate",
#'                ylab = "True Positive Rate", main = "ROC Curve", col = 3,
#'                lty = 2)
#' graphics::lines(mfpr, mtpr, col = 4, lty = 1)
#' graphics::abline(0, 1, lty = 2, col = 1)
#' graphics::legend("bottomright", legend = c("MOUTHWASH", "ASH"), col = c(4, 3),
#'                  lty = c(1, 2))
#'
mouthwash <- function(Y, X, k = NULL, cov_of_interest = ncol(X),
                      include_intercept = TRUE, limmashrink = TRUE,
                      fa_func = pca_naive, fa_args = list(),
                      likelihood = c("normal", "t"),
                      mixing_dist = c("normal", "uniform", "+uniform", "sym_uniform"),
                      lambda_type = c("zero_conc", "uniform"),
                      pi_init_type = c("zero_conc", "uniform", "random"),
                      degrees_freedom = NULL, pi_init = NULL,
                      grid_seq = NULL, lambda_seq = NULL,
                      lambda0 = 10, scale_var = TRUE,
                      plot_update = FALSE,
                      sprop = 0, var_inflate_pen = 0,
                      subsample = FALSE,
                      num_sub = min(1000, ncol(Y)),
                      same_grid = FALSE,
                      use_t_adjust = FALSE,
                      detailed_output = FALSE,
                      verbose = TRUE) {

    ## Make sure input is correct -------------------------------------------
    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(is.numeric(cov_of_interest))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.logical(limmashrink))
    assertthat::assert_that(is.function(fa_func))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(lambda0 >= 1)
    assertthat::assert_that(sprop >= 0)
    assertthat::assert_that(var_inflate_pen >= 0)
    assertthat::assert_that(num_sub >= 1, num_sub <= ncol(Y))
    check_same <- apply(Y, 2, stats::sd) > 0
    if (!all(check_same)) {
      stop(paste0("Columns [", paste(which(!check_same, arr.ind = TRUE), collapse = ", "), "] of Y have sample SD of 0. This is prohibited."))
    }
    if (length(cov_of_interest) > 1) {
      stop("We do not currently support more than one covariate of interest.")
    }

    likelihood   <- match.arg(likelihood)
    mixing_dist  <- match.arg(mixing_dist)
    pi_init_type <- match.arg(pi_init_type)
    lambda_type  <- match.arg(lambda_type)

    if (use_t_adjust & likelihood == "normal") {
      stop("to use the use_t_adjust option, please set likelihood = 't'")
    }

    if (likelihood == "t" & mixing_dist == "normal" & !use_t_adjust) {
      stop("normal mixtures not implemented for t-likelihood unless use_t_adjust = TRUE.")
    }
    if (sprop == 1 & scale_var & var_inflate_pen < 10 ^ -6) {
      stop("if sprop is 1 and scale_var = TRUE, then var_inflate_pen should be > 0")
    }

    if (!subsample & ncol(Y) >= 20000) {
      message("This will take awhile. To speed things up, try setting `subsample = TRUE`")
    }

    if (verbose)
      cat(sprintf(paste("Running mouthwash on %d x %d matrix X and",
                        "%d x %d matrix Y.\n"),
                  nrow(X),ncol(X),nrow(Y),ncol(Y)))

    ## Rotate -------------------------------------------------------------
    if (verbose)
      cat(" - Computing independent basis using QR decomposition.\n")
    timing <- system.time(
      rotate_out <- rotate_model(Y = Y, X = X, k = k,
                                 cov_of_interest = cov_of_interest,
                                 include_intercept = include_intercept,
                                 limmashrink = limmashrink, fa_func = fa_func,
                                 fa_args = fa_args, do_factor = TRUE))
    if (verbose)
      cat(" - Computation took",timing["elapsed"],"seconds.\n")
    if (rotate_out$k == 0) {
      stop("k estimated to be 0. You might not need mouthwash")
    }

    ## Deal with degrees of freedom -----------------------------------------
    if (verbose)
        cat(" - Running additional preprocessing steps.\n")
    timing <- system.time({
    if (likelihood == "normal") {
        if (!is.null(degrees_freedom)) {
            message("likelihood = \"normal\" but degrees_freedom not NULL. Setting degrees_freedom to Inf")
        }
        degrees_freedom <- Inf
    }
    if (!is.null(rotate_out$prior_df) & is.null(degrees_freedom) & likelihood == "t") {
        degrees_freedom <- rotate_out$prior_df + nrow(X) - ncol(X) - k
        if (degrees_freedom == Inf) {
            message("limma estimated df = Inf . Changing likelihood to \"normal\".")
            likelihood <- "normal"
        }
    } else if (is.null(degrees_freedom) & likelihood == "t") {
        degrees_freedom <- nrow(X) - ncol(X) - k
    }
    assertthat::assert_that(length(degrees_freedom) == 1 | length(degrees_freedom) == ncol(Y))
    assertthat::assert_that(all(degrees_freedom > 0))

    ## rescale alpha and sig_diag by R22 to get data for second step --------
    alpha_tilde <- rotate_out$alpha / c(rotate_out$R22)
    S_diag      <- c(rotate_out$sig_diag / c(rotate_out$R22 ^ 2))
    betahat_ols <- matrix(rotate_out$betahat_ols, ncol = 1)

    ## use adjust_by_t to use normal ----------------------------------------
    if (use_t_adjust) {
      S_diag <- adjust_by_t(betahat = betahat_ols, sebetahat = sqrt(S_diag),
                            df = degrees_freedom) ^ 2
      likelihood <- "normal" ## should switch from t to normal likelihood
      degrees_freedom <- Inf ## because dealt with degrees_freedom
                             ## earlier, so need to do it here after
                             ## changing the likelihood.
    }

    ## Exchangeable versions of the models ----------------------------------
    if (sprop > 0) {
        sgamma           <- S_diag ^ (-1 * sprop / 2)
        alpha_tilde_star <- alpha_tilde * sgamma
        betahat_ols_star <- betahat_ols * sgamma
        S_diag_star      <- S_diag ^ (1 - sprop)
    } else {
        alpha_tilde_star <- alpha_tilde
        betahat_ols_star <- betahat_ols
        S_diag_star      <- S_diag
    }

    ## Set grid and penalties ----------------------------------------------
    if (!is.null(lambda_seq) & is.null(grid_seq)) {
        stop("lambda_seq specified but grid_seq is NULL")
    }

    if (is.null(grid_seq)) {
        grid_vals <- get_grid_var(betahat_ols = betahat_ols_star, S_diag = S_diag_star)
        if (mixing_dist == "normal") {
            grid_seq <- sign(grid_vals$tau2_seq) * sqrt(abs(grid_vals$tau2_seq))
        } else {
            grid_seq <- grid_vals$tau2_seq
        }
    }

    if (mixing_dist == "normal") {
        tau2_seq <- grid_seq
        M <- length(tau2_seq)
        a_seq <- NULL
        b_seq <- NULL
        zero_spot <- which(abs(tau2_seq) < 10 ^ -14)
    } else if (mixing_dist == "uniform") {
        a_seq <- c(-1 * grid_seq[length(grid_seq):2], rep(0, length(grid_seq)))
        b_seq <- c(rep(0, length(grid_seq)), grid_seq[2:length(grid_seq)])
        M <- length(a_seq)
        tau2_seq <- NULL
        zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    } else if (mixing_dist == "+uniform") {
        a_seq <- rep(0, length(grid_seq))
        b_seq <- grid_seq
        M <- length(a_seq)
        tau2_seq <- NULL
        zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    } else if (mixing_dist == "sym_uniform") {
        a_seq <- -1 * grid_seq
        b_seq <- grid_seq
        M <- length(a_seq)
        tau2_seq <- NULL
        zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    }
    assertthat::are_equal(length(zero_spot), 1)

    if (is.null(lambda_seq)) {
        if (lambda_type == "uniform") {
            lambda_seq <- rep(1, M)
        } else if (lambda_type == "zero_conc") {
            lambda_seq <- rep(1, M)
            lambda_seq[zero_spot] <- lambda0
        }
    }})
    if (verbose)
      cat(" - Computation took",timing["elapsed"],"seconds.\n")

    ## Run MOUTHWASH --------------------------------------------------------
    if (!subsample) {
      if (verbose)
        cat(" - Running second step of mouthwash:\n")
      timing <- system.time(
        val <- mouthwash_second_step(betahat_ols = betahat_ols_star,
                                     S_diag = S_diag_star,
                                     alpha_tilde = alpha_tilde_star,
                                     tau2_seq = tau2_seq, a_seq = a_seq,
                                     b_seq = b_seq,
                                     degrees_freedom = degrees_freedom,
                                     lambda_seq = lambda_seq,
                                     mixing_dist = mixing_dist,
                                     likelihood = likelihood,
                                     pi_init_type = pi_init_type,
                                     scale_var = scale_var,
                                     plot_update = plot_update,
                                     sprop = sprop,
                                     var_inflate_pen = var_inflate_pen,
                                     verbose = verbose))
      if (verbose)
        cat(" - Second step took",timing["elapsed"],"seconds.\n")
    } else {
      cat("Running second step of mouthwash:\n")
      timing <- system.time({
        col_keep <- sort(sample(x = 1:ncol(Y), size = num_sub))
        betahat_ols_star <- betahat_ols_star[col_keep]
        S_diag_star <- S_diag_star[col_keep]
        alpha_tilde_star <- alpha_tilde_star[col_keep, , drop = FALSE]

        val2 <- mouthwash_second_step(betahat_ols = betahat_ols_star,
                                      S_diag = S_diag_star,
                                      alpha_tilde = alpha_tilde_star,
                                      tau2_seq = tau2_seq, a_seq = a_seq,
                                      b_seq = b_seq,
                                      degrees_freedom = degrees_freedom,
                                      lambda_seq = lambda_seq,
                                      mixing_dist = mixing_dist,
                                      likelihood = likelihood,
                                      pi_init_type = pi_init_type,
                                      scale_var = scale_var,
                                      plot_update = plot_update,
                                      sprop = sprop,
                                      var_inflate_pen = var_inflate_pen,
                                      verbose = verbose)
      })
      if (verbose)
        cat(" - Second step took",timing["elapsed"],"seconds.\n")
      if (verbose)
        cat(" - Running adaptive shrinkage method.\n")
      timing <- system.time({
          az <- alpha_tilde %*% val2$z2

          mixcompdist <- mixing_dist
          if (mixcompdist == "uniform") {
              mixcompdist <- "halfuniform"
          } else if (mixcompdist == "sym_uniform") {
              mixcompdist <- "uniform"
          }

          ashr_df <- degrees_freedom
          if (likelihood == "normal") {
              ashr_df <- NULL
          }

          if (same_grid) {
              ash_g <- val2$fitted_g
          } else {
              ash_g <- NULL
          }
          val <- ashr::ash.workhorse(betahat = c(betahat_ols - az),
                                     sebetahat = c(sqrt(val2$xi * S_diag)),
                                     df = ashr_df,
                                     prior = "nullbiased",
                                     nullweight = lambda_seq[zero_spot],
                                     mixcompdist = mixcompdist,
                                     g = ash_g,
                                     alpha = sprop)
          val$pi0 <- val2$pi0
          val$xi  <- val2$xi
          val$z2  <- val2$z2
      })
      if (verbose)
        cat(" - Computation took",timing["elapsed"],"seconds.\n")
   }

    ## Estimate rest of the hidden confounders ------------------------------
    if (verbose)
      cat(" - Estimating additional hidden confounders.\n")
    timing <- system.time({
    Y1  <- rotate_out$Y1
    Z2 <- val$z2
    Z3 <- rotate_out$Z3
    if (!is.null(Y1)) {
        R12 <- rotate_out$R12
        R11 <- rotate_out$R11
        Q   <- rotate_out$Q
        beta1_ols <- solve(R11) %*% (Y1 - R12 %*% t(betahat_ols))
        resid_top <- Y1 - R12 %*% t(val$result$PosteriorMean) - R11 %*% beta1_ols
        Z1  <- solve(t(alpha_tilde) %*% diag(1 / rotate_out$sig_diag) %*% alpha_tilde) %*%
            t(alpha_tilde) %*% diag(1 / rotate_out$sig_diag) %*% t(resid_top)
        Zhat <- Q %*% rbind(t(Z1), t(Z2), Z3)
    } else {
        Q   <- rotate_out$Q
        Zhat <- Q %*% rbind(t(Z2), Z3)
    }

    val$Zhat <- Zhat
    val$alphahat <- t(rotate_out$alpha)
    val$sig_diag <- rotate_out$sig_diag

    class(val) <- "mouthwash"

    if (detailed_output) {
      val$extra                  <- list()
      val$extra$az               <- alpha_tilde_star %*% val$z2
      val$extra$alpha_tilde      <- alpha_tilde
      val$extra$alpha_tilde_star <- alpha_tilde_star
      val$extra$rotate_out       <- rotate_out
    } else {
      val$z2 <- NULL
    }})
    if (verbose)
      cat(" - Computation took",timing["elapsed"],"seconds.\n")

    return(val)
}

#' The second step of MOUTHWASH.
#'
#' @inheritParams mouthwash
#' @param betahat_ols A vector of numerics. The OLS estimates of the
#'     coefficients of interest.
#' @param S_diag A vector of positive numerics. The estimated standard
#'     errors.
#' @param alpha_tilde A matrix. The number of rows should be equal the
#'     length of betahat_ols. The number of columns should equal the
#'     number of hidden confounders.
#' @param tau2_seq The grid of variances of the mixing distributions
#'     if \code{mixing_dist = "normal"}. Only one of \code{tau2_seq}
#'     or \code{a_seq} and \code{b_seq} need be specified.
#' @param a_seq The grid of lower bounds for the uniforms if
#'     \code{mixing_dist} is one of the uniforms.
#' @param b_seq The grid of upper bounds for the uniforms if
#'     \code{mixing_dist} is one of the uniforms.
#' @param degrees_freedom The degrees of freedom of the t-distribution
#'     if \code{likelihood = "t"}.
#' @param plot_update A logical. Should I plot the the path of the
#'     log-likelihood (\code{TRUE}) or not (\code{FALSE})? Only
#'     applicable when \code{mixing_dist} is not \code{"normal"}.
#' @param lambda_seq A numeric vector with elements all greater than
#'     or equal to 1. These are the tuning parameters for the mixing
#'     proportions.
#'
#' @author David Gerard
#'
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., and Stephens, M. 2020. "Empirical Bayes shrinkage and false discovery rate estimation, allowing for unwanted variation", \emph{Biostatistics}, 21(1), 15-32 \doi{10.1093/biostatistics/kxy029}}
#' }
#'
#' @export
#'
mouthwash_second_step <-
  function(betahat_ols, S_diag, alpha_tilde,
           lambda_seq, tau2_seq = NULL,
           a_seq = NULL, b_seq = NULL,
           mixing_dist = c("normal", "uniform", "+uniform", "sym_uniform"),
           likelihood = c("normal", "t"),
           pi_init_type = c("zero_conc", "uniform", "random"),
           scale_var = TRUE,
           degrees_freedom = NULL,
           plot_update = FALSE,
           sprop = 0, var_inflate_pen = 0,
           verbose = TRUE) {

    ## Make sure input is correct ------------------------------------------
    mixing_dist  <- match.arg(mixing_dist)
    likelihood   <- match.arg(likelihood)
    pi_init_type <- match.arg(pi_init_type)

    ## Check df
    if (likelihood == "normal") {
      if (!is.null(degrees_freedom)) {
        if (degrees_freedom != Inf) {
         message('degrees_freedom provided but likelihood = "normal". Ignoring degrees_freedom.')
        }
      }
      degrees_freedom <- Inf
    } else if (likelihood == "t") {
      if (is.null(degrees_freedom)) {
        stop('likelihood = "t" but degrees_freedom not provided')
      }
    }

    if (mixing_dist == "normal") {
        assertthat::assert_that(!is.null(tau2_seq))
        M <- length(tau2_seq)
        zero_spot <- which(abs(tau2_seq) < 10 ^ -14)
    } else if (mixing_dist == "uniform" | mixing_dist == "+uniform" | mixing_dist == "sym_uniform") {
        assertthat::assert_that(!is.null(a_seq))
        assertthat::assert_that(!is.null(b_seq))
        M <- length(a_seq)
        assertthat::are_equal(length(b_seq), M)
        zero_spot <- which(abs(a_seq) < 10 ^ -14 & abs(b_seq) < 10 ^ -14)
    }
    assertthat::are_equal(length(zero_spot), 1)

    k <- ncol(alpha_tilde)
    assertthat::are_equal(length(betahat_ols), nrow(alpha_tilde))
    assertthat::are_equal(length(S_diag), length(betahat_ols))
    assertthat::are_equal(length(lambda_seq), M)
    assertthat::assert_that(sprop >= 0)
    assertthat::assert_that(var_inflate_pen >= 0)

    if (sprop == 1 & scale_var & var_inflate_pen < 10 ^ -6) {
        stop("if sprop is 1 and scale_var = TRUE, then var_inflate_pen should be > 0")
    }

    ## initialize parameters and run EM ------------------------------------
    if (verbose)
      cat("    + Estimating model parameters using EM.\n")
    timing <- system.time({
    z2_init <- matrix(stats::rnorm(k), ncol = 1)
    pi_init <- initialize_mixing_prop(pi_init_type = pi_init_type, zero_spot = zero_spot, M = M)

    if (likelihood == "normal" & mixing_dist == "normal") {
        pizxi_init <- c(pi_init, z2_init, 1)
        sqout <- SQUAREM::squarem(par = pizxi_init, fixptfn = normal_mix_fix_wrapper,
                                  objfn = normal_mix_llike_wrapper, betahat_ols = betahat_ols,
                                  S_diag = S_diag, alpha_tilde = alpha_tilde, tau2_seq = tau2_seq,
                                  lambda_seq = lambda_seq, scale_var = scale_var,
                                  control = list(tol = 10 ^ -4),
                                  var_inflate_pen = var_inflate_pen)
        pi_vals  <- sqout$par[1:M]
        z2_final <- sqout$par[(M + 1):(M + k)]
        xi_final <- sqout$par[M + k + 1]
    } else if (mixing_dist == "uniform" | mixing_dist == "+uniform" | mixing_dist == "sym_uniform") {
        opt_out <- mouthwash_coordinate(pi_init = pi_init, z_init = z2_init, xi_init = 1,
                                        betahat_ols = betahat_ols, S_diag = S_diag,
                                        alpha_tilde = alpha_tilde, a_seq = a_seq,
                                        b_seq = b_seq, lambda_seq = lambda_seq,
                                        degrees_freedom = degrees_freedom, scale_var = scale_var,
                                        plot_update = plot_update,
                                        var_inflate_pen = var_inflate_pen)
        pi_vals  <- opt_out$pi_vals
        z2_final <- opt_out$z2
        xi_final <- opt_out$xi
    }})
    if (verbose)
      cat("    + Computation took",timing["elapsed"],"seconds.\n")

    ## make mix object  ----------------------------------------------------
    if (verbose)
        cat("    + Generating adaptive shrinkage (ash) output.\n")
    timing <- system.time({
    if (mixing_dist == "uniform" | mixing_dist == "+uniform" | mixing_dist == "sym_uniform") {
        ghat <- ashr::unimix(pi = pi_vals, a = a_seq, b = b_seq)
    } else if (mixing_dist == "normal") {
        ghat <- ashr::normalmix(pi = pi_vals, mean = rep(0, M), sd = sqrt(tau2_seq))
    }

    ## For ashr compatibility -----------------------------------------------
    mixcompdist <- mixing_dist
    if (mixcompdist == "uniform") {
        mixcompdist <- "halfuniform"
    } else if (mixcompdist == "sym_uniform") {
        mixcompdist <- "uniform"
    }

    ashr_df <- degrees_freedom
    if (likelihood == "normal") {
        ashr_df <- NULL
    }

    ## deal with non-zero sprop before returning ash output -----------------
    ## Recall that betahat_ols, alpha_tilde_ols, and S_diag are
    ## actually modified based on sprop before being sent to
    ## mouthwash_second_step. The following udoes this modification
    ## before sending these values to ashr::ash.workhorse to obtain
    ## summary values. Note that sprop for me = alpha for ashr.
    if (sprop > 0) {
      S_diag_real      <- S_diag ^ (1 / (1 - sprop))
      sgamma           <- S_diag_real ^ (sprop / 2)
      betahat_ols_real <- betahat_ols * sgamma
      alpha_tilde_real <- alpha_tilde * sgamma
    } else {
        betahat_ols_real <- betahat_ols
        alpha_tilde_real <- alpha_tilde
        S_diag_real      <- S_diag
    }

    az <- alpha_tilde_real %*% z2_final

    ## Call ashr for summaries ---------------------------------------------
    val <- ashr::ash.workhorse(betahat = c(betahat_ols_real - az),
                               sebetahat = c(sqrt(xi_final * S_diag_real)),
                               df = ashr_df,
                               prior = "nullbiased",
                               nullweight = lambda_seq[zero_spot],
                               g = ghat,
                               fixg = TRUE,
                               mixcompdist = mixcompdist,
                               alpha = sprop) ## really need this
    val <- c(val, list(pi0 = pi_vals[zero_spot]))
    val <- c(val, list(z2 = z2_final))
    val <- c(val, list(xi = xi_final))
    })
    if (verbose)
      cat("    + Computation took",timing["elapsed"],"seconds.\n")
    return(val)
}

#' Function for initializing mixing proportions.
#'
#' @inheritParams mouthwash
#' @param M The number of mixing distributions.
#' @param zero_spot The location of the zero.
#'
#' @author David Gerard
initialize_mixing_prop <- function(pi_init_type = c("zero_conc", "uniform", "random"),
                                   zero_spot, M) {
    pi_init_type <- match.arg(pi_init_type)
    assertthat::assert_that(M > 0)

    if (pi_init_type == "zero_conc") {
        pi_init <- rep(0.1 / (M - 1), length = M)
        pi_init[zero_spot] <- 0.9
    } else if (pi_init_type == "uniform") {
        pi_init <- rep(1 / M, length = M)
    } else if (pi_init_type == "random") {
        pi_init <- stats::runif(n = M)
        pi_init <- pi_init / sum(pi_init)
    }
    return(pi_init)
}



#' Default way to set grid of variances.
#'
#' @param betahat_ols The OLS estimates of the coefficients of interest.
#' @param S_diag The standard error estimates.
#'
#' @author David Gerard
get_grid_var <- function(betahat_ols, S_diag) {

    ## Check input ----------------------------------------------------
    assertthat::are_equal(length(betahat_ols), length(S_diag))
    assertthat::assert_that(all(S_diag > 0))

    ## default grid to be same as in ASH ------------------------------
    tau2_min <- min(S_diag) / 100
    tau2_max <- 16 * max(betahat_ols ^ 2 - S_diag) ## used to be 4 * max...
    if (tau2_max < 0) {
        tau2_max <- 64 * tau2_min
    }
    tau2_current <- tau2_min
    tau2_seq <- c(0, tau2_current)
    mult_fact <- 2
    while (tau2_current <= tau2_max) {
        tau2_current <- tau2_current * mult_fact
        tau2_seq <- c(tau2_seq, tau2_current)
    }
    M <- length(tau2_seq)
    return(list(tau2_seq = tau2_seq, M = M))
}


#' Wrapper for dt with a non-zero mean and non-1 scale parameter.
#'
#' @param x A numeric. Where to evaluate the density function.
#' @param df A positive numeric. The degrees of freedom.
#' @param mean A numeric. The mean of the t. Defaults to 0.
#' @param sd A positive numeric. The standard deviation of the
#'     t. Defaults to 1.
#' @param log A logical. Should we return the log-density
#'     (\code{TRUE}) or the density (\code{FALSE})?
dt_wrap <- function(x, df, mean = 0, sd = 1, log = FALSE) {
    x_new <- (x - mean) / sd
    if (!log) {
        dval <- stats::dt(x_new, df = df) / sd
    } else {
        dval <- stats::dt(x_new, df = df, log = TRUE) - log(sd)
    }
    return(dval)
}

#' Wrapper for pt with a non-zero mena and non-1 scale parameter.
#'
#' @param x A numeric. Where to evaluate the cdf.
#' @inheritParams dt_wrap
#'
pt_wrap <- function(x, df, mean = 0, sd = 1) {
    x_new <- (x - mean) / sd
    pval <- stats::pt(x_new, df = df)
    return(pval)
}

#' Returns adjusted sebetahat's based on t likelihood so that we can use a normal likelihood.
#'
#' @param betahat The estimates of the effects.
#' @param sebetahat The estimates of the standard errors of \code{betahat}.
#' @param df The degrees of freedom of the t. Can either be of length 1 or the same length of
#'     \code{betahat}.
#'
#' @author David Gerard
#'
adjust_by_t <- function(betahat, sebetahat, df) {
  ## Check input -----------------------------------------
  assertthat::are_equal(length(betahat), length(sebetahat))
  assertthat::assert_that(length(betahat) == length(df) | length(df) == 1)
  assertthat::assert_that(df > 0)

  ## Convert t to z --------------------------------------
  zstats <- stats::qnorm(stats::pt(q = betahat / sebetahat, df = df))
  snew <- c(betahat / zstats)

  return(snew)
}
