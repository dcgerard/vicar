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
#' The assumed mode is \deqn{Y = X\beta + Z\alpha + E.} \eqn{Y} is a
#' \eqn{n} by \code{p} matrix of response varaibles. For example, each
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
#' part only contains nuisance parameters, the second part contains
#' the coefficients of interest, and the third part contains the
#' confounders. \code{mouthwash} applies a user-provided factor
#' analysis to the third part to estimate the confounding factors,
#' then runs an EM (or coordinate-ascent) algorithm on the second part
#' to estimate the coefficients of interest.
#'
#' Many forms of factor analyses are avaiable in this package. The
#' default is PCA with the column-wise residual mean-squares as the
#' estimates of the column-wise variances.
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
#'     corresponding the covariates of interest. Mostly output for
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
#'     }
#'
#' @export
#'
#' @author David Gerard
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
                      sprop = 0) {


    ## Make sure input is correct --------------------------------------------------------
    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(is.numeric(cov_of_interest))
    assertthat::are_equal(length(cov_of_interest), 1)
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.logical(limmashrink))
    assertthat::assert_that(is.function(fa_func))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(lambda0 >= 1)
    assertthat::assert_that(sprop >= 0)

    likelihood   <- match.arg(likelihood)
    mixing_dist  <- match.arg(mixing_dist)
    pi_init_type <- match.arg(pi_init_type)
    lambda_type  <- match.arg(lambda_type)

    if (likelihood == "t" & mixing_dist == "normal") {
        stop("normal mixtures not implemented for t-likelihood yet (or likely ever).")
    }
    if (scale_var & sprop == 1) {
        stop("sprop cannot be 1 when scale_var is TRUE")
    }


    ## Rotate ---------------------------------------------------------------------------
    rotate_out <- rotate_model(Y = Y, X = X, k = k,
                               cov_of_interest = cov_of_interest,
                               include_intercept = include_intercept,
                               limmashrink = limmashrink, fa_func = fa_func,
                               fa_args = fa_args, do_factor = TRUE)


    ## Deal with degrees of freedom -----------------------------------------------------
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

    ## rescale alpha and sig_diag by R22 to get data for second step ------------------
    alpha_tilde <- rotate_out$alpha / c(rotate_out$R22)
    S_diag      <- c(rotate_out$sig_diag / c(rotate_out$R22 ^ 2))
    betahat_ols <- matrix(rotate_out$betahat_ols, ncol = 1)

    ## Exchangeable versions of the models ---------------------------------------------
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

    ## Set grid and penalties ---------------------------------------------------------
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
    }

    ## Run MOUTHWASH -------------------------------------------------------------
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
                                 sprop = sprop)

    ## Estimate rest of the hidden confounders -----------------------------------
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
#'
#'
#' @author David Gerard
#'
#' @export
#'
mouthwash_second_step <- function(betahat_ols, S_diag, alpha_tilde,
                                  lambda_seq, tau2_seq = NULL,
                                  a_seq = NULL, b_seq = NULL,
                                  mixing_dist = c("normal", "uniform", "+uniform", "sym_uniform"),
                                  likelihood = c("normal", "t"),
                                  pi_init_type = c("zero_conc", "uniform", "random"),
                                  scale_var = TRUE, degrees_freedom,
                                  plot_update = FALSE,
                                  sprop = 0) {

    ## Make sure input is correct -------------------------------------------------
    mixing_dist  <- match.arg(mixing_dist)
    likelihood   <- match.arg(likelihood)
    pi_init_type <- match.arg(pi_init_type)

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

    if (sprop == 1 & scale_var) {
        stop("if sprop is 1, then scale_var cannot be TRUE")
    }

    ## initialize parameters and run EM --------------------------------------------
    z2_init <- matrix(stats::rnorm(k), ncol = 1)
    pi_init <- initialize_mixing_prop(pi_init_type = pi_init_type, zero_spot = zero_spot, M = M)


    if (likelihood == "normal" & mixing_dist == "normal") {
        pizxi_init <- c(pi_init, z2_init, 1)
        sqout <- SQUAREM::squarem(par = pizxi_init, fixptfn = normal_mix_fix_wrapper,
                                  objfn = normal_mix_llike_wrapper, betahat_ols = betahat_ols,
                                  S_diag = S_diag, alpha_tilde = alpha_tilde, tau2_seq = tau2_seq,
                                  lambda_seq = lambda_seq, scale_var = scale_var,
                                  control = list(tol = 10 ^ -4))
        pi_vals  <- sqout$par[1:M]
        z2_final <- sqout$par[(M + 1):(M + k)]
        xi_final <- sqout$par[M + k + 1]
    } else if (mixing_dist == "uniform" | mixing_dist == "+uniform" | mixing_dist == "sym_uniform") {
        ## sqout <- SQUAREM::squarem(par = pizxi_init, fixptfn = uniform_mix_fix_wrapper,
        ##                           objfn = uniform_mix_llike_wrapper, betahat_ols = betahat_ols,
        ##                           S_diag = S_diag, alpha_tilde = alpha_tilde, a_seq = a_seq,
        ##                           b_seq = b_seq, lambda_seq = lambda_seq,
        ##                           degrees_freedom = degrees_freedom, scale_var = scale_var,
        ##                           control = list(tol = 10 ^ -4))

        opt_out <- mouthwash_coordinate(pi_init = pi_init, z_init = z2_init, xi_init = 1,
                                        betahat_ols = betahat_ols, S_diag = S_diag,
                                        alpha_tilde = alpha_tilde, a_seq = a_seq,
                                        b_seq = b_seq, lambda_seq = lambda_seq,
                                        degrees_freedom = degrees_freedom, scale_var = scale_var,
                                        plot_update = plot_update)
        pi_vals  <- opt_out$pi_vals
        z2_final <- opt_out$z2
        xi_final <- opt_out$xi
    }




    ## make mix object  --------------------------------------------------------------
    if (mixing_dist == "uniform" | mixing_dist == "+uniform" | mixing_dist == "sym_uniform") {
        ghat <- ashr::unimix(pi = pi_vals, a = a_seq, b = b_seq)
    } else if (mixing_dist == "normal") {
        ghat <- ashr::normalmix(pi = pi_vals, mean = rep(0, M), sd = sqrt(tau2_seq))
    }

    ## For ashr compatibility -------------------------------------------------------
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


    ## deal with non-zero sprop before returning ash output -------------------------
    ## Recall that betahat_ols, alpha_tilde_ols, and S_diag are
    ## actually modified based on sprop before being sent to
    ## mouthwash_second_step. The following udoes this modification
    ## before sending these values to ashr::ash.workhorse to obtain
    ## summary values. Note that sprop for me = alpha for ashr.
    if (sprop > 0) {
        sgamma <- S_diag ^ (sprop / 2)
        betahat_ols_real <- betahat_ols * sgamma
        alpha_tilde_real <- alpha_tilde * sgamma
        S_diag_real      <- S_diag ^ (1 + sprop)
    } else {
        betahat_ols_real <- betahat_ols
        alpha_tilde_real <- alpha_tilde
        S_diag_real      <- S_diag
    }

    az <- alpha_tilde_real %*% z2_final

    ## Call ashr for summaries ------------------------------------------------------
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
    tau2_max <- 4 * max(betahat_ols ^ 2 - S_diag)
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


#' Wrapper for dt with a non-zero mena and non-1 scale parameter.
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
