## Functions for performing MOUTHWASH. This is a re-coding from
## succotashr because I'm slightly better at programming now.


#' MOUTHWASH: Maximize Over Unobservables To Help With Adaptive SHrinkage.
#'
#'
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
#'
#' @export
#'
#' @author David Gerard
mouthwash <- function(Y, X, k = NULL, cov_of_interest = ncol(X),
                      include_intercept = TRUE, limmashrink = TRUE,
                      fa_func = pca_naive, fa_args = list(),
                      likelihood = c("normal", "t"),
                      mixing_dist = c("normal", "uniform", "+uniform", "sym_uniform", "normal"),
                      lambda_type = c("zero_conc", "uniform"),
                      pi_init_type = c("zero_conc", "uniform", "random"),
                      degrees_freedom = NULL,
                      pi_init = NULL, grid_seq = NULL, lambda_seq = NULL,
                      lambda0 = 10,
                      scale_var = TRUE) {


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

    likelihood   <- match.arg(likelihood)
    mixing_dist  <- match.arg(mixing_dist)
    pi_init_type <- match.arg(pi_init_type)
    lambda_type  <- match.arg(lambda_type)

    if (likelihood == "t" & mixing_dist == "normal") {
        stop("normal mixtures not implemented for t-likelihood yet (or likely ever).")
    }

    if (ncol(Y) > 100 & mixing_dist != "normal") {
        warning("uniform mixtures somewhat janky for moderate to large p")
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

    ## Set grid and penalties ---------------------------------------------------------
    if (!is.null(lambda_seq) & is.null(grid_seq)) {
        stop("lambda_seq specified but grid_seq is NULL")
    }

    if (is.null(grid_seq)) {
        grid_vals <- get_grid_var(betahat_ols = betahat_ols, S_diag = S_diag)
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

    val <- mouthwash_second_step(betahat_ols = betahat_ols, S_diag = S_diag,
                                 alpha_tilde = alpha_tilde, tau2_seq = tau2_seq,
                                 a_seq = a_seq, b_seq = b_seq,
                                 degrees_freedom = degrees_freedom,
                                 lambda_seq = lambda_seq, mixing_dist = mixing_dist,
                                 likelihood = likelihood, pi_init_type = pi_init_type,
                                 scale_var = scale_var)
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
                                  scale_var = TRUE, degrees_freedom) {

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

    ## initialize parameters and run EM --------------------------------------------
    z2_init <- matrix(stats::rnorm(k), ncol = 1)
    pi_init <- initialize_mixing_prop(pi_init_type = pi_init_type, zero_spot = zero_spot, M = M)
    pizxi_init <- c(pi_init, z2_init, 1)

    if (likelihood == "normal" & mixing_dist == "normal") {
        sqout <- SQUAREM::squarem(par = pizxi_init, fixptfn = normal_mix_fix_wrapper,
                                  objfn = normal_mix_llike_wrapper, betahat_ols = betahat_ols,
                                  S_diag = S_diag, alpha_tilde = alpha_tilde, tau2_seq = tau2_seq,
                                  lambda_seq = lambda_seq, scale_var = scale_var,
                                  control = list(tol = 10 ^ -4))

        normal_mix_fix_wrapper(pizxi_vec = pizxi_init, betahat_ols = betahat_ols,
                               S_diag = S_diag, alpha_tilde = alpha_tilde, tau2_seq = tau2_seq,
                               lambda_seq = lambda_seq, scale_var = scale_var)

    } else if (mixing_dist == "uniform" | mixing_dist == "+uniform" | mixing_dist == "sym_uniform") {
        sqout <- SQUAREM::squarem(par = pizxi_init, fixptfn = uniform_mix_fix_wrapper,
                                  objfn = uniform_mix_llike_wrapper, betahat_ols = betahat_ols,
                                  S_diag = S_diag, alpha_tilde = alpha_tilde, a_seq = a_seq,
                                  b_seq = b_seq, lambda_seq = lambda_seq,
                                  degrees_freedom = degrees_freedom, scale_var = scale_var,
                                  control = list(tol = 10 ^ -4))
    }

    pi_vals  <- sqout$par[1:M]
    z2_final <- sqout$par[(M + 1):(M + k)]
    xi_final <- sqout$par[M + k + 1]

    az <- alpha_tilde %*% z2_final

    ## ash summaries --------------------------------------------------------------
    if (mixing_dist == "uniform" | mixing_dist == "+uniform" | mixing_dist == "sym_uniform") {
        ghat <- ashr::unimix(pi = pi_vals, a = a_seq, b = b_seq)
    } else if (mixing_dist == "normal") {
        ghat <- ashr::normalmix(pi = pi_vals, mean = rep(0, M), sd = sqrt(tau2_seq))
    }

    if (likelihood == "normal" & mixing_dist == "normal") {
        data <- ashr::set_data(betahat = c(betahat_ols - az),
                               sebetahat = c(sqrt(xi_final * S_diag)),
                               lik = ashr::normal_lik())
    } else {
        data <- ashr::set_data(betahat = c(betahat_ols - az),
                               sebetahat = c(sqrt(xi_final * S_diag)),
                               lik = ashr::t_lik(degrees_freedom))
    }

    val <- list()
    val <- c(val, list(fitted_g = ghat))
    val <- c(val, list(loglik = ashr::calc_loglik(ghat, data)))
    val <- c(val, list(logLR = ashr::calc_logLR(ghat, data)))
    val <- c(val, list(data = data))
    val <- c(val, list(pi0 = pi_vals[zero_spot]))
    val <- c(val, list(z2 = z2_final))
    val <- c(val, list(xi = xi_final))
    NegativeProb  <- ashr:::calc_np(g = ghat, data = data)
    PositiveProb  <- ashr:::calc_pp(g = ghat, data = data)
    lfsr          <- ashr:::calc_lfsr(g = ghat, data = data)
    svalue        <- ashr:::calc_svalue(g = ghat, data = data)
    lfdr          <- ashr:::calc_lfdr(g = ghat, data = data)
    qvalue        <- ashr:::calc_qvalue(g = ghat, data = data)
    PosteriorMean <- ashr:::calc_pm(g = ghat, data = data)
    PosteriorSD   <- ashr:::calc_psd(g = ghat, data = data)
    result <- data.frame(NegativeProb, PositiveProb, lfsr, svalue, lfdr, qvalue,
                         PosteriorMean, PosteriorSD)
    val <- c(val, list(result = result))

    class(val) <- "ash"
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
