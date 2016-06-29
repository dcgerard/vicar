#' Calibrated LEAPP where the tuning parameter is chosen by a modified
#' cross-validation.
#'
#' Doesn't work too well.
#'
#' @inheritParams vruv4
#'
#' @author David Gerard
#'
vicarius_leapp <- function(Y, X, k = NULL,
                          cov_of_interest = ncol(X),
                          limmashrink = FALSE,
                          include_intercept = TRUE,
                          fa_func = pca_naive, fa_args = list()) {

    if(!requireNamespace("leapp", quietly = TRUE)) {
        stop(paste("LEAPP needs to be installed to use", match.call()))
    }

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(cov_of_interest >= 1 & cov_of_interest <= ncol(X))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.logical(limmashrink))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(is.null(fa_args$Y))
    assertthat::assert_that(is.null(fa_args$r))
    assertthat::assert_that(is.function(fa_func))

    rotate_out <- rotate_model(Y = Y, X = X, k = k,
                               cov_of_interest = cov_of_interest,
                               include_intercept = include_intercept,
                               limmashrink = limmashrink, fa_func = fa_func,
                               fa_args = fa_args)

    betahat_ols     <- rotate_out$betahat_ols
    alpha_scaled    <- rotate_out$alpha_scaled
    sig_diag_scaled <- rotate_out$sig_diag_scaled
    k               <- rotate_out$k

    betaSigma <- betahat_ols / sqrt(sig_diag_scaled)

    alphaSigma <- diag(1 / sqrt(sig_diag_scaled)) %*% alpha_scaled

    H <- alphaSigma %*% solve(t(alphaSigma) %*% alphaSigma) %*% t(alphaSigma)

    result <- cleapp(X = alphaSigma, Y = betaSigma, H = H)
    return(result)
}


#' Cleaned version of LEAPP.
#'
#' Doesn't work too well.
#'
#' @inheritParams vicarius_leapp
#' @param method Either "hard" or "soft"
#' @param H Equal to X(X'X)^{-1}X'
#'
cleapp <- function(X, Y, H, method = "hard") {

    if(!requireNamespace("leapp", quietly = TRUE)) {
        stop(paste("LEAPP needs to be installed to use", match.call()))
    }

    assertthat::assert_that(is.matrix(X))
    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(H))

    r        <- ncol(X)
    N        <- nrow(Y)
    betaInit <- MASS::rlm(Y ~ X - 1)$coefficients
    tmp      <- t(diag(ncol(H)) - H) %*% Y/sqrt(1 - diag(H))
    lambdas  <- seq(round(norm(tmp, "I")/1 + 1), 0, length = 100)

    lambdas <- seq(30, 5, length = 100)
    err_vec <- rep(NA, length = length(lambdas))
    for (index in 1:length(lambdas)) {
        sigma <- lambdas[index] / sqrt(2 * log(N))
        err_vec[index] <- leapp_cross(X = X, Y = Y, sigma = sigma, method = method)
        cat(err_vec[index], "\n")
    }

    spout <- stats::smooth.spline(x = lambdas / sqrt(2 * log(N)), y = err_vec)
    sigma_final <- spout$x[which.min(spout$y)]

    vout <- vicarius_principis_IPOD(X = X, Y = Y, H = H, sigma = sigma_final)
    return(vout)
}

#' Modified version of leapp::IPOD to run a version of
#' cross-validation to choose the tuning parameter.
#'
#' Doesn't work too well.
#'
#' @inheritParams cleapp
#' @param sigma The variance.
#' @param TOL The stopping criterion tolerance.
#'
vicarius_principis_IPOD <- function (X, Y, H, sigma, method = "hard",
                                     TOL = 1e-04)
{
    if (is.null(X)) {
        stop("X cannot be null")
    }

    if(!requireNamespace("leapp", quietly = TRUE)) {
        stop(paste("LEAPP needs to be installed to use", match.call()))
    }


    assertthat::assert_that(is.matrix(X))
    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(H))
    assertthat::assert_that(sigma >= 0)
    assertthat::are_equal(length(sigma), 1)

    p     <- nrow(Y)
    ress  <- NULL
    gamma <- NULL

    betaInit <- MASS::rlm(Y ~ X - 1)$coefficients

    sigma  <- sigma/sqrt(2 * log(p))
    result <- leapp::IPODFUN(X, Y, H, sigma, betaInit, method = method,
                            TOL = TOL)
    gamma <- result$gamma
    ress <- result$ress

    beta <- solve(t(X) %*% X) %*% t(X) %*% (Y - gamma)

    return(list(gamma = gamma, beta = beta))
}


#' 10-fold cross-validation-ish method.
#'
#' Doesn't work too well.
#'
#' @inheritParams vicarius_principis_IPOD
#' @param num_sets The fold of the cross-validation.
#' @param betaInit The initial value of beta.
#'
leapp_cross <- function(X, Y, sigma, betaInit, method = "hard",
                        TOL = 1e-04, num_sets = 10) {

    if(!requireNamespace("leapp", quietly = TRUE)) {
        stop(paste("LEAPP needs to be installed to use", match.call()))
    }

    p <- nrow(Y)
    rand_perm <- sample(1:p)

    locs <- round(seq(0, p, length = num_sets + 1))
    start_vals <- locs[1:(length(locs) - 1)] + 1
    end_vals <- locs[2:length(locs)]

    obs_pos_list <- list()
    for (index in 1:num_sets) {
        obs_pos_list[[index]] <- rand_perm[start_vals[index]:end_vals[index]]
    }


    err <- 0
    for (index in 1:num_sets) {
        current_obs <- (1:p)[-obs_pos_list[[index]]]
        nv <- length(obs_pos_list[[index]])
        Xsub <- X[current_obs, , drop = FALSE]
        Ysub <- Y[current_obs, , drop = FALSE]
        Hsub <- Xsub %*% solve(t(Xsub) %*% Xsub) %*% t(Xsub)
        vout <- vicarius_principis_IPOD(X = Xsub, Y = Ysub, H = Hsub,
                                        sigma = sigma, method = method)

        egamma  <- mean(vout$gamma)
        e2gamma <- mean(vout$gamma ^ 2)

        Xval <- X[obs_pos_list[[index]], , drop = FALSE]
        Yval <- Y[obs_pos_list[[index]], , drop = FALSE]
        rval <- Yval - Xval %*% vout$beta
        current_err <- sum(rval ^ 2) -  2 * sum(rval) * egamma + nv * e2gamma
        err <- err + current_err
    }

    return(err)
}
