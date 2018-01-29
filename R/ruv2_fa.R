#' PCA when first \code{vr} rows have a variance multiplicatively
#' different from the rest of the rows.
#'
#' Modified truncated SVD. The variance estimates are just the
#' column-wise mean-squares of the last \eqn{n - vr} rows. This form
#' of factor analysis is mostly for variance inflation with RUV2.
#'
#' This doesn't work too well. I think the Z's and sigmas need to be
#' estimated jointly.
#'
#'
#' @param Y A matrix of numerics. The data.
#' @param r the rank.
#' @param vr The number of the first few rows whose variances differ
#'     by a multiplicative factor.
#' @param mle A logical. Should we run an MLE on the residuals of PCA
#'     (\code{TRUE}) or just use a two-step estimator (\code{FALSE}).
#'
#' @export
#'
#' @author David Gerard
pca_ruv2 <- function(Y, r, vr, mle = FALSE) {
    assertthat::assert_that(is.matrix(Y))
    assertthat::are_equal(length(r), 1)
    assertthat::assert_that(r >= 0 & r < min(dim(Y)))
    assertthat::assert_that(vr > 0 & vr < nrow(Y))

    p <- ncol(Y)
    n <- nrow(Y)

    if (r == 0) {
        Gamma <- NULL
        Z <- NULL
        resids <- Y
    } else {
        svd_Y <- svd(Y)
        Gamma <- svd_Y$v[, 1:r, drop = FALSE] %*% diag(svd_Y$d[1:r], r, r) /
            sqrt(n)
        Z <- sqrt(n) * svd_Y$u[, 1:r, drop = FALSE]
        resids <- Y - Z %*% t(Gamma)
    }

    ## MLE to find variances
    R1 <- resids[1:vr, ]
    R2 <- resids[(vr + 1): n, ]

    r1 <- colSums(R1 ^ 2)
    r2 <- colSums(R2 ^ 2)

    Sigma_init <- r2 / (n - vr)
    lambda_init <- mean(r1 / Sigma_init) / vr

    if (mle) {
        sqout <- SQUAREM::squarem(par = c(Sigma_init, lambda_init),
                                  fixptfn = pcaruv2_fix,
                                  objfn = pcaruv2_obj, r1 = r1,
                                  r2 = r2, n = n, vr = vr)
        sig_diag <- sqout$par[-length(sqout$par)]
        lambda   <- sqout$par[length(sqout$par)]
    } else {
        sig_diag <- Sigma_init
        lambda   <- lambda_init
    }

    return(list(alpha = Gamma, Z = Z, sig_diag = sig_diag, lambda = lambda))
}

#' Fix point for mle in pca_ruv2.
#'
#' This is mostly for use in SQUAREM.
#'
#' @param sig_lambda A vector. All but the first elements are the
#'     current values of Sigma. The last element is the current valur
#'     of lambda.
#' @param r1 The sum of squares of the first vr columns of the
#'     residuals.
#' @param r2 The sum of squares of the last n - vr columns of the
#'     residuals.
#' @param n The number of samples.
#' @param vr The number of columns where variance inflation is in
#'     effect.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{pca_ruv2}}, \code{\link{pcaruv2_obj}} for the
#'     objective function that this maximizes.
pcaruv2_fix <- function(sig_lambda, r1, r2, n, vr) {
    Sigma  <- sig_lambda[-length(sig_lambda)]
    lambda <- sig_lambda[length(sig_lambda)]

    Sigma_new <- (r1 / lambda + r2) / n
    lambda_new <- mean(r1 / Sigma_new) / vr

    return(c(Sigma_new, lambda_new))
}

#' The objective function for mle in pca_ruv2.
#'
#' This is mostly for use in SQUAREM.
#'
#' @inheritParams pcaruv2_fix
#'
#' @author David Gerard
#'
#' @seealso \code{\link{pca_ruv2}} \code{\link{pcaruv2_fix}} for the
#'     fixed point iteration that maximizes this objective function.
pcaruv2_obj <- function(sig_lambda, r1, r2, n, vr) {
    Sigma  <- sig_lambda[-length(sig_lambda)]
    lambda <- sig_lambda[length(sig_lambda)]

    tr1 <- sum(r1 / Sigma) / lambda
    tr2 <- sum(r2 / Sigma)

    p <- length(Sigma)

    llike <- -vr * p * log(lambda) / 2 -
        n * sum(log(Sigma)) / 2 -
        (tr1 + tr2) / 2
    return(llike)
}


#' Quasi-mle when first \code{vr} rows of \code{Y} have variances
#' which differ by a scale factor from the rest of the rows.
#'
#' @inheritParams pca_ruv2
#' @param tol The function tolerance for stopping. A positive numeric.
#' @param maxit The maximum number of iterations to run.
#'
#' @return \code{alpha} A matrix of the estimated factor loadings.
#'
#'         \code{Z} A matrix of the estimated factors.
#'
#'         \code{sig_diag} A vector of the estimated column-wise
#'         variances.
#'
#'         \code{lambda} A numeric scalar of the estimated variance
#'         inflation parameter.
#'
#' @export
#'
#' @author David Gerard
qmle_ruv2 <- function(Y, r, vr, tol = 10 ^ -3, maxit = 1000) {

    assertthat::assert_that(tol > 0)

    ## Get starting values using CATE's fa.em
    cate_faout <- cate::fa.em(Y = Y, r = r)

    alpha    <- t(cate_faout$Gamma)
    sig_diag <- cate_faout$Sigma
    lambda   <- 2

    Y1 <- Y[1:vr, , drop = FALSE]
    Y2 <- Y[(vr + 1):nrow(Y), , drop = FALSE]

    ## Run through a few iterations with lambda = 2
    for (index in 1:20) {
        uout <- update_sig_alpha(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                                 Y1 = Y1, Y2 = Y2)
        alpha <- uout$alpha
        sig_diag <- uout$sig_diag
    }

    llike_vec <- c()
    llike <- qmle_obj(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                      Y1 = Y1, Y2 = Y2)

    for (index in 1:maxit) {
        lambda_old <- lambda
        alpha_old <- alpha
        sig_diag_old <- sig_diag
        llike_old <- llike
        uout <- update_sig_alpha(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                                 Y1 = Y1, Y2 = Y2)
        lambda <- update_lambda(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                                    Y1 = Y1, Y2 = Y2)
        alpha <- uout$alpha
        sig_diag <- uout$sig_diag

        llike <- qmle_obj(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                          Y1 = Y1, Y2 = Y2)
        llike_vec <- c(llike_vec, llike)
        if (abs(llike / llike_old - 1) < tol) {
            break
        }
    }


    ## Slow old way to get Z
    ## Z <- Y %*% diag(1 / sig_diag) %*% t(alpha) %*%
    ##     solve(alpha %*% diag(1 / sig_diag) %*% t(alpha))

    ## New fast way to get G
    sig_inv <- 1 / sig_diag
    asig <- sweep(alpha, MARGIN = 2, sig_inv, `*`)
    asiga <- tcrossprod(alpha, asig)
    eigenAsiga <- eigen(asiga, symmetric = TRUE)
    YSA <- tcrossprod(Y, asig)
    Z <- tcrossprod(YSA, tcrossprod(eigenAsiga$vectors,
                                    sweep(eigenAsiga$vectors,
                                          MARGIN = 2,
                                          1 / eigenAsiga$values,
                                          FUN = `*`)))

    return(list(alpha = alpha, sig_diag = sig_diag, lambda = lambda, Z = Z))
}


#' Quasi-mle when first \code{vr} rows of \code{Y} have variances
#' which differ by a scale factor from the rest of the rows.
#'
#' This one will iterate through a grid of lambdas and choose the one
#' with the highest likelihood.
#'
#' @inheritParams pca_ruv2
#' @param tol The function tolerance for stopping. A positive numeric.
#' @param maxit The maximum number of iterations to run.
#' @param gridsize A positive integer. The size of the grid for the
#'     variance inflation parameter.
#' @param gridlow A positive numeric. The minimum variance inflation
#'     parameter to try out.
#' @param gridhigh A positive numeric. The maximum variance inflation
#'     parameter to try out.
#' @param plot_llike A logical. Should I plot the profiled
#'     log-likelihood?
#'
#' @return \code{alpha} A matrix of the estimated factor loadings.
#'
#'         \code{Z} A matrix of the estimated factors.
#'
#'         \code{sig_diag} A vector of the estimated column-wise
#'         variances.
#'
#'         \code{lambda} A numeric scalar of the estimated variance
#'         inflation parameter.
#'
#' @export
#'
#' @author David Gerard
qmle_ruv2_lambda_grid <- function(Y, r, vr, gridsize = 20, gridlow = 1,
                                  gridhigh = 10, tol = 10 ^ -3, maxit = 1000,
                                  plot_llike = FALSE) {

    assertthat::assert_that(tol > 0)
    assertthat::assert_that(gridlow < gridhigh)
    assertthat::assert_that(gridlow > 0)

    Y1 <- Y[1:vr, , drop = FALSE]
    Y2 <- Y[(vr + 1):nrow(Y), , drop = FALSE]

    lambda_vec <- seq(gridlow, gridhigh, length = gridsize)

    ## Get starting values using CATE's fa.em
    cate_faout <- cate::fa.em(Y = Y, r = r)

    ## run through each lambda
    llike_vec <- rep(NA, length = gridsize)
    for(lambda_index in 1:gridsize) {
        lambda   <- lambda_vec[lambda_index]
        alpha    <- t(cate_faout$Gamma)
        sig_diag <- cate_faout$Sigma

        llike <- qmle_obj(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                          Y1 = Y1, Y2 = Y2)

        for (index in 1:maxit) {
            alpha_old <- alpha
            sig_diag_old <- sig_diag
            llike_old <- llike
            uout <- update_sig_alpha(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                                     Y1 = Y1, Y2 = Y2)
            alpha <- uout$alpha
            sig_diag <- uout$sig_diag

            llike <- qmle_obj(alpha = alpha, sig_diag = sig_diag, lambda = lambda,
                              Y1 = Y1, Y2 = Y2)
            if (abs(llike / llike_old - 1) < tol) {
                break
            }
        }

        if (lambda_index == 1) {
            alpha_final <- alpha
            sig_diag_final <- sig_diag
            lambda_final <- lambda
        } else if (llike > max(llike_vec, na.rm = TRUE)) {
            alpha_final <- alpha
            sig_diag_final <- sig_diag
            lambda_final <- lambda
        }
        llike_vec[lambda_index] <- llike
    }

    which_lambda <- match(lambda_final, lambda_vec)

    if (which_lambda == 1) {
        warning("smallest lambda chosen. Maybe you should choose a smaller gridlow.")
    } else if (which_lambda == gridsize) {
        warning("largest lambda chosen. Maybe you should choose a larger gridhigh.")
    }

    ## Slow old way to get Z
    ## Z <- Y %*% diag(1 / sig_diag) %*% t(alpha) %*%
    ##     solve(alpha %*% diag(1 / sig_diag) %*% t(alpha))

    ## New fast way to get G
    sig_inv <- 1 / sig_diag_final
    asig <- sweep(alpha_final, MARGIN = 2, sig_inv, `*`)
    asiga <- tcrossprod(alpha_final, asig)
    eigenAsiga <- eigen(asiga, symmetric = TRUE)
    YSA <- tcrossprod(Y, asig)
    Z <- tcrossprod(YSA, tcrossprod(eigenAsiga$vectors,
                                    sweep(eigenAsiga$vectors,
                                          MARGIN = 2,
                                          1 / eigenAsiga$values,
                                          FUN = `*`)))

    if (plot_llike) {
        graphics::plot(lambda_vec, llike_vec,
                       xlab = "variance inflation parameter",
                       ylab = "log-likelihood",
                       main = "Profiled log-likelihood for lambda",
                       type = "l")
    }

    return(list(alpha = alpha_final, sig_diag = sig_diag_final, lambda = lambda_final, Z = Z))
}

#' Wrapper for \code{\link[stats]{optim}} and \code{\link{qmle_obj}}.
#'
#' Uses Brent's method to obtimize the variance inflation parameter
#' given sig_diag and alpha.
#'
#' @inheritParams qmle_obj
#' @param maxit The number of iterations to run. Shouldn't be too
#'     large since this is only an intermediate step.
#'
#' @author David Gerard
update_lambda <- function(alpha, sig_diag, lambda, Y1, Y2, maxit = 10) {
    oout <- stats::optim(par = lambda, fn = qmle_obj, method = "Brent",
                         lower = 0, upper = 10, control = list(fnscale = -1, maxit = maxit),
                         alpha = alpha, sig_diag = sig_diag, Y1 = Y1, Y2 = Y2)
    return(oout$par)
}


#' Basic no optimized code update for sig_diag and alpha.
#'
#' @inheritParams qmle_obj
#'
#' @author David Gerard
update_sig_alpha_basic <- function(alpha, sig_diag, lambda, Y1, Y2) {
    r <- nrow(alpha)
    m <- nrow(Y1)
    n <- nrow(Y2)
    p <- ncol(Y1)
    assertthat::are_equal(ncol(Y1), ncol(Y2))
    assertthat::are_equal(length(sig_diag), ncol(Y2))
    assertthat::are_equal(ncol(alpha), ncol(Y2))

    if(p > 200) {
        stop("this won't work for large p")
    }

    Sigma1 <- diag(sig_diag * lambda)

    Sigma2 <- diag(sig_diag)

    del1 <- solve(Sigma1 + t(alpha) %*% alpha) %*% t(alpha)
    capdel1 <- diag(r) - alpha %*% del1
    ezz1 <- t(del1) %*% t(Y1) %*% Y1 %*% del1 + capdel1 * m

    del2 <- solve(Sigma2 + t(alpha) %*% alpha) %*% t(alpha)
    capdel2 <- diag(r) - alpha %*% del2
    ezz2 <- t(del2) %*% t(Y2) %*% Y2 %*% del2 + capdel2 * n

    A <- t(Y1) %*% Y1 %*% del1 / lambda + t(Y2) %*% Y2 %*% del2
    B <- ezz1 / lambda + ezz2

    alpha_new <- solve(B) %*% t(A)

    t1 <- colSums(Y1 ^ 2) / lambda + colSums(Y2 ^ 2)
    t2 <- colSums(t(A) * alpha)
    t3 <- rowSums(crossprod(alpha, B) * t(alpha))
    sig_diag_new <- (t1 - 2 * t2 + t3) / (n + m)
    return(list(alpha = alpha_new, sig_diag = sig_diag_new))
}

#' Fast update for sig_diag and alpha.
#'
#' @inheritParams qmle_obj
#'
#' @author David Gerard
update_sig_alpha <- function(alpha, sig_diag, lambda, Y1, Y2) {
    r <- nrow(alpha)
    m <- nrow(Y1)
    n <- nrow(Y2)
    p <- ncol(Y1)
    assertthat::are_equal(ncol(Y1), ncol(Y2))
    assertthat::are_equal(length(sig_diag), ncol(Y2))
    assertthat::are_equal(ncol(alpha), ncol(Y2))

    sig_inv <- 1 / sig_diag
    asig <- sweep(alpha, MARGIN = 2, sig_inv, `*`)
    asiga <- tcrossprod(alpha, asig)
    M2 <- diag(r) + asiga
    eigenM2 <- eigen(M2, symmetric = TRUE)
    capdel2 <- tcrossprod(sweep(eigenM2$vectors, MARGIN = 2, 1 / eigenM2$values, `*`),
                          eigenM2$vectors)
    del2 <- crossprod(asig, capdel2)
    Ydel2 <- Y2 %*% del2
    ezz2 <- crossprod(Ydel2, Ydel2) + capdel2 * n

    M1 <- diag(r) + asiga / lambda
    eigenM1 <- eigen(M1, symmetric = TRUE)
    capdel1 <- tcrossprod(sweep(eigenM1$vectors, MARGIN = 2, 1 / eigenM1$values, `*`),
                          eigenM1$vectors)
    del1 <- crossprod(asig, capdel1) / lambda
    Ydel1 <- Y1 %*% del1
    ezz1 <- crossprod(Ydel1, Ydel1) + capdel1 * m

    A <- crossprod(Y1, Ydel1) / lambda + crossprod(Y2, Ydel2)
    B <- ezz1 / lambda + ezz2
    eigenB <- eigen(B, symmetric = TRUE)
    Binv <- tcrossprod(sweep(eigenB$vectors, MARGIN = 2, 1 / eigenB$values, `*`),
                       eigenB$vectors)
    alpha_new <- tcrossprod(Binv, A)

    t1 <- colSums(Y1 ^ 2) / lambda + colSums(Y2 ^ 2)
    t2 <- colSums(t(A) * alpha)
    t3 <- rowSums(crossprod(alpha, B) * t(alpha))

    ## implementation checks
    ## t1p <- diag(t(Y1) %*% Y1 / lambda + t(Y2) %*% Y2)
    ## plot(t1p, t1); abline(0, 1)
    ## t2p <- diag(A %*% alpha)
    ## plot(t2p, t2); abline(0, 1)
    ## t3p <- diag(t(alpha) %*% B %*% alpha)
    ## plot(t3, t3p); abline(0, 1)

    sig_diag_new <- (t1 - 2 * t2 + t3) / (n + m)

    return(list(alpha = alpha_new, sig_diag = sig_diag_new))
}

#' Objective function for quasi-mle approach.
#'
#' This isn't really the log-likelihood, but the log-likelihood plus
#' some constant.
#'
#' For separate parts of hte likelihood, this code borrows heavily
#' from the function \code{\link[cate]{fa.em}}.
#'
#' @param alpha The current factor loadings
#' @param sig_diag A vector of the current variances
#' @param lambda The current inflation parameter for Y1
#' @param Y1 The first part of the observations.
#' @param Y2 The second part of the observations.
qmle_obj <- function(alpha, sig_diag, lambda, Y1, Y2) {
    r <- nrow(alpha)
    m <- nrow(Y1)
    n <- nrow(Y2)
    p <- ncol(Y1)
    assertthat::are_equal(ncol(Y1), ncol(Y2))
    assertthat::are_equal(length(sig_diag), ncol(Y2))
    assertthat::are_equal(ncol(alpha), ncol(Y2))

    ## Fast way, courtasy of CATE ------------------------------------
    sample_var2 <- colMeans(Y2 ^ 2)
    sig_inv <- 1 / sig_diag
    asig <- sweep(alpha, MARGIN = 2, sig_inv, `*`)
    asiga <- tcrossprod(alpha, asig)
    M2 <- diag(r) + asiga
    eigenM2 <- eigen(M2, symmetric = TRUE)
    Y2SG <- tcrossprod(Y2, asig)
    B2 <- (1 / sqrt(eigenM2$values)) * tcrossprod(t(eigenM2$vectors), Y2SG)
    tr2 <- -sum(B2 ^ 2) / n + sum(sample_var2 * sig_inv)
    logdet2 <- sum(log(eigenM2$values)) + sum(log(sig_diag))
    llike2 <- -tr2 - logdet2

    sample_var1 <- colMeans(Y1 ^ 2)
    M1 <- diag(r) + asiga / lambda
    eigenM1 <- eigen(M1, symmetric = TRUE)
    Y1SG <- tcrossprod(Y1, asig) / lambda
    B1 <- (1 / sqrt(eigenM1$values)) * tcrossprod(t(eigenM1$vectors), Y1SG)
    tr1 <- -sum(B1 ^ 2) / m + sum(sample_var1 * sig_inv / lambda)
    logdet1 <- sum(log(eigenM1$values)) + sum(log(sig_diag * lambda))
    llike1 <- -tr1 - logdet1

    llike <- (llike1 * m + llike2 * n) / (n + m)

    return(llike)
}


#' Basic no optimized code objective function.
#'
#' @inheritParams qmle_obj
#'
#' @author David Gerard
qmle_obj_basic <- function(alpha, sig_diag, lambda, Y1, Y2) {
    r <- nrow(alpha)
    m <- nrow(Y1)
    n <- nrow(Y2)
    p <- ncol(Y1)
    assertthat::are_equal(ncol(Y1), ncol(Y2))
    assertthat::are_equal(length(sig_diag), ncol(Y2))
    assertthat::are_equal(ncol(alpha), ncol(Y2))

    if(p > 200) {
        stop("this won't work for large p")
    }

    Sigma_inv <- diag(1 / sig_diag)
    ## cov_inv1b <- Sigma_inv / lambda - Sigma_inv %*% t(alpha) %*%
    ##     solve(diag(r) + alpha %*% Sigma_inv %*% t(alpha) / lambda) %*%
    ##     alpha %*% Sigma_inv / lambda ^ 2
    cov_inv1 <- solve(t(alpha) %*% alpha + diag(sig_diag) * lambda)


    ## cov_inv2b <- Sigma_inv - Sigma_inv %*% t(alpha) %*%
    ##     solve(diag(r) + alpha %*% Sigma_inv %*% t(alpha)) %*%
    ##     alpha %*% Sigma_inv
    cov_inv2 <-solve(t(alpha) %*% alpha + diag(sig_diag))

    siginv1_evals <- eigen(cov_inv1, symmetric = TRUE, only.values = TRUE)$values

    siginv2_evals <- eigen(cov_inv2, symmetric = TRUE, only.values = TRUE)$values


    llike1 <- sum(log(siginv1_evals)) * m - sum(diag(Y1 %*% cov_inv1 %*% t(Y1)))
    llike2 <- sum(log(siginv2_evals)) * n - sum(diag(Y2 %*% cov_inv2 %*% t(Y2)))

    llike <- (llike1 + llike2) / (n + m)

    return(llike)
}


trim <- function(x, tol = 10 ^ -13) {
    x[abs(x) < tol] <- 0
    return(x)
}



#' PCA when first \code{vr} rows have a variance multiplicatively
#' different from the rest of the rows.
#'
#' Modified truncated SVD. The variance estimates are just the
#' column-wise mean-squares of the last \eqn{n - vr} rows. This form
#' of factor analysis is mostly for variance inflation with RUV2.
#'
#' This doesn't work too well. I think the Z's and sigmas need to be
#' estimated jointly.
#'
#'
#' @param Y A matrix of numerics. The data.
#' @param r the rank.
#' @param vr The number of the first few rows whose variances differ
#'     by a multiplicative factor.
#' @param limmashrink A logical. Should we shrink the variance
#'     estimates prior to performing gls to get Z1 (\code{TRUE}) or
#'     not (\code{FALSE})?
#' @param likelihood What should the likelihood be when we estimate
#'     Z1? Options are \code{"t"} and \code{"normal"}.
#'
#' @export
#'
#' @author David Gerard
pca_2step <- function(Y, r, vr, limmashrink = TRUE, likelihood = c("t", "normal")) {
    assertthat::assert_that(is.matrix(Y))
    assertthat::are_equal(length(r), 1)
    assertthat::assert_that(r > 0 & r < min(dim(Y)))
    assertthat::assert_that(vr > 0 & vr < nrow(Y))

    likelihood <- match.arg(likelihood)

    p <- ncol(Y)
    n <- nrow(Y)

    Y1 <- Y[1:vr, , drop = FALSE]
    Y2 <- Y[(vr + 1):n, , drop = FALSE]

    svd_Y2 <- svd(Y2)
    alpha <- t(svd_Y2$v[, 1:r, drop = FALSE] %*% diag(svd_Y2$d[1:r], r, r)) /
        sqrt(n - vr)
    Z2 <- sqrt(n - vr) * svd_Y2$u[, 1:r, drop = FALSE]
    sig_diag <- colSums((Y2 - Z2 %*% alpha) ^ 2) / (n - vr - r)

    if (limmashrink) {
        limma_out <- limma::squeezeVar(var = sig_diag,
                                       df = n - vr - r)
        sig_diag <- limma_out$var.post
        prior_df <- limma_out$df.prior
        degrees_freedom <- prior_df + n - vr - r
    } else {
        degrees_freedom <- n - vr - r
    }

    ## estimate Z1 either using a normal or a t likelihood.
    if (likelihood == "normal") {
        Z1     <- Y1 %*% diag(1 / sig_diag) %*% t(alpha) %*%
            solve(alpha %*% diag(1 / sig_diag) %*% t(alpha))
        r1     <- colMeans((Y1 - Z1 %*% alpha) ^ 2)
        Z1 <- t(Z1)
        lambda <- mean(r1 / sig_diag)
    } else if (likelihood == "t") {
        tout   <- tregress_em(Y = t(Y1), alpha = t(alpha),
                              sig_diag = sig_diag,
                              nu = degrees_freedom)
        Z1     <- tout$Z
        lambda <- tout$lambda
    }

    Z <- rbind(t(Z1), Z2)

    return(list(alpha = alpha, Z = Z, sig_diag = sig_diag, lambda = lambda))
}
