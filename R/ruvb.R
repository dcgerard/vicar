#' Bayesian version of Removing Unwanted Variation.
#'
#' @inheritParams vruv4
#' @param fa_func A function that takes as input a matrix names
#'     \code{Y} that has missing values and returns a list of matrices
#'     called \code{Yhat} where each matrix is the same dimension of
#'     \code{Y} with the missing values filled in.
#' @param fa_args A list of additional parameters to pass to
#'     \code{fa_func}.
#'
#' @return \code{beta2hat} The estimates of the coefficients of the
#'     covariates of interest that do not correspond to control genes.
#'
#'     \code{betahat_long} The estimates of the coefficients. Those
#'     corresponding to control genes are set to 0.
#'
#'     \code{sebetahat} If \code{do_variance = TRUE}, then these are
#'     the "standard errors" of \code{beta2hat} (but not really).
#'
#'     \code{tstats} If \code{do_variance = TRUE}, then these are
#'     the "t-statistics" of \code{beta2hat} (but not really).
#'
#'     \code{pvalues} If \code{do_variance = TRUE}, then these are
#'     the "p-values" of \code{tstats} (but not really).
#'
#'
#' @author David Gerard
#'
#' @export
ruvb <- function(Y, X, ctl, k = NULL, fa_func, fa_args = list(),
                 cov_of_interest = ncol(X), include_intercept = TRUE) {

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

    Ytilde <- matrix(NA, nrow = nrow(Y21) + nrow(Y31), ncol = ncol(Y21) + ncol(Y22))
    Ytilde[1:ncovariates, 1:ncontrols] <- Y21
    Ytilde[(ncovariates + 1):nrow(Ytilde), 1:ncontrols] <- Y31
    Ytilde[(ncovariates + 1):nrow(Ytilde), (ncontrols + 1):ncol(Ytilde)] <- Y32

    R22inv <- backsolve(R22, diag(nrow(R22)))



    return()
}


#' Gibbs sampler for Bayesian SVD.
#'
#' This is a modification of the Bayesian approach from Hoff (2007) to
#' allow for heteroscedastic columns. We start the missing values from
#' the RUV4 solution.
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
#' @export
#'
#' @references Hoff, P. D. (2012). Model averaging and dimension
#'     selection for the singular value decomposition. Journal of the
#'     American Statistical Association.
bsvd <- function(Y21, Y31, Y32, k, nsamp = 10000,
                 burnin = round(nsamp / 4), keep = 20,
                 print_update = TRUE,
                 plot_update = FALSE) {

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

    Y22_array <- array(NA, dim = c(ncovs, p - ncontrols, round(nsamp / keep) + 1))
    Y22_array[, , 1] <- Y22_current
    keep_index <- 2
    mu_psi_phi <- matrix(NA, nrow = round(nsamp / keep) + 1, ncol = 3)
    mu_psi_phi[1, ] <- c(mu_current, psi_current, phi_current)
    plot_iters <- round(seq(0, 1, by = 0.01) * nsamp)
    if (print_update) {
        cat("Progress:\n")
        pb <- utils::txtProgressBar(style = 3)
    }
    for (gindex in 1:nsamp) {
        if (print_update) {
            utils::setTxtProgressBar(pb = pb, value = gindex / nsamp)
        }
        ## Update U --------------------------------------
        Cmat <- Yxi %*% sweep(v_current, 2, delta_current, `*`)
        u_current <- rstiefel::rmf.matrix.gibbs(M = Cmat, X = u_current)

        ## Update V --------------------------------------
        Cmat <- crossprod(Yxi, sweep(u_current, 2, delta_current, `*`))
        v_current <- rstiefel::rbmf.matrix.gibbs(A = diag(xi_current),
                                                 B = diag(-delta_current / 2),
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
        psi_rate    <- (k + eta_0) / 2
        psi_shape   <- (eta_0 * tau_0 + sum((delta_current - mu_current) ^ 2)) / 2
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

        if (gindex %% keep == 0) {
            Y22_array[, , keep_index] <- Y22_current
            mu_psi_phi[keep_index, 1] <- mu_current
            mu_psi_phi[keep_index, 2] <- psi_current
            mu_psi_phi[keep_index, 3] <- phi_current
            keep_index <- keep_index + 1
        }

        if (gindex %in% plot_iters) {
            if (plot_update) {
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
