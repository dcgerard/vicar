#' Bayesian version of Removing Unwanted Variation.
#'
#' Right now, this only allows you to use a uniform prior on the
#' beta's. I plan on the future to allow you to input a function
#' that's the prior.
#'
#' I have three versions of Bayesian factor analyses that I
#' recommend. The first, \code{\link{bfa_gs}} is the Bayesian factor
#' analysis used in Gerard and Stephens (2016). This is the default
#' version. The second is \code{\link{bfl}}. This version links the
#' variances between the factors and observations. Operationally,
#' there is very little practicle difference between these two but
#' \code{\link{bfa_gs}} would probably sit better with some
#' people. The last is \code{bfa_wrapper}, which is just a wrapper for
#' the R package bfa. The main thing about this version is that they
#' do not use a hierarchical prior on the variances.
#'
#' @inheritParams vruv4
#' @param fa_func A function that takes as input matrices named
#'     \code{Y21}, \code{Y31}, \code{Y32}, and \code{k} and returns a
#'     list, one of whose elements is called \code{Y22_array}. See
#'     \code{\link{bsvd}} for an example function.
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
#' @seealso \code{\link{bfa_gs}}, \code{\link{bfl}}, and
#'     \code{bfa_wrapper} for implemented Bayesian factor analyeses.
#'
#' @export
ruvb <- function(Y, X, ctl, k = NULL, fa_func = bfa_gs, fa_args = list(),
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

    return_list <- list()
    return_list$posterior_means   <- apply(betahat_post, c(1, 2), mean)
    return_list$posterior_medians <- apply(betahat_post, c(1, 2), stats::median)
    return_list$posterior_upper   <- apply(betahat_post, c(1, 2), stats::quantile, c(0.975))
    return_list$posterior_lower   <- apply(betahat_post, c(1, 2), stats::quantile, c(0.025))
    return_list$lfsr              <- apply(betahat_post, c(1, 2), calc_lfsr)

    lfsr_order <- order(c(return_list$lfsr))
    svalues <- matrix((cumsum(c(return_list$lfsr)[lfsr_order]) /
                       (1:(prod(dim(return_list$lfsr)))))[order(lfsr_order)],
                      nrow = nrow(return_list$lfsr), ncol = ncol(return_list$lfsr))
    return_list$svalues <- svalues

    return(return_list)
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
    bfout <- bfa::bfa_gauss(form1, data = Y, num.factor = k,
                            keep.scores = TRUE, thin = keep,
                            nburn = burnin, nsim = nsamp,
                            center.data = FALSE, scale.data = FALSE,
                            factor.scales = TRUE, loading.prior = "normal",
                            print.status = print_status)

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
