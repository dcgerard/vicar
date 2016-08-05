
#' EM algorithm for factor analysis with missing block matrix.
#'
#' @author David Gerard
#'
#' @param Y21 Top left of matrix.
#' @param Y31 Bottom left of matrix.
#' @param Y32 Top right of matrix.
#' @param k The rank of the mean matrix.
#' @param gls A logical. Should we estimate Z by generalized least
#'     squares (\code{TRUE}) or by a multivariate normality assumption
#'     (\code{FALSE})?
#' @param init_type A string. Should we initialize using PCA on
#'     \code{cbind(Y31, Y32)} (\code{"pca"}) or using maximum
#'     likelihood on \code{cbind(Y31, Y32)} (\code{"ml"})?
#' @return The top right of the matrix.
#'
#' @export
em_miss <- function(Y21, Y31, Y32, k, gls = TRUE, init_type = c("ml", "pca")) {

    p <- ncol(Y31) + ncol(Y32)
    
    init_type <- match.arg(init_type)

    if (is.null(Y21)) {
        cout <- cate::fa.em(Y = cbind(Y31, Y32), r = k)
        alpha_final <- t(cout$Gamma)
        sig_diag_final <- cout$Sigma
        return(list(alpha = alpha_final, sig_diag = sig_diag_final))
    }

    assertthat::are_equal(ncol(Y21), ncol(Y31))
    assertthat::are_equal(nrow(Y31), nrow(Y32))

    ## Get initial values -------------------------------------------------
    if (init_type == "pca") {
        pcout <- pca_naive(cbind(Y31, Y32), r = k)
        alpha_init <- t(pcout$alpha)
        sig_diag_init <- pcout$sig_diag
    } else if (init_type == "ml") {
        mlout <- cate::fa.em(Y = cbind(Y31, Y32), r = k)
        alpha_init <- t(mlout$Gamma)
        sig_diag_init <- mlout$Sigma
    }
    alpha_sigma_init <- c(c(alpha_init), sig_diag_init)

    ## Run SQUAREM --------------------------------------------------------
    sqout <- SQUAREM::squarem(par = alpha_sigma_init, fixptfn = em_miss_fix_fast,
                              objfn = em_miss_obj_fast,
                              Y21 = Y21, Y31 = Y31,
                              Y32 = Y32, k = k, control = list(tol = 10 ^ -6))

    ## em_miss_obj_fast(alpha_sigma = sqout$par, Y21 = Y21, Y31 = Y31, Y32 = Y32, k = k)

    sig_diag_final <- sqout$par[((k * p) + 1):length(sqout$par)]
    alpha_final <- matrix(sqout$par[1:(k * p)], nrow = k, ncol = p)

    ## Estimate Z --------------------------------------------------------
    alpha_c <- alpha_final[, 1:ncol(Y21)]
    alpha_nc <- alpha_final[, (ncol(Y21) + 1):ncol(alpha_final)]
    sig_diag_c <- sig_diag_final[1:ncol(Y21)]
    if (gls) {
        ## Zhat <- Y21 %*% diag(1 / sig_diag_c) %*% t(alpha_c) %*%
        ##     solve(alpha_c %*% diag(1 / sig_diag_c) %*% t(alpha_c))

        Zhat <- tcrossprod(sweep(Y21, 2, 1 / sig_diag_c, `*`), alpha_c) %*%
            solve(tcrossprod(sweep(alpha_c, 2, 1 / sig_diag_c, `*`), alpha_c))
    } else {
        Zhat <- tcrossprod(sweep(Y21, 2, 1 / sig_diag_c, `*`), alpha_c) %*%
            solve(tcrossprod(sweep(alpha_c, 2, 1 / sig_diag_c, `*`), alpha_c) + diag(k))
    }

    Y22hat <- Zhat %*% alpha_nc

    return(list(alpha = alpha_final, sig_diag = sig_diag_final, Z = Zhat, Y22hat = Y22hat))
}


#' Faster version of \code{\link{em_miss_fix}}.
#'
#'
#' @author David Gerard
#'
#' @param alpha_sigma A vector. The first k * p elements are the
#'     vectorization of alpha. The last elements are sig_diag.
#' @param k The rank of the mean.
#' @param Y21 The top left matrix.
#' @param Y31 The bottom left matrix.
#' @param Y32 The bottom right matrix.
em_miss_fix_fast <- function(alpha_sigma, Y21, Y31, Y32, k) {
    assertthat::assert_that(is.matrix(Y21))
    assertthat::assert_that(is.matrix(Y31))
    assertthat::assert_that(is.matrix(Y32))
    assertthat::are_equal(ncol(Y21), ncol(Y31))
    assertthat::are_equal(nrow(Y31), nrow(Y32))

    p <- ncol(Y31) + ncol(Y32)
    n <- nrow(Y21) + nrow(Y31)
    ncovs <- nrow(Y21)
    ncontrols <- ncol(Y21)

    assertthat::are_equal(length(alpha_sigma), p + k * p)

    alphavec <- alpha_sigma[1:(k * p)]
    alphamat <- matrix(alphavec, nrow = k, ncol = p)
    alphamat_c <- alphamat[, 1:ncontrols, drop = FALSE]
    alphamat_nc <- alphamat[, (ncontrols + 1):p, drop = FALSE]
    sig_diag <- alpha_sigma[((k * p) + 1):length(alpha_sigma)]
    sig_diag_c <- sig_diag[1:ncontrols]
    sig_diag_nc <- sig_diag[(ncontrols + 1):p]

    Y3 <- cbind(Y31, Y32)

    alpha_times_sigma <- sweep(alphamat, 2, 1 / sig_diag, `*`)
    covpreinv <-  tcrossprod(alpha_times_sigma, alphamat)+ diag(k)
    eigen_covpreinv <- eigen(covpreinv, symmetric = TRUE)
    delta_small <- eigen_covpreinv$vectors %*%
        sweep(crossprod(eigen_covpreinv$vectors, alpha_times_sigma), 1,
              1 / eigen_covpreinv$values, `*`)


    alpha_times_sigma_c <- sweep(alphamat_c, 2, 1 / sig_diag_c, `*`)
    covpreinv_c <-  tcrossprod(alpha_times_sigma_c, alphamat_c) + diag(k)
    eigen_covpreinv_c <- eigen(covpreinv_c, symmetric = TRUE)
    delta_small_c <- eigen_covpreinv_c$vectors %*%
        sweep(crossprod(eigen_covpreinv_c$vectors, alpha_times_sigma_c), 1,
              1 / eigen_covpreinv_c$values, `*`)


    Delta <- tcrossprod(delta_small, alphamat)
    Delta_c <- tcrossprod(delta_small_c, alphamat_c)


    S_diag <- colSums(Y3 ^ 2)
    S2_diag <- colSums(Y21 ^ 2)

    dY3 <- tcrossprod(delta_small, Y3)
    dSd <- tcrossprod(dY3)
    dY2 <- tcrossprod(delta_small_c, Y21)
    dS2d_c <- tcrossprod(tcrossprod(delta_small_c, Y21))

    F <- (n - ncovs) * (diag(k) - Delta) + dSd
    G <- ncovs * (diag(k) - Delta_c) + dS2d_c


    C <- crossprod(Y3, t(dY3))
    S2d <- crossprod(Y21, t(dY2))

    D <- rbind(S2d, crossprod(alphamat_nc, ncovs * (diag(k) - Delta_c) + dS2d_c))

    CD <- C + D

    FG_inv <- solve(F + G)

    alpha_new <- tcrossprod(FG_inv, CD)

    third_tr <- colSums(t(CD) * alpha_new)

    second_tr1 <- ncovs * sig_diag_nc
    second_tr2 <- rowSums(crossprod(alphamat_nc, dS2d_c + ncovs * (diag(k) - Delta_c)) *
                          t(alphamat_nc))
    second_tr <- c(S2_diag, second_tr1 + second_tr2)


    sig_diag_new <- (S_diag + second_tr - third_tr) / n

    alpha_sigma_new <- c(c(alpha_new), sig_diag_new)
    return(alpha_sigma_new)
}

#' A faster version of em_miss_obj.
#'
#' @author David Gerard
#'
#' @inheritParams em_miss_fix_fast
#'
#' @seealso \code{\link{em_miss_fix_fast}} for a fixed point iteration
#'     that will maximize this function. \code{\link{em_miss_obj}} for
#'     the slower version of this function.
em_miss_obj_fast <- function(alpha_sigma, Y21, Y31, Y32, k) {
    assertthat::assert_that(is.matrix(Y21))
    assertthat::assert_that(is.matrix(Y31))
    assertthat::assert_that(is.matrix(Y32))
    assertthat::are_equal(ncol(Y21), ncol(Y31))
    assertthat::are_equal(nrow(Y31), nrow(Y32))

    p <- ncol(Y31) + ncol(Y32)
    n <- nrow(Y21) + nrow(Y31)
    ncovs <- nrow(Y21)
    ncontrols <- ncol(Y21)

    assertthat::are_equal(length(alpha_sigma), p + k * p)

    alphavec <- alpha_sigma[1:(k * p)]
    alphamat <- matrix(alphavec, nrow = k, ncol = p)
    alphamat_c <- alphamat[, 1:ncontrols, drop = FALSE]
    alphamat_nc <- alphamat[, (ncontrols + 1):p, drop = FALSE]
    sig_diag <- alpha_sigma[((k * p) + 1):length(alpha_sigma)]
    sig_diag_c <- sig_diag[1:ncontrols]
    sig_diag_nc <- sig_diag[(ncontrols + 1):p]

    ## this is mostly so that SQUAREM doesn't overshoot variances
    if (any(sig_diag <= 0)) {
        return(-Inf)
    }
    
    Y3 <- cbind(Y31, Y32)


    Sig_inv <- 1 / sig_diag
    asiginv <- sweep(alphamat, MARGIN = 2, Sig_inv, `*`)
    asiginv_half <- sweep(alphamat, MARGIN = 2, sqrt(Sig_inv), `*`)
    asiga <- tcrossprod(alphamat, asiginv)
    ecent <- eigen(diag(k) + asiga)
    right_mat <- sweep(crossprod(ecent$vectors, asiginv_half), 1, 1 / sqrt(ecent$values), `*`)
    svout <- svd(right_mat)
    tr_part1 <- sum(sweep(Y3, 2, sqrt(Sig_inv), `*`) ^ 2) -
        sum(sweep(sweep(Y3, 2, sqrt(Sig_inv), `*`) %*% svout$v, 2, svout$d, `*`) ^ 2)

    det_part1 <- sum(log((1 - svout$d ^ 2))) + sum(log(Sig_inv))


    ## tr_part1 <- sum(sweep(Y3, 2, sqrt(Sig_inv), `*`) ^ 2) -
    ##     sum(asiginvY3 * (solve(diag(k) + asiga) %*% asiginvY3))
    ## sum(diag(Y3 %*% solve(t(alphamat) %*% alphamat + diag(sig_diag)) %*% t(Y3)))

    ## invtemp <- diag(Sig_inv) - diag(Sig_inv) %*% t(alphamat) %*%
    ##     solve(diag(k) + alphamat %*% diag(Sig_inv) %*% t(alphamat)) %*%
    ##     alphamat %*% diag(Sig_inv)
    ## invtrue <- solve(t(alphamat) %*% alphamat + diag(sig_diag))

    Sig_inv_c <- 1 / sig_diag_c
    asiginv_c <- sweep(alphamat_c, MARGIN = 2, Sig_inv_c, `*`)
    asiginv_half_c <- sweep(alphamat_c, MARGIN = 2, sqrt(Sig_inv_c), `*`)
    asiga_c <- tcrossprod(alphamat_c, asiginv_c)
    ecent_c <- eigen(diag(k) + asiga_c)
    right_mat_c <- sweep(crossprod(ecent_c$vectors, asiginv_half_c), 1,
                         1 / sqrt(ecent_c$values), `*`)
    svout_c <- svd(right_mat_c)
    tr_part2 <- sum(sweep(Y21, 2, sqrt(Sig_inv_c), `*`) ^ 2) -
        sum(sweep(sweep(Y21, 2, sqrt(Sig_inv_c), `*`) %*% svout_c$v, 2, svout_c$d, `*`) ^ 2)

    det_part2 <- sum(log((1 - svout_c$d ^ 2))) + sum(log(Sig_inv_c))


    llike <- (n - ncovs) * det_part1 - tr_part1 + ncovs * det_part2 - tr_part2
    return(llike)
}

#' Fixed point iteration for em algorithm with missing block.
#'
#' @author David Gerard
#'
#' @param alpha_sigma A vector. The first k * p elements are the
#'     vectorization of alpha. The last elements are sig_diag.
#' @param S A matrix. This is the sample covariance matrix of the
#'     samples that don't have a missing block.
#' @param S2 A matrix. This is the sample covariance matrix of the
#'     samples that do have a missing block.
#' @param ncovs The number of covariates.
#' @param ncontrols The number of controls
#' @param n The sample size
#' @param p The number of genes.
#' @param k The rank of the mean.
#'
#'
em_miss_fix <- function(alpha_sigma, S, S2, ncovs, ncontrols, n, p, k) {

    if (is.null(S2) & ncovs != 0) {
        stop("If S2 is NULL then ncovs must equal 0")
    } else if (ncovs == 0 & !is.null(S2)) {
        stop("If ncovs is 0 then S2 must be NULL")
    }


    assertthat::are_equal(length(alpha_sigma), k * p + p)
    alphavec <- alpha_sigma[1:(k * p)]
    alphamat <- matrix(alphavec, nrow = k, ncol = p)
    alphamat_c <- alphamat[, 1:ncontrols, drop = FALSE]
    alphamat_nc <- alphamat[, (ncontrols + 1):p, drop = FALSE]
    sig_diag <- alpha_sigma[((k * p) + 1):length(alpha_sigma)]
    sig_diag_c <- sig_diag[1:ncontrols]
    sig_diag_nc <- sig_diag[(ncontrols + 1):p]

    alpha_times_sigma <- sweep(alphamat, 2, 1 / sig_diag, `*`)
    covpreinv <-  tcrossprod(alpha_times_sigma, alphamat)+ diag(k)
    eigen_covpreinv <- eigen(covpreinv, symmetric = TRUE)
    delta_small <- eigen_covpreinv$vectors %*%
        sweep(crossprod(eigen_covpreinv$vectors, alpha_times_sigma), 1,
              1 / eigen_covpreinv$values, `*`)

    ## delta_small should equal temp
    ## temp <- solve(alphamat %*% diag(1 / sig_diag) %*% t(alphamat) + diag(k)) %*%
    ##     alphamat %*% diag(1 / sig_diag)
    ## summary(c(temp - delta_small))

    alpha_times_sigma_c <- sweep(alphamat_c, 2, 1 / sig_diag_c, `*`)
    covpreinv_c <-  tcrossprod(alpha_times_sigma_c, alphamat_c) + diag(k)
    eigen_covpreinv_c <- eigen(covpreinv_c, symmetric = TRUE)
    delta_small_c <- eigen_covpreinv_c$vectors %*%
        sweep(crossprod(eigen_covpreinv_c$vectors, alpha_times_sigma_c), 1,
              1 / eigen_covpreinv_c$values, `*`)

    ## delta_small_c should equal temp_c
    ## temp_c <- solve(alphamat_c %*% diag(1 / sig_diag_c) %*%
    ##                        t(alphamat_c) + diag(k)) %*%
    ##     alphamat_c %*% diag(1 / sig_diag_c)
    ## summary(c(temp_c - delta_small_c))


    Delta <- tcrossprod(delta_small, alphamat)
    Delta_c <- tcrossprod(delta_small_c, alphamat_c)

    C <- tcrossprod(S, delta_small)

    F <- (n - ncovs) * diag(k) - (n - ncovs) * Delta + delta_small %*% C

    if (ncovs != 0) {
        S2_delta_c <- tcrossprod(S2, delta_small_c)
        dS2d_c <- delta_small_c %*% S2_delta_c
        Dbottom <- crossprod(alphamat_nc, ncovs * (diag(k) - Delta_c) + dS2d_c)
        D <- rbind(S2_delta_c, Dbottom)

        G <- ncovs * (diag(k) - Delta_c) + dS2d_c


        eFG <- eigen(F + G, symmetric = TRUE)
        CD <- C + D
        alpha_new <- tcrossprod(sweep(eFG$vectors, 2, 1 / eFG$values, `*`), CD %*% eFG$vectors)

        ## alpha_new should equal temp below
        ## temp <- solve(F + G) %*% t(C + D)
        ## summary(c(temp - alpha_new))
    } else {
        alpha_new <- solve(F) %*% t(C)
    }


    first_tr <- diag(S)
    if (ncovs != 0) {
        third_tr <- rowSums(CD * t(alpha_new))
        ## diag(CD %*% alpha_new)
        second_tr1 <- ncovs * sig_diag_nc
        second_tr2 <- rowSums(crossprod(alphamat_nc, dS2d_c + ncovs * (diag(k) - Delta_c)) *
                              t(alphamat_nc))
        second_tr3 <- diag(S2)
        second_tr <- c(second_tr3, second_tr1 + second_tr2)

        sig_diag_new <- (first_tr + second_tr - third_tr) / n
    } else {
        third_tr <- rowSums(C * t(alpha_new))
        sig_diag_new <- (first_tr - third_tr) / n
    }

    ## see if update of alpha is correct
    ## sig_diag_new <- sig_diag


    ## sig_diag_new should equal temp below
    ## temp2 <- ncovs * sig_diag_nc +
    ## diag(t(alphamat_nc) %*% (ncovs * diag(k) + ncovs * Delta_c + delta_small_c %*% S2 %*%
    ##                          t(delta_small_c)) %*% alphamat_nc)
    ## temp <- (diag(S) + c(temp2, diag(S2)) - diag((C + D) %*% solve(F + G) %*% t(C + D))) /
    ##     (n - ncovs)
    ## summary(temp - sig_diag_new)

    alpha_sigma_new <- c(c(alpha_new), sig_diag_new)
    return(alpha_sigma_new)
}

#' The objective function for em with a missing block.
#'
#' @inheritParams em_miss_fix
#'
#' @author David Gerard
em_miss_obj <- function(alpha_sigma, S, S2, ncovs, ncontrols, n, p, k) {

    if (is.null(S2) & ncovs != 0) {
        stop("If S2 is NULL then ncovs must equal 0")
    } else if (ncovs == 0 & !is.null(S2)) {
        stop("If ncovs is 0 then S2 must be NULL")
    }

    assertthat::are_equal(length(alpha_sigma), k * p + p)
    alphavec <- alpha_sigma[1:(k * p)]
    alphamat <- matrix(alphavec, nrow = k, ncol = p)
    alphamat_c <- alphamat[, 1:ncontrols, drop = FALSE]
    alphamat_nc <- alphamat[, (ncontrols + 1):p, drop = FALSE]
    sig_diag <- alpha_sigma[((k * p) + 1):length(alpha_sigma)]
    sig_diag_c <- sig_diag[1:ncontrols]
    sig_diag_nc <- sig_diag[(ncontrols + 1):p]


    alpha_times_sigma <- sweep(alphamat, 2, 1 / sig_diag, `*`)
    aSa <- tcrossprod(alpha_times_sigma, alphamat)
    eprecov <- eigen(diag(k) + aSa, symmetric = TRUE)
    right_prod <- sweep(crossprod(eprecov$vectors, alpha_times_sigma), 1,
                        1 / sqrt(eprecov$values), FUN = `*`)
    tot_prod <- crossprod(right_prod)
    aasig_inv <- diag(1 / sig_diag) - tot_prod

    ## aasig_inv should be the same as temp below
    ## temp <- solve(t(alphamat) %*% alphamat + diag(sig_diag))
    ## plot(temp, aasig_inv)
    ## abline(0, 1)

    if (ncovs != 0) {
        alpha_times_sigma_c <- sweep(alphamat_c, 2, 1 / sig_diag_c, `*`)
        aSa_c <- tcrossprod(alpha_times_sigma_c, alphamat_c)
        eprecov_c <- eigen(diag(k) + aSa_c, symmetric = TRUE)
        right_prod_c <- sweep(crossprod(eprecov_c$vectors, alpha_times_sigma_c), 1,
                              1 / sqrt(eprecov_c$values), FUN = `*`)
        tot_prod_c <- crossprod(right_prod_c)
        aasig_inv_c <- diag(1 / sig_diag_c) - tot_prod_c
    }
    ## aasig_inv_c should be the same as temp_c below
    ## temp_c <- solve(t(alphamat_c) %*% alphamat_c + diag(sig_diag_c))
    ## plot(temp_c, aasig_inv_c)
    ## abline(0, 1)

    exp_part1 <- sum(aasig_inv * S)
    det_part1 <- sum(log(eigen(aasig_inv, symmetric = TRUE, only.values = TRUE)$values))

    if (ncovs != 0) {
        exp_part2 <- sum(aasig_inv_c * S2)
        det_part2 <- sum(log(eigen(aasig_inv_c, symmetric = TRUE, only.values = TRUE)$values))
        llike <- (n - ncovs) * det_part1 - exp_part1 + ncovs * det_part2 - exp_part2
    } else {
        llike <- n * det_part1 - exp_part1
    }
    return(llike)
}


#' Very inneficient copy of rubin and thayer iteration mostly meant for debugging.
#'
#' @author David Gerard
#'
#' @inheritParams em_miss_fix
rubin_copy <- function(alpha_sigma, S, n, p, k) {
    assertthat::are_equal(length(alpha_sigma), k * p + p)
    alphavec <- alpha_sigma[1:(k * p)]
    alphamat <- matrix(alphavec, nrow = k, ncol = p)
    sig_diag <- alpha_sigma[((k * p) + 1):length(alpha_sigma)]

    delta_small <- solve(diag(sig_diag) + t(alphamat) %*% alphamat) %*% t(alphamat)
    Delta <- diag(k) - alphamat %*% delta_small

    Cyy <- S / n

    alphamat_new <- solve(t(delta_small) %*% Cyy %*% delta_small + Delta) %*%
        t(delta_small) %*% t(Cyy)
    sig_diag_new <- diag(Cyy - Cyy %*% delta_small %*%
                         solve((t(delta_small) %*% Cyy %*% delta_small + Delta)) %*%
                         t(delta_small) %*% t(Cyy))

    alpha_sigma_new <- c(c(alphamat_new), sig_diag_new)
    return(alpha_sigma_new)
}
