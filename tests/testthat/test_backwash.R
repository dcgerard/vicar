context("Backwash")

test_that("backwash works", {

  set.seed(887)
    n <- 11
    p <- 73
    q <- 3
    k <- 5

    X <- matrix(rnorm(n * q), nrow = n)
    beta <- matrix(rnorm(q * p), nrow = q)
    Y <- X %*% beta + matrix(rnorm(n * p), nrow = n)
    cov_of_interest = ncol(X)
    include_intercept = TRUE
    limmashrink = TRUE
    fa_func = vicar::pca_naive
    fa_args = list()
    lambda_type = "zero_conc"
    pi_init_type = "zero_conc"
    grid_seq = NULL
    lambda_seq = NULL
    lambda0 = 10
    scale_var = TRUE
    sprop = 0

}
)

test_that("BACKWASH ELBO increases for each update", {

    set.seed(16)
    p <- 101
    k <- 11

    beta <- stats::rnorm(p, sd = 10)
    beta[1:77] <- 0
    alpha_tilde <- matrix(stats::rnorm(p * k), nrow = p)
    z <- stats::rnorm(k)
    S_diag <- stats::rchisq(p, df = 1)
    err <- stats::rnorm(p, sd = sqrt(S_diag))

    eigen_alpha <- eigen(crossprod(alpha_tilde, alpha_tilde), symmetric = TRUE)
    a2_half_inv <- eigen_alpha$vectors %*% diag(1 / sqrt(eigen_alpha$values)) %*% t(eigen_alpha$vectors)
    Amat <- alpha_tilde %*% a2_half_inv

    ## m1 <- Amat %*% t(Amat)
    ## m2 <- alpha_tilde %*% solve(t(alpha_tilde) %*% alpha_tilde) %*% t(alpha_tilde)
    ## all(abs(m1 - m2) < 10 ^ -14)

    betahat_ols <- beta + alpha_tilde %*% a2_half_inv  %*% z + err

    M <- 13

    tau2_seq <- seq(0, 3, length = M)
    pivec <- stats::runif(M)
    pivec <- pivec / sum(pivec)
    lambda_seq <- c(10, rep(1, length = M - 1))

    mubeta <- rnorm(p)
    muv <- matrix(rnorm(k), ncol = 1)
    xi <- 1
    phi <- 1


    qbout <- back_update_qbeta(betahat_ols = betahat_ols,
                               S_diag = S_diag, Amat = Amat,
                               pivec = pivec, tau2_seq = tau2_seq,
                               muv = muv, xi = xi, phi = phi)
    mubeta <- qbout$mubeta
    mubeta_matrix <- qbout$mubeta_matrix
    sig2beta_matrix <- qbout$sig2beta_matrix
    gamma_mat <- qbout$gamma_mat

    ## test that pi is updated correctly ---------------------------------------
    pivec1 <- pivec
    pivec <- back_update_pi(gamma_mat = gamma_mat, lambda_seq = lambda_seq)
    pivec2 <- pivec

    a1 <- sum(colSums(gamma_mat) * log(pivec1)) + sum((lambda_seq - 1) * log(pivec1))
    a2 <- sum(colSums(gamma_mat) * log(pivec2)) + sum((lambda_seq - 1) * log(pivec2))
    expect_true(a1 <= a2)

    ## test that v is updated correctly --------------------------------------
    ## This sequence will decrease the elbo but it should increase it
    qvout <- back_update_v(betahat_ols = betahat_ols, S_diag = S_diag,
                           Amat = Amat, mubeta = mubeta, xi = xi,
                           phi = phi)
    muv <- qvout$muv
    Sigma_v <- qvout$Sigma_v

    phi <- back_update_phi(betahat_ols = betahat_ols, S_diag = S_diag,
                               Amat = Amat, mubeta = mubeta, muv = muv,
                               Sigma_v= Sigma_v)
    xi <- back_update_xi(betahat_ols = betahat_ols, S_diag = S_diag,
                             Amat = Amat, mubeta = mubeta,
                             mubeta_matrix = mubeta_matrix,
                             sig2beta_matrix = sig2beta_matrix,
                             gamma_mat = gamma_mat, muv = muv,
                             Sigma_v = Sigma_v, phi = phi)

    qvout <- back_update_v(betahat_ols = betahat_ols, S_diag = S_diag,
                           Amat = Amat, mubeta = mubeta, xi = xi,
                           phi = phi)
    muv2 <- qvout$muv
    Sigma_v2 <- qvout$Sigma_v

    ASA <- t(Amat) %*% diag(1 / S_diag) %*% Amat
    b11 <- - sum(diag(ASA %*% (muv %*% t(muv) + Sigma_v))) * (phi ^ 2) / (2 * xi)
    b12 <- 2 * phi * t(betahat_ols) %*% diag(1 / S_diag) %*% Amat %*% muv / (2 * xi)
    b13 <- - 2 * phi * t(mubeta) %*% diag(1 / S_diag) %*% Amat %*% muv / (2 * xi)
    b14 <- - t(muv) %*% muv / 2
    b15 <-  determinant(Sigma_v, logarithm = TRUE)$modulus / 2
    b16 <- - sum(diag(Sigma_v)) / 2
    b1 <- b11 + b12 + b13 + b14 + b15 + b16

    b21 <- - sum(diag(ASA %*% (muv2 %*% t(muv2) + Sigma_v2))) * (phi ^ 2) / (2 * xi)
    b22 <- 2 * phi * t(betahat_ols) %*% diag(1 / S_diag) %*% Amat %*% muv2 / (2 * xi)
    b23 <- - 2 * phi * t(mubeta) %*% diag(1 / S_diag) %*% Amat %*% muv2 / (2 * xi)
    b24 <- - t(muv2) %*% muv2 / 2
    b25 <-  determinant(Sigma_v2, log = TRUE)$modulus / 2
    b26 <- - sum(diag(Sigma_v2)) / 2
    b2 <- b21 + b22 + b23 + b24 + b25 + b26
    expect_true(b2 > b1)

    elbo1 <- back_elbo(betahat_ols = betahat_ols, S_diag = S_diag,
                      Amat = Amat, tau2_seq = tau2_seq,
                      pivec = pivec, lambda_seq = lambda_seq,
                      mubeta = mubeta,
                      mubeta_matrix = mubeta_matrix,
                      sig2beta_matrix = sig2beta_matrix,
                      gamma_mat = gamma_mat, muv = muv,
                      Sigma_v = Sigma_v, phi = phi, xi = xi)

    elbo2 <- back_elbo(betahat_ols = betahat_ols, S_diag = S_diag,
                      Amat = Amat, tau2_seq = tau2_seq,
                      pivec = pivec, lambda_seq = lambda_seq,
                      mubeta = mubeta,
                      mubeta_matrix = mubeta_matrix,
                      sig2beta_matrix = sig2beta_matrix,
                      gamma_mat = gamma_mat, muv = muv2,
                      Sigma_v = Sigma_v2, phi = phi, xi = xi)

    expect_true(elbo2 > elbo1)

    ## test that elbo always increases ---------------------------------------
    itermax <- 20
    elbo_mat <- matrix(NA, nrow = 5, ncol = itermax - 1)
    rownames(elbo_mat) <- c("beta", "pi", "v", "phi", "xi")
    for (iterindex in 1:itermax) {

        qbout <- back_update_qbeta(betahat_ols = betahat_ols,
                                   S_diag = S_diag, Amat = Amat,
                                   pivec = pivec, tau2_seq = tau2_seq,
                                   muv = muv, xi = xi, phi = phi)
        mubeta <- qbout$mubeta
        mubeta_matrix <- qbout$mubeta_matrix
        sig2beta_matrix <- qbout$sig2beta_matrix
        gamma_mat <- qbout$gamma_mat
        if (iterindex > 1) {
            elbo <- back_elbo(betahat_ols = betahat_ols, S_diag = S_diag,
                              Amat = Amat, tau2_seq = tau2_seq,
                              pivec = pivec, lambda_seq = lambda_seq,
                              mubeta = mubeta,
                              mubeta_matrix = mubeta_matrix,
                              sig2beta_matrix = sig2beta_matrix,
                              gamma_mat = gamma_mat, muv = muv,
                              Sigma_v = Sigma_v, phi = phi, xi = xi)
            elbo_mat[1, iterindex - 1] <- elbo
        }


        pivec <- back_update_pi(gamma_mat = gamma_mat, lambda_seq = lambda_seq)
        if (iterindex > 1) {
            elbo <- back_elbo(betahat_ols = betahat_ols, S_diag = S_diag,
                              Amat = Amat, tau2_seq = tau2_seq,
                              pivec = pivec, lambda_seq = lambda_seq,
                              mubeta = mubeta,
                              mubeta_matrix = mubeta_matrix,
                              sig2beta_matrix = sig2beta_matrix,
                              gamma_mat = gamma_mat, muv = muv,
                              Sigma_v = Sigma_v, phi = phi, xi = xi)
            elbo_mat[2, iterindex - 1] <- elbo
        }

        qvout <- back_update_v(betahat_ols = betahat_ols, S_diag = S_diag,
                               Amat = Amat, mubeta = mubeta, xi = xi,
                               phi = phi)
        muv <- qvout$muv
        Sigma_v <- qvout$Sigma_v
        if (iterindex > 1) {
            elbo <- back_elbo(betahat_ols = betahat_ols, S_diag = S_diag,
                              Amat = Amat, tau2_seq = tau2_seq,
                              pivec = pivec, lambda_seq = lambda_seq,
                              mubeta = mubeta,
                              mubeta_matrix = mubeta_matrix,
                              sig2beta_matrix = sig2beta_matrix,
                              gamma_mat = gamma_mat, muv = muv,
                              Sigma_v = Sigma_v, phi = phi, xi = xi)
            elbo_mat[3, iterindex - 1] <- elbo
        }

        phi <- back_update_phi(betahat_ols = betahat_ols, S_diag = S_diag,
                               Amat = Amat, mubeta = mubeta, muv = muv,
                               Sigma_v= Sigma_v)
        if (iterindex > 1) {
            elbo <- back_elbo(betahat_ols = betahat_ols, S_diag = S_diag,
                              Amat = Amat, tau2_seq = tau2_seq,
                              pivec = pivec, lambda_seq = lambda_seq,
                              mubeta = mubeta,
                              mubeta_matrix = mubeta_matrix,
                              sig2beta_matrix = sig2beta_matrix,
                              gamma_mat = gamma_mat, muv = muv,
                              Sigma_v = Sigma_v, phi = phi, xi = xi)
            elbo_mat[4, iterindex - 1] <- elbo
        }

        xi <- back_update_xi(betahat_ols = betahat_ols, S_diag = S_diag,
                             Amat = Amat, mubeta = mubeta,
                             mubeta_matrix = mubeta_matrix,
                             sig2beta_matrix = sig2beta_matrix,
                             gamma_mat = gamma_mat, muv = muv,
                             Sigma_v = Sigma_v, phi = phi)
        if (iterindex > 1) {
            elbo <- back_elbo(betahat_ols = betahat_ols, S_diag = S_diag,
                              Amat = Amat, tau2_seq = tau2_seq,
                              pivec = pivec, lambda_seq = lambda_seq,
                              mubeta = mubeta,
                              mubeta_matrix = mubeta_matrix,
                              sig2beta_matrix = sig2beta_matrix,
                              gamma_mat = gamma_mat, muv = muv,
                              Sigma_v = Sigma_v, phi = phi, xi = xi)
            elbo_mat[5, iterindex - 1] <- elbo
        }
    }

    did_increase <- matrix(c(NA, c(elbo_mat)[2:length(c(elbo_mat))] - c(elbo_mat)[1:(length(c(elbo_mat)) - 1)]) >= -10^-10, nrow = 5)
    rownames(did_increase) <- rownames(elbo_mat)
    # did_increase

    expect_true(all(did_increase[-1]))


    # mout <- mouthwash_second_step(betahat_ols = betahat_ols,
    #                               S_diag = S_diag,
    #                               alpha_tilde = alpha_tilde,
    #                               lambda_seq = lambda_seq,
    #                               tau2_seq = tau2_seq, a_seq = NULL,
    #                               b_seq = NULL,
    #                               mixing_dist = "normal",
    #                               likelihood = "normal",
    #                               pi_init_type = "zero_conc",
    #                               scale_var = TRUE,
    #                               degrees_freedom = Inf)
    #
    # mout$pi0
    # pivec[1]
    # plot(mout$result$lfdr)
    # abline(v = 77)
    # plot(gamma_mat[, 1])
    # abline(v = 77)
    # plot(mout$result$lfdr, gamma_mat[, 1], col = (1:p <= 77) + 2)
    # abline(0, 1)

  # muv <- matrix(rnorm(k), ncol = 1)
  # xi <- 1
  # phi <- 1

  ## see that qv works
  # qvout <- back_update_v(betahat_ols = betahat_ols, S_diag = S_diag, Amat = Amat, mubeta = mubeta,
  #                        xi = xi, phi = phi)
  #
  # musig <- c(qvout$muv, chol(qvout$Sigma_v))
  # qvoptim_wrapper <- function(musig, betahat_ols, S_diag, Amat, tau2_seq, pivec,
  #                             lambda_seq, mubeta, mubeta_matrix,
  #                             sig2beta_matrix, gamma_mat, phi, xi) {
  #   q <- ncol(Amat)
  #   muv <- musig[1:q]
  #   Sigma_half <- matrix(musig[-(1:q)], nrow = q)
  #   Sigma_v <- t(Sigma_half) %*% Sigma_half
  #
  #   elbo <- back_elbo(betahat_ols = betahat_ols, S_diag = S_diag, Amat = Amat, tau2_seq = tau2_seq, pivec = pivec,
  #                     lambda_seq = lambda_seq, mubeta = mubeta, mubeta_matrix = mubeta_matrix,
  #                     sig2beta_matrix = sig2beta_matrix, gamma_mat = gamma_mat, muv = muv,
  #                     Sigma_v = Sigma_v, phi = phi, xi = xi)
  #   return(elbo)
  # }
  # oout <- stats::optim(par = musig, fn = qvoptim_wrapper, betahat_ols = betahat_ols, S_diag = S_diag,
  #                      Amat = Amat, tau2_seq = tau2_seq, pivec = pivec,
  #                      lambda_seq = lambda_seq, mubeta = mubeta, mubeta_matrix = mubeta_matrix,
  #                      sig2beta_matrix = sig2beta_matrix, gamma_mat = gamma_mat, phi = phi, xi = xi,
  #                      control = list(maxit = 10000, fnscale = -1))
  #
  # qvoptim_wrapper(musig, betahat_ols = betahat_ols, S_diag = S_diag,
  #                 Amat = Amat, tau2_seq = tau2_seq, pivec = pivec,
  #                 lambda_seq = lambda_seq, mubeta = mubeta, mubeta_matrix = mubeta_matrix,
  #                 sig2beta_matrix = sig2beta_matrix, gamma_mat = gamma_mat, phi = phi, xi = xi)
  # qvoptim_wrapper(oout$par, betahat_ols = betahat_ols, S_diag = S_diag,
  #                 Amat = Amat, tau2_seq = tau2_seq, pivec = pivec,
  #                 lambda_seq = lambda_seq, mubeta = mubeta, mubeta_matrix = mubeta_matrix,
  #                 sig2beta_matrix = sig2beta_matrix, gamma_mat = gamma_mat, phi = phi, xi = xi)
  # nfac <- length(muv)
  # plot(oout$par[1:nfac], qvout$muv)
  # abline(0, 1)
  # Sigma_half <- matrix(oout$par[-(1:nfac)], ncol = nfac)
  # Sigma_v <- t(Sigma_half) %*% Sigma_half
  # plot(qvout$Sigma_v, Sigma_v)
  # abline(0, 1)
  #
  # summary(c(qvout$Sigma_v - Sigma_v))
  # summary(c(qvout$muv - oout$par[1:nfac]))

}
)

