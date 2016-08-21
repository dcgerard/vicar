#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
using namespace Rcpp;
using namespace arma;



//' Fast Gibbs sampler for Bayesian factor analysis.
//'
//'
//'
//' @param Y22init A matrix of numerics. The initial value of Y22.
//' @param Y21 A matrix of numerics.
//' @param Y31 A matrix of numerics.
//' @param Y32 A matrix of numerics.
//' @param nsamp The number of iterations to run in the Gibbs sampler, not including the burnin.
//' @param burnin The number of iterations to burnin.
//' @param thin We only collect samples every \code{thin} iterations.
//' @param Linit A numeric matrix. The initial values of the loadings.
//' @param Finit A numeric matrix. The initial values for the factors.
//' @param xi_init A numeric vector. The initial values of the precisions.
//' @param phi_init A numeric scalar. The initial value of the mean of the precitions.
//' @param zeta_init A numeric vector. The initial values of the augmented row precisions.
//' @param theta_init A numeric vector. The initial values the the factor precisions.
//' @param kappa_init A numeric scalar. The initial value of the mean of the factor precisions.
//' @param rho_0 The prior sample size for column-specific the
//'     precisions.
//' @param alpha_0 The prior sample size for the mean of the
//'     column-specific precisions.
//' @param beta_0 The prior mean of the mean of the column-specific
//'     precisions.
//' @param eta_0 The prior sample size of the expanded parameters.
//' @param tau_0 The prior mean of of the expanded parameters.
//' @param delta_0 The prior sample size of the column-specific
//'     precisions of the factors.
//' @param lambda_0 The prior sample size of the mean of the
//'     column-specific precisions of the factors.
//' @param nu_0 The prior mean of the mean of the column-specific
//'     precisions of the factors.
//' @param hetero_factors A logical. Should we also update the precisions of the
//'     factors (\code{TRUE}), or not (\code{FALSE})?
//' @param display_progress A logical. If \code{TRUE}, then a progress bar will
//'     be displayed and you'll be able to interupt the C++ code. If \code{FALSE},
//'     then neither of these capabilities will be provided.
//'
//' @author David Gerard
//' @export
//'
//'
// [[Rcpp::export]]
List bfa_gd_gibbs(NumericMatrix Linit, NumericMatrix Finit, NumericVector xi_init,
		    double phi_init, NumericVector zeta_init, NumericVector theta_init,
		    double kappa_init, NumericMatrix Y22init, NumericMatrix Y21,
		    NumericMatrix Y31, NumericMatrix Y32, int nsamp, int burnin,
		    int thin, double rho_0, double alpha_0, double delta_0,
		    double lambda_0, double nu_0, double beta_0, double eta_0,
		    double tau_0, bool hetero_factors, bool display_progress) {

  // Get dimensions of all matrices -------------------------------------------
  int n         = Linit.nrow(); // sample size
  int p         = Finit.ncol(); // number of genes
  int nfac      = Linit.ncol(); // number of factors
  int ncov      = Y21.nrow();   // number of covariates
  int ncontrols = Y21.ncol();   // number of controls

  // Initialize arma matrices and vectors -------------------------------------
  // mat, colvec, join_rows, and join_cols are from arma
  mat Lcurrent(Linit.begin(), n, nfac, false);
  mat Fcurrent(Finit.begin(), nfac, p, false);
  colvec xi_current(xi_init.begin(), p, false);
  colvec theta_current(theta_init.begin(), p, false);
  colvec zeta_current(zeta_init.begin(), nfac, false);
  double phi_current = phi_init;
  double kappa_current = kappa_init;
  mat Y22current(Y22init.begin(), ncov, p - ncontrols, false);
  mat Y21a(Y21.begin(), ncov, ncontrols, false);
  mat Y31a(Y31.begin(), n - ncov, ncontrols, false);
  mat Y32a(Y32.begin(), n - ncov, p - ncontrols, false);
  mat Y3a = join_rows(Y31a, Y32a);
  mat Ycurrent = join_cols(join_rows(Y21a, Y22current), Y3a);

  // calculate which indices I keep samples -----------------------------------
  int nkeeps = floor(nsamp / thin);
  IntegerVector keep_indices(nkeeps);
  for (int keep_index = 0; keep_index < nkeeps; keep_index++) {
    keep_indices[keep_index] = burnin + thin * keep_index;
  }
  cube Y22array(ncov, p - ncontrols, nkeeps);
  mat xi_mat(p, nkeeps);
  vec phi_vec(nkeeps);
  int thin_index = 0;
  int current_thin_index = keep_indices[0];

  // calculate which indices I update progress message ------------------------
  int nchecks = 100;
  int check_thin = floor((nsamp + burnin) / nchecks);
  Progress prog(nchecks, display_progress); // instance for progress bar
  IntegerVector check_indices(nchecks);
  for (int init_ind = 0; init_ind < nchecks; init_ind++) {
    check_indices[init_ind] = init_ind * check_thin;
  }
  int current_check_index = check_indices[0];
  int check_index = 0;

  for (int gindex = 0; gindex < nsamp + burnin; gindex++) {

    // Update progress bar and check for aborted job --------------------------
    if (gindex == current_check_index) {
      if (Progress::check_abort()) {
	return List::create("Aborted");
      }
      prog.increment();
      check_index++;
      current_check_index = check_indices[check_index];
    }


    // Update Lcurrent --------------------------------------------------------
    vec Leigval;
    mat Leigvec;
    eig_sym(Leigval, Leigvec, Fcurrent * diagmat(xi_current) * trans(Fcurrent) +
	    diagmat(zeta_current));
    mat Lcovhalf = Leigvec * inv(diagmat(sqrt(Leigval))) * trans(Leigvec);
    mat Lmean    = Ycurrent * diagmat(xi_current) * trans(Fcurrent) * Leigvec *
      inv(diagmat(Leigval)) * trans(Leigvec);
    Lcurrent = Lmean + randn(n, nfac) * Lcovhalf;

    // Update Fcurrent --------------------------------------------------------
    for (int fcol = 0; fcol < p; fcol++) {
      colvec thetaj(nfac);
      thetaj.fill(theta_current[fcol]);
      vec Feigval;
      mat Feigvec;
      eig_sym(Feigval, Feigvec, trans(Lcurrent) * Lcurrent * xi_current[fcol] + diagmat(thetaj));
      mat Fcovhalf = Feigvec * inv(diagmat(sqrt(Feigval))) * trans(Feigvec);
      mat Fmean = Feigvec * inv(diagmat(Feigval)) * trans(Feigvec) * trans(Lcurrent) *
	Ycurrent.col(fcol) * xi_current[fcol];
      Fcurrent.col(fcol) = Fmean + Fcovhalf * randn(nfac);
    }

    // Get current mean estimate ----------------------------------------------
    mat mean_current = Lcurrent * Fcurrent;

    // Update xi_current ------------------------------------------------------
    // They only allow for scale, not rate
    colvec resvec = trans(sum(square(Ycurrent - mean_current)));
    double xi_shape = (n + rho_0) / 2;
    for (int xindex = 0; xindex < p; xindex++) {
      double xi_scale = 2 / (resvec[xindex] + rho_0 * phi_current);
      xi_current[xindex] = R::rgamma(xi_shape, xi_scale);
    }

    // Update phi_current -----------------------------------------------------
    double phi_shape = (p * rho_0 + alpha_0) / 2;
    double phi_scale = 2 / (alpha_0 * beta_0 + rho_0 * sum(xi_current));
    phi_current = R::rgamma(phi_shape, phi_scale);

    // Update zeta_current ----------------------------------------------------
    colvec Lsumvec = trans(sum(square(Lcurrent), 0));
    double zeta_shape = (n + eta_0) / 2;
    for (int zindex = 0; zindex < nfac; zindex++) {
      double zeta_scale = 2 / (eta_0 * tau_0 + Lsumvec[zindex]);
      zeta_current[zindex] = R::rgamma(zeta_shape, zeta_scale);
    }

    if (hetero_factors) {
      // Update theta_current -------------------------------------------------
      colvec Fsumvec = trans(sum(square(Fcurrent), 0));
      double theta_shape = (nfac + delta_0) / 2;
      for (int theta_index = 0; theta_index < p; theta_index++) {
	double theta_scale = 2 / (delta_0 * kappa_current + Fsumvec[theta_index]);
	theta_current[theta_index] = R::rgamma(theta_shape, theta_scale);
      }

      // Update kappa_current -------------------------------------------------
      double kappa_shape = (p * delta_0 + lambda_0) / 2;
      double kappa_scale = 2 / (lambda_0 * nu_0 + delta_0 * sum(theta_current));
      kappa_current = R::rgamma(kappa_shape, kappa_scale);
    }

    // Update Y22current and Ycurrent -----------------------------------------
    mat Y22mean = mean_current.submat(0, ncontrols, ncov - 1, p - 1);
    Y22current = Y22mean + randn(ncov, p - ncontrols) *
      inv(diagmat(sqrt(xi_current.subvec(ncontrols, p - 1))));
    Ycurrent.submat(0, ncontrols, ncov - 1, p - 1) = Y22current;

    // Keep samples if correct index -------------------------------------------
    if (gindex == current_thin_index) {
      Y22array.slice(thin_index) = Y22current;
      xi_mat.col(thin_index) = xi_current;
      phi_vec[thin_index] = phi_current;
      thin_index++;
      current_thin_index = keep_indices[thin_index];
    }
  }

  return List::create(_["Y22_array"] = Y22array,
		      _["xi"] = xi_mat,
		      _["phi"] = phi_vec);
}





//' Fast Gibbs sampler for Bayesian factor analysis.
//'
//' This is very similar to \code{\link{bfa_gd_gibbs}} except that we link the precisions
//' of the observations with the precisions of the factors. For some reason, this works
//' very well in practice.
//'
//'
//'
//'
//'
//'
//' @param Y22init A matrix of numerics. The initial value of Y22.
//' @param Y21 A matrix of numerics.
//' @param Y31 A matrix of numerics.
//' @param Y32 A matrix of numerics.
//' @param nsamp The number of iterations to run in the Gibbs sampler, not including the burnin.
//' @param burnin The number of iterations to burnin.
//' @param thin We only collect samples every \code{thin} iterations.
//' @param Linit A numeric matrix. The initial values of the loadings.
//' @param Finit A numeric matrix. The initial values for the factors.
//' @param xi_init A numeric vector. The initial values of the precisions.
//' @param phi_init A numeric scalar. The initial value of the mean of the precitions.
//' @param zeta_init A numeric vector. The initial values of the augmented row precisions.
//' @param rho_0 The prior sample size for column-specific the
//'     precisions.
//' @param alpha_0 The prior sample size for the mean of the
//'     column-specific precisions.
//' @param beta_0 The prior mean of the mean of the column-specific
//'     precisions.
//' @param eta_0 The prior sample size of the expanded parameters.
//' @param tau_0 The prior mean of of the expanded parameters.
//' @param display_progress A logical. If \code{TRUE}, then a progress bar will
//'     be displayed and you'll be able to interupt the C++ code. If \code{FALSE},
//'     then neither of these capabilities will be provided.
//'
//' @author David Gerard
//' @export
//'
//'
// [[Rcpp::export]]
List bfa_gs_linked_gibbs(NumericMatrix Linit, NumericMatrix Finit, NumericVector xi_init,
			 double phi_init, NumericVector zeta_init,
			 NumericMatrix Y22init, NumericMatrix Y21,
			 NumericMatrix Y31, NumericMatrix Y32, int nsamp, int burnin,
			 int thin, double rho_0, double alpha_0,
			 double beta_0, double eta_0,
			 double tau_0, bool display_progress) {

  // Get dimensions of all matrices -------------------------------------------
  int n         = Linit.nrow(); // sample size
  int p         = Finit.ncol(); // number of genes
  int nfac      = Linit.ncol(); // number of factors
  int ncov      = Y21.nrow();   // number of covariates
  int ncontrols = Y21.ncol();   // number of controls

  // Initialize arma matrices and vectors -------------------------------------
  // mat, colvec, join_rows, and join_cols are from arma
  mat Lcurrent(Linit.begin(), n, nfac, false);
  mat Fcurrent(Finit.begin(), nfac, p, false);
  colvec xi_current(xi_init.begin(), p, false);
  colvec zeta_current(zeta_init.begin(), nfac, false);
  double phi_current = phi_init;
  mat Y22current(Y22init.begin(), ncov, p - ncontrols, false);
  mat Y21a(Y21.begin(), ncov, ncontrols, false);
  mat Y31a(Y31.begin(), n - ncov, ncontrols, false);
  mat Y32a(Y32.begin(), n - ncov, p - ncontrols, false);
  mat Y3a = join_rows(Y31a, Y32a);
  mat Ycurrent = join_cols(join_rows(Y21a, Y22current), Y3a);

  // calculate which indices I keep samples -----------------------------------
  int nkeeps = floor(nsamp / thin);
  IntegerVector keep_indices(nkeeps);
  for (int keep_index = 0; keep_index < nkeeps; keep_index++) {
    keep_indices[keep_index] = burnin + thin * keep_index;
  }
  cube Y22array(ncov, p - ncontrols, nkeeps);
  mat xi_mat(p, nkeeps);
  vec phi_vec(nkeeps);
  int thin_index = 0;
  int current_thin_index = keep_indices[0];

  // calculate which indices I update progress message ------------------------
  int nchecks = 100;
  int check_thin = floor((nsamp + burnin) / nchecks);
  Progress prog(nchecks, display_progress); // instance for progress bar
  IntegerVector check_indices(nchecks);
  for (int init_ind = 0; init_ind < nchecks; init_ind++) {
    check_indices[init_ind] = init_ind * check_thin;
  }
  int current_check_index = check_indices[0];
  int check_index = 0;

  for (int gindex = 0; gindex < nsamp + burnin; gindex++) {

    // Update progress bar and check for aborted job --------------------------
    if (gindex == current_check_index) {
      if (Progress::check_abort()) {
	return List::create("Aborted");
      }
      prog.increment();
      check_index++;
      current_check_index = check_indices[check_index];
    }


    // Update Lcurrent --------------------------------------------------------
    vec Leigval;
    mat Leigvec;
    eig_sym(Leigval, Leigvec, Fcurrent * diagmat(xi_current) * trans(Fcurrent) +
	    diagmat(zeta_current));
    mat Lcovhalf = Leigvec * inv(diagmat(sqrt(Leigval))) * trans(Leigvec);
    mat Lmean    = Ycurrent * diagmat(xi_current) * trans(Fcurrent) * Leigvec *
      inv(diagmat(Leigval)) * trans(Leigvec);
    Lcurrent = Lmean + randn(n, nfac) * Lcovhalf;

    // Update Fcurrent --------------------------------------------------------
    vec Feigval;
    mat Feigvec;
    eig_sym(Feigval, Feigvec, trans(Lcurrent) * Lcurrent + eye(nfac, nfac));
    mat Frowcovhalf = Feigvec * inv(diagmat(sqrt(Feigval))) * trans(Feigvec);
    mat Fmean       = Feigvec * inv(diagmat(Feigval)) * trans(Feigvec) * trans(Lcurrent) *
      Ycurrent;
    Fcurrent = Fmean + Frowcovhalf * randn(nfac, p) * inv(diagmat(sqrt(xi_current)));

    // Get current mean estimate ----------------------------------------------
    mat mean_current = Lcurrent * Fcurrent;

    // Update xi_current ------------------------------------------------------
    // They only allow for scale, not rate
    colvec resvec = trans(sum(square(Ycurrent - mean_current)));
    colvec uvec  = trans(sum(square(Fcurrent)));
    double xi_shape = (n + rho_0) / 2;
    for (int xindex = 0; xindex < p; xindex++) {
      double xi_scale = 2 / (resvec[xindex] + uvec[xindex] + rho_0 * phi_current);
      xi_current[xindex] = R::rgamma(xi_shape, xi_scale);
    }

    // Update phi_current -----------------------------------------------------
    double phi_shape = (p * rho_0 + alpha_0) / 2;
    double phi_scale = 2 / (alpha_0 * beta_0 + rho_0 * sum(xi_current));
    phi_current = R::rgamma(phi_shape, phi_scale);

    // Update zeta_current ----------------------------------------------------
    colvec Lsumvec = trans(sum(square(Lcurrent), 0));
    double zeta_shape = (n + eta_0) / 2;
    for (int zindex = 0; zindex < nfac; zindex++) {
      double zeta_scale = 2 / (eta_0 * tau_0 + Lsumvec[zindex]);
      zeta_current[zindex] = R::rgamma(zeta_shape, zeta_scale);
    }


    // Update Y22current and Ycurrent -----------------------------------------
    mat Y22mean = mean_current.submat(0, ncontrols, ncov - 1, p - 1);
    Y22current = Y22mean + randn(ncov, p - ncontrols) *
      inv(diagmat(sqrt(xi_current.subvec(ncontrols, p - 1))));
    Ycurrent.submat(0, ncontrols, ncov - 1, p - 1) = Y22current;

    // Keep samples if correct index -------------------------------------------
    if (gindex == current_thin_index) {
      Y22array.slice(thin_index) = Y22current;
      xi_mat.col(thin_index) = xi_current;
      phi_vec[thin_index] = phi_current;
      thin_index++;
      current_thin_index = keep_indices[thin_index];
    }
  }

  return List::create(_["Y22_array"] = Y22array,
		      _["xi"] = xi_mat,
		      _["phi"] = phi_vec);
}
