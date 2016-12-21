# Plotting methods for mouthwash, backwash, vruv4, and ruv3

#' Plotting method for \code{moutwash}.
#'
#' This will produce
#'
#' @param x The output of \code{\link{mouthwash}}, of class \code{mouthwash}.
#' @param ... Not used.
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso \code{\link{mouthwash}}.
#'
plot.mouthwash <- function(x, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needs to be installed to use the plot method.")
  }

  if (class(x$fitted_g) == "normalmix") {
    ## Plot unimodal distribution ------------------------------------------------------
    dfdat <- data.frame(pi = x$fitted_g$pi, sd = x$fitted_g$sd)
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "sd", xend = "sd", y = 0, yend = "pi")) +
      ggplot2::geom_segment() +
      ggplot2::theme_bw() +
      ggplot2::xlab("Mixing Standard Deviations") +
      ggplot2::ylab("Mixing Proportions") +
      ggplot2::ggtitle("Mixing SD's")
    print(pl)

    cat ("Press [enter] to continue")
    line <- readline()

    xlim <- max(3 * max(x$fitted_g$sd[x$fitted_g$pi > 0.01]), 1)

    quant_seq <- seq(-xlim, xlim, length = 100)
    dense_seq <- dnormalmix(x = quant_seq, mixdense = x$fitted_g)
    dfdat <- data.frame(quantile = quant_seq, density = dense_seq)
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "quantile", y = "density")) +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::ggtitle("Effects Mixture Density") +
      ggplot2::geom_segment(x = 0, xend = 0, y = 0, yend = x$pi0, lty = 2) +
      ggplot2::ylim(0, max(x$pi0, dense_seq))
    print(pl)

  } else if (class(x$fitted_g) == "unimix") {
    dfdat <- data.frame(pi = c(x$fitted_g$pi, x$fitted_g$pi), bounds = c(x$fitted_g$a, x$fitted_g$b))
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "bounds", xend = "bounds", y = 0, yend = "pi")) +
      ggplot2::geom_segment() +
      ggplot2::theme_bw() +
      ggplot2::xlab("Mixing Uniform Bounds (Lower and Upper)") +
      ggplot2::ylab("Mixing Proportions")
    print(pl)
  }

  cat ("Press [enter] to continue")
  line <- readline()

  pl <- ggplot2::ggplot(data = x$result, mapping = ggplot2::aes_string(x = "betahat", y = "PosteriorMean")) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::xlab("OLS Estimates") +
    ggplot2::ylab("MOUTHWASH Estimates") +
    ggplot2::ggtitle("Shrinkage of OLS Estimates") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = 2, col = "gray50") +
    ggplot2::geom_hline(yintercept = 0, lty = 2, col = "gray50")
  print(pl)

  line <- ""
  index <- 1
  while (line != "Y" & line != "y" & line != "yes" & line != "YES" & line != "Yes" &
         line != "N" & line != "n" & line != "no" & line != "NO" & line != "No" & index < 4) {
    cat ("Compare to ASH? Y/N")
    line <- readline()
    index <- index + 1
  }

  if (index == 4) {
    cat("OK. Try again later.\n")
  }

  if (line == "Y" | line == "y" | line == "yes" | line == "YES" | line == "Yes") {
    ash_args <- list()
    ash_args$betahat <- x$result$betahat
    ash_args$sebetahat <- x$result$sebetahat
    ash_out <- do.call(what = ashr::ash, args = ash_args)
    dfdat <- data.frame(MOUTHWASH = x$result$lfdr, ASH = ashr::get_lfdr(ash_out))
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "ASH", y = "MOUTHWASH")) +
      ggplot2::geom_point() +
      ggplot2::ggtitle("LFDR's") +
      ggplot2::theme_bw() +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = 2, col = "gray50")
    print(pl)

    cat ("Press [enter] to continue")
    line <- readline()

    dfdat <- data.frame(MOUTHWASH = x$result$PosteriorMean, ASH = ashr::get_pm(ash_out))
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "ASH", y = "MOUTHWASH")) +
      ggplot2::geom_point() +
      ggplot2::ggtitle("Posterior Means") +
      ggplot2::theme_bw() +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = 2, col = "gray50")
    print(pl)
  }
}

#' Plotting method for \code{\link{backwash}}.
#'
#' @param x An object of class \code{backwash}, returned from \code{\link{backwash}}.
#' @param ... Not used.
#'
#' @author David Gerard
#'
#' @export
#'
plot.backwash <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needs to be installed to use the plot method.")
  }

  ## Plot unimodal distribution ------------------------------------------------------
  fitted_g <- list()
  fitted_g$pi <- x$fitted_g$pivec
  fitted_g$sd <- sqrt(x$fitted_g$tau2_seq)
  fitted_g$mean <- rep(0, length = length(x$fitted_g$tau2_seq))
  class(fitted_g) <- "normalmix"

  dfdat <- data.frame(pi = fitted_g$pi, sd = fitted_g$sd)
  pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "sd", xend = "sd", y = 0, yend = "pi")) +
    ggplot2::geom_segment() +
    ggplot2::theme_bw() +
    ggplot2::xlab("Mixing Standard Deviations") +
    ggplot2::ylab("Mixing Proportions") +
    ggplot2::ggtitle("Mixing SD's")
  print(pl)

  cat ("Press [enter] to continue")
  line <- readline()

  xlim <- max(3 * max(fitted_g$sd[fitted_g$pi > 0.01]), 1)

  quant_seq <- seq(-xlim, xlim, length = 100)
  dense_seq <- dnormalmix(x = quant_seq, mixdense = fitted_g)
  dfdat <- data.frame(quantile = quant_seq, density = dense_seq)
  pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "quantile", y = "density")) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Effects Mixture Density") +
    ggplot2::geom_segment(x = 0, xend = 0, y = 0, yend = x$pi0, lty = 2) +
    ggplot2::ylim(0, max(x$pi0, dense_seq))
  print(pl)

  cat ("Press [enter] to continue")
  line <- readline()

  pl <- ggplot2::ggplot(data = x$result, mapping = ggplot2::aes_string(x = "betahat", y = "PosteriorMean")) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::xlab("OLS Estimates") +
    ggplot2::ylab("BACKWASH Estimates") +
    ggplot2::ggtitle("Shrinkage of OLS Estimates") +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = 2, col = "gray50") +
    ggplot2::geom_hline(yintercept = 0, lty = 2, col = "gray50")
  print(pl)

  line <- ""
  index <- 1
  while (line != "Y" & line != "y" & line != "yes" & line != "YES" & line != "Yes" &
         line != "N" & line != "n" & line != "no" & line != "NO" & line != "No" & index < 4) {
    cat ("Compare to ASH? Y/N")
    line <- readline()
    index <- index + 1
  }

  if (index == 4) {
    cat("OK. Try again later.\n")
  }

  if (line == "Y" | line == "y" | line == "yes" | line == "YES" | line == "Yes") {
    ash_args <- list()
    ash_args$betahat <- x$result$betahat
    ash_args$sebetahat <- x$result$sebetahat
    ash_out <- do.call(what = ashr::ash, args = ash_args)
    dfdat <- data.frame(BACKWASH = x$result$lfdr, ASH = ashr::get_lfdr(ash_out))
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "ASH", y = "BACKWASH")) +
      ggplot2::geom_point() +
      ggplot2::ggtitle("LFDR's") +
      ggplot2::theme_bw() +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = 2, col = "gray50")
    print(pl)

    cat ("Press [enter] to continue")
    line <- readline()

    dfdat <- data.frame(BACKWASH = x$result$PosteriorMean, ASH = ashr::get_pm(ash_out))
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "ASH", y = "BACKWASH")) +
      ggplot2::geom_point() +
      ggplot2::ggtitle("Posterior Means") +
      ggplot2::theme_bw() +
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = 2, col = "gray50")
    print(pl)
  }

}


#' Plotting method for \code{ruvb}.
#'
#' For the traceplots, the blue dashed line is a non-parameteric smoother.
#' This should let you see if there are any major trends in the traceplot.
#'
#' @param x The output of \code{\link{ruvb}}, of class \code{ruvb}.
#' @param ... Not used.
#'
#' @author David Gerard
#'
#'
#' @seealso \code{\link{ruvb}}
#'
plot.ruvb <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("To use this plotting method, you need to install ggplot2.")
  }

  if (is.null(x$betahat_post)) {
    message("If you want traceplots, rerun ruvb with return_mcmc = TRUE")
  } else {

    effective_sizes <- apply(X = x$betahat_post, MARGIN = c(1, 2), FUN = coda::effectiveSize)
    dfdat <- data.frame(effective_size = c(effective_sizes))
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "effective_size")) +
      ggplot2::geom_histogram(bins = 30, col = "black", fill = "white") +
      ggplot2::theme_bw() +
      ggplot2::xlab("Effective Sample Size") +
      ggplot2::ggtitle("Effective Sample Sizes")
    print(pl)

    cat("Would you care for a traceplot? Y/N")
    line <- readline()
    while (line == "Y" | line == "y" | line == "YES" | line == "Yes" | line == "yes") {
      rownum <- sample(1:dim(x$betahat_post)[1], size = 1)
      colnum <- sample(1:dim(x$betahat_post)[2], size = 1)

      betahat_post <- x$betahat_post[rownum, colnum, ]
      dfdat <- data.frame(index = 1:length(betahat_post), betahat = betahat_post)
      pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "index", y = "betahat")) +
        ggplot2::geom_line() +
        ggplot2::geom_smooth(method = "loess", se = FALSE, lty = 2) +
        ggplot2::theme_bw() +
        ggplot2::xlab("index") +
        ggplot2::ylab("betahat") +
        ggplot2::ggtitle(paste0("Traceplot of (", rownum, ", ", colnum, ")th Element"))
      print(pl)

      cat("Would you care for another traceplot? Y/N")
      line <- readline()
    }

    cat("Would you care for a histogram? Y/N")
    line <- readline()
    while (line == "Y" | line == "y" | line == "YES" | line == "Yes" | line == "yes") {
      rownum <- sample(1:dim(x$betahat_post)[1], size = 1)
      colnum <- sample(1:dim(x$betahat_post)[2], size = 1)

      betahat_post <- x$betahat_post[rownum, colnum, ]
      dfdat <- data.frame(betahat = betahat_post)
      pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "betahat")) +
        ggplot2::geom_histogram(fill = "white", color = "black", bins = 30) +
        ggplot2::theme_bw() +
        ggplot2::xlab("index") +
        ggplot2::ylab("betahat") +
        ggplot2::ggtitle(paste0("Histogram of (", rownum, ", ", colnum, ")th Element")) +
        ggplot2::geom_vline(xintercept = mean(betahat_post), lty = 2, col = "red")
      print(pl)

      cat("Would you care for another histogram? Y/N")
      line <- readline()
    }
  }

}



