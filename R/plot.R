#' Make a plot of the output of cv_trac
#'
#' @param cvfit_trac an object as returned by \code{\link{cv_trac}}.
#' @param iw vector of indices specifying which weight sequence solutions to
#' include
#' @param superimpose whether to superimpose
#' @export
plot_cv_trac <- function(cvfit_trac, iw = NULL, superimpose = TRUE) {
  num_w <- length(cvfit_trac$cv)
  if (is.null(iw)) iw <- 1:num_w
  if (!superimpose) {
    r <- floor(sqrt(num_w))
    graphics::par(mfrow = c(r, ceiling(length(iw) / r)))
    lapply(1:num_w, function(iw) plot_cv_trac_single_w(cvfit_trac$cv[[iw]]))
  } else {
    x <- cvfit_trac$cv
    graphics::par(mar = c(5, 5, 5, 1))
    xrang <- range(lapply(x, function(xx) xx$nonzeros))
    yrang <- range(lapply(x, function(xx) c(xx$m - xx$se, xx$m + xx$se)),
                   na.rm = TRUE)
    graphics::plot(0, 0, xlab = "Number of nonzero gamma",
                   ylab = "Cross-validation Error",
                   type = "n", xlim = xrang, ylim = yrang)
    for (iw in seq_along(x)) {
      ggb:::error_bars(x[[iw]]$nonzeros,
                       x[[iw]]$m - x[[iw]]$se,
                       x[[iw]]$m + x[[iw]]$se, width = 0.01,
                       col = "darkgrey")
      graphics::lines(x[[iw]]$nonzeros, x[[iw]]$m, col = "darkgrey", pch = 19)
      graphics::points(x[[iw]]$nonzeros, x[[iw]]$m, pch = 19, col = iw)
      graphics::abline(v = x[[iw]]$nonzeros[x[[iw]]$ibest], lty = 3, col = iw)
      graphics::abline(v = x[[iw]]$nonzeros[x[[iw]]$i1se], lty = 3, col = iw)
    }
  }
  invisible()
}

plot_cv_trac_single_w <- function(cvfit_trac_single_w) {
  x <- cvfit_trac_single_w
  graphics::par(mar = c(5, 5, 5, 1))
  yrang <- range(c(x$m - x$se, x$m + x$se))
  graphics::plot(log(x$fraclist), x$m, xlab = "log(lambda)",
                 ylab = "Cross-validation Error",
                 type = "n", ylim = c(0,1))
  graphics::axis(3, at = log(x$fraclist), labels = paste(x$nonzero), srt = 90,
                 adj = 0)
  graphics::mtext("Number of nonzero gamma", 3, 4, cex = 1.2)
  ggb:::error_bars(log(x$fraclist), x$m - x$se, x$m + x$se, width = 0.01,
                   col = "darkgrey")
  graphics::points(log(x$fraclist), x$m, col = 2, pch = 19)
  graphics::abline(v = log(x$lambda_best), lty = 3)
  graphics::abline(v = log(x$lambda_1se), lty = 3)
  invisible()
}

#' Plot trac coefficient path
#'
#' @param fit an object as returned by \code{\link{trac}}.
#' @param iw index specifying which weight sequence solution to plot
#' @param coef which type of parameter to show in the plot
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
plot_trac_path <- function(fit, iw = 1, coef = c("alpha", "beta", "gamma")) {
  coef <- match.arg(coef)
  fit[[iw]][[coef]] %>%
    as.matrix() %>%
    t() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(frac = fit[[iw]]$fraclist) %>%
    tidyr::pivot_longer(cols = -.data$frac) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$frac,
                                 y = .data$value,
                                 group = .data$name,
                                 color = .data$name)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_log10() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(y = coef, x = "Fraction of lambda_max")
}
