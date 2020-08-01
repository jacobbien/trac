#' Make a plot of the output of cv_trac
#' @export
plot_cv_trac <- function(cvfit_trac, iw = NULL, superimpose = TRUE) {
  num_w <- length(cvfit_trac$cv)
  if (is.null(iw)) iw <- 1:num_w
  if (!superimpose) {
    r <- floor(sqrt(num_w))
    par(mfrow = c(r, ceiling(length(iw) / r)))
    lapply(1:num_w, function(iw) plot_cv_trac_single_w(cvfit_trac$cv[[iw]]))
  } else {
    x <- cvfit_trac$cv
    graphics::par(mar = c(5, 5, 5, 1))
    xrang = range(lapply(x, function(xx) xx$nonzeros))
    yrang = range(lapply(x, function(xx) c(xx$m - xx$se, xx$m + xx$se)),
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
  yrang = range(c(x$m - x$se, x$m + x$se))
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

plot_cv_lrtrac <- function(cvlrfit, iw) {
  nlam1 <- length(cvlrfit$cv[[iw]])
  m <- lapply(cvlrfit$cv[[iw]], function(cc) cc$m)
  flam2 <- lapply(cvlrfit$cv[[iw]],
                  function(cc) cc$lam2list / max(cc$lam2list))
  plot(1, 1, type = "n", xlim = range(unlist(flam2)),
       ylim = range(unlist(m)), xlab = "Fraction of lam2max",
       ylab = "CV error", log = "x")
  for (ii in seq(nlam1)) {
    if (length(m[[ii]]) == 0) next
    lines(flam2[[ii]], m[[ii]], col = rgb(0, 0, 0, 0.5))
  }
  ilam1 <- cvlrfit$ilam1_best[[iw]]
  ilam2 <- cvlrfit$ilam2_best[[iw]]
  lines(flam2[[ilam1]], m[[ilam1]], col = 2, lwd = 2)
  abline(v = flam2[[ilam1]][ilam2], col = 2)
  legend("topleft", legend = c(paste0("CV best")), fill = c(2))
}
