#' Refit subject to sparsity constraints
#'
#' Given output of \code{\link{trac}}, solves the least squares problem with
#' compositional constraint on the features selected by trac.
#'
#' minimize_{beta, beta0, gamma} 1/(2n) || y - beta0 1_n - Z beta ||^2
#' subject to gamma_{trac nonselected} = 0, beta = A gamma, 1_p^T beta = 0
#'
#' @param fit output of trac
#' @param Z,y,A same arguments as passed to \code{\link{trac}}
#' @param tol tolerance for deciding whether a gamma value is zero
#' @export
refit_trac <- function(fit, Z, y, A, tol = 1e-5) {
  n <- nrow(Z)
  p <- ncol(Z)
  t_size <- ncol(A) + 1
  nlam <- lapply(fit, function(ft) length(ft$fraclist))

  for (iw in seq_along(fit)) {
    for (i in seq(nlam[[iw]])) {
      # let nz be nonzero elements from trac and z be the zero ones:
      z <- which(abs(fit[[iw]]$gamma[, i]) <= tol)
      if (length(z) == t_size - 1) {
        # all zero
        fit[[iw]]$beta0[i] = fit[[iw]]$beta0[i]
        fit[[iw]]$beta[, i] = fit[[iw]]$beta[, i]
        fit[[iw]]$gamma[, i] = fit[[iw]]$gamma[, i]
        next
      }

      # minimize_{beta0, gamma} 1/(2n) || y - beta0 1_n - Z A gamma ||^2
      # subject to gamma_z = 0, 1_p^T A gamma = 0
      # can be written as
      # minimize_{beta0, gamma} 1/(2n) || y - beta0 1_n - Z A_nz gamma_nz ||^2
      # subject to gamma_z = 0, v^T gamma_nz = 0
      # with constraint vector v = A_nz^T 1_p
      v <- Matrix::colSums(A[, -z]) # has length of nz

      # the above is equiv to solving
      # minimize_{beta0, gamma_nz} 1/(2n) || y - beta0 1_n - Z A_nz (I - P)gamma_nz ||^2
      # where P = vv^T / ||v||^2
      # or
      # minimize_{beta0, gamma_nz} 1/(2n) || y - M [beta0; gamma_nz] ||^2
      # where M = [1_n, Z A_nz (I - P)]
      # or [beta0; gamma_nz] = M^{+} y
      ZA <- Z %*% A[, -z]
      ZAP <- (ZA %*% v) %*% matrix(v, nrow = 1) / sum(v^2)
      M <- cbind(1, ZA - ZAP)
      sv <- svd(M)
      # M^{+} y = V D^{+} U^T y
      solution <- sv$v %*% (c(1 / sv$d[-length(sv$d)], 0) * crossprod(sv$u, y))
      fit[[iw]]$gamma[z, i] = 0
      fit[[iw]]$gamma[-z, i] = solution[-1]
      fit[[iw]]$alpha[, i] <- fit[[iw]]$gamma[, i]
      fit[[iw]]$alpha[-z, i] <- fit[[iw]]$gamma[-z, i] * v
      fit[[iw]]$beta0[i] = solution[1]
      fit[[iw]]$beta[, i] = A[, -z] %*% fit[[iw]]$gamma[-z, i]
    }
    fit[[iw]]$tol <- tol
    fit[[iw]]$refit <- TRUE
  }
  fit
}
