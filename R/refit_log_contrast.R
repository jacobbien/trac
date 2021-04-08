#' Refit subject to sparsity constraints
#'
#' Given output of \code{\link{sparse_log_contrast}}, solves the least squares problem with
#' compositional constraint on the features selected by sparse_log_contrast.
#'
#' minimize_{beta, beta0} 1/(2n) || y - beta0 1_n - Z beta ||^2
#' subject to beta_{nonselected} = 0, 1_p^T beta = 0
#'
#' @param fit output of trac
#' @param Z,y same arguments as passed to \code{\link{sparse_log_contrast}}
#' @param tol tolerance for deciding whether a beta value is zero
#' @export
refit_sparse_log_contrast <- function(fit, Z, y, tol = 1e-5) {
  n <- nrow(Z)
  p <- ncol(Z)
  nlam <- length(fit$fraclist)

  for (i in seq(nlam)) {
    # let nz be nonzero elements from sparse log contrast and z be the zero ones:
    z <- which(abs(fit$beta[, i]) <= tol)
    if (length(z) == p) {
      # all zero
      fit$beta0[i] <- fit$beta0[i]
      fit$beta[, i] <- fit$beta[, i]
      next
    }

    # minimize_{beta0, beta} 1/(2n) || y - beta0 1_n - Z beta ||^2
    # subject to beta_z = 0, 1_p^T beta = 0
    # can be written as
    # minimize_{beta0, beta} 1/(2n) || y - beta0 1_n - Z_nz beta_nz ||^2
    # subject to beta_z = 0, 1_nz^T beta_nz = 0
    v <- rep(1, p - length(z))

    # the above is equiv to solving
    # minimize_{beta0, beta_nz} 1/(2n) || y - beta0 1_n - Z_nz (I - P)beta_nz ||^2
    # where P = vv^T / ||v||^2
    # or
    # minimize_{beta0, beta_nz} 1/(2n) || y - M [beta0; beta_nz] ||^2
    # where M = [1_n, Z_nz (I - P)]
    # or [beta0; beta_nz] = M^{+} y
    if (length(z) == 0)
      Zz <- Z
    else
      Zz <- Z[, -z]
    ZzP <- (Zz %*% v) %*% matrix(v, nrow = 1) / sum(v^2)
    M <- cbind(1, Zz - ZzP)
    sv <- svd(M)
    # M^{+} y = V D^{+} U^T y
    solution <- sv$v %*% (c(1 / sv$d[-length(sv$d)], 0) * crossprod(sv$u, y))
    fit$beta[z, i] <-  0
    if (length(z) == 0)
      fit$beta[, i] <- solution[-1]
    else
      fit$beta[-z, i] <- solution[-1]
    fit$beta0[i] <- solution[1]
  }
  fit$tol <- tol
  fit$refit <- TRUE
  fit
}
