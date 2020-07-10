#' Make predictions based on a trac fit
#'
#' @param fit output of the function \code{\link{trac}}
#' @param new_Z a new data matrix (see \code{Z} from \code{\link{trac}})
#'
#' @return a vector of \code{nrow(new_Z)} predictions.
#' @export
predict_trac <- function(fit, new_Z) {
  # fit: output of wag
  # new_Z: n_new by p matrix
  yhat <- list()
  for (iw in seq_along(fit)) {
    yhat[[iw]] <- t(as.numeric(fit[[iw]]$beta0) + t(as.matrix(new_Z %*% fit[[iw]]$beta)))
    rownames(yhat[[iw]]) <- rownames(new_Z)
  }
  return(yhat)
}


