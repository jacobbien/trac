#' Make predictions based on a trac fit
#'
#' @param fit output of the function \code{\link{trac}}
#' @param new_Z a new data matrix (see \code{Z} from \code{\link{trac}})
#' @param new_X a new data matrix (see \code{X} from \code{\link{trac}})
#' @return a vector of \code{nrow(new_Z) + nrow(new_X)} predictions.
#' @export
predict_trac <- function(fit, new_Z, new_X = NULL) {
  # fit: output of wag
  # new_Z: n_new by p matrix
  classification <- fit[[1]]$method %in% c("classificiation",
                                           "classification_huber")
  intercept_classif <- fit[[1]]$intercept_classif
  yhat <- list()
  for (iw in seq_along(fit)) {
    if (!is.null(new_X)) {
      new_Z <- cbind(new_Z, new_X)
    }
    if (intercept_classif) {
      yhat[[iw]] <- t(as.numeric(fit[[iw]]$beta0) +
        t(as.matrix(new_Z %*% fit[[iw]]$beta)))
    } else {
      yhat[[iw]] <- t(t(as.matrix(new_Z %*% fit[[iw]]$beta)))
    }
#    if (classification) {
#      yhat[[iw]] <- yhat[[iw]] >= 0
#      yhat[[iw]] <- yhat[[iw]] * 2 - 1
#    }
    rownames(yhat[[iw]]) <- rownames(new_Z)
  }
  return(yhat)
}
