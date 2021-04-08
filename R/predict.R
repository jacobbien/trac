#' Make predictions based on a trac fit
#'
#' @param fit output of the function \code{\link{trac}}
#' @param new_Z a new data matrix (see \code{Z} from \code{\link{trac}})
#' @param new_X a new data matrix (see \code{X} from \code{\link{trac}})
#' @param output string  either "raw", "probability" or "class" only relevant
#'   classification tasks
#' @return a vector of \code{nrow(new_Z) + nrow(new_X)} predictions.
#' @export
predict_trac <- function(fit, new_Z, new_X = NULL,
                         output = c("raw", "probability", "class")) {
  # fit: output of wag
  # new_Z: n_new by p matrix
  # add to make make code backwards compatible
  if(is.null(fit[[1]]$method)) {
    fit[[1]]$method <- "regr"
    fit[[1]]$intercept <- TRUE
  }

  output <- match.arg(output)
  classification <- fit[[1]]$method %in% c("classif",
                                           "classif_huber")
  intercept <- fit[[1]]$intercept

  yhat <- list()
  for (iw in seq_along(fit)) {
    if (!is.null(new_X)) {
      categorical_list <- get_categorical_variables(new_X)
      categorical <- categorical_list[["categorical"]]
      n_categorical <- categorical_list[["n_categorical"]]
      if (n_categorical > 0) {
        new_X[, categorical] <- sapply(new_X[, categorical], as.numeric)
        new_X[, categorical] <- sapply(new_X[, categorical], function(x) x - 1)
      }
      new_Z <- cbind(new_Z, new_X)
      new_Z <- as.matrix(new_Z)
    }
    if (intercept) {
      yhat[[iw]] <- t(as.numeric(fit[[iw]]$beta0) +
        t(as.matrix(new_Z %*% fit[[iw]]$beta)))
    } else {
      yhat[[iw]] <- t(t(as.matrix(new_Z %*% fit[[iw]]$beta)))
    }
    if(classification) {
      if (output == "class") {
        yhat[[iw]] <- yhat[[iw]] >= 0
        yhat[[iw]] <- yhat[[iw]] * 2 - 1
      }
      if (output == "probability") {
        yhat[[iw]] <- probability_transform(yhat = yhat[[iw]],
                                            A = fit[[iw]]$hyper_prob[1, ],
                                            B = fit[[iw]]$hyper_prob[2, ])

      }
    }

    rownames(yhat[[iw]]) <- rownames(new_Z)
  }
  return(yhat)
}


probability_transform <- function(yhat, A, B) {
  # following the idea of Chapter 3.2 of
  # Lin, H. T., Lin, C. J., & Weng, R. C. (2007). A note on Platt’s
  # probabilistic outputs for support vector machines.
  # Machine learning, 68(3), 267-276.

  fApB <- t(A * t(yhat) + B)
  fApB_index <- fApB >= 0
  prob <- matrix(nrow = nrow(yhat), ncol = ncol(yhat))
  prob[fApB_index] <- exp(-fApB[fApB_index]) / (1 + exp(-fApB[fApB_index]))
  prob[!fApB_index] <- 1 / (1 + exp(fApB[!fApB_index]))
  prob
}
