#' Make predictions based on a sparse log contrast fit
#'
#' @param fit output of the function \code{\link{sparse_log_contrast}}
#' @param new_Z a new data matrix (see \code{Z} from
#'     \code{\link{sparse_log_contrast}})
#' @param new_additional_covariates a new data matrix
#'    (see \code{additional_covariates}
#'    from \code{\link{sparse_log_contrast}})
#' @param output string  either "raw", "probability" or "class" only relevant
#'   classification tasks
#' @return a vector of \code{nrow(new_Z)} predictions.
#' @export
predict_sparse_log_contrast <- function(fit, new_Z,
                                        new_additional_covariates = NULL,
                                        output = c("raw", "probability",
                                                   "class")) {
  # fit: output of wag
  # new_Z: n_new by p matrix
  # new_additional_covariates: n_new by p' matrix



  # add to make make code backwards compatible, regression
  # and with intercept as default
  if (is.null(fit$method)) {
    fit$method <- "regr"
    fit$intercept <- TRUE
  }
  # check which method to use
  output <- match.arg(output)
  classification <- fit$method %in% c("classif",
                                      "classif_huber")
  intercept <- fit$intercept

  # Transform additional covariates

  if (!is.null(new_additional_covariates)) {
    categorical_list <-
      get_categorical_variables(new_additional_covariates)
    categorical <- categorical_list[["categorical"]]
    n_categorical <- categorical_list[["n_categorical"]]
    if (n_categorical > 0) {
      new_additional_covariates[, categorical] <-
        transform_categorical_variables(new_additional_covariates,
                                        categorical)
    }
    new_Z <- cbind(as.matrix(new_Z), new_additional_covariates)
    new_Z <- as.matrix(new_Z)
  }

  # save output

  if (intercept) {
    yhat <- t(as.numeric(fit$beta0) +
                t(as.matrix(new_Z %*% fit$beta)))
  } else {
    yhat <- t(t(as.matrix(new_Z %*% fit$beta)))
  }
  # classification: transform output if not raw score
  if (classification) {
    if (output == "class") {
      yhat <- yhat >= 0
      yhat <- yhat * 2 - 1
    }
    if (output == "probability") {
      yhat <- probability_transform(
        yhat = yhat,
        A = fit$hyper_prob[1, ],
        B = fit$hyper_prob[2, ]
      )

    }
  }

  rownames(yhat) <- rownames(new_Z)

  yhat
}
