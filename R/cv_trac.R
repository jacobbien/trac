#' Perform cross validation for tuning parameter selection
#'
#' This function is to be called after calling \code{\link{trac}}.  It performs
#' \code{nfold}-fold cross validation.
#'
#' @param fit output of \code{\link{trac}} function.
#' @param Z,y,A,X same arguments as passed to \code{\link{trac}}
#' @param folds a partition of \code{1:nrow(Z)}.
#' @param nfolds number of folds for cross-validation
#' @param summary_function how to combine the errors calculated on each
#' observation within a fold (e.g. mean or median)
#' @export
cv_trac <- function(fit, Z, y, A, X = NULL, folds = NULL, nfolds = 5,
                    summary_function = stats::median) {
  n <- nrow(Z)
  p <- ncol(Z)
  stopifnot(length(y) == n)
  if (is.null(folds)) {
    folds <- ggb:::make_folds(n, nfolds)
  } else {
    nfolds <- length(folds)
  }
  cv <- list()
  fit_folds <- list() # save this to reuse by log-ratio's cv function
  for (iw in seq_along(fit)) {
    if (length(fit) > 1) cat("CV for weight sequence #", iw, fill = TRUE)
    errs <- matrix(NA, ncol(fit[[iw]]$beta), nfolds)
    for (i in seq(nfolds)) {
      cat("fold", i, fill = TRUE)
      if(is.null(fit[[iw]]$method)) fit[[iw]]$method <- "regression"
      # train on all but i-th fold (and use settings from fit):
      fit_folds[[i]] <- trac(Z[-folds[[i]], ],
        y[-folds[[i]]],
        A,
        X[-folds[[i]]],
        fraclist = fit[[iw]]$fraclist,
        w = fit[[iw]]$w,
        method = fit[[iw]]$method
      )
      if (fit[[iw]]$refit) {
        fit_folds[[i]] <- refit_trac(
          fit_folds[[i]], Z[-folds[[i]], ],
          y[-folds[[i]]], A
        )
      }
      if (fit[[iw]]$method == "regression" | is.null(fit[[iw]]$method)) {
        errs[, i] <- apply(
          (predict_trac(
            fit_folds[[i]],
            Z[folds[[i]], ],
            X[folds[[i]], ])[[1]] - y[folds[[i]]])^2,
          2, summary_function
        )
      }
      if (fit[[iw]]$method == "classification") {
        # loss: max(0, 1 - y_hat * y)^2
        r <- 1 - predict_trac(
          fit_folds[[i]], Z[folds[[i]], ], X[folds[[i]], ])[[1]] *
          c(y[folds[[i]]])
        r_max <- r > 0
        r <- (r * r_max)^2
        errs[, i] <- apply(r, 2, summary_function)
      }
      if (fit[[iw]]$method == "classification_huber") {
        rho_classification <- fit[[iw]]$rho_classification
        r <- predict_trac(
          fit_folds[[i]], Z[folds[[i]], ], X[folds[[i]], ])[[1]] *
          c(y[folds[[i]]])
        r_case_3 <- !(r >= 1)
        r_case_2 <- r <= rho_classification
        r_case_1 <- (r > rho_classification) & (r < 1)

        r <- r_case_3 * r
        r[r_case_2] <- (1 - rho_classification) *
          (1 + rho_classification - 2 * r[r_case_2])
        r[r_case_1] <- (1 - r[r_case_1])^2

        errs[, i] <- apply(r, 2, summary_function)
      }
    }
    m <- rowMeans(errs)
    se <- apply(errs, 1, stats::sd) / sqrt(nfolds)
    ibest <- which.min(m)
    i1se <- min(which(m < m[ibest] + se[ibest]))
    cv[[iw]] <- list(
      errs = errs, m = m, se = se,
      lambda_best = fit[[iw]]$fraclist[ibest], ibest = ibest,
      lambda_1se = fit[[iw]]$fraclist[i1se], i1se = i1se,
      fraclist = fit[[iw]]$fraclist, w = fit[[iw]]$w,
      nonzeros = colSums(abs(fit[[iw]]$gamma) > 1e-5),
      fit_folds = fit_folds
    )
  }
  list(
    cv = cv,
    iw_best = which.min(lapply(cv, function(cvv) cvv$m[cvv$ibest])),
    iw_1se = which.min(lapply(cv, function(cvv) cvv$m[cvv$i1se])),
    folds = folds
  )
}
