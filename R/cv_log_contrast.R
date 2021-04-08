#' Perform cross validation for tuning parameter selection for sparse log contrast
#'
#' This function is to be called after calling \code{\link{sparse_log_contrast}}.
#' It performs \code{nfold}-fold cross validation.
#'
#' @param fit output of \code{\link{sparse_log_contrast}} function.
#' @param Z,y same arguments as passed to \code{\link{sparse_log_contrast}}. C
#' that is used will be taken from fit object.
#' @param folds a partition of \code{1:nrow(Z)}.
#' @param nfolds number of folds for cross-validation
#' @param summary_function how to combine the errors calculated on each
#' observation within a fold (e.g. mean or median)
#' @export
cv_sparse_log_contrast <- function(fit, Z, y, folds = NULL, nfolds = 5, summary_function = stats::median) {
  n <- nrow(Z)
  p <- ncol(Z)
  stopifnot(length(y) == n)
  if(is.null(folds)) folds <- ggb:::make_folds(n, nfolds)
  else
    nfolds <- length(folds)
  cv <- list()
  fit_folds <- list() # save this to reuse by log-ratio's cv function
  errs <- matrix(NA, ncol(fit$beta), nfolds)
  for (i in seq(nfolds)) {
    cat("fold", i, fill = TRUE)
    # train on all but i-th fold (and use settings from fit):
    fit_folds[[i]] <- sparse_log_contrast(Z[-folds[[i]], ],
                                          y[-folds[[i]]],
                                          fit$C,
                                          fraclist = fit$fraclist)
    if (fit$refit) stop("Not yet supported.")
    errs[, i] <- apply((predict_trac(list(fit_folds[[i]]), Z[folds[[i]], ])[[1]] - y[folds[[i]]])^2,
                       2, summary_function)
  }
  m <- rowMeans(errs)
  se <- apply(errs, 1, stats::sd) / sqrt(nfolds)
  ibest <- which.min(m)
  i1se <- min(which(m < m[ibest] + se[ibest]))
  cv <- list(errs = errs, m = m, se = se,
             lambda_best = fit$fraclist[ibest], ibest = ibest,
             lambda_1se = fit$fraclist[i1se], i1se = i1se,
             fraclist = fit$fraclist,
             nonzeros = colSums(abs(fit$beta) > 1e-5),
             fit_folds = fit_folds)
  list(cv = cv, folds = folds)
}
