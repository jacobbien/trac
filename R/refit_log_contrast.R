#' Refit subject to sparsity constraints
#'
#' Given output of \code{\link{sparse_log_contrast}}, solves the least squares problem with
#' compositional constraint on the features selected by sparse_log_contrast.
#'
#' minimize_{beta, beta0} 1/(2n) || y - beta0 1_n - Z beta ||^2
#' subject to beta_{nonselected} = 0, 1_p^T beta = 0
#'
#' @param fit output of sparse_log_contrast
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



#' Refit log-contrast for classification to sparsity constraint
#'
#' Given output of \code{\link{sparse_log_contrast}}, solves the classification
#' problem with compositional constraint on the features selected by
#' sparse_log_contrast. In contrast to \code{\link{refit_sparse_log_contrast}}
#' this function only does the refit for one model specified by i_selected
#' or component_selected.
#' This i_selected usually comes from the i1se from the cross-validation output.
#'
#'
#' @param fit output of sparse_log_contrast
#' @param i_selected indicator which lambda is selected based on for example the
#'   cross-validation procedure
#' @param Z,y,additional_covariates same arguments as passed to
#'   \code{\link{sparse_log_contrast}}
#' @param tol tolerance for deciding whether a beta value is zero
#' @param component_selected vector with indices which component to include
#' @export
refit_sparse_log_contrast_classif <- function(fit, i_selected = NULL, Z, y,
                                              additional_covariates = NULL,
                                              tol = 1e-5,
                                              component_selected = NULL) {
  # check which components have non-zero coefficients
  if (is.null(component_selected)) {
    selected_variables <- abs(fit$beta[, i_selected]) > tol
    if (sum(selected_variables) < 1) stop("There are no selected components")
  } else {
    p_test <- ifelse(is.null(additional_covariates),  ncol(Z),
           (ncol(Z) + ncol(additional_covariates)))
    if (max(component_selected) > p_test) {
      stop("The values of component_selected must be smaller than ncol(Z)
           + ncol(additional_covariates")
    }
    if (min(component_selected) < 1) {
      stop("The values of component_selected must be greater than 1")
    }
    selected_variables <- component_selected
  }
  # for covariates check which of them need to satisfy the zero sum constraint
  C <- fit$C[, selected_variables]
  C <- matrix(C, nrow = 1)
  # extract meta information from the fit object
  rho <- fit$rho
  normalized <- fit$normalized
  intercept <- fit$intercept

  # Transformation of the additional covariates
  if (!is.null(additional_covariates)) {
    # transform additional covariates to data.frame --> easier to work with
    # different types of input
    if (!is.data.frame(additional_covariates)) {
      additional_covariates <- data.frame(additional_covariates)
    }
    # define the number of additional covariates
    p_x <- ncol(additional_covariates)
  }

  # normalize the non-compositional data if wanted
  if (!is.null(additional_covariates)) {
    if (normalized) {
      # call the normalization helper function
      normalized_values <-
        normalization_additional_covariates(additional_covariates =
                                              additional_covariates,
                                            p_x = p_x,
                                            intercept = intercept)
      additional_covariates <- normalized_values$X
    } else {
      # get the number of categorical variables if no normalization is applied
      categorical_list <- get_categorical_variables(additional_covariates)
      categorical <- categorical_list[["categorical"]]
      n_categorical <- categorical_list[["n_categorical"]]
      if (n_categorical > 0) {
        additional_covariates[, categorical] <-
          transform_categorical_variables(additional_covariates, categorical)
      }
    }
  }


  # Concenate the non-compositional data if available
  if (!is.null(additional_covariates)) {
    # add the non-compositional covariates
    X_classo <- as.matrix(cbind(Z, additional_covariates))
  } else {
    X_classo <- as.matrix(Z)
  }
  #
  X_classo <- X_classo[, selected_variables]

  n <- nrow(X_classo)
  p <- ncol(X_classo)

  # set up CLASSO problem:
  prob <- classo$classo_problem(X = X_classo,
                                C = C,
                                y = array(y))
  prob$formulation$classification <- TRUE
  prob$formulation$concomitant <- FALSE
  if (intercept) {
    prob$formulation$intercept <- TRUE
  }
  if (fit$method == "classif_huber") {
    prob$formulation$huber <- TRUE
    prob$formulation$rho_classification <- rho
  } else {
    prob$formulation$huber <- FALSE
  }
  prob$model_selection$PATH <- FALSE
  prob$model_selection$CV <- FALSE
  prob$model_selection$StabSel <- FALSE
  prob$model_selection$LAMfixed <- TRUE
  prob$model_selection$LAMfixedparameters$rescaled_lam <- TRUE
  prob$model_selection$LAMfixedparameters$lam <- 0.0



  # solve it
  prob$solve()
  # extract outputs
  beta <- (prob$solution$LAMfixed$beta)
  if (intercept) {
    # c-lasso can estimate beta0 --> select first column (estimated beta0)
    # delete the first column afterwards
    beta0 <- beta[1]
    beta <- beta[-1]
  }
  beta <- t(beta)

  list(beta0 = beta0,
       beta = beta,
       C = C,
       refit = TRUE,
       method = fit$method,
       intercept = intercept,
       rho = rho,
       normalized = normalized)
}
