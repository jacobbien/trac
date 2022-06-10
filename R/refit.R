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

#' #' Refit trac for classification to sparsity constraint
#' #'
#' #' Given output of \code{\link{trac}}, solves the classification
#' #' problem with compositional constraint on the features selected by
#' #' trac. In contrast to \code{\link{refit_trac}}
#' #' this function only does the refit for one model specified by i_selected or
#' #' component_selected.
#' #' This i_selected usually comes from the i1se from the cross-validation output.
#' #'
#' #'
#' #' @param fit output of sparse_log_contrast
#' #' @param i_selected indicator which lambda is selected based on for example the
#' #'   cross-validation procedure
#' #' @param Z,y,additional_covariates same arguments as passed to
#' #'   \code{\link{sparse_log_contrast}}
#' #' @param tol tolerance for deciding whether a beta value is zero
#' #' @param component_selected vector with indices which column of A to include
#' #' @export
#' refit_trac_classif <- function(fit, i_selected = NULL, Z, y, A,
#'                                additional_covariates = NULL,
#'                                tol = 1e-5, component_selected = NULL) {
#'   n <- nrow(Z)
#'   p <- ncol(Z)
#'   t_size <- ncol(A) + 1
#'   browser()
#'   rho <- fit[[1]]$rho
#'   intercept <- fit[[1]]$intercept
#'   normalized <- fit[[1]]$normalized
#'   method <- fit[[1]]$method
#'   num_w <- length(fit)
#'   if (!is.null(additional_covariates)) {
#'     # transform additional covariates to data.frame --> easier to work with
#'     # different types of input
#'     if (!is.data.frame(additional_covariates)) {
#'       additional_covariates <- data.frame(additional_covariates)
#'     }
#'     # define the number of additional covariates
#'     p_x <- ncol(additional_covariates)
#'   }
#'
#'   # normalize the non-compositional data if wanted
#'   if (!is.null(additional_covariates)) {
#'     if (normalized) {
#'       # call the normalization helper function
#'       normalized_values <-
#'         normalization_additional_covariates(additional_covariates =
#'                                               additional_covariates,
#'                                             p_x = p_x,
#'                                             intercept = intercept)
#'       additional_covariates <- normalized_values$X
#'     } else {
#'       # get the number of categorical variables if no normalization is applied
#'       categorical_list <- get_categorical_variables(additional_covariates)
#'       categorical <- categorical_list[["categorical"]]
#'       n_categorical <- categorical_list[["n_categorical"]]
#'       if (n_categorical > 0) {
#'         additional_covariates[, categorical] <-
#'           transform_categorical_variables(additional_covariates, categorical)
#'       }
#'     }
#'   }
#'   yt <- y
#'   # clr transformation on Z
#'   Zbar <- Matrix::rowMeans(Z)
#'   Z_clr <- Z - Zbar
#'   # add the additional covariates
#'   Z_clrA <- as.matrix(Z_clr %*% A)
#'
#'   # define number of nodes and leafs under the node in order to
#'   # calculate the geom mean for the compositional data only
#'   v <- Matrix::colMeans(Z_clrA)
#'   M <- Matrix::t(Matrix::t(Z_clrA) - v)
#'
#'
#'   refit <- list()
#'   for (iw in seq(num_w)) {
#'
#'     if (is.null(component_selected)) {
#'       selected_variables <- abs(fit[[iw]]$gamma[, i_selected]) > tol
#'       if (sum(selected_variables) < 1) stop("There are no selected components")
#'     } else {
#'       p_test <- ifelse(is.null(additional_covariates),  ncol(A),
#'                        (ncol(A) + ncol(additional_covariates)))
#'       if (max(component_selected) > p_test) {
#'         stop("The values of component_selected must be smaller than ncol(A)
#'            + ncol(additional_covariates")
#'       }
#'       if (min(component_selected) < 1) {
#'         stop("The values of component_selected must be greater than 1")
#'       }
#'       selected_variables <- component_selected
#'     }
#'     w <- fit[[iw]]$w
#'     C <- matrix(Matrix::colSums(A %*% diag(1 / w)), 1, t_size - 1)
#'     X_classo <- M %*% diag(1 / w)
#'     if (!is.null(additional_covariates)) {
#'       # add the non-compositional covariates
#'       X_classo <- as.matrix(cbind(X_classo, additional_covariates))
#'       # set C = 0 for non compositional data, since
#'       # no zero sum constraint necessary
#'       C <- cbind(C, matrix(rep(0, p_x), nrow = 1))
#'     }
#'     p_X_classo <- ncol(X_classo)
#'     X_classo <- X_classo[, selected_variables]
#'     C <- C[, selected_variables]
#'     C <- as.matrix(C, nrow = 1)
#'     # set up CLASSO problem:
#'     prob <- classo$classo_problem(
#'       X = X_classo,
#'       C = C,
#'       y = array(yt))
#'
#'
#'     prob$formulation$classification <- TRUE
#'     prob$formulation$concomitant <- FALSE
#'     if (intercept) {
#'       prob$formulation$intercept <- TRUE
#'     }
#'     if (method == "classif_huber") {
#'       prob$formulation$huber <- TRUE
#'       prob$formulation$rho_classification <- rho
#'     } else {
#'       prob$formulation$huber <- FALSE
#'     }
#'     prob$model_selection$PATH <- FALSE
#'     prob$model_selection$CV <- FALSE
#'     prob$model_selection$StabSel <- FALSE
#'     prob$model_selection$LAMfixed <- TRUE
#'     prob$model_selection$LAMfixedparameters$rescaled_lam <- TRUE
#'     prob$model_selection$LAMfixedparameters$lam <- 0.0
#'
#'     # solve  it
#'     prob$solve()
#'     # extract outputs
#'
#'     delta <- (prob$solution$LAMfixed$beta)
#'     # Warn the user if the solver could not solve the problem
#'     if (sum(is.nan(delta)) > 0) {
#'       warning("There is a problem with the estimation,
#'               try a different method or without interception")
#'       delta[is.nan(delta)] <- 0
#'     }
#'     if (intercept) {
#'       # c-lasso can estimate beta0 --> select first column (estimated beta0)
#'       # delete the first column afterwards
#'       beta0 <- delta[1]
#'       delta <- delta[-1]
#'     }
#'     delta <- t(delta)
#'     delta_tmp <- rep(0, p_X_classo)
#'     delta_tmp[selected_variables] <- delta
#'     delta <- delta_tmp
#'     # delta <- as.numeric(delta)
#'     # gammahat = W^-1 deltahat and betahat = A gammahat
#'     gamma <- 1 / w * delta[1:(t_size - 1)]
#'     gamma <- as.matrix(gamma)
#'
#'     beta <- A %*% gamma
#'     # alphahat = diag(1^T A) %*% gammahat:
#'     nleaves <- Matrix::colSums(A)
#'     alpha <- nleaves * gamma
#'     lambda_classo <- fit[[iw]]$fraclist
#'     if (!intercept) beta0 <- 0
#'     rownames(beta) <- rownames(A)
#'     if (!is.null(additional_covariates)) {
#'       beta <- rbind(beta, delta[(t_size):(t_size + p_x - 1), ])
#'       alpha <- rbind(alpha, delta[(t_size):(t_size + p_x - 1), ])
#'       gamma <- rbind(gamma, delta[(t_size):(t_size + p_x - 1), ])
#'       if (!is.null(rownames(A)) & !is.null(colnames(A)) &
#'           !is.null( colnames(additional_covariates))) {
#'         rownames(beta) <- c(rownames(A), colnames(additional_covariates))
#'         rownames(gamma) <- rownames(alpha) <- c(colnames(A),
#'                                                 colnames(additional_covariates))
#'       }
#'       if (normalized && (normalized_values$n_numeric != 0)) {
#'         # rescale betas for numerical values
#'         # rescale only if beta not 0
#'         beta <- rescale_betas(
#'           beta = beta,
#'           p_x = p_x,
#'           p = p,
#'           n_numeric = normalized_values$n_numeric,
#'           categorical = normalized_values$categorical,
#'           xs = normalized_values$xs,
#'           xm = normalized_values$xm
#'         )
#'       }
#'     } else {
#'       rownames(gamma) <- rownames(alpha) <- colnames(A)
#'     }
#'     refit[[iw]] <- list(
#'       beta0 = beta0,
#'       beta = beta,
#'       gamma = gamma,
#'       alpha = alpha,
#'       fraclist = lambda_classo, # / (2 * n),
#'       w = w,
#'       w_additional_covariates = fit[[1]]$w_additional_covariates,
#'       fit_classo = prob,
#'       refit = TRUE,
#'       method = method,
#'       intercept = intercept,
#'       rho = rho,
#'       normalized = normalized
#'     )
#'   }
#'   refit
#' }
#'
