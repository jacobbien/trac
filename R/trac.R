#' Perform tree-based aggregation
#'
#' Solves the weighted aggregation problem using the CLASSO module
#' in Python.  The optimization problems are:
#'
#' Regression:
#' minimize_{beta, beta0, gamma} 1/(2n) || y - beta0 1_n - Z_clr beta ||^2
#'                                     + lamda_max * frac || W * gamma ||_1
#' subject to beta = A gamma, 1_p^T beta = 0
#' where W = diag(w) with w_u > 0 for all u
#'
#' Classification:
#' minimize_{beta, beta0, gamma} max(1 - y_i(beta0 + Z_clr_i * beta), 0)^2 +
#'                               lambda_max * frac || W * gamma ||_1
#' subject to beta = A gamma, 1_p^T beta = 0
#' where W = diag(w) with w_u > 0 for all u
#'
#' Classification Huberized:
#' see C2 formulation in c-lasso
#'
#'
#' Observe that the tuning parameter is specified through "frac", the fraction
#' of lamda_max (which is the smallest value for which gamma is nonzero).
#'
#' @param Z n by p matrix containing log(X)
#' @param y n vector (response)
#' @param A p by (t_size-1) binary matrix giving tree structure (t_size is the
#'   total number of nodes and the -1 is because we do not include the root)
#' @param additional_covariates n by p' matrix containing additional covariates
#'    / features
#' @param fraclist (optional) vector of tuning parameter multipliers.  Or a list
#'   of length num_w of such vectors. Should be in (0, 1].
#' @param nlam number of tuning parameters (ignored if fraclist non-NULL)
#' @param min_frac smallest value of tuning parameter multiplier (ignored if
#'   fraclist non-NULL)
#' @param w vector of positive weights of length t_size - 1 (default: all equal
#'   to 1). Or a list of num_w such vectors.
#' @param w_additional_covariates vector of positive weights of
#'   length ncol(additional_covariates) (default: all equal to 1).
#' @param method string which estimation method to use should be in
#'   ("regr", "classif", "classif_huber")
#' @param intercept only works for classification! Should the intercept be
#'   fitted. Default is TRUE, set to FALSE if the intercept should not be
#'   included
#' @param normalized if `TRUE` normalize the additional covariates.
#'   In this case the calculation for each covariate / feature:
#'   (X-X_mean) / ||X||_2
#'   The weights will be transformed back to the original scale.
#' @param rho value for huberized classification loss.
#'   Default = -0.0.
#' @param output only relevant for classification. String indicating whether
#'   the raw score output or probability for class 1 should be used.
#'   The probability is estimated with Platt’s probibalistic output
#'
#' @return a list of length num_w, where each list element corresponds to the
#'   solution for that choice of w.  Note that the fraclist depends on the
#'   choice of w. beta0 is the intercept; beta is the coefficient vector on the
#'   scale of leaves of the tree; gamma is the coefficient vector on nodes of
#'   the tree where the features are sums of logs of the leaf-features within
#'   that node's subtree; alpha is the coefficient vector on nodes of the tree
#'   where the features are log of the geometric mean of leaf-features within
#'   that node's subtree.
#' @export
trac <- function(Z, y, A, additional_covariates = NULL, fraclist = NULL,
                 nlam = 20, min_frac = 1e-4, w = NULL,
                 w_additional_covariates = NULL,
                 method = c("regr", "classif", "classif_huber"),
                 intercept = TRUE, normalized = TRUE,
                 rho = 0.0,
                 output = c("raw", "probability")) {
  # input check
  n <- length(y)
  stopifnot(nrow(Z) == n)
  p <- ncol(Z)
  t_size <- ncol(A) + 1
  # create the weights and transform to list
  if (is.null(w)) w <- list(rep(1, t_size - 1))
  # if ((!is.null(w_additional_covariates)) & (!is.null(w)) & is.numeric(w) &
  #     (length(w) != t_size - 1)) {
  #   w <- c(w, rep(1, p_x))
  # }
  if (is.numeric(w)) w <- list(w)
  if (!is.list(w)) stop("w must be a list.")
  if (any(lapply(w, length) != t_size - 1)) {
    stop("every element of w must be of length (t_size - 1).")
  }
  if (any(unlist(lapply(w, function(ww) any(ww <= 0))))) {
    # for simplicity, for now we require this.
    stop("every element of w must be positive.")
  }
  num_w <- length(w)
  # partial matching for the output and method
  output <- match.arg(output)
  method <- match.arg(method)

  if (!is.null(additional_covariates)) {
    # transform additional covariates to data.frame --> easier to work with
    # different types of input
    if (!is.data.frame(additional_covariates)) {
      additional_covariates <- data.frame(additional_covariates)
    }
    # define the number of additional covariates
    p_x <- ncol(additional_covariates)
    # create new weights
    if (is.null(w_additional_covariates)) w_additional_covariates <- rep(1, p_x)
    # check input additional_covariates
    check_additional_covariates(additional_covariates, n,
                                w_additional_covariates, p_x)
  }


  # define supported methods and call helper function
  # to check input
  method_check <- check_method(method = method, y = y,
                               rho = rho)
  classification <- method_check$classification
  y <- method_check$y

  # CLASSO can solve a problem of the form
  #
  # (based on https://github.com/muellsen/classo)
  # [R1] Standard constrained Lasso regression:
  # minimize_{beta} || y - X beta ||^2 + lam ||beta||_1
  # subject to C beta = 0
  #
  # Define
  # yt = y - ybar 1_n
  # M = Z_clr * A - 1_n v^T, where v = colMeans(Z_clr * A)
  #
  # The problem we want to solve can be written equivalently as
  # (see classo_wag_derivation.pdf)
  #
  # minimize_{delta} || yt - M W^-1 delta ||^2 + 2*n*lam || delta ||_1
  # subject to 1_p^T A W^-1 delta = 0
  #
  # Given a solution deltahat, we can get
  #   gammahat = W^-1 deltahat and betahat = A gammahat
  # create faclist if not given
  if (!is.null(fraclist)) {
    if (num_w == 1 & !is.list(fraclist)) fraclist <- list(fraclist)
    stopifnot(unlist(fraclist) >= 0 & unlist(fraclist) <= 1)
  }
  if (is.null(fraclist)) {
    fraclist <- lapply(1:num_w,
                       function(x) exp(seq(0, log(min_frac), length = nlam)))
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
    # define the weights as 1 for every non compositional covariates
    # since c-lasso can now handle weights and we do not need to transform
    # them anymore and pass them directly to the solver
    # the weights for the compositional effects is taken into account by
    # modifying C and the weights for the non-compositional effects through
    # w in c-lasso
    w_not_additional_covariates <- rep(1, (t_size - 1))
    w_x <- c(w_not_additional_covariates, w_additional_covariates)
  }
  if (classification) {
    # for classification we do not need to scale the outcome
    yt <- y
  } else {
    # scale y
    ybar <- mean(y)
    yt <- y - ybar
  }
  # clr transformation on Z
  Zbar <- Matrix::rowMeans(Z)
  Z_clr <- Z - Zbar
  # add the additional covariates
  Z_clrA <- as.matrix(Z_clr %*% A)

  # define number of nodes and leafs under the node in order to
  # calculate the geom mean for the compositional data only
  v <- Matrix::colMeans(Z_clrA)
  M <- Matrix::t(Matrix::t(Z_clrA) - v)


  # always use a intercept when not classification due to the nature of the
  # estimation
  if (!classification) intercept <- TRUE
  fit <- list()
  # define the weighted l1 norm penalization
  for (iw in seq(num_w)) {
    # for each w...
    #  C is 1_p^T A W^-1
    C <- matrix(Matrix::colSums(A %*% diag(1 / w[[iw]])), 1, t_size - 1)
    X_classo <- M %*% diag(1 / w[[iw]])
    if (!is.null(additional_covariates)) {
      # add the non-compositional covariates
      X_classo <- as.matrix(cbind(X_classo, additional_covariates))
      # set C = 0 for non compositional data, since
      # no zero sum constraint necessary
      C <- cbind(C, matrix(rep(0, p_x), nrow = 1))
    }
    # set up CLASSO problem:
    prob <- classo$classo_problem(
      X = X_classo,
      C = C,
      y = array(yt))
    prob$formulation$classification <- classification
    prob$formulation$concomitant <- FALSE
    prob$model_selection$PATH <- TRUE
    prob$model_selection$CV <- FALSE
    prob$model_selection$LAMfixed <- FALSE
    prob$model_selection$StabSel <- FALSE
    prob$model_selection$PATHparameters$lambdas <- fraclist[[iw]]
    if (classification & intercept) {
      prob$formulation$intercept <- TRUE
    }
    if (method == "classif_huber") {
      prob$formulation$huber <- TRUE
      prob$formulation$rho_classification <- rho
    } else {
      prob$formulation$huber <- FALSE
    }
    if (!is.null(additional_covariates)) prob$formulation$w <- w_x
    # solve  it
    prob$solve()
    # extract outputs

    delta <- as.matrix(prob$solution$PATH$BETAS)
    # Warn the user if the solver could not solve the problem
    if (sum(is.nan(delta)) > 0) {
      warning("There is a problem with the estimation,
              try a different method or without interception")
      delta[is.nan(delta)] <- 0
    }
    if (classification & intercept) {
      # c-lasso can estimate beta0 --> select first column (estimated beta0)
      # delete the first column afterwards
      beta0 <- delta[, 1]
      delta <- delta[, -1]
    }
    delta <- t(delta)
    # delta <- as.numeric(delta)
    # gammahat = W^-1 deltahat and betahat = A gammahat
    gamma <- diag(1 / w[[iw]]) %*% delta[1:(t_size - 1), ]
    beta <- A %*% gamma
    # alphahat = diag(1^T A) %*% gammahat:
    nleaves <- Matrix::colSums(A)
    alpha <- nleaves * gamma
    lambda_classo <- prob$model_selection$PATHparameters$lambdas
    if (!classification) beta0 <- ybar - crossprod(gamma, v)
    if (!intercept) beta0 <- rep(0, times = length(lambda_classo))
    rownames(beta) <- rownames(A)
    if (!is.null(additional_covariates)) {
      beta <- rbind(beta, delta[(t_size):(t_size + p_x - 1), ])
      alpha <- rbind(alpha, delta[(t_size):(t_size + p_x - 1), ])
      gamma <- rbind(gamma, delta[(t_size):(t_size + p_x - 1), ])
      if (!is.null(rownames(A)) & !is.null(colnames(A)) &
          !is.null( colnames(additional_covariates))) {
      rownames(beta) <- c(rownames(A), colnames(additional_covariates))
      rownames(gamma) <- rownames(alpha) <- c(colnames(A),
                                              colnames(additional_covariates))
      }
      if (normalized && (normalized_values$n_numeric != 0)) {
        # rescale betas for numerical values
        # rescale only if beta not 0
        beta <- rescale_betas(
          beta = beta,
          p_x = p_x,
          p = p,
          n_numeric = normalized_values$n_numeric,
          categorical = normalized_values$categorical,
          xs = normalized_values$xs,
          xm = normalized_values$xm
        )
      }
    } else {
      rownames(gamma) <- rownames(alpha) <- colnames(A)
    }
    if (output == "probability") {
      eps <- 1e-3
      if (!is.null(A) & !is.null(additional_covariates)) {
        A <- A[1:p, 1:(ncol(A) - p_x)]
      }
      hyper_prob <- get_probability_cv(Z = Z,
                                       additional_covariates =
                                         additional_covariates, A = A, y = y,
                                       method = method,
                                       w = w[[iw]],
                                       w_additional_covariates =
                                         w_additional_covariates,
                                       fraclist = lambda_classo, nfolds = 5,
                                       eps = eps,
                                       n_lambda = length(lambda_classo))
    } else {
      hyper_prob <- NULL
    }
    fit[[iw]] <- list(
      beta0 = beta0,
      beta = beta,
      gamma = gamma,
      alpha = alpha,
      fraclist = lambda_classo, # / (2 * n),
      w = w[[iw]],
      w_additional_covariates = w_additional_covariates,
      fit_classo = prob,
      refit = FALSE,
      method = method,
      intercept = intercept,
      rho = rho,
      hyper_prob = hyper_prob,
      normalized = normalized
    )
  }
  fit
}
