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
#' minimize_{beta, beta0, gamma} max(1 - y(beta0 1_n * Z_clr beta), 0)^2 +
#'                               lambda_max * frac || W * gamma ||_1
#' subject to beta = A gamma, 1_p^T beta = 0
#' where W = diag(w) with w_u > 0 for all u
#'
#' Classification Huber:
#' minimize_{beta, beta0, gamma} l_p(y * (beta0 +Z_clr^T* beta)) +
#'                               lambda_max * frac || W * gamma ||_1
#' subject to beta = A gamma, 1_p^T beta = 0
#' where W = diag(w) with w_u > 0 for all u
#'
#'
#' Observe that the tuning parameter is specified through "frac", the fraction
#' of lamda_max (which is the smallest value for which gamma is nonzero).
#'
#' @param Z n by p matrix containing log(X)
#' @param y n vector (response)
#' @param A p by (t_size-1) binary matrix giving tree structure (t_size is the
#'   total number of nodes and the -1 is because we do not include the root)
#' @param X n by p' matrix containing metadata
#' @param fraclist (optional) vector of tuning parameter multipliers.  Or a list
#'   of length num_w of such vectors. Should be in (0, 1].
#' @param nlam number of tuning parameters (ignored if fraclist non-NULL)
#' @param min_frac smallest value of tuning parameter multiplier (ignored if
#'   fraclist non-NULL)
#' @param w vector of positive weights of length t_size - 1 (default: all equal
#'   to 1). Or a list of num_w such vectors.
#' @param method string which estimation method to use should be in
#'   ("regression", "classification", "classification_huber")
#' @param intercept_classif boolean indicating if the intercept should be
#'   included for classification
#' @param normalized logical: indicate if metadata should be normalized.
#'   In this case the calculation for each covariate / feature:
#'   (X-X_mean) / ||X||_2
#' @param rho_classification value for huberized classification loss.
#'   Default = -0.0
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
trac <- function(Z, y, A, X = NULL, fraclist = NULL, eps = 1e-3, nlam = 20,
                 min_frac = 1e-4, w = NULL, method = "regression",
                 intercept_classif = TRUE, normalized = TRUE,
                 rho_classification = 0.0) {
  n <- length(y)
  stopifnot(nrow(Z) == n)
  p <- ncol(Z)
  t_size <- ncol(A) + 1


  if (!is.null(X)) {
    # basic check of input metadata X
    stopifnot(nrow(X) == n)
    if (!is.data.frame(X)) X <- data.frame(X)
    if (any(is.na(X))) {
      stop(paste(
        "missing data is currently not supported.",
        "There seems to be missing values in the metadata (X)",
        sep = " "
      ))
    }
    p_x <- ncol(X)
    # add the covariates from X to tree, the variables should not be rooted to
    # the same tree as the compositional data
    # create A as:
    # (A 0)
    # (0 1)
    # where 1 is the diagonal matrix with the metadata's covariates dimensions
    A_rows <- matrix(data = rep(0, time = p_x * (t_size - 1)), nrow = p_x)
    A_rows <- Matrix::Matrix(A_rows, sparse = TRUE)
    rownames(A_rows) <- colnames(X)
    # add names for later
    A <- rbind(A, A_rows)
    A_column <- matrix(data = rep(0, time = p_x * p), ncol = p_x)
    A_column <- Matrix::Matrix(A_column, sparse = TRUE)
    A_diagonal <- Matrix::Diagonal(p_x, x = 1)
    A_column <- rbind(A_column, A_diagonal)
    A <- cbind(A, A_column)
    t_size <- ncol(A) + 1
  }


  if (is.null(w)) w <- list(rep(1, t_size - 1))
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
  # define supported methods .. maybe add partial matching?
  supported_methods <- c("regression", "classification", "classification_huber")
  if (!(method %in% supported_methods)) {
    stop(paste("trac currently supports the following methods: ",
      paste(supported_methods, collapse = ", "),
      sep = " "
    ))
  }
  if (length(unique(y)) == 2 & method == "regression") {
    warning("this looks like a classification task, check the method argument")
  }
  # create a dummy variable with information about weather it is a classifcation
  # task or not for later
  if (method %in% c("classification", "classification_huber")) {
    classification <- TRUE
  } else {
    classification <- FALSE
  }
  if (length(unique(y)) != 2 & classification) {
    stop("the response variable should be binary for the classification task")
  }
  if (!all(unique(y) %in% c(-1, 1)) & classification) {
    warning(paste(
      "values of y must be -1 and 1 for classification.",
      "The fitted model is based on a transformation of y to",
      "(-1,1). The first value of y is coded as 1."
    ))
    y <- y == y[0]
    y <- y * 2 - 1
  }
  # If not a classification task --> intercept_classification is TRUE for
  # prediction
  if (!classification) intercept_classif <- TRUE

  stopifnot(rho_classification < 1)
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

  if (!is.null(fraclist)) {
    if (num_w == 1 & !is.list(fraclist)) fraclist <- list(fraclist)
    stopifnot(unlist(fraclist) >= 0 & unlist(fraclist) <= 1)
  }

  if (is.null(fraclist)) {
    fraclist <- lapply(1:num_w,
                       function(x) exp(seq(0, log(min_frac), length = nlam)))
  }

  if (classification) {
    yt <- y
    Zbar <- rowMeans(Z)
    Z_clr <- Z - Zbar
    M <- as.matrix(Z_clr %*% A)
  } else {
    Zbar <- rowMeans(Z)
    Z_clr <- Z - Zbar
    ybar <- mean(y)
    yt <- y - ybar

    Z_clrA <- as.matrix(Z_clr %*% A)
    v <- Matrix::colMeans(Z_clrA)
    M <- Matrix::t(Matrix::t(Z_clrA) - v)
  }

  if (!is.null(X)) {
    if (normalized) {
      classes_x <- sapply(X, class)
      factors <- c("factor")
      categorical <- classes_x %in% factors
      n_categorical <- sum(categorical)
      n_numeric <- p_x - n_categorical
      if (n_categorical == 1) {
        # apply transformation to categorical variable
        # smth like
        # first value as reference
        # X[, categorical] <- X[, categorical] == X[, categorical][0]
        # lives on the same scale as normalized numerical variables ranging
        # mainly in (-1, 1)
        # X[, categorical] <- X[, categorical] * 2 - 1
      }
      if (n_categorical > 1) {
        # apply transformation to categorical variables
      }
      if (n_numeric == 1) {
        if (!intercept_classif) {
          xm <- 0
        } else {
          xm <- mean(X[, !categorical])
        }
        xs <- sqrt(mean(X[, !categorical]^2) - xm^2)
        X[, !categorical] <- (X[, !categorical] - xm) / xs
      }
      if (n_numeric > 1) {
        if (!intercept_classif) {
          xm <- rep(0, length = nrow(X[, !categorical]))
        } else {
          xm <- apply(X[, !categorical], 2, function(x) mean(x))
        }
        xs <- apply(X[, !categorical], 2, function(x) norm(x, type = "2"))
        X[, !categorical] <- t((t(X[, !categorical]) - xm) / xs)
      }
    }
    M <- cbind(M, X)
  }

  fit <- list()
  for (iw in seq(num_w)) {
    # for each w...
    #  C is 1_p^T A W^-1
    C <- matrix(Matrix::colSums(A %*% diag(1 / w[[iw]])), 1, t_size - 1)
    if (!is.null(X)) {
      # set C to 0 for not compositional data
      C[length(C) - ((p_x - 1):0)] <- rep(0, p_x)
    }
    X_classo <- M %*% diag(1 / w[[iw]])

    # set up CLASSO problem:
    prob <- classo$classo_problem(
      X = X_classo,
      C = C,
      y = array(yt))
    if (classification) {
      prob$formulation$classification <- TRUE
    } else {
      prob$formulation$classification <- FALSE
    }
    prob$formulation$concomitant <- FALSE
    prob$formulation$huber <- FALSE
    prob$model_selection$PATH <- TRUE
    prob$model_selection$CV <- FALSE
    prob$model_selection$LAMfixed <- FALSE
    prob$model_selection$StabSel <- FALSE
    prob$model_selection$PATHparameters$lambdas <- fraclist[[iw]]
    if (classification & intercept_classif) {
      prob$formulation$intercept <- TRUE
    }
    if (method == "classification_huber") {
      prob$formulation$huber <- TRUE
      prob$formulation$rho_classification <- rho_classification
    }

    # solve  it
    prob$solve()

    # extract outputs
    delta <- as.matrix(prob$solution$PATH$BETAS)
    if (classification & intercept_classif) {
      # c-lasso can estimate beta0 --> select first column (estimated beta0)
      # delete the first column afterwards
      beta0 <- delta[, 1]
      delta <- delta[, -1]
    }
    delta <- t(delta)
    # gammahat = W^-1 deltahat and betahat = A gammahat
    gamma <- diag(1 / w[[iw]]) %*% delta
    beta <- A %*% gamma
    # rescale betas for numerical values
    if (!is.null(X) & normalized) {
      beta[length(C) - ((p_x - 1):0)][, !categorical] <-
        beta[length(C) - ((p_x - 1):0)][, !categorical] * xs + xm
    }
    # alphahat = diag(1^T A) %*% gammahat:
    nleaves <- Matrix::colSums(A)
    alpha <- nleaves * gamma
    lambda_classo <- prob$model_selection$PATHparameters$lambdas
    if (!classification) beta0 <- ybar - crossprod(gamma, v)
    rownames(beta) <- rownames(A)
    rownames(gamma) <- rownames(alpha) <- colnames(A)
    fit[[iw]] <- list(
      beta0 = beta0,
      beta = beta,
      gamma = gamma,
      alpha = alpha,
      fraclist = lambda_classo, # / (2 * n),
      w = w[[iw]],
      fit_classo = prob,
      refit = FALSE,
      method = method,
      intercept_classif = intercept_classif,
      rho_classification = rho_classification
    )
  }
  fit
}

