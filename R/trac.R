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
#' @param X n by p' matrix containing additional covariates
#' @param fraclist (optional) vector of tuning parameter multipliers.  Or a list
#'   of length num_w of such vectors. Should be in (0, 1].
#' @param nlam number of tuning parameters (ignored if fraclist non-NULL)
#' @param min_frac smallest value of tuning parameter multiplier (ignored if
#'   fraclist non-NULL)
#' @param w vector of positive weights of length t_size - 1 (default: all equal
#'   to 1). Or a list of num_w such vectors.
#' @param w_meta vector of positive weights of length ncol(X) (default:
#'   all equal to 1).
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
trac <- function(Z, y, A, X = NULL, fraclist = NULL, nlam = 20,
                 min_frac = 1e-4, w = NULL, w_meta = NULL,
                 method = c("regr", "classif", "classif_huber"),
                 intercept = TRUE, normalized = TRUE,
                 rho = 0.0,
                 output = c("raw", "probability")) {

  n <- length(y)
  stopifnot(nrow(Z) == n)
  p <- ncol(Z)
  t_size <- ncol(A) + 1

  output <- match.arg(output)
  method <- match.arg(method)

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
    if (is.null(w_meta)) w_meta <- rep(1, p_x)
    if (!is.vector(w_meta)) stop("w_meta must be a matrix or a vector.")
    if (length(w_meta) != p_x) stop("w_meta must be of length ncol(X)")
    A <- A_add_X(X = X, A = A, p = p, p_x = p_x)
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

  if (!is.null(fraclist)) {
    if (num_w == 1 & !is.list(fraclist)) fraclist <- list(fraclist)
    stopifnot(unlist(fraclist) >= 0 & unlist(fraclist) <= 1)
  }

  if (is.null(fraclist)) {
    fraclist <- lapply(1:num_w,
                       function(x) exp(seq(0, log(min_frac), length = nlam)))
  }

  if (!is.null(X)) {
    if (normalized) {
      normalized_values <-
        normalization_x(X = X, p_x = p_x, intercept = intercept)
      X <- normalized_values$X
    } else {
      categorical_list <- get_categorical_variables(X)
      categorical <- categorical_list[["categorical"]]
      n_categorical <- categorical_list[["n_categorical"]]
      if (n_categorical > 0) {
        X[, categorical] <- sapply(X[, categorical], as.numeric)
        X[, categorical] <- sapply(X[, categorical], function(x) x - 1)
      }
    }
    w_not_meta <- rep(1, (t_size - 1 - p_x))
    w_x <- c(w_not_meta, w_meta)
  }
  if (classification) {
    yt <- y
    Zbar <- Matrix::rowMeans(Z)
    Z_clr <- Z - Zbar
    if (!is.null(X)) Z_clr <- as.matrix(cbind(Z_clr, X))
    M <- as.matrix(Z_clr %*% A)
  } else {
    Zbar <- Matrix::rowMeans(Z)
    Z_clr <- Z - Zbar
    ybar <- mean(y)
    yt <- y - ybar

    if (!is.null(X)) Z_clr <- as.matrix(cbind(Z_clr, X))
    Z_clrA <- as.matrix(Z_clr %*% A)
    v <- Matrix::colMeans(Z_clrA)
    if (!is.null(X)) {
      M <- Z_clrA
      M[, 1:(t_size - 1 - ncol(X))] <-
        Matrix::t(Matrix::t(Z_clrA[, 1:(t_size - 1 - ncol(X))]) -
                    v[1:(t_size - 1 - ncol(X))])
    } else {
    M <- Matrix::t(Matrix::t(Z_clrA) - v)
    }
  }

  if (!classification) intercept <- TRUE
  fit <- list()
  for (iw in seq(num_w)) {
    # for each w...
    #  C is 1_p^T A W^-1
    C <- matrix(Matrix::colSums(A %*% diag(1 / w[[iw]])), 1, t_size - 1)
    if (!is.null(X)) {
      # set C = 0 for non compositional data
      C[length(C) - ((p_x - 1):0)] <- rep(0, p_x)
    }
    X_classo <- M %*% diag(1 / w[[iw]])

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
    if (!is.null(X)) prob$formulation$w <- w_x
    # solve  it
    prob$solve()
    # extract outputs
    delta <- as.matrix(prob$solution$PATH$BETAS)
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
    gamma <- diag(1 / w[[iw]]) %*% delta
    beta <- A %*% gamma
    # rescale betas for numerical values
    if ((!is.null(X)) & normalized) {
      if(normalized_values$n_numeric != 0) {
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
    }
    # alphahat = diag(1^T A) %*% gammahat:
    nleaves <- Matrix::colSums(A)
    alpha <- nleaves * gamma
    lambda_classo <- prob$model_selection$PATHparameters$lambdas
    if (!classification) beta0 <- ybar - crossprod(gamma, v)
    if (!intercept) beta0 <- rep(0, times = length(lambda_classo))
    rownames(beta) <- rownames(A)
    rownames(gamma) <- rownames(alpha) <- colnames(A)
    if (output == "probability") {
      eps <- 1e-3
      if (!is.null(A)) A <- A[1:p, 1:(ncol(A) - p_x)]
      hyper_prob <-
        get_probability_cv(
          Z = Z,
          X = X,
          A = A,
          y = y,
          method = method,
          w = w[[iw]],
          w_meta = w_meta,
          fraclist = lambda_classo,
          nfolds = 3,
          eps = eps,
          n_lambda = length(lambda_classo)
          )
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
      w_meta = w_meta,
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


A_add_X <- function(X, A, p, p_x) {
  # add the covariates from X to tree, the variables should not be rooted to
  # the same tree as the compositional data
  # create A as:
  # (A_old 0)
  # (0     1)
  # where 1 is the diagonal matrix with the metadata's covariates dimensions
  # X is the dataframe with the additional covariates and A is the tree
  # for the compositional data
  # return: updated A
  A_diagonal <- Matrix::Diagonal(p_x, x = 1)
  A_rownames <- rownames(A)
  A_colnames <- colnames(A)
  X_rownames <- colnames(X)
  A <- Matrix::bdiag(A, A_diagonal)
  # add rownames
  if (!is.null(A_rownames)) {
    if (!is.null(X_rownames)) {
      if (!is.null(A_colnames)) {
        dimnames(A) <- list(c(A_rownames, X_rownames),
                            c(A_colnames, X_rownames))
      } else {
        rownames(A) <- c(A_rownames, X_rownames)
      }
    } else {
      dimnames(A) <- list(c(A_rownames, rep("", times = p_x)),
                          c(A_colnames, rep("", times = p_x)))
    }
  } else {
    if (!is.null(X_rownames)) {
      rownames(A) <- c(rep("", times = p), X_rownames)
    }
  }
  A
}


check_method <- function(method, y, rho = 0.0) {
  # check the inputs for classification tasks
  supported_methods <- c("regr", "classif", "classif_huber")
  if (!(method %in% supported_methods)) {
    stop(paste("trac currently supports the following methods: ",
               paste(supported_methods, collapse = ", "),
               sep = " "
    ))
  }
  if (length(unique(y)) == 2 & method == "regr") {
    warning("this looks like a classification task, check the method argument")
  }
  # create a dummy variable with information about weather it is a
  # classification task or not for later
  if (method %in% c("classif", "classif_huber")) {
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
    y <- y == y[1]
    y <- y * 2 - 1
  }
  # If not a classification task --> intercept is TRUE for
  # prediction
  # only accept rho for huberized loss smaller 1
  stopifnot(rho < 1)
  # return classification indicator and transformed y
  list(classification = classification,
       y = y)
}

get_categorical_variables <- function(X) {
  categorical <- sapply(X, is.factor)
  n_categorical <- sum(categorical)
  list(categorical = categorical,
       n_categorical = n_categorical)
}

normalization_x <- function(X, p_x, intercept) {
  # normalize metadata. (x-mean(x)) / norm(x)
  categorical_list <- get_categorical_variables(X)
  categorical <- categorical_list[["categorical"]]
  n_categorical <- categorical_list[["n_categorical"]]
  n_numeric <- p_x - n_categorical
  if (n_numeric == 1) {
    if (!intercept) {
      xm <- 0
    } else {
      xm <- mean(X[, !categorical])
    }
    xs <- norm(X[, !categorical], type = "2")
    X[, !categorical] <- (X[, !categorical] - xm) / xs
  }
  if (n_numeric > 1) {
    if (!intercept) {
      xm <- rep(0, length = ncol(X[, !categorical]))
    } else {
      xm <- apply(X[, !categorical], 2, function(x) mean(x))
    }
    xs <- apply(X[, !categorical], 2, function(x) norm(x, type = "2"))
    X[, !categorical] <- t((t(X[, !categorical]) - xm) / xs)
  }
  if(n_numeric == 0) {
    xm <- NULL
    xs <- NULL
  }
  if (n_categorical > 0) {
    X[, categorical] <- sapply(X[, categorical], as.numeric)
    X[, categorical] <- sapply(X[, categorical], function(x) x - 1)
  }
  list(categorical = categorical,
       n_numeric = n_numeric,
       xm = xm,
       xs = xs,
       X = X)
}

rescale_betas <- function(beta, p_x, p, n_numeric, categorical, xs, xm) {
  # rescales the betas to the original scale
  if (n_numeric > 0) {
    if ((p_x == 1) & (n_numeric == 1)) {
      normalized_zero <- beta[(p + p_x), ] == 0
      beta[(p + p_x), ] <-
        beta[(p + p_x), ] * xs + xm
      beta[(p + p_x), ][normalized_zero] <- 0
    } else if (p_x != 1) {
      normalized_zero <-
        beta[(p + p_x) - ((p_x - 1):0), ][!categorical, ] == 0
      beta[(p + p_x) - ((p_x - 1):0), ][!categorical, ] <-
        beta[(p + p_x) - ((p_x - 1):0), ][!categorical, ] * xs + xm
      beta[(p + p_x) - ((p_x - 1):0), ][!categorical, ][normalized_zero] <- 0
    }
  }
  beta
}


get_probability_cv <- function(Z, X, A, y, method, w, w_meta, fraclist,
                               nfolds = 3, eps, n_lambda) {
  # calculate the hyperparameter for platts method. Do three fold cv
  # to obtain the decision values and pass to the platt algorithm for each
  # lambda
  n <- nrow(Z)
  if (!is.null(X)) {
    if (!is.data.frame(X)) X <- data.frame(X)
  }
  folds <- ggb:::make_folds(n, nfolds)
  decision_values <- matrix(ncol = n_lambda)
  label <- c()
  for (i in seq(nfolds)) {
    fit_folds <- trac(Z[-folds[[i]], ],
                           y[-folds[[i]]],
                           A,
                           X[-folds[[i]], ],
                           fraclist = fraclist,
                           w = w,
                      w_meta = w_meta,
                           method = method,
                           output = "raw")
    decision_value <- predict_trac(fit_folds, Z[folds[[i]], ],
                                   X[folds[[i]], ])[[1]]
    decision_values <- rbind(decision_values, decision_value)
    label <- c(label, c(y[folds[[i]]]))
  }
  decision_values <- decision_values[-1, ]
  hyper_prob <- matrix(nrow = 2, ncol = n_lambda)
  for (i in seq_len(ncol(decision_values))) {
    hyper_tmp <- get_probability_platt(decision_values = decision_values[, i],
                                       label = label, eps = eps)
    hyper_prob[1, i] <- hyper_tmp$A
    hyper_prob[2, i] <- hyper_tmp$B
  }
  hyper_prob
}


get_probability_platt <- function(decision_values, label, eps) {
  # This algorithm is based on the pseudo-code of
  # Lin, H. T., Lin, C. J., & Weng, R. C. (2007). A note on Platt’s
  # probabilistic outputs for support vector machines. Machine learning,
  # 68(3), 267-276.
  # Appendix 3 Pseudo code of Algorithm 1 of the paper

  # Parameter setting
  max_iter <- 100
  min_step <- 1e-10
  sigma <- 1e-12

  n_label <- length(label)
  label <- label > 0
  prior1 <- sum(label)
  prior0 <- n_label - prior1

  hi_target <- (prior1 + 1) / (prior1 + 2)
  lo_target <- 1 / (prior0 + 2)
  t_target <- c()
  t_target[label] <- hi_target
  t_target[!label] <- lo_target

  A <- 0
  B <- log((prior0 + 1) / (prior1 + 1))

  fApB <- decision_values * A + B
  fApB_index <- fApB >= 0
  fval <- sum(t_target[fApB_index] * fApB[fApB_index] +
                log(1 + exp(-fApB[fApB_index])))
  fval <- fval + sum((t_target[!fApB_index] - 1) * fApB[!fApB_index] +
                       log(1 + exp(fApB[!fApB_index])))

  p <- c()
  q <- c()


  for (i in 1:max_iter) {
    h_11 <- sigma
    h_22 <- sigma
    h_21 <- 0
    g_1 <- 0
    g_2 <- 0

    fApB <- decision_values * A + B
    fApB_index <- fApB >= 0

    p[fApB_index] <- exp(-fApB[fApB_index]) / (1 + exp(-fApB[fApB_index]))
    q[fApB_index] <- 1 / (1 + exp(-fApB[fApB_index]))

    p[!fApB_index] <- 1 / (1 + exp(fApB[!fApB_index]))
    q[!fApB_index] <- exp(fApB[!fApB_index]) / (1 + exp(fApB[!fApB_index]))

    d_2 <- p * q
    h_11 <- h_11 + sum(decision_values * decision_values * d_2)
    h_22 <- h_22 + sum(d_2)
    h_21 <- h_21 + sum(decision_values * d_2)

    d_1 <- t_target - p
    g_1 <- g_1 + sum(decision_values * d_1)
    g_2 <- g_2 + sum(d_1)

    if (abs(g_1) < eps && abs(g_2) < eps) {
      break
    }

    deter <- h_11 * h_22 - h_21^2
    d_a <- - (h_22 * g_1 - h_21 * g_2) / deter
    d_b <- - (- h_21 * g_1 + h_11 * g_2) / deter
    gd <- g_1 * d_a + g_2 * d_b
    step_size <- 1
    while (step_size >= min_step) {
      new_A <- A + step_size * d_a
      new_B <- B + step_size * d_b
      new_f <- 0

      fApB <- decision_values * new_A + new_B
      fApB_index <- fApB >= 0


      new_f <- sum(t_target[fApB_index] * fApB[fApB_index] +
                   log(1 + exp(-fApB[fApB_index])))
      new_f <- new_f + sum((t_target[!fApB_index] - 1) * fApB[!fApB_index] +
                             log(1 + exp(fApB[!fApB_index])))

      if (new_f < fval + 0.0001 * step_size * gd) {
        A <- new_A
        B <- new_B
        fval <- new_f
        break
      } else {
        step_size <- step_size / 2
      }
      if (step_size < min_step) {
        warning("The optimisation step (to find the probabilities) did
                not converge")
        break
      }
    }
  }
  list(A = A,
       B = B)
}
