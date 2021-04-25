#' Perform tree-based aggregation
#'
#' Solves the weighted aggregation problem using the CLASSO module
#' in Python.  The optimization problem is
#'
#' minimize_{beta, beta0, gamma} 1/(2n) || y - beta0 1_n - Z_clr beta ||^2
#'                                     + lamda_max * frac || W * gamma ||_1
#' subject to beta = A gamma, 1_p^T beta = 0
#' where W = diag(w) with w_u > 0 for all u
#'
#' Observe that the tuning parameter is specified through "frac", the fraction
#' of lamda_max (which is the smallest value for which gamma is nonzero).
#'
#' @param Z n by p matrix containing log(X)
#' @param y n vector (response)
#' @param A p by (t_size-1) binary matrix giving tree structure (t_size is the
#'   total number of nodes and the -1 is because we do not include the root)
#' @param fraclist (optional) vector of tuning parameter multipliers.  Or a list
#'   of length num_w of such vectors. Should be in (0, 1].
#' @param nlam number of tuning parameters (ignored if fraclist non-NULL)
#' @param min_frac smallest value of tuning parameter multiplier (ignored if
#'   fraclist non-NULL)
#' @param w vector of positive weights of length t_size - 1 (default: all equal
#'   to 1). Or a list of num_w such vectors.
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
trac <- function(Z, y, A, fraclist = NULL, nlam = 20, min_frac = 1e-4, w = NULL) {
  n <- length(y)
  stopifnot(nrow(Z) == n)
  p <- ncol(Z)
  t_size <- ncol(A) + 1
  if (is.null(w)) w <- list(rep(1, t_size - 1))
  if (is.numeric(w)) w <- list(w)
  if (!is.list(w)) stop("w must be a list.")
  if (any(lapply(w, length) != t_size - 1))
    stop("every element of w must be of length (t_size - 1).")
  if (any(unlist(lapply(w, function(ww) any(ww <= 0)))))
    stop("every element of w must be positive.") # for simplicity, for now we require this.
  num_w <- length(w)


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

  if(!is.null(fraclist)) {
    if (num_w == 1 & !is.list(fraclist)) fraclist <- list(fraclist)
    stopifnot(unlist(fraclist) >= 0 & unlist(fraclist) <= 1)
  }

  if(is.null(fraclist))
    fraclist <- lapply(1:num_w,
                      function(x) exp(seq(0, log(min_frac), length = nlam)))

  Zbar <- Matrix::rowMeans(Z)
  Z_clr <- Z - Zbar
  ybar <- mean(y)
  yt <- y - ybar

  Z_clrA <- as.matrix(Z_clr %*% A)
  v <- Matrix::colMeans(Z_clrA)
  M <- Matrix::t(Matrix::t(Z_clrA) - v)

  fit <- list()
  for (iw in seq(num_w)) {
    # for each w...
    #  C is 1_p^T A W^-1
    C <- matrix(Matrix::colSums(A %*% diag(1 / w[[iw]])), 1, t_size - 1)
    X_classo <- M %*% diag(1 / w[[iw]])

    # set up CLASSO problem:
    prob <- classo$classo_problem(X = X_classo,
                                   C = C,
                                   y = array(yt))
    prob$formulation$classification <- FALSE
    prob$formulation$concomitant <- FALSE
    prob$formulation$huber <- FALSE
    prob$model_selection$PATH <- TRUE
    prob$model_selection$CV <- FALSE
    prob$model_selection$LAMfixed <- FALSE
    prob$model_selection$StabSel <- FALSE
    prob$model_selection$PATHparameters$lambdas <- fraclist[[iw]]

    # solve  it
    prob$solve()

    # extract outputs
    delta <- purrr::map(prob$solution$PATH$BETAS, as.numeric) %>%
      rlang::set_names(paste0("V", 1:length(prob$solution$PATH$BETAS))) %>%
      dplyr::bind_cols() %>%
      as.matrix()
    # gammahat = W^-1 deltahat and betahat = A gammahat
    gamma <- diag(1 / w[[iw]]) %*% delta
    beta <- A %*% gamma
    # alphahat = diag(1^T A) %*% gammahat:
    nleaves <- Matrix::colSums(A)
    alpha <- nleaves * gamma
    lambda_classo <- prob$model_selection$PATHparameters$lambdas
    beta0 <- ybar - crossprod(gamma, v)
    rownames(beta) <- rownames(A)
    rownames(gamma) <- rownames(alpha) <- colnames(A)
    fit[[iw]] <- list(beta0 = beta0,
                      beta = beta,
                      gamma = gamma,
                      alpha = alpha,
                      fraclist = lambda_classo, # / (2 * n),
                      w = w[[iw]],
                      fit_classo = prob,
                      refit = FALSE)
  }
  fit
}


