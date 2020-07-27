#' Perform tree-based aggregation
#'
#' Solves the weighted aggregation problem using the CLASSO module
#' in Python.  The optimization problem is
#'
#' minimize_{beta, beta0, gamma} 1/(2n) || y - beta0 1_n - Z_clr beta ||^2
#'                                     + lam || W * gamma ||_1
#' subject to beta = A gamma, 1_p^T beta = 0
#' where W = diag(w) with w_u > 0 for all u
#'
#' @param Z n by p matrix containing log(X)
#' @param y n vector (response)
#' @param A p by t_size binary matrix giving tree structure (t_size is the total
#'    number of nodes)
#' @param lamlist (optional) vector of tuning parameters.  Or a list of length num_w
#'          of such vectors.
#' @param nlam number of tuning parameters (ignored if lamlist non-NULL)
#' @param flmin ratio of smallest to largest tuning parameter (ignored if lamlist
#'        non-NULL)
#' @param w vector of positive weights of length t_size (default: all equal to 1). Or
#'    a list of num_w such vectors.
#'
#' @return a list of length num_w, where each list element corresponds to the solution
#' for that choice of w.  Note that the lamlist depends on the choice of w.
#' @export
trac <- function(Z, y, A, lamlist = NULL, eps = 1e-3, nlam = 20, flmin = 1e-4, w = NULL) {
  n <- length(y)
  stopifnot(nrow(Z) == n)
  p <- ncol(Z)
  t_size <- ncol(A)
  if (is.null(w)) w <- list(rep(1, t_size))
  if (is.numeric(w)) w <- list(w)
  if (!is.list(w)) stop("w must be a list.")
  if (any(lapply(w, length) != t_size))
    stop("every element of w must be of length t_size.")
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

  if(!is.null(lamlist)) stopifnot(lamlist >= 0)

  Zbar <- rowMeans(Z)
  Z_clr <- Z - Zbar
  ybar <- mean(y)
  yt <- y - ybar

  Z_clrA <- as.matrix(Z_clr %*% A)
  v <- colMeans(Z_clrA)
  M <- t(t(Z_clrA) - v)

  fit <- list()
  for (iw in seq(num_w)) {
    # for each w...
    #  C is 1_p^T A W^-1
    C <- matrix(Matrix::colSums(A %*% diag(1 / w[[iw]])), 1, t_size)
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
    prob$model_selection$PATHparameters$lamin <- flmin

    # solve  it
    prob$solve()

    # extract outputs
    delta <- as.matrix(map_dfc(prob$solution$PATH$BETAS, as.numeric))
    # gammahat = W^-1 deltahat and betahat = A gammahat
    gamma <- diag(1 / w[[iw]]) %*% delta
    beta <- A %*% gamma

    lambda_classo <- prob$model_selection$PATHparameters$lambdas
    beta0 <- ybar - crossprod(gamma, v)
    rownames(beta) <- rownames(A)
    rownames(gamma) <- colnames(A)
    fit[[iw]] <- list(beta0 = beta0,
                      beta = beta,
                      gamma = gamma,
                      lamlist = lambda_classo / (2 * n),
                      w = w[[iw]],
                      fit_classo = prob,
                      refit = FALSE)
  }
  if (!is.null(lamlist)) {
    # user specified some lambda values of interest.
    # Here we do linear interpolation to go from the path to the
    # solutions at specifically the lambda values in lamlist:
    fit2 <- list()
    for (iw in seq(num_w)) {
      fit2[[iw]] <- fit[[iw]]
      lam_knots <- fit[[iw]]$lamlist
      lin_interp <- function(b) approx(lam_knots,
                                       b,
                                       lamlist,
                                       rule = 1:2)$y
      fit2[[iw]]$beta <- t(apply(fit[[iw]]$beta, 1, lin_interp))
      fit2[[iw]]$gamma <- t(apply(fit[[iw]]$gamma, 1, lin_interp))
      fit2[[iw]]$beta0 <- lin_interp(fit[[iw]]$beta0)
      fit2[[iw]]$lamlist <- lamlist
    }
    return(fit2)
  }
  fit
}


