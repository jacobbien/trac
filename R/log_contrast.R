#' Perform sparse log-contrast regression
#'
#' Solves the constrained lasso problem using the CLASSO module
#' in Python.  The optimization problem is
#'
#' minimize_{beta, beta0} 1/(2n) || y - beta0 1_n - Zagg_clr beta ||^2
#'                                     + lamda_max * frac || beta ||_1
#' subject to C beta = 0
#'
#' Default is C = 1_p^T, but C can be a general matrix.
#'
#' Observe that the tuning parameter is specified through "frac", the fraction
#' of lamda_max (which is the smallest value for which beta is nonzero).
#'
#' @param Z n by p matrix containing log(X)
#' @param y n vector (response)
#' @param C m by p matrix. Default is a row vector of ones.
#' @param fraclist (optional) vector of tuning parameter multipliers.  Should be in (0, 1].
#' @param nlam number of tuning parameters (ignored if fraclist non-NULL)
#' @param min_frac smallest value of tuning parameter multiplier (ignored if
#'   fraclist non-NULL)
#'
#' @export
sparse_log_contrast <- function(Z, y, C = NULL, fraclist = NULL, nlam = 20, min_frac = 1e-4) {
  n <- length(y)
  stopifnot(nrow(Z) == n)
  p <- ncol(Z)

  if (is.null(C)) C <- matrix(1, nrow = 1, ncol = p)

  # CLASSO can solve a problem of the form
  #
  # (based on https://github.com/Leo-Simpson/c-lasso)
  # [R1] Standard constrained Lasso regression:
  # minimize_{beta} || y - X beta ||^2 + lam ||beta||_1
  # subject to C beta = 0
  #
  # Define
  # yt = y - ybar 1_n
  # M = Z_clr - 1_n v^T, where v = colMeans(Z_clr)
  #
  if(is.null(fraclist))
    fraclist <- exp(seq(0, log(min_frac), length = nlam))

  Zbar <- Matrix::rowMeans(Z)
  Z_clr <- Z - Zbar
  ybar <- mean(y)
  yt <- y - ybar

  v <- Matrix::colMeans(Z_clr)
  M <- Matrix::t(Matrix::t(Z_clr) - v)

  fit <- list()
  X_classo <- as.matrix(M)

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
  prob$model_selection$PATHparameters$lambdas <- fraclist

  # solve  it
  prob$solve()

  # extract outputs
  beta <- purrr::map(prob$solution$PATH$BETAS, as.numeric) %>%
    rlang::set_names(paste0("V", 1:length(prob$solution$PATH$BETAS))) %>%
    dplyr::bind_cols() %>%
    as.matrix()
  lambda_classo <- prob$model_selection$PATHparameters$lambdas
  beta0 <- ybar - crossprod(beta, v)
  rownames(beta) <- colnames(Z)
  list(beta0 = beta0,
       beta = beta,
       C = C,
       fraclist = lambda_classo, # / (2 * n),
       fit_classo = prob,
       refit = FALSE)
}
