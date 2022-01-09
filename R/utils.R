# Transform non-compositional input --------------------------------------------

#' Functions transform non-compositional input
#'
#' normalize metadata. (x-mean(x)) / norm(x) for numerical variables and
#' 0/1 for categorical variables provided as a factor
#'
#' @keywords internal
#' @param additional_covariates new data matrix
#'    (see \code{additional_covariates} from \code{\link{trac}})
#' @param p_x number of additional covariates
#' @param intercept Boolean if normalization should be applied
#' @return list with vector indicating if categorical or not, number of
#'    categorical variables, vector of means and l2 norm and normalized
#'    variables
#'

normalization_additional_covariates <- function(additional_covariates,
                                                p_x, intercept) {
  # check which variables are categorical and which are not categorical
  categorical_list <- get_categorical_variables(additional_covariates)
  categorical <- categorical_list[["categorical"]]
  n_categorical <- categorical_list[["n_categorical"]]
  n_numeric <- p_x - n_categorical
  # calculate the mean, mean 0 for models without intercept
  if (!intercept) {
    xm <- rep(0, length = ncol(additional_covariates[, !categorical,
                                                     drop = FALSE]))
  } else {
    xm <- apply(additional_covariates[, !categorical, drop = FALSE], 2,
                function(x) mean(x))
  }
  # calculate l2 norm
  xs <- apply(additional_covariates[, !categorical, drop = FALSE], 2,
              function(x) norm(x, type = "2"))
  # (x-mean(x)) / norm(x)
  additional_covariates[, !categorical] <-
    t((t(additional_covariates[, !categorical, drop = FALSE]) - xm) / xs)
  if (n_numeric == 0) {
    xm <- NULL
    xs <- NULL
  }
  # one hot encoding for categorical input
  if (n_categorical > 0) {
    additional_covariates[, categorical] <-
      sapply(additional_covariates[, categorical, drop = FALSE], as.numeric)
    additional_covariates[, categorical] <-
      sapply(additional_covariates[, categorical, drop = FALSE],
             function(x) x - 1)
  }
  # return Boolean if categorical , number of numerical inputs, mean vector,
  # l2 norm vector, transformed additional covariates
  list(categorical = categorical,
       n_numeric = n_numeric,
       xm = xm,
       xs = xs,
       X = additional_covariates)
}

#' Check if additional variables are categorical
#'
#' Check if the additional non-compositional covariates are categorical or not
#' based on the type of the column. If the column is a binary factor then it is
#' assumed to be a categorical variable. Useful to transform the input.
#'
#' @keywords internal
#' @param X see \code{additional_covariates} from \code{\link{trac}}
#' @return list with vector indicating if categorical or not, number of
#'    categorical variables
#'

get_categorical_variables <- function(X) {
# fct to get categorical covariates
categorical <- sapply(X, is.factor)
n_categorical <- sum(categorical)
list(categorical = categorical,
     n_categorical = n_categorical)
}

#' Transform non-compositional categorical input
#'
#' Transform factors to numerical values of 0 and 1
#' @keywords internal
#' @param X see \code{additional_covariates} from \code{\link{trac}}
#' @param categorical vector indicating if categorical or not
#' @return data frame with transformed categorical variables from factor to
#'    numerical
#'

transform_categorical_variables <- function(X, categorical) {
  X[, categorical] <- sapply(X[, categorical], as.numeric)
  sapply(X[, categorical], function(x) x - 1)
}

#' Re-scale betas of normalized variables
#'
#' Re-scales betas of the normalized covariates to the original scale for
#' non-zero betas
#' @keywords internal
#' @param beta estimated coefficients
#' @param p_x count of additional non-compositional covariates
#' @param p count of compositional covariates
#' @param n_numeric count of numerical non-compositional covariates
#' @param categorical boolean vector indicating if the additional
#'    non-compositional covariates are categorical or not (see
#'    \code{get_categorical_variables})
#' @param xs vector l2 norm (see\code{normalization_additional_covariates})
#' @param xm vector of means (see\code{normalization_additional_covariates})
#' @return matrix of re-scaled betas
#' @export
#'

rescale_betas <- function(beta, p_x, p, n_numeric, categorical, xs, xm) {
  # re-scales the betas to the original scale
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

#' Check non-compositional inputs
#'
#' Check if the additional non-compositional inputs have NAs and have the same
#' number of observations as the compositional inputs
#' @keywords internal
#' @param additional_covariates new data matrix
#'    (see \code{additional_covariates} from \code{\link{trac}})
#' @param n number of observations
#' @param w_additional_covariates weights for the estimation of the coefficients
#' @param p_x vector with number of additional non-compositional covariates
#' @return errors if the requirements are not met
#' @export
#'


check_additional_covariates <-
  function(additional_covariates, n, w_additional_covariates, p_x) {
    # basic check input additional covariates
    stopifnot(nrow(additional_covariates) == n)
    # No missing data allowed
    if (any(is.na(additional_covariates))) {
      stop(paste(
        "missing data is currently not supported.",
        "There seems to be missing values in the non compositional covariates ",
        "(additional_covariates).",
        sep = " "
      ))
    }
    # Weight vector needs to be the same length as the number of covariates
    if (!is.vector(w_additional_covariates)) {
      stop("w_additional_covariates must be a matrix or a vector.")
    }
    if (length(w_additional_covariates) != p_x) {
      stop("w_additional_covariates must be
           of length ncol(additional_covariates)")
    }
}


# Check Input for classification -----------------------------------------------

#' Check method input and other hyperparameter
#'
#' Check if the method input and other hyperparameter for classification
#' are correctly specified
#'
#' @keywords internal
#' @param method The method (see \code{method}
#'     from \code{\link{trac}})
#' @param y The outcome (see \code{y}
#'     from \code{\link{trac}})
#' @param rho The hyperparameter for huberized loss for classification
#'     (see \code{rho}from \code{\link{trac}})
#' @return list with vector indicating which method is used and the outcome
#' @export
#'

check_method <- function(method, y, rho = 0.0) {
  # check the inputs for classification tasks
  supported_methods <- c("regr", "classif", "classif_huber")
  if (!(method %in% supported_methods)) {
    stop(paste("trac currently supports the following methods: ",
               paste(supported_methods, collapse = ", "),
               sep = " "
    ))
  }
  # If the outcome is binary warn the user
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
  # The response variable should be binary for classification
  if (length(unique(y)) != 2 & classification) {
    stop("the response variable should be binary for the classification task")
  }
  # Transform the output for classification to 1, -1 if not done yet
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


# Get probabilistic output of raw score ----------------------------------------

#' 5-fold cross-validation to find the optimal parameters A and B
#' for probabilistic output. Obtain the decision values and pass to the platt
#' algorithm for each lambda
#'
#' @keywords internal
#' @param Z same arguments as passed to \code{\link{trac}}
#' @param additional_covariates same arguments as passed to \code{\link{trac}}
#' @param A same arguments as passed to \code{\link{trac}}
#' @param y same arguments as passed to \code{\link{trac}}
#' @param method same arguments as passed to \code{\link{trac}}
#' @param w same arguments as passed to \code{\link{trac}}
#' @param w_additional_covariates same arguments as passed to \code{\link{trac}}
#' @param fraclist same arguments as passed to \code{\link{trac}}
#' @param nfolds number of folds
#' @param eps convergence criterion
#' @param n_lambda number of lambdas
#' @return matrix with hyperparameter A and B as columns
#' @export
#'

get_probability_cv <- function(Z, additional_covariates, A, y, method, w,
                               w_additional_covariates, fraclist,
                               nfolds = 5, eps, n_lambda) {
  # check additional covariates and number of observations
  n <- nrow(Z)
  if (!is.null(additional_covariates)) {
    if (!is.data.frame(additional_covariates)) {
      additional_covariates <- data.frame(additional_covariates)
    }
  }
  # create folds
  folds <- ggb:::make_folds(n, nfolds)
  # pre-create matrix with decision values
  decision_values <- matrix(ncol = n_lambda)
  label <- c()
  # run cv
  for (i in seq(nfolds)) {
    # fit model
    fit_folds <- trac(Z[-folds[[i]], ],
                      y[-folds[[i]]],
                      A,
                      additional_covariates[-folds[[i]], ],
                      fraclist = fraclist,
                      w = w,
                      w_additional_covariates = w_additional_covariates,
                      method = method,
                      output = "raw")
    # run prediction
    decision_value <- predict_trac(fit_folds, Z[folds[[i]], ],
                                   additional_covariates[folds[[i]], ])[[1]]
    # add the prediction value to the matrix
    decision_values <- rbind(decision_values, decision_value)
    # extract labels
    label <- c(label, c(y[folds[[i]]]))
  }
  decision_values <- decision_values[-1, ]
  hyper_prob <- matrix(nrow = 2, ncol = n_lambda)
  for (i in seq_len(ncol(decision_values))) {
    # run platt probability algorithm to get A and B
    hyper_tmp <- get_probability_platt(decision_values = decision_values[, i],
                                       label = label, eps = eps)
    hyper_prob[1, i] <- hyper_tmp$A
    hyper_prob[2, i] <- hyper_tmp$B
  }
  hyper_prob
}

#' Calculate parameters A and B based on Platt
#'
#' This algorithm is based on the pseudo-code of
#  Lin, H. T., Lin, C. J., & Weng, R. C. (2007). A note on Platt’s
#  probabilistic outputs for support vector machines. Machine learning,
#  68(3), 267-276.
#  Appendix 3 Pseudo code of Algorithm 1 of the paper
#'
#' @keywords internal
#' @param decision_values output score of trac for classification
#'    (see \code{y} from \code{\link{trac}})
#' @param label labels see \code{y} from \code{\link{trac}}
#' @param eps threshold for convergence
#' @return A and B
#' @export
#'


get_probability_platt <- function(decision_values, label, eps) {

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
    d_a <- -(h_22 * g_1 - h_21 * g_2) / deter
    d_b <- -(-h_21 * g_1 + h_11 * g_2) / deter
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

#' Transform the score output to probabilistic output
#'
#' following the idea of Chapter 3.2 of
#  Lin, H. T., Lin, C. J., & Weng, R. C. (2007). A note on Platt’s
#  probabilistic outputs for support vector machines.
#  Machine learning, 68(3), 267-276.
#'
#' @keywords internal
#' @param yhat score based on the output of trac
#'    (see \code{y} from \code{\link{trac}})
#' @param A Hyperparameter for probabilistic transformation
#' @param B Hyperparameter for probabilistic transformation
#' @return Probabilistic output
#' @export
#'

probability_transform <- function(yhat, A, B) {


  fApB <- t(A * t(yhat) + B)
  fApB_index <- fApB >= 0
  prob <- matrix(nrow = nrow(yhat), ncol = ncol(yhat))
  prob[fApB_index] <- exp(-fApB[fApB_index]) / (1 + exp(-fApB[fApB_index]))
  prob[!fApB_index] <- 1 / (1 + exp(fApB[!fApB_index]))
  prob
}
