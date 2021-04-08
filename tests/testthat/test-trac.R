test_that("Test TRAC functionalities", {
  # set seed for reproducibility
  set.seed(1)

  ## Simple simulated data
  # set number of observations and number of otus/asvs
  n <- 100
  # p <- 10
  p <- 3

  # construct data matrix
  Z <- matrix(rnorm(n = n * p), nrow = n)
  # no negative values allowed
  Z <- abs(Z)
  # take log
  Z <- log(Z)
  # construct a taxonomic tree
  A <- diag(p)
  # A <- cbind(A, c(rep(x = 1, times = p/2), rep(x = 0, times = p/2)))
  # A <- cbind(A, c(rep(x = 1, times = p/5), rep(x = 0, times = 4*p/5)))
  # A <- cbind(A, c(rep(x = 1, times = p)))
  A <- cbind(A, c(1, 1, 0))
  A <- cbind(A, rep(x = 1, times = p))

  # add additional covariates
  X <- data.frame(numeric_feature = rnorm(n = n),
                  categorical_feature = sample(c(0, 1), replace = TRUE,
                                               size = n))

  X$categorical_feature <- as.factor(X$categorical_feature)


  # construct an augmented dataset with the aggregrated levels
  Z_mod <- Z %*% A
  # build y as the combination of column 12 and 10
  # y <- Z_mod[, 12] - Z_mod[, 10] + rnorm(n) #response
  y <- Z_mod[, 4] - Z_mod[, 3] + rnorm(n)
  # for classifcation take the sign --> values in c(-1, 1)
  y_classif <- sign(y)
  # table(y_classif)


  # regression
  fit <- trac(Z, y, A)
  expect_true(all(fit[[1]]$alpha[, 2][c(3,4)] != 0))
  # plot_cv_trac(cvfit_trac = cvfit)


  # classification

  # without intercept
  fit_classif <- trac(Z, y_classif, A, method = "classif", intercept = FALSE)
  expect_true(all(fit_classif[[1]]$alpha[, 2][c(3,4)] != 0))
  expect_true(all((fit_classif[[1]]$beta0 == 0)))



  fit_huber <- trac(Z, y_classif, A, method = "classif_huber",
                    intercept = FALSE)
  expect_true(all(fit_huber[[1]]$alpha[, 2][c(3,4)] != 0))
  expect_true(all(fit_huber[[1]]$beta0 == 0))

  # test nlam parameter
  fit_classif <- trac(Z, y_classif, A, method = "classif", intercept = FALSE,
                      nlam = 10)
  expect_equal(dim(fit_classif[[1]]$beta)[2], 10)

  fit_huber <- trac(Z, y_classif, A, method = "classif_huber",
                    intercept = FALSE,
                    nlam = 10)
  expect_equal(dim(fit_huber[[1]]$beta)[2], 10)


  # with intercept
  y <- 0.5 + Z_mod[, 4] - Z_mod[, 3] + rnorm(n)
  # for classifcation take the sign --> values in c(-1, 1)
  y_classif <- sign(y)
  # table(y_classif)



  fit_classif <- trac(Z, y_classif, A, method = "classif", intercept = TRUE)
  expect_true(all(fit_classif[[1]]$alpha[, 2][c(3,4)] != 0))

  fit_huber <- trac(Z, y_classif, A, method = "classif_huber", intercept = TRUE)
  expect_true(all(fit_huber[[1]]$alpha[, 2][c(3,4)] != 0))

  # try rho

  fit_huber <- trac(Z, y_classif, A, method = "classif_huber", intercept = TRUE,
                    rho = -1)
  expect_true(all(fit_huber[[1]]$alpha[, 2][c(3,4)] != 0))

  # additional covariates

  y <- Z_mod[, 4] - Z_mod[, 3] + X[, 1] + (as.numeric(X[, 2]) - 1) + rnorm(n)
  # for classifcation take the sign --> values in c(-1, 1)
  y_classif <- sign(y)
  # table(y_classif)


  fit_classif <- trac(Z, y_classif, A, X = X, method = "classif",
                      intercept = FALSE)
  expect_true(all(fit_classif[[1]]$alpha[, 10][c(3, 4, 6, 7)] != 0))

  fit_huber <- trac(Z, y_classif, A, X = X, method = "classif_huber",
                    intercept = FALSE)
  expect_true(all(fit_huber[[1]]$alpha[, 10][c(3, 4, 6, 7)] != 0))




  y <- 0.5 + Z_mod[, 4] - Z_mod[, 3] + X[, 1] + (as.numeric(X[, 2]) - 1) +
    rnorm(n)
  # for classifcation take the sign --> values in c(-1, 1)
  y_classif <- sign(y)
  # table(y_classif)


  fit <- trac(Z, y, A, X = X, normalized = TRUE)
  expect_true(all(fit[[1]]$alpha[, 2][c(3,4)] != 0))

  fit_classif <- trac(Z, y_classif, A, X = X, method = "classif",
                      intercept = TRUE,
                      normalized = TRUE)
  expect_true(all(fit_classif[[1]]$alpha[, 2][c(3,4)] != 0))

  fit_huber <- trac(Z, y_classif, A, X = X, method = "classif_huber",
                    intercept = TRUE, normalized = TRUE)
  expect_true(all(fit_huber[[1]]$alpha[, 2][c(3,4)] != 0))





  X$numeric_feature1 <- rnorm(n = n, mean = 5, sd = 5)

  fit_classif <- trac(Z, y_classif, A, X = X, method = "classif",
                      intercept = TRUE,
                      normalized = TRUE, output = "probability")
  expect_true(all(fit_classif[[1]]$alpha[, 2][c(3,4)] != 0))

  fit_classif <- trac(Z, y_classif, A, X = X, method = "classif",
                      intercept = TRUE,
                      normalized = FALSE)
  expect_true(all(fit_classif[[1]]$alpha[, 2][c(3,4)] != 0))

  fit_huber <- trac(Z, y_classif, A = A, X = X, method = "classif_huber",
                    intercept = TRUE, normalized = TRUE, output = "probability")
  expect_true(all(fit_huber[[1]]$alpha[, 2][c(3,4)] != 0))



  p_x <- 3
  x_transform <- trac:::normalization_x(X, p_x = p_x, intercept = TRUE)
  norm_x <- norm(x_transform$X[, 1], type = "2")
  expect_equal(norm_x, 1, tolerance = 5e-3)

  expect_equal(sum(x_transform$categorical == TRUE), 1)
  expect_equal(length(x_transform$xm), 2)
  expect_equal(length(x_transform$xs), 2)


  expect_equal(x_transform$X$numeric_feature * x_transform$xs[1] +
              x_transform$xm[1], X$numeric_feature)

  expect_equal(x_transform$X$numeric_feature1 * x_transform$xs[2] +
              x_transform$xm[2], X$numeric_feature1)


  beta <- c(0, 0, 1, 1, 1,
            0, 0, 0, 0, 0)
  beta <- matrix(beta, ncol = 2, byrow = FALSE)
  p <- 2

  rescale_beta <- trac:::rescale_betas(beta = beta, p_x = p_x, p = p,
                                       n_numeric = x_transform$n_numeric,
                                       categorical = x_transform$categorical,
                                       xs = x_transform$xs,
                                       xm = x_transform$xm)


  expect_equal(rescale_beta[, 2], c(0, 0, 0, 0, 0))
  expect_equal(rescale_beta[3, 1],
               unname(beta[3, 1] * x_transform$xs[1] + x_transform$xm[1]))
  expect_equal(rescale_beta[5, 1],
               unname(beta[5, 1] * x_transform$xs[2] + x_transform$xm[2]))
  expect_equal(dim(rescale_beta), dim(beta))

  expect_error(trac(Z, y, A = A, X = X, method = "classification"))
  expect_error(trac(Z, y, A = A, X = X, method = "classif"))
  expect_error(trac(Z, y_classif, A = A, X = X, method = "classif_huber",
                    rho = 1))
  X_mod_na <- X
  X_mod_na <- X_mod_na[1,1] <- NA
  expect_error(trac(Z, y_classif, A = A, X = X_mod_na, method = "classif"))
  expect_error(trac(Z, y_classif, A = A, X = X, method = "classifc"))
  expect_warning(trac(Z, y_classif + 1, A = A, X = X, method = "classif"))

  y <- 0.5 + Z_mod[, 4] - Z_mod[, 3] + rnorm(n)
  # for classifcation take the sign --> values in c(-1, 1)
  y_classif <- sign(y)
  # table(y_classif)


  fit <- trac(Z, y, A)
  cvfit <- cv_trac(fit, Z = Z, y = y, A = A)
  expect_true(all(fit[[1]]$alpha[, cvfit$cv[[1]]$i1se][c(3,4)] != 0))

  fit_classif <- trac(Z, y_classif, A, method = "classif")
  cvfit_classif <- cv_trac(fit_classif, Z = Z, y = y_classif, A = A)
  expect_true(all(
    fit_classif[[1]]$alpha[, cvfit_classif$cv[[1]]$i1se][c(3,4)] != 0))



  fit_classif <- trac(Z, y_classif, A, method = "classif", intercept = FALSE)
  cvfit_classif <- cv_trac(fit_classif, Z = Z, y = y_classif, A = A)
  expect_true(all(
    fit_classif[[1]]$alpha[, cvfit_classif$cv[[1]]$i1se][c(3,4)] != 0))

  fit_huber <- trac(Z, y_classif, A, method = "classif_huber",
                    intercept = TRUE)
  cvfit_huber <- cv_trac(fit_huber, Z = Z, y = y_classif, A = A)
  expect_true(all(
    fit_huber[[1]]$alpha[, cvfit_huber$cv[[1]]$i1se][c(3,4)] != 0))

  fit_huber <- trac(Z, y_classif, A, method = "classif_huber",
                    intercept = FALSE)
  cvfit_huber <- cv_trac(fit_huber, Z = Z, y = y_classif, A = A)

  expect_true(all(
    fit_huber[[1]]$alpha[, cvfit_huber$cv[[1]]$i1se][c(3,4)] != 0))



  y <- 0.5 + Z_mod[, 4] - Z_mod[, 3]
  # for classifcation take the sign --> values in c(-1, 1)
  y_classif <- sign(y)
  # table(y_classif)

  fit <- trac(Z, y, A, X = X)
  cvfit <- cv_trac(fit, Z = Z, y = y, A = A, X = X)
  expect_true(all(fit[[1]]$alpha[, cvfit$cv[[1]]$i1se][c(3,4)] != 0))


  fit_classif <- trac(Z, y_classif, A, X = X, method = "classif",
                      intercept = TRUE,
                      normalized = TRUE)
  cvfit_classif <- cv_trac(fit_classif, Z = Z, y = y_classif, A = A, X = X)

  expect_true(all(
    fit_classif[[1]]$alpha[, cvfit_classif$cv[[1]]$i1se][c(3,4)] != 0))

  fit_huber <- trac(Z, y_classif, A, X = X, method = "classif_huber",
                    intercept = TRUE, normalized = TRUE, rho = -1)
  cvfit_huber <- cv_trac(fit_huber, Z = Z, y = y_classif, A = A, X = X)

  expect_true(all(
    fit_huber[[1]]$alpha[, cvfit_huber$cv[[1]]$i1se][c(3,4)] != 0))


  fit_classif <- trac(Z, y_classif, A, X = X, method = "classif",
                    normalized = TRUE,
                    output = "probability")
  probability <- predict_trac(fit_classif, new_Z = Z, new_X = X,
                              output = "probability")
  expect_lt(max(probability[[1]][, 10]), 1)
  expect_gt(min(probability[[1]][, 10]), 0)
  class <- predict_trac(fit_classif, new_Z = Z, new_X = X,
                        output = "class")
  expect_equal(length(unique(class[[1]][, 10])), 2)


  fit_classif <- trac(Z, y_classif, A, X = X, method = "classif",
                      normalized = TRUE, intercept = FALSE,
                      output = "probability")
  probability <- predict_trac(fit_classif, new_Z = Z, new_X = X,
                              output = "probability")
  expect_lt(max(probability[[1]][, 10]), 1)
  expect_gt(min(probability[[1]][, 10]), 0)
  cvfit_classif <- cv_trac(fit_classif, Z = Z, y = y_classif, A = A, X = X,
                           stratified = TRUE)
})
