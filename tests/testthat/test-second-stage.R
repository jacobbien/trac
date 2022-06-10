library(reticulate)

skip_if_no_classo <- function() {
  # helper function to skip tests if we don't have the 'classo' module
  # got this from here:
  # https://rstudio.github.io/reticulate/articles/package.html
  if (!py_module_available("classo"))
    skip("classo not available for testing")
}

set.seed(123)
log_pseudo <- function(x, pseudo_count = 1) log(x + pseudo_count)
get_nz <- function(x) x[x != 0]

## Simple simulated data -------------------------------------------------------
# set number of observations and number of otus/asvs
n <- 100
p <- 300

# construct data matrix
binom_sim <- matrix(rbinom(n = n * p, size = 1, prob = 1/3), nrow = n)
otu_sim <- matrix(ifelse(binom_sim == 0, 0,
                         rnbinom(n = n * p,
                                 mu = exp(10),
                                 size = 2)), nrow = n, ncol = p)

# Z <- matrix(rnorm(n = n * p), nrow = n)
# no negative values allowed
# Z <- abs(Z)
# take log
Z <- log_pseudo(otu_sim)
# construct a taxonomic tree
A <- diag(p)
A <- cbind(A, c(rep(x = 1, times = p/3), rep(x = 0, times = (2*p/3))))
A <- cbind(A, c(rep(x = 0, times = p/3), rep(x = 1, times = (p/3)),
                rep(x = 0, times = (p/3))))
A <- cbind(A, c(rep(x = 0, times = 2*p/3), rep(x = 1, times = (p/3))))
A <- cbind(A, c(rep(x = 1, times = p/3), rep(x = 0, times = (2*p/3))))
A <- cbind(A, c(rep(x = 0, times = 2*p/3), rep(x = 1, times = (p/3))))

# rooted
A <- cbind(A, rep(x = 1, times = p))
A_n <- colSums(A)

# add additional covariates: one numerical and one categorical
X <- data.frame(numeric_feature = rnorm(n = n),
                categorical_feature = sample(c(0, 1), replace = TRUE,
                                             size = n))

X$categorical_feature <- as.factor(X$categorical_feature)

y <- Z[, 1] - Z[, p] + Z[, 2] - Z[, (p - 1)] + rnorm(n)
y_classif <- sign(y)

test_that("two-stage sparse log contrast regression on simulated data", {
  skip_if_no_classo()
  set.seed(123)
  # regression
  fit_sparse_log_contrast <-
    sparse_log_contrast(Z = Z, y = y, min_frac = 1e-2, nlam = 30)
  reg_cv_sparse_log_contrast <-
    cv_sparse_log_contrast(fit = fit_sparse_log_contrast, Z = Z, y = y)
  pre_selected <-
    fit_sparse_log_contrast$beta[,
                                 reg_cv_sparse_log_contrast[[1]]$i1se]
  reg_two_stage <- second_stage(Z = Z, y = y, betas = pre_selected)
  expect_true(all(reg_two_stage$log_ratios[c(
    "comp_1///comp_300", "comp_2///comp_299")] != 0))
  # classification
  fit_sparse_log_contrast <-
    sparse_log_contrast(Z = Z, y = y_classif, min_frac = 1e-2, nlam = 30,
                        method = "classif")
  clasif_cv_sparse_log_contrast <-
    cv_sparse_log_contrast(fit = fit_sparse_log_contrast, Z = Z, y = y_classif)
  pre_selected <-
    fit_sparse_log_contrast$beta[,
                                 clasif_cv_sparse_log_contrast[[1]]$i1se]
  classif_two_stage <- second_stage(Z = Z, y = y_classif, betas = pre_selected,
                                method = "classif", criterion = "min")
  expect_true(all(classif_two_stage$log_ratios[c(
    "comp_1///comp_300", "comp_2///comp_299")] != 0))
})

y <- 2*Z[, 1] - 2*Z[, p] + 3*Z[, 2] - 3*Z[, (p - 1)] + X[, 1] +
  (as.numeric(X[, 2]) - 1) + rnorm(n)
y_classif <- sign(y)

test_that("two-stage sparse log contrast regression on simulated data with
          additional covariates", {
            skip_if_no_classo()
            set.seed(123)
            # regression
            fit_sparse_log_contrast <-
              sparse_log_contrast(Z = Z, y = y, additional_covariates = X,
                                  min_frac = 1e-2, nlam = 30)
            reg_cv_sparse_log_contrast <-
              cv_sparse_log_contrast(fit = fit_sparse_log_contrast, Z = Z,
                                     y = y, additional_covariates = X)
            pre_selected <-
              fit_sparse_log_contrast$beta[,
                                           reg_cv_sparse_log_contrast[[1]]$i1se]
            reg_two_stage <- second_stage(Z = Z, y = y,
                                          additional_covariates = X,
                                          betas = pre_selected)
            expect_true(all(reg_two_stage$log_ratios[c(
              "comp_1///comp_300", "comp_2///comp_299",
              "numeric_feature", "categorical_feature")] != 0))

            # classification
            fit_sparse_log_contrast <-
              sparse_log_contrast(Z = Z, y = y_classif,
                                  additional_covariates = X,
                                  min_frac = 1e-2,
                                  method = "classif")
            classif_cv_sparse_log_contrast <-
              cv_sparse_log_contrast(fit = fit_sparse_log_contrast, Z = Z,
                                     y = y_classif, additional_covariates = X)
            pre_selected <-
              fit_sparse_log_contrast$beta[,
                                           classif_cv_sparse_log_contrast[[1]]$
                                             i1se]
            classif_two_stage <- second_stage(Z = Z, y = y_classif,
                                          additional_covariates = X,
                                          betas = pre_selected,
                                          method = "classif")
            expect_true(all(classif_two_stage$log_ratios[c(
              "comp_1///comp_300", "comp_2///comp_299",
              "numeric_feature", "categorical_feature")] != 0))
})



# construct an augmented dataset with the aggregrated levels as geom means
Z_mod <- Z %*% A
Z_mod <- t(t(Z_mod) - A_n)

# build y as the combination of column p + 2, p + 3, 1 and 2
y <- 3*Z_mod[, (p + 2)] - 3*Z_mod[, (p + 3)] +
  10*Z_mod[, 1] - 10*Z_mod[, 2] + rnorm(n)
# for classifcation take the sign --> values in c(-1, 1)
y_classif <- sign(y)

test_that("two-stage trac on simulated data", {
  skip_if_no_classo()
  set.seed(123)
  # regression
  fit_trac <-
    trac(Z = Z, A = A, y = y, min_frac = 1e-2, nlam = 30)
  reg_cv_trac <-
    cv_trac(fit = fit_trac, Z = Z, A = A, y = y)
  pre_selected <- fit_trac[[1]]$gamma[, reg_cv_trac$cv[[1]]$ibest]
  reg_two_stage <- second_stage(Z = Z, y = y, A = A, betas = pre_selected)
  expect_true(all(reg_two_stage$log_ratios[c(
    "comp_1///comp_2", "comp_302///comp_303")] != 0))
  # classification
  fit_trac <-
    trac(Z = Z, A = A, y = y_classif, min_frac = 1e-2, nlam = 30,
                        method = "classif")
  clasif_cv_trac <-
    cv_trac(fit = fit_trac, Z = Z, A = A, y = y_classif)
  pre_selected <-
    fit_trac[[1]]$gamma[, clasif_cv_trac$cv[[1]]$ibest]
  classif_two_stage <- second_stage(Z = Z, A = A, y = y_classif,
                                    betas = pre_selected,
                                    method = "classif", criterion = "min")
  expect_true(all(classif_two_stage$log_ratios[c(
    "comp_1///comp_2", "comp_302///comp_303")] != 0))
})



# build y as the combination of column p + 2, p + 3, 1 and 2
y <- 3*Z_mod[, (p + 2)] - 3*Z_mod[, (p + 3)] +
  15*Z_mod[, 1] - 15*Z_mod[, 2] + X[, 1] + (as.numeric(X[, 2]) - 1) + rnorm(n)
# for classifcation take the sign --> values in c(-1, 1)
y_classif <- sign(y)

test_that("two-stage trac on simulated data", {
  skip_if_no_classo()
  set.seed(123)
  # regression
  fit_trac <-
    trac(Z = Z, A = A, y = y, min_frac = 1e-2, nlam = 30,
         additional_covariates = X)
  reg_cv_trac <-
    cv_trac(fit = fit_trac, Z = Z, A = A, y = y, additional_covariates = X)
  pre_selected <- fit_trac[[1]]$gamma[, reg_cv_trac$cv[[1]]$ibest]
  reg_two_stage <- second_stage(Z = Z, y = y, A = A, betas = pre_selected,
                                additional_covariates = X)
  expect_true(all(reg_two_stage$log_ratios[c(
    "comp_1///comp_2", "comp_302///comp_303",
    "numeric_feature", "categorical_feature")] != 0))
  # classification
  fit_trac <-
    trac(Z = Z, A = A, y = y_classif, min_frac = 1e-2, nlam = 30,
         method = "classif",
         additional_covariates = X, w_additional_covariates = c(0.02, 0.02))
  clasif_cv_trac <-
    cv_trac(fit = fit_trac, Z = Z, A = A, y = y_classif,
            additional_covariates = X)
  pre_selected <-
    fit_trac[[1]]$gamma[, clasif_cv_trac$cv[[1]]$ibest]
  classif_two_stage <- second_stage(Z = Z, A = A, y = y_classif,
                                    betas = pre_selected,
                                    method = "classif", criterion = "min",
                                    additional_covariates = X,
                                    alpha = 0.05)
  expect_true(all(classif_two_stage$log_ratios[c(
    "comp_1///comp_2", "comp_302///comp_303")] != 0))
})
