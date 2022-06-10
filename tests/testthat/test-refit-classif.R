library(reticulate)

# run devtools::test() to run these tests.

skip_if_no_classo <- function() {
  # helper function to skip tests if we don't have the 'classo' module
  # got this from here: https://rstudio.github.io/reticulate/articles/package.html
  if (!py_module_available("classo"))
    skip("classo not available for testing")
}

set.seed(123)
log_pseudo <- function(x, pseudo_count = 1)
  log(x + pseudo_count)

## Simple simulated data -------------------------------------------------------
# set number of observations and number of otus/asvs
n <- 100
p <- 300

# construct data matrix
binom_sim <-
  matrix(rbinom(n = n * p, size = 1, prob = 1 / 4), nrow = n)
otu_sim <- matrix(ifelse(binom_sim == 0, 0,
                         rnbinom(
                           n = n * p,
                           mu = exp(10),
                           size = 2
                         )),
                  nrow = n,
                  ncol = p)

# Z <- matrix(rnorm(n = n * p), nrow = n)
# no negative values allowed
# Z <- abs(Z)
# take log
Z <- log_pseudo(otu_sim)

# add additional covariates: one numerical and one categorical
X <- data.frame(
  numeric_feature = rnorm(n = n, mean = 0, sd = 2),
  categorical_feature = sample(c(0, 1), replace = TRUE,
                               size = n)
)

X$categorical_feature <- as.factor(X$categorical_feature)

y <- +Z[, 1] - Z[, 2] + Z[, 3] - Z[, 4] + X[, 1] +
  (as.numeric(X[, 2]) - 1) + rnorm(n)
# for classifcation take the sign --> values in c(-1, 1)
y_classif <- sign(y)

test_that("refit classification log contrast on simulated data", {
  skip_if_no_classo()
  fit_sparse_log_contrast <-
    sparse_log_contrast(
      Z = Z,
      y = y_classif,
      min_frac = 1e-2,
      nlam = 30,
      method = "classif"
    )
  clasif_cv_sparse_log_contrast <-
    cv_sparse_log_contrast(fit = fit_sparse_log_contrast, Z = Z, y = y_classif)
  i_selected <- clasif_cv_sparse_log_contrast$cv$i1se
  refit_sparse_log_contrast <- refit_sparse_log_contrast_classif(
    fit = fit_sparse_log_contrast,
    i_selected = i_selected,
    Z = Z,
    y = y_classif
  )
  # check zero-sum constraint
  expect_true(sum(refit_sparse_log_contrast$beta) <= 1e-10)
  # check the number of selected betas (should be equal to the )
  expect_true(
    length(refit_sparse_log_contrast$beta) ==
      sum(fit_sparse_log_contrast$beta[, i_selected] != 0)
  )
})

test_that("refit classification log contrast on simulated data with add covar",
          {
            skip_if_no_classo()
            fit_sparse_log_contrast <-
              sparse_log_contrast(
                Z = Z,
                y = y_classif,
                additional_covariates = X,
                min_frac = 1e-2,
                nlam = 30,
                method = "classif",
                w_additional_covariates = c(0.02, 0.02)
              )
            clasif_cv_sparse_log_contrast <-
              cv_sparse_log_contrast(
                fit = fit_sparse_log_contrast,
                Z = Z,
                y = y_classif,
                additional_covariates = X
              )
            i_selected <- clasif_cv_sparse_log_contrast$cv$i1se
            refit_sparse_log_contrast <-
              refit_sparse_log_contrast_classif(
                fit = fit_sparse_log_contrast,
                i_selected = i_selected,
                Z = Z,
                y = y_classif,
                additional_covariates = X
              )
            # check zero-sum constraint for compositional covariates
            p_compositional <-
              length(refit_sparse_log_contrast$beta) - 2
            sum_compositional <- sum(refit_sparse_log_contrast$beta[, 1:p_compositional])
            expect_true(sum_compositional <= 1e-10)
            # check the number of selected betas (should be equal to the )
            expect_true(
              length(refit_sparse_log_contrast$beta) ==
                sum(fit_sparse_log_contrast$beta[, i_selected] != 0)
            )
          })
