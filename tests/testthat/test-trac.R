library(reticulate)

# run devtools::test() to run these tests.

skip_if_no_classo <- function() {
  # helper function to skip tests if we don't have the 'classo' module
  # got this from here: https://rstudio.github.io/reticulate/articles/package.html
  if (!py_module_available("classo"))
    skip("classo not available for testing")
}

set.seed(123)
ntot <- length(sCD14$y)
n <- round(2/3 * ntot)
tr <- sample(ntot, n)
log_pseudo <- function(x, pseudo_count = 1) log(x + pseudo_count)
ytr <- sCD14$y[tr]
yte <- sCD14$y[-tr]
ztr <- log_pseudo(sCD14$x[tr, ])
zte <- log_pseudo(sCD14$x[-tr, ])

get_nz <- function(x) x[x!=0]

test_that("trac (a=1) on sCD14", {
  skip_if_no_classo()
  fit <- trac(ztr, ytr, A = sCD14$A, min_frac = 1e-2, nlam = 30)
  alpha_nz <- get_nz(fit[[1]]$alpha[, 2])
  expect_equal(as.numeric(alpha_nz), c(306.4128 * c(1, -1)), tolerance = 1e-4)
  expect_equal(names(alpha_nz),
               c("Life::Bacteria::Firmicutes::Clostridia::Clostridiales::Ruminococcaceae",
                 "Life::Bacteria::Firmicutes::Clostridia::Clostridiales::Lachnospiraceae")
               )
  expect_equal(length(unique(fit[[1]]$beta[, 2])), 3)
})


test_that("trac (a=1/2) on sCD14", {
  skip_if_no_classo()
  fit2 <- trac(ztr, ytr, A = sCD14$A, min_frac = 1e-2,
               nlam = 30, w = Matrix::colSums(sCD14$A)^0.5)
  alpha_nz2 <- get_nz(fit2[[1]]$alpha[, 2])
  expect_equal(as.numeric(alpha_nz2[1:3]),
               c(-43.64246, -94.29725,  56.53847),
               tolerance = 1e-4)
  expect_equal(sum(alpha_nz2), 0, tolerance = 1e-4)
  expect_equal(names(alpha_nz2)[1:3],
               c("Life::Bacteria::Firmicutes::Negativicutes::Selenomonadales::Veillonellaceae::Mitsuokella::Otu000070",
                "Life::Bacteria::Actinobacteria",
                "Life::Bacteria::Bacteroidetes::Bacteroidia::Bacteroidales::Bacteroidaceae::Bacteroides"))
  expect_equal(length(unique(fit2[[1]]$beta[, 2])), 6)
})

test_that("sparse log contrast on sCD14", {
  skip_if_no_classo()
  fit3 <- sparse_log_contrast(ztr, ytr, min_frac = 1e-2, nlam = 30)
  beta_nz <- get_nz(fit3$beta[, 2])
  expect_equal(as.numeric(beta_nz), c(-1.029745, -27.359612, 28.389357),
               tolerance = 1e-4)
  expect_equal(names(beta_nz), c("Otu000038", "Otu000070", "Otu000014"))
})


## Simple simulated data
# set number of observations and number of otus/asvs
n <- 150
p <- 3

# construct data matrix
Z <- matrix(rnorm(n = n * p), nrow = n)
# no negative values allowed
Z <- abs(Z)
# take log
Z <- log(Z)
# construct a taxonomic tree
A <- diag(p)
A <- cbind(A, c(1, 1, 0))
# rooted
A <- cbind(A, rep(x = 1, times = p))

# add additional covariates
X <- data.frame(numeric_feature = rnorm(n = n),
                categorical_feature = sample(c(0, 1), replace = TRUE,
                                             size = n))

X$categorical_feature <- as.factor(X$categorical_feature)


# construct an augmented dataset with the aggregrated levels
Z_mod <- Z %*% A

# without intercept
# build y as the combination of column 4 and 3
y <- Z_mod[, 4] - Z_mod[, 3] + rnorm(n)
# for classifcation take the sign --> values in c(-1, 1)
y_classif <- sign(y)

test_that("trac classification on simulated data", {
  skip_if_no_classo()
  fit_classif <- trac(Z, y_classif, A, method = "classif", intercept = FALSE)
  expect_true(all(fit_classif[[1]]$alpha[, 2][c(3,4)] != 0))
  expect_true(all((fit_classif[[1]]$beta0 == 0)))

  fit_huber <- trac(Z, y_classif, A, method = "classif_huber",
                    intercept = FALSE)
  expect_true(all(fit_huber[[1]]$alpha[, 2][c(3,4)] != 0))
  expect_true(all(fit_huber[[1]]$beta0 == 0))
})

y <- 0.1 + Z_mod[, 4] - Z_mod[, 3] + X[, 1] + (as.numeric(X[, 2]) - 1) +
  rnorm(n)
# for classifcation take the sign --> values in c(-1, 1)
y_classif <- sign(y)

test_that("trac classification on simulated data with additional covariates", {
  fit_classif_1 <- trac(Z, y_classif, A, X = X, method = "classif",
                      intercept = TRUE,
                      normalized = TRUE)
  expect_true(all(fit_classif_1[[1]]$alpha[, 2][c(3,4)] != 0))

  fit_huber_1 <- trac(Z, y_classif, A, X = X, method = "classif_huber",
                    intercept = TRUE, normalized = TRUE)
  expect_true(all(fit_huber_1[[1]]$alpha[, 2][c(3,4)] != 0))

  fit_huber_2 <- trac(Z, y_classif, A, X = X, method = "classif_huber",
                      intercept = TRUE, normalized = FALSE, rho = -1,
                      w_meta = rep(0.0001, ncol(X)))
  expect_true(all(fit_huber_2[[1]]$alpha[, 20][c(3, 4, 6, 7)] != 0))
})




