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
