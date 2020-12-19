context("conv")
library("Matrix")

test_that("Construct A for non compositional data works", {
  # construct A with 2 compositional hiracical features
  A <- diag(x = 1, nrow = 2, ncol = 2)
  A <- cbind(A, c(1, 1))
  # name two compositional features compositional1, compositional2
  rownames(A) <- c("compositional1", "compositional1")
  # create a non compositional and non hiracical feature named p1
  X <- data.frame(p1 = c(0, 0, 0))

    t_size <- ncol(A) + 1
  p <- 2
  p_x <- 1

  # Expect new A as (A 0
  #                  0 1)
  test_A <- matrix(data = c(1,0,1,0,0,1,1,0,0,0,0,1), nrow = 3, ncol = 4,
                   byrow = TRUE)
  # Expect the names to be the names of compositional data + non compositional
  # data
  test_A_names <- c(rownames(A), "p1")
  expect_equal(unname(as.matrix(A_add_X(X, A, t_size, p, p_x))),
               test_A)
  expect_equal(rownames(A_add_X(X, A, t_size, p, p_x)),
               test_A_names)
  })
