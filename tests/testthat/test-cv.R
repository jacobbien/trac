test_that("Check stratified fold creation", {
  y <- c(1,1,1,1,1,-1,-1,-1,1,1)
  n <- length(y)
  nfolds <- 3
  folds <- trac:::make_folds_stratified(n, nfolds, y)
  # check if every fold contains 1 observation of label -1
  class_per_group <- matrix(ncol = 2, nrow = 0)
  colnames(class_per_group) <- c("-1", "1")
  for (i in 1:nfolds) {
    class_per_group <- rbind(class_per_group, c(table(y[folds[[i]]])))
  }
  expect_equal(class_per_group[, 1], c(1, 1, 1))
  # check if there are duplicates in folds
  index <- c()
  for (i in 1:nfolds) {
    index <- c(index, folds[[i]])
  }
  expect_equal(sort(index), 1:n)
})


