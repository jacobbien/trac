# This file is adapted from here:
# https://rstudio.github.io/reticulate/articles/package.html

# global reference to classo (will be initialized in .onLoad)
classo <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  classo <<- reticulate::import("classo", delay_load = TRUE)
}

install_classo <- function(method = "auto", conda = "auto") {
  reticulate::py_install("c-lasso", method = method, conda = conda, pip = TRUE)
}
