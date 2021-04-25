# trac: Tree-based Aggregation of Compositional Data

This R package, called `trac`, performs tree-based aggregation of compositional data with a particular focus on microbiome data.  It implements the method proposed in [Tree-Aggregated Predictive Modeling of Microbiome Data](https://www.biorxiv.org/content/10.1101/2020.09.01.277632v1).  It uses internally the [c-lasso](https://github.com/Leo-Simpson/c-lasso) solver in Python (described in [this paper](https://joss.theoj.org/papers/10.21105/joss.02844)).

The easiest way to install `trac` is by using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) R package (if not already installed, open R and run `install.packages("devtools")`).

To install `trac`, run

``` r
devtools::install_github("jacobbien/trac")
```

in R.
