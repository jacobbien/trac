#' Aggregate data to a fixed level
#'
#' This function aggregates a data set to a fixed level of a taxonomy.  For
#' example, if aggregating to the Phylum level, it would sum up the columns of x
#' corresponding to OTUs/ASVs in that Phylum.  It also adjusts the tree information
#' (in particular A, tree, and tax) for this new smaller problem.
#'
#' @param x n by p data matrix of counts, where n is the number of samples, p is
#' the number of OTUs/ASVs.
#' @param y response vector of length n
#' @param A p by (t_size-1) binary matrix giving tree structure (t_size is the
#' total number of nodes and the -1 is because we do not include the root).  This
#' is created by the function \code{\link{phylo_to_A}}.
#' @param tax a tax table, which is a p by number-of-levels matrix, where cell ij gives
#' the name of the level j ancestor of OTU/ASV i.
#' @param level which level are we aggregating to?
#' @param collapse see \code{\link{tax_table_to_phylo}}
#' @export
aggregate_to_level <- function(x, y, A, tax, level = 7, collapse = FALSE) {
  # x: n by p (p = num OTUs)
  # y: n vector
  # tax: tax table, which is p by number-of-levels
  # level: column of tax (corresponding to taxonomic rank)
  tax <- unique(tax[, 1:level])
  levels <- colnames(tax)[1:level]
  formula_str <- paste0("~", paste(levels, collapse = "/"))
  tree1 <- tax_table_to_phylo(eval(parse(text = formula_str)),
                              data = tax, collapse = collapse)
  A_agg <- phylo_to_A(tree1)
  taxa_of_otus <- rownames(A) %>%
    stringr::str_split("::") %>%
    purrr::map_chr(~ paste(.x[1:level], collapse = "::"))
  unique_taxa <- unique(taxa_of_otus)
  otu_to_taxa <- Matrix::sparseMatrix(i = 1:length(taxa_of_otus),
                              j = match(taxa_of_otus, unique_taxa))
  x_agg <- x %*% otu_to_taxa
  colnames(x_agg) <- unique_taxa

  list(x = x_agg, y = y, A = A_agg, tree = tree1, tax = tax)
}
