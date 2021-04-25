#' sCD14 data
#'
#' This is the Gut (HIV) data used in [Tree-Aggregated Predictive Modeling of
#' Microbiome
#' Data](https://www.biorxiv.org/content/10.1101/2020.09.01.277632v1.full). This
#' was processed based on a phyloseq object provided by Javier Rivera-Pinto. It
#' is the BCN0 (Bacelona test dataset) from [Gut Microbiota Linked to Sexual
#' Preference and HIV
#' Infection](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4816837/).
#'
#' @format A named list:
#' \describe{
#' \item{y}{Vector of n = 152 soluble CD14
#'   levels in units of pg/ml}
#' \item{x}{Matrix of fecal 16S rRNA amplicon data,
#'   with n = 152 rows corresponding to people and p = 539 columns corresponding
#'   to OTUs.}
#' \item{tree}{Taxonomic tree of class `phylo`}
#' \item{tax}{A data frame containing the taxonomic information for each OTU}
#' \item{A}{A binary matrix encoding the tree structure with p = 539 rows,
#' corresponding to leaves, and 626 columns, corresponding to all non-root nodes
#' in the tree.  See [Tree-Aggregated Predictive Modeling of Microbiome
#' Data](https://www.biorxiv.org/content/10.1101/2020.09.01.277632v1.full) for
#' the definition of A.} }
"sCD14"
