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


#' malawi vs venezuela, adults only
#'
#' This is a subset of the "malawi vs venezuela, adults only" task from the
#' [Microbiome Learning
#' Repo](
#' https://knights-lab.github.io/MLRepo/docs/yatsunenko_malawi_venezuela.html).
#' The task is to predict whether the individuals live in Malawi or Venezuela
#' based on the OTU count table.
#' OTUs that occur in less than 10% of the samples are excluded. Only
#' Bacteria are included. The
#' taxonomic table for the OTUs was obtained by the [greengenes reference
#' database](https://www.nature.com/articles/ismej2011139?report=reader).
#' The data was originally published in [Yatsunenko, T.,
#' Rey, F. E., Manary, M. J., Trehan, I., Dominguez-Bello, M. G., Contreras, M.,
#'  ... & Gordon, J. I. (2012).
#'  Human gut microbiome viewed across age and geography.
#'  nature, 486(7402),
#'  222-227.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3376388/).
#'
#'
#' @format A phyloseq object:
#' \describe{
#'   \item{otu_table}{OTU count table with 54 rows corresponding to subjects
#'     and 5008 columns corresponding to OTUs}
#'   \item{tax_table}{Taxonomic table with the taxonomic assignment on different
#'     levels for the OTUs. The names of the OTUs are saved in the rownames.
#'     The 5008 rows correspond to the different OTUS and the 7 columns to
#'     the different levels of the taxonomic tree (Kingdom, Phylum, Class,
#'     Order, Family, Genus, Species)}
#'   \item{sam_data}{Meta data about the different observations. Contains
#'     the labels Vars and additional covariates sex (binary) and age (numeric)}
#' }
"malawi"
