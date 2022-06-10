library(tidyverse) # data wrangling
library(phyloseq) # microbiome data wrangling

### Get metadata with labels and additional covarariates sex and age
# get labels malawi vs venezuela
sampledata_link <-
  "https://knights-lab.github.io/MLRepo/datasets/yatsunenko/task-malawi-venezuela.txt"
sampledata <-
  read.table(
    url(sampledata_link),
    sep = "\t",
    quote = "",
    row.names = 1,
    comment.char = "",
    header = TRUE
  )
# replace GAZ from the labels
sampledata <- sampledata %>%
  mutate(Var = str_replace_all(Var, "GAZ:", ""))
# add additional covariates like sex and age
additional_covariates_link <-
  "https://knights-lab.github.io/MLRepo/datasets/yatsunenko/mapping-orig.txt"
mapping_original <- read.delim(additional_covariates_link) %>%
  filter(X.SampleID %in% rownames(sampledata)) %>%
  select(X.SampleID, AGE, SEX) %>%
  rename(rowname = X.SampleID) %>%
  rename(sex = SEX) %>%
  rename(age = AGE) %>%
  mutate(age = as.numeric(age))

sampledata <- sampledata %>%
  rownames_to_column() %>%
  left_join(., mapping_original, by = "rowname") %>%
  column_to_rownames()

### get otu count table and filter for the observations with known labels
otu_link <-
  "http://metagenome.cs.umn.edu/public/MLRepo/yatsunenko2012.gg.otutable.txt"
otu <- read.table(
  otu_link,
  sep = "\t",
  comment = "",
  row.names = 1,
  head = TRUE,
  check.names = FALSE,
  quote = ""
)
otu <- t(otu)
# keep only samples with label
otu <- otu[rownames(otu) %in% rownames(sampledata),]

### get taxonomic table
# get greengenes db
gg_13_5_link <-
  "https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5_taxonomy.txt.gz"
gg_13_5_taxonomy <-
  read_delim(gg_13_5_link, delim = "\t", col_names = FALSE)
colnames(gg_13_5_taxonomy) <- c("value", "tax")
# transform the otu names to numeric and join the database with the
# otu names from the dataset
tax_preprocessing <- enframe(colnames(otu)) %>%
  mutate(value = as.numeric(value))
tax_preprocessing <- left_join(tax_preprocessing,
                               gg_13_5_taxonomy, by = "value")
# define the different taxonomic ranks
rank_names <- c("Kingdom",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species")
# create the taxonomic table
tax_preprocessing <- tax_preprocessing %>%
  select(!name) %>%
  column_to_rownames(var = "value") %>%
  separate(
    tax,
    into = rank_names,
    sep = ";",
    extra = "drop",
    fill = "right"
  ) %>%
  mutate(across(everything(), str_trim))
tax_preprocessing <- as.matrix(tax_preprocessing)

### create phyloseq object
otu <- otu_table(otu, taxa_are_rows = FALSE)
tax <- tax_table(tax_preprocessing)
sample_variable <- sample_data(sampledata)
phylo <- phyloseq(otu, tax, sample_variable)
# take only bacteria into consideration
phylo <- subset_taxa(phylo , Kingdom == "k__Bacteria")

# filter otus with less then 10% prevalence
freq <- colSums(sign(phylo@otu_table@.Data))
malawi <- prune_taxa(freq > 0.1 * nsamples(phylo), phylo)

usethis::use_data(malawi, overwrite = TRUE)

