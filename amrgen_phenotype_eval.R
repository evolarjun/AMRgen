#!/usr/bin/env Rscript

# amrgen_phenotype_eval.R

# Possibly need the following to get it to run:
# Install google cloud CLI and authenticate
# (e.g., for maximum danger and permissions: gcloud auth login --no-launch-browser)
# install.packages("remotes")
# remotes::install_github("AMRverse/AMRgen")
# install.packages("logistf")
# install.packages("bigrquery")
# install.packages("dplyr")
# install.packages("optparse")

library(optparse)

# Option Parsing
option_list <- list(
  make_option(c("--drug"),
    type = "character", default = NULL,
    help = "Antibiotic name (e.g., 'Ciprofloxacin')", metavar = "character"
  ),
  make_option(c("--taxgroup"),
    type = "character", default = NULL,
    help = "Taxonomic group (e.g., 'Pseudomonas aeruginosa')", metavar = "character"
  ),
  make_option(c("--drug_class_list"),
    type = "character", default = NULL,
    help = "Comma-separated list of drug classes (e.g., 'Quinolones,Fluoroquinolones')", metavar = "character"
  ),
  make_option(c("--standard"),
    type = "character", default = "clsi",
    help = "SIR interpretation standard ('clsi' or 'eucast') [default= %default]", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$drug) || is.null(opt$taxgroup) || is.null(opt$drug_class_list)) {
  print_help(opt_parser)
  stop("Missing required arguments. Please provide --drug, --taxgroup, and --drug_class_list", call. = FALSE)
}

if (!(opt$standard %in% c("clsi", "eucast"))) {
  stop("Invalid --standard argument. Must be 'clsi' or 'eucast'.", call. = FALSE)
}

# Analysis config
#################
# moved to after option parsing to fail and print usage message quicker
library(bigrquery)
library(dplyr)
library(AMRgen)
# library(logistf)

drug <- opt$drug
taxgroup <- opt$taxgroup
drug_class_list <- unlist(strsplit(opt$drug_class_list, ","))
# need to use AMRgen classes which
# I think are derived from the R AMR package, but I'm not sure:
# https://amr-for-r.org/reference/antimicrobial_selectors.html
# Actually it looks like they're set in AMRgen/data-raw/amrfp_drug_classes_agents.tsv
# https://github.com/AMRverse/AMRgen/blob/main/data-raw/amrfp_drug_classes_agents.tsv

# GCP project to run the queries under
project_id <- "ncbi-pd-amr"

# Adjust the following queries if needed

# Select which results to grab from MicroBIGG-E Can be more restrictive to
# reduce the number of elements you're examining (could be misleading in that
# case though)
mb_query <- "
WITH ast as (
  SELECT distinct(biosample_acc)
  FROM `ncbi-pathogen-detect.pdbrowser.ast`
  WHERE taxgroup_name = @taxgroup
)
SELECT mb.biosample_acc, element_symbol, subclass, type, subtype, amr_method,
  hierarchy_node, scientific_name
FROM
`ncbi-pathogen-detect.pdbrowser.microbigge` mb
INNER JOIN ast on ast.biosample_acc = mb.biosample_acc
WHERE type = 'AMR'
-- AND subclass IN UNNEST(@subclasses)"

# Grab AST data, limiting to drug and taxgroup of interest
ast_query <- "
select biosample_acc, platform, standard, disk_diffusion, mic, measurement_sign,
antibiotic, scientific_name, bioproject_acc, phenotype
from `ncbi-pathogen-detect.pdbrowser.ast`
where taxgroup_name = @taxgroup
and antibiotic = LOWER(@drug)
"

#####################
# End analysis config

message("Running BigQuery query to get phenotype data")
ast_query_job <- bq_project_query(project_id,
  query = ast_query,
  parameters = list(taxgroup = taxgroup, drug = drug)
)
pd_ast <- bq_table_download(ast_query_job) %>%
  rename(
    "Laboratory typing platform" = platform,
    "Testing standard" = standard,
    "Disk diffusion (mm)" = disk_diffusion,
    "MIC (mg/L)" = mic,
    "Measurement sign" = measurement_sign,
    "Antibiotic" = antibiotic,
    "BioProject" = bioproject_acc,
    "Resistance phenotype" = phenotype,
    "Scientific name" = scientific_name
  )

message("Importing phenotypes and calling SIR")
interpret_eucast_flag <- FALSE
interpret_clsi_flag <- FALSE
sir_col_name <- ""
message("Interpreting SIR using ", opt$standard)
if (opt$standard == "clsi") {
  interpret_clsi_flag <- TRUE
  sir_col_name <- "pheno_clsi"
} else if (opt$standard == "eucast") {
  interpret_eucast_flag <- TRUE
  sir_col_name <- "pheno_eucast"
}

ast <- import_ncbi_ast(
  input = pd_ast, sample_col = "biosample_acc",
  interpret_eucast = interpret_eucast_flag, interpret_clsi = interpret_clsi_flag
)

# grab the relevant entries from MicroBIGG-E
message("Running BigQuery query to get genotype data")
mb_query_job <- bq_project_query(project_id,
  query = mb_query,
  parameters = list(taxgroup = taxgroup)
)
pd_amrfp <- bq_table_download(mb_query_job) %>%
  rename(
    "Gene symbol" = element_symbol,
    "Element type" = type,
    "Method" = amr_method,
    "Hierarchy_node" = hierarchy_node,
    "Element subtype" = subtype,
    "Subclass" = subclass,
  )

geno <- import_amrfp(pd_amrfp, sample_col = "biosample_acc")

# now produce the graphs
message("Running analysis")
filename <- paste(
  sep = ".", drug, gsub(" ", "_", taxgroup),
  paste(collapse = ".", drug_class_list), "pdf"
)
pdf(file = filename)
# subclasses = sub("(.)", "\\U\\1", tolower(subclasses), perl = TRUE)
soloPPV <- solo_ppv_analysis(geno, ast,
  antibiotic = drug,
  drug_class_list = drug_class_list,
  sir_col = sir_col_name
)
upset <- amr_upset(soloPPV$amr_binary,
  min_set_size = 5, assay = "mic",
  order = "value", print_category_counts = T
)
dev.off()
message("Wrote analysis output to file ", filename)
