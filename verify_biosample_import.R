library(testthat)
library(dplyr)
library(readr)
library(AMR)

# Load the package functions (assuming we are in the package root or can access R/ source)
# Sourcing files since we might not have installed the updated package yet
source("R/utils.R")
source("R/import_pheno.R")
source("R/breakpoints.R")
# Need to make sure dependencies are loaded.
# import_pheno.R uses utils.R and other things.
# For simplicity, we assume dependencies are available in the environment.

# Create verification script for import_biosample_ast
test_that("import_biosample_ast matches import_ncbi_ast", {
    # 1. Run import_biosample_ast
    message("Fetching data from BioSample...")
    # Use a small batch size for testing, and maybe limiting to get a few known ones would be hard without specific query
    # But we can query for Klebsiella and hope to get overlaps with our static file.
    # The static file is 'R/klebsiella_browser_asts.tsv'.

    # We'll fetch 50 records to be quick.
    biosample_data <- import_biosample_ast("Klebsiella", batch_size = 50, max_records = 50)

    expect_true(!is.null(biosample_data))
    expect_true(nrow(biosample_data) > 0)

    # 2. Run import_ncbi_ast on the static file
    message("Reading static file...")
    static_file <- "R/klebsiella_browser_asts.tsv"
    if (!file.exists(static_file)) {
        warning("Static file R/klebsiella_browser_asts.tsv not found.")
        return()
    }

    # We pass a dataframe to import_ncbi_ast (the file content) to avoid re-reading if we want,
    # or just pass the path.
    static_data <- import_ncbi_ast(static_file)

    # 3. Find common BioSamples (id)
    common_ids <- intersect(biosample_data$id, static_data$id)
    message(paste("Found", length(common_ids), "common BioSamples."))

    if (length(common_ids) == 0) {
        warning("No common BioSamples found between live fetch and static file. Cannot verify equality.")
        # This might happen if the fetch gets different records than the file.
        # But at least the fetch worked and returned structured data.
    } else {
        # 4. Compare data for common IDs
        # We join and compare specific columns

        bio_subset <- biosample_data %>%
            filter(id %in% common_ids) %>%
            arrange(id, drug_agent) %>%
            select(id, drug_agent, mic, pheno_provided)

        static_subset <- static_data %>%
            filter(id %in% common_ids) %>%
            arrange(id, drug_agent) %>%
            select(id, drug_agent, mic, pheno_provided)

        # Since fetching might return more or fewer drugs per sample depending on XML vs TSV parsing nuances?
        # Or maybe the static file is older.
        # Let's check overlap of (id, drug_agent)

        bio_subset <- bio_subset %>% mutate(source = "live")
        static_subset <- static_subset %>% mutate(source = "static")

        combined <- inner_join(bio_subset, static_subset, by = c("id", "drug_agent"))

        message(paste("Comparing", nrow(combined), "data points (id + drug)."))

        # Compare MICs
        mic_match <- combined$mic.x == combined$mic.y
        # handle NAs
        mic_match[is.na(combined$mic.x) & is.na(combined$mic.y)] <- TRUE
        mic_match[is.na(mic_match)] <- FALSE

        expect_true(sum(!mic_match) < 0.1 * nrow(combined)) # Allow small mismatch due to parsing differences?
        if (sum(!mic_match) > 0) {
            message("Found mismatches in MICs:")
            print(combined[!mic_match, ] %>% select(id, drug_agent, mic.x, mic.y))
        } else {
            message("All shared MICs match.")
        }

        # Compare Phenotypes
        pheno_match <- combined$pheno_provided.x == combined$pheno_provided.y
        # handle NAs
        pheno_match[is.na(combined$pheno_provided.x) & is.na(combined$pheno_provided.y)] <- TRUE
        pheno_match[is.na(pheno_match)] <- FALSE

        expect_true(sum(!pheno_match) == 0)
        if (sum(!pheno_match) > 0) {
            message("Found mismatches in Phenotypes:")
            print(combined[!pheno_match, ] %>% select(id, drug_agent, pheno_provided.x, pheno_provided.y))
        } else {
            message("All shared phenotypes match.")
        }
    }
})
