# ===================================================================== #
#  Licensed as GPL-v3.0.                                                #
#                                                                       #
#  Developed as part of the AMRverse (https://github.com/AMRverse):     #
#  https://github.com/AMRverse/AMRgen                                   #
#                                                                       #
#  We created this package for both routine data analysis and academic  #
#  research and it was publicly released in the hope that it will be    #
#  useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                       #
#  This R package is free software; you can freely use and distribute   #
#  it for both personal and commercial purposes under the terms of the  #
#  GNU General Public License version 3.0 (GNU GPL-3), as published by  #
#  the Free Software Foundation.                                        #
# ===================================================================== #

#' AMR Logistic Regression Analysis
#'
#' Performs logistic regression to analyze the relationship between genetic markers and phenotype (R, and NWT) for a specified antibiotic.
#' @param geno_table A data frame containing the genotype data.
#' @param pheno_table A data frame containing the phenotypic data.
#' @param antibiotic A character string specifying the antibiotic to model using logistic regression.
#' @param drug_class_list A vector of drug class names. Used to subset the relevant markers for analysis.
#' @param geno_sample_col (Optional) A character string specifying the column in `geno_table` that identifies the sample IDs. Defaults to `NULL`.
#' @param pheno_sample_col (Optional) A character string specifying the column in `pheno_table` that identifies the sample IDs. Defaults to `NULL`.
#' @param sir_col (Optional) A character string specifying the column in `pheno_table` that contains the phenotype values (e.g., resistance/susceptibility). Defaults to `"pheno"`.
#' @param ecoff_col (Optional) A character string specifying the column in `pheno_table` containing the ECOFF (epidemiological cutoff) values. Defaults to `"ecoff"`.
#' @param maf (Optional) An integer specifying the minimum allele frequency (MAF) threshold. Markers with a MAF lower than this value will be excluded. Defaults to 10.
#' @param fit_glm (Optional) Change to TRUE to fit model with glm. Otherwise fit model with logistf (default).
#' @param single_plot (Optional) A logical value. If `TRUE`, a single plot is produced comparing the estimates for resistance (`R`) and non-resistance (`NWT`). Otherwise, two plots are printed side-by-side. Defaults to `TRUE`.
#' @param colors (Optional) A vector of two colors, to use for R and NWT models in the plots. Defaults to `c("maroon", "blue4")`.
#' @param axis_label_size (Optional) A numeric value controlling the size of axis labels in the plot. Defaults to 9.
#' @param marker_col (Optional) Name of the column containing the marker identifiers, whose unique values will be treate as predictors in the regression. Defaults to `"marker"`.
#' @importFrom dplyr any_of select where
#' @importFrom ggplot2 ggtitle
#' @importFrom stats glm
#' @importFrom logistf logistf
#' @return A list with the following components:
#' - `bin_mat`: The binary matrix of genetic data and phenotypic resistance information.
#' - `modelR`: The fitted logistic regression model for resistance (`R`).
#' - `modelNWT`: The fitted logistic regression model for non-resistance (`NWT`).
#' - `plot`: A ggplot object comparing the estimates for resistance and non-resistance with corresponding statistical significance indicators.
#' @export
#' @examples
#' # Example usage of the amr_logistic function
#' result <- amr_logistic(
#'   geno_table = import_amrfp(ecoli_geno_raw, "Name"),
#'   pheno_table = ecoli_ast,
#'   sir_col = "pheno_clsi",
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   maf = 10
#' )
#' # To access the plot:
#' print(result$plot)
amr_logistic <- function(geno_table, pheno_table, antibiotic, drug_class_list,
                         geno_sample_col = NULL, pheno_sample_col = NULL,
                         sir_col = "pheno", ecoff_col = "ecoff",
                         maf = 10, fit_glm = FALSE, single_plot = TRUE,
                         colors = c("maroon", "blue4"),
                         axis_label_size = 9, marker_col = "marker.label") {

  bin_mat <- get_binary_matrix(
    geno_table = geno_table,
    pheno_table = pheno_table,
    antibiotic = antibiotic,
    drug_class_list = drug_class_list,
    geno_sample_col = geno_sample_col,
    pheno_sample_col = pheno_sample_col,
    sir_col = sir_col,
    ecoff_col = ecoff_col,
    marker_col = marker_col
  )

  if (fit_glm) {
    cat("Fitting logistic regression models using glm\n")
    if (sum(!is.na(bin_mat$R)) > 0) {
      modelR <- glm(R ~ ., data = bin_mat %>% select(-any_of(c("id", "pheno", "ecoff", "mic", "disk", "NWT"))) %>% select(where(~ sum(., na.rm = TRUE) >= maf)), family = stats::binomial(link = "logit"))
      modelR <- glm_details(modelR) %>%
        mutate(marker = gsub("\\.\\.", ":", marker)) %>%
        mutate(marker = gsub("`", "", marker))
    }
    if (sum(!is.na(bin_mat$NWT)) > 0) {
      modelNWT <- glm(NWT ~ ., data = bin_mat %>% select(-any_of(c("id", "pheno", "ecoff", "mic", "disk", "R"))) %>% select(where(~ sum(., na.rm = TRUE) >= maf)), family = stats::binomial(link = "logit"))
      modelNWT <- glm_details(modelNWT) %>%
        mutate(marker = gsub("\\.\\.", ":", marker)) %>%
        mutate(marker = gsub("`", "", marker))
    }
  } else {
    cat("...Fitting logistic regression model to R using logistf\n")
    if (sum(!is.na(bin_mat$R))>0) {
      to_fit <- bin_mat %>% filter(!is.na(R)) %>% select(-any_of(c("id", "pheno", "ecoff", "mic", "disk", "NWT"))) %>% select(where(~ sum(., na.rm = TRUE) >= maf))
      modelR <- logistf::logistf(R ~ ., data = to_fit, pl = FALSE)
      modelR <- logistf_details(modelR) %>%
        mutate(marker = gsub("\\.\\.", ":", marker)) %>%
        mutate(marker = gsub("`", "", marker))
    }
    cat("...Fitting logistic regression model to NWT using logistf\n")
    if (sum(!is.na(bin_mat$NWT))>0) {
      to_fit <- bin_mat %>% filter(!is.na(NWT)) %>% select(-any_of(c("id", "pheno", "ecoff", "mic", "disk", "R"))) %>% select(where(~ sum(., na.rm = TRUE) >= maf))
      modelNWT <- logistf::logistf(NWT ~ ., data = to_fit, pl = FALSE)
      modelNWT <- logistf_details(modelNWT) %>%
        mutate(marker = gsub("\\.\\.", ":", marker)) %>%
        mutate(marker = gsub("`", "", marker))
    }
  }

  cat("Generating plots\n")
  if (exists("modelR") & exists("modelNWT")) { # if we have 2 models, plot them together
    plot <- compare_estimates(modelR, modelNWT,
                              single_plot = single_plot,
                              title1 = "R", title2 = "NWT",
                              colors = colors, axis_label_size = axis_label_size
    )
    if (single_plot) {
      ggtitle(
        label = paste("R and NWT for", antibiotic),
        subtitle = paste("for", paste(drug_class_list, collapse = ","), "markers present in at least", maf, "samples")
      )
    }
  } else if (exists("modelR")) {
    plot <- plot_estimates(modelR)
  } else if (exists("modelNWT")) {
    plot <- plot_estimates(modelNWT)
  }
  if (!exists("modelNWT")) {modelNWT <- NULL} # need an object to return, set to null
  if (!exists("modelR")) {modelR <- NULL} # need an object to return, set to null

  print(plot)

  return(list(
    bin_mat = bin_mat,
    modelR = modelR,
    modelNWT = modelNWT,
    plot = plot
  ))
}
