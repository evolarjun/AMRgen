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

#' Import and Process AST Data from an NCBI File
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset, processes the data, and optionally interprets the results based on MIC or disk diffusion data. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be commpressed) and parses relevant columns (antibiotic names, species names, MIC or disk data) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF). If expected columns are not found warnings will be given, and interpretation may not be possible.
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in NCBI antibiogram format. These files can be downloaded from NCBI AST browser, e.g. https://www.ncbi.nlm.nih.gov/pathogens/ast#Pseudomonas%20aeruginosa
#' @param sample_col A string indicating the name of the column with sample identifiers. If `NULL`, assume this is 'BioSample'.
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against ECOFF human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the field 'Scientific name' will be assumed to specify the species for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the field 'Antibiotic' will be assumed to specify the antibiotic for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @param source (optional) A single value to record as the source of these data points, e.g. "NCBI_browser". By default, the field 'BioProject' will be used to indicate the source for each row in the input file, but if this is missing or you want to override it with a single value for all samples, you may provide a source name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of coalesce if_else mutate relocate rename
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The biological sample identifier (renamed from `BioSample` or specified column).
#' - `spp_pheno`: The species phenotype, formatted using the `as.mo` function.
#' - `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function.
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the `as.mic` function.
#' - `disk`: The disk diffusion measurement (in mm), formatted using the `as.disk` function.
#' - `method`: The AST platform recorded in the input file as the source of the measurement.
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay.
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R), based on the MIC or disk diffusion data.
#' - `pheno_provided`: The original phenotype interpretation provided in the input file.
#' - `source`: The source of each data point (renamed from `BioProject` in the input file, or replaced with a single value passed in as the 'source' parameter).
#' @export
#' @examples
#' # small example E. coli AST data from NCBI
#' ecoli_ast_raw
#'
#' # import without re-interpreting resistance
#' pheno <- import_ncbi_ast(ecoli_ast_raw)
#' head(pheno)
#'
#' # import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
#' pheno <- import_ncbi_ast(ecoli_ast_raw, interpret_eucast = TRUE, interpret_ecoff = TRUE)
#' head(pheno)
import_ncbi_ast <- function(input, sample_col = "BioSample", source = NULL, species = NULL, ab = NULL,
                            interpret_eucast = FALSE, interpret_clsi = FALSE, interpret_ecoff = FALSE) {
  ast <- process_input(input)

  # find id column
  if (!is.null(sample_col)) {
    if (sample_col %in% colnames(ast)) {
      ast <- ast %>% rename(id = any_of(sample_col))
    } else {
      stop(paste("Invalid column name:", sample_col))
    }
  } else {
    stop("Please specify the column containing sample identifiers, via parameter 'sample_col'")
  }

  if ("Laboratory typing platform" %in% colnames(ast)) {
    ast <- ast %>% mutate(method = `Laboratory typing platform`)
  } else {
    cat("Warning: Expected AST platform column 'Laboratory typing platform' not found in input\n")
  }

  # parse guideline column
  if ("Testing standard" %in% colnames(ast)) {
    ast <- ast %>% mutate(guideline = `Testing standard`)
  } else {
    cat("Warning: Expected AST standard column 'Testing standard' not found in input\n")
  }

  # parse disk column
  if ("Disk diffusion (mm)" %in% colnames(ast)) {
    ast <- ast %>% mutate(disk = as.disk(`Disk diffusion (mm)`))
  } else {
    cat("Warning: Expected column 'Disk diffusion (mm)' not found in input\n")
  }

  # parse mic column
  if ("MIC (mg/L)" %in% colnames(ast)) {
    if ("Measurement sign" %in% colnames(ast)) {
      ast <- ast %>%
        mutate(mic = paste0(`Measurement sign`, `MIC (mg/L)`)) %>%
        mutate(mic = gsub("==", "", mic))
    } else {
      ast <- ast %>% mutate(mic = `MIC (mg/L)`)
      cat("Warning: Expected column 'Measurement sign' not found in input, be careful interpreting MIC\n")
    }
    ast <- ast %>% mutate(mic = as.mic(mic))
  } else {
    cat("Warning: Expected column 'MIC (mg/L)' not found in input\n")
  }

  # parse antibiotic column
  if ("Antibiotic" %in% colnames(ast)) {
    ast <- ast %>% mutate(drug_agent = as.ab(Antibiotic))
  } else {
    stop("Expected column 'Antibiotic' not found in input.")
  }

  # parse species column
  if ("Scientific name" %in% colnames(ast)) {
    ast <- ast %>% mutate(spp_pheno = as.mo(`Scientific name`))
  } else {
    cat("Warning: Expected species column 'Scientific name' not found in input\n")
  }

  # parse bioproject column as source
  if ("BioProject" %in% colnames(ast)) {
    ast <- ast %>% mutate(source = BioProject)
  } else {
    if (!is.null(source)) {
      ast <- ast %>% mutate(source = source)
      cat(paste0("Setting source to user-provided value: ", source, "\n"))
    } else {
      cat("Warning: Expected column 'BioProject' not found in input\n")
    }
  }

  # parse phenotype SIR column
  if ("Resistance phenotype" %in% colnames(ast)) {
    ast <- ast %>% mutate(pheno_provided = as.sir(`Resistance phenotype`))
  } else {
    cat("Warning: Expected phenotype SIR column 'Resistance phenotype' not found in input\n")
  }

  ast <- interpret_ast(ast, interpret_ecoff = interpret_ecoff, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, species = species, ab = ab)

  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", "pheno_eucast", "pheno_clsi", "ecoff", "guideline", "method", "source", "pheno_provided", "spp_pheno")))

  return(ast)
}


#' Import and Process AST Data from an EBI File
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset, processes the data, and optionally interprets the results based on MIC or disk diffusion data. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be commpressed) and parses relevant columns (antibiotic names, species names, MIC or disk data) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF). If expected columns are not found warnings will be given, and interpretation may not be possible.
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in EBI antibiogram format. These files can be downloaded from the EBI AMR browser, e.g. https://www.ebi.ac.uk/amr/data/?view=experiments
#' @param sample_col A string indicating the name of the column with sample identifiers. If `NULL`, assume this is 'phenotype-BioSample_ID'.
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against ECOFF human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the field 'phenotype-organism' will be assumed to specify the species for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the field 'phenotype-antibiotic_name' will be assumed to specify the antibiotic for each row in the input file, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @param source (optional) A single value to record as the source of these data points, e.g. "EBI_browser". By default, the field 'phenotype-AMR_associated_publications' will be used to indicate the source for each row in the input file, but if this is missing or you want to override it with a single value for all samples, you may provide a source name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of coalesce if_else mutate relocate rename
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The biological sample identifier (renamed from `phenotype-BioSample_ID` or specified column).
#' - `spp_pheno`: The species phenotype, formatted using the `as.mo` function.
#' - `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function.
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the `as.mic` function.
#' - `disk`: The disk diffusion measurement (in mm), formatted using the `as.disk` function.
#' - `method`: The AST platform recorded in the input file as the source of the measurement.
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay.
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R), based on the MIC or disk diffusion data.
#' - `pheno_provided`: The original phenotype interpretation provided in the input file.
#' - `source`: The source of each data point (renamed from the publications field in the input file, or replaced with a single value passed in as the 'source' parameter).
#' @export
#' @examples
#' \dontrun{
#' # import without re-interpreting resistance
#' pheno <- import_ebi_ast("EBI_AMR_data.csv.gz")
#' head(pheno)
#'
#' # import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
#' pheno <- import_ebi_ast("EBI_AMR_data.csv.gz", interpret_eucast = TRUE, interpret_ecoff = TRUE)
#' head(pheno)
#' }
import_ebi_ast <- function(input, sample_col = "phenotype-BioSample_ID", source = NULL, species = NULL, ab = NULL,
                           interpret_eucast = FALSE, interpret_clsi = FALSE, interpret_ecoff = FALSE) {
  ast <- process_input(input)

  # find id column
  if (!is.null(sample_col)) {
    if (sample_col %in% colnames(ast)) {
      ast <- ast %>% rename(id = any_of(sample_col))
    } else {
      stop(paste("Invalid column name:", sample_col))
    }
  } else {
    stop("Please specify the column containing sample identifiers, via parameter 'sample_col'")
  }

  # parse disk data
  if ("phenotype-gen_measurement" %in% colnames(ast)) {
    ast <- ast %>%
      mutate(disk = if_else(grepl("mm", `phenotype-gen_measurement`),
        as.disk(`phenotype-gen_measurement`), NA
      )) %>%
      mutate(disk = as.disk(disk))
  } else {
    ast <- ast %>% mutate(disk = as.disk(NA))
    cat("No disk data (units 'mm') found in input\n")
  }

  # parse mic data
  if ("phenotype-gen_measurement" %in% colnames(ast)) {
    ast <- ast %>%
      mutate(mic = if_else(grepl("mg/L", `phenotype-gen_measurement`),
        as.mic(`phenotype-gen_measurement`), NA
      )) %>%
      mutate(mic = as.mic(mic))
  } else {
    ast <- ast %>% mutate(mic = as.mic(NA))
    cat("No MIC data (units 'mg/L') found in input\n")
  }

  # parse antibiotic column
  if ("phenotype-antibiotic_name" %in% colnames(ast)) {
    ast <- ast %>% mutate(drug_agent = as.ab(`phenotype-antibiotic_name`))
  } else {
    stop("Expected drug name column 'phenotype-antibiotic_name' not found in input.")
  }

  # parse species column
  if ("phenotype-organism" %in% colnames(ast)) {
    ast <- ast %>% mutate(spp_pheno = as.mo(`phenotype-organism`))
  } else {
    cat("Warning: Expected species column 'phenotype-organism' not found in input\n")
  }

  if ("phenotype-platform" %in% colnames(ast)) {
    ast <- ast %>% mutate(method = `phenotype-platform`)
  } else {
    cat("Warning: Expected AST platform column 'phenotype-platform' not found in input\n")
  }

  if ("phenotype-ast_standard" %in% colnames(ast)) {
    ast <- ast %>% mutate(guideline = `phenotype-ast_standard`)
  } else {
    cat("Warning: Expected AST standard column 'phenotype-ast_standard' not found in input\n")
  }

  if ("phenotype-AMR_associated_publications" %in% colnames(ast)) {
    ast <- ast %>% mutate(source = `phenotype-AMR_associated_publications`)
  } else {
    if (!is.null(source)) {
      ast <- ast %>% mutate(source = source)
    }
    cat("Warning: Expected pubmed ID column 'phenotype-AMR_associated_publications' not found in input\n")
  }

  if ("phenotype-resistance_phenotype" %in% colnames(ast)) {
    ast <- ast %>% mutate(pheno_provided = as.sir(`phenotype-resistance_phenotype`))
  } else {
    cat("Warning: Expected pubmed ID column 'phenotype-resistance_phenotype' not found in input\n")
  }

  ast <- interpret_ast(ast, interpret_ecoff = interpret_ecoff, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, species = species, ab = ab)

  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", "pheno_eucast", "pheno_clsi", "ecoff", "guideline", "method", "source", "pheno_provided", "spp_pheno")))

  return(ast)
}


#' Interpret AST data in a standard format tibble
#'
#' This function applies human EUCAST or CLSI breakpoints, and/or ECOFF, to interpret AST data.
#' @param ast A tibble containing the AST measures in standard AMRgen format, as output by `import_ast`. It must contain assay measurements in columns 'mic' (class mic) and/or 'disk'. Interpretation requires an organism (column 'spp_pheno' of class 'mo', or a single value passed via the 'species' parameter) and an antibiotic (column 'drug_agent' of class 'ab', or a single value passed via the 'ab' parameter).
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against ECOFF human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the spp_pheno field in the input file will be assumed to specify the species for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the drug_agent field in the input file will be assumed to specify the antibiotic for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr across coalesce mutate
#' @return A copy of the input table, with additional columns:
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R), based on the MIC or disk diffusion data.
#' - `spp_pheno`: The species phenotype, formatted using the `as.mo` function (either taken from the input table, or the single value specified by 'species' parameter).
#' - `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function (either taken from the input table, or the single value specified by 'ab' parameter)..
#' @export
#' @examples
#' # import without re-interpreting resistance
#' pheno <- import_ncbi_ast(ecoli_ast_raw)
#' head(pheno)
#'
#' # interpret phenotypes
#' pheno <- interpret_ast(pheno)
#' 
#' \dontrun{
#' pheno <- read_csv("AST.csv") %>% 
#'   mutate(drug_agent=as.ab(antibiotic)) %>% # convert antibiotic field to 'drug_agent' of class 'ab'
#'   mutate(mic=paste0(sign,MIC)) %>% 
#'   mutate(mic=as.mic(mic)) # create a single 'mic' column of class 'mic'
#'   
#' pheno <- interpret_ast(pheno, species="Escherichia coli")
#' }
interpret_ast <- function(ast, interpret_ecoff = TRUE, interpret_eucast = TRUE, interpret_clsi = TRUE, species = NULL, ab = NULL) {
  if (interpret_ecoff | interpret_eucast | interpret_clsi) {
    # check we have species
    if (!is.null(species)) {
      cat(paste0("Interpreting all data as species: ", mo_name(as.mo(species)), "\n"))
      if ("spp_pheno" %in% colnames(ast)) {
        if (length(unique(ast$spp_pheno)) > 1) {
          cat("Warning: ignoring 'spp_pheno' column in input table, which contains multiple species:\n")
          cat(paste(unique(ast$spp_pheno), collapse = ", "))
          cat("\n")
        } else if (as.mo(unique(ast$spp_pheno)) != as.mo(species)) {
          cat(paste0("Warning: ignoring 'spp_pheno' column in input table, which indicates a different species: ", unique(ast$spp_pheno), "\n"))
        }
      }
      ast <- ast %>% mutate(spp_pheno = as.mo(species))
    } else if (!("spp_pheno" %in% colnames(ast))) {
      stop("Warning: Could not interpret data, need to provide 'species' parameter or 'spp_pheno' column in input table")
    }

    # check we have antibiotic
    if (!is.null(ab)) {
      cat(paste0("Interpreting all data as drug: ", ab_name(as.ab(ab)), "\n"))
      if ("drug_agent" %in% colnames(ast)) {
        if (length(unique(ast$drug_agent)) > 1) {
          cat("Warning: ignoring 'drug_agent' column in input table, which contains multiple drugs:\n")
          cat(paste(unique(ast$drug_agent), collapse = ", "))
          cat("\n")
        } else if (as.ab(unique(ast$drug_agent)) != as.ab(ab)) {
          cat(paste0("Warning: ignoring 'drug_agent' column in input table, which indicates a different drug: ", unique(ast$drug_agent), "\n"))
        }
      }
      ast <- ast %>% mutate(drug_agent = as.ab(ab))
    } else if (!("drug_agent" %in% colnames(ast))) {
      stop("Warning: Could not interpret data, need to provide 'species' parameter or 'drug_agent' column in input table")
    }

    # interpret data
    if (interpret_eucast) {
      if ("mic" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.mic), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "EUCAST", .names = "pheno_eucast_mic", capped_mic_handling="conservative"))
      }
      if ("disk" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.disk), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "EUCAST", .names = "pheno_eucast_disk"))
      }
      if (("pheno_eucast_mic" %in% colnames(ast)) & ("pheno_eucast_disk" %in% colnames(ast))) {
        ast <- ast %>%
          mutate(pheno_eucast = coalesce(pheno_eucast_mic, pheno_eucast_disk))
      } else if ("pheno_eucast_mic" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_eucast=pheno_eucast_mic)
      } else if ("pheno_eucast_disk" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_eucast=pheno_eucast_disk)
      }
    }
    if (interpret_clsi) {
      if ("mic" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.mic), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "CLSI", .names = "pheno_clsi_mic", capped_mic_handling="conservative"))
      }
      if ("disk" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.disk), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "CLSI", .names = "pheno_clsi_disk"))
      }
      if (("pheno_clsi_mic" %in% colnames(ast)) & ("pheno_clsi_disk" %in% colnames(ast))) {
        ast <- ast %>%
          mutate(pheno_clsi = coalesce(pheno_clsi_mic, pheno_clsi_disk))
      } else if ("pheno_clsi_mic" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_clsi=pheno_clsi_mic)
      } else if ("pheno_eucast_disk" %in% colnames(ast)) {
        ast <- ast %>% rename(pheno_clsi_disk=pheno_clsi_disk)
      }
    }
    if (interpret_ecoff) {
      if ("mic" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.mic), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "EUCAST", breakpoint_type = "ECOFF", .names = "ecoff_mic", capped_mic_handling="conservative"))
      }
      if ("disk" %in% colnames(ast)) {
        ast <- ast %>%
          mutate(across(where(is.disk), as.sir, mo = "spp_pheno", ab = "drug_agent", guideline = "EUCAST", breakpoint_type = "ECOFF", .names = "ecoff_disk"))
      }
      if (("ecoff_mic" %in% colnames(ast)) & ("ecoff_disk" %in% colnames(ast))) {
        ast <- ast %>%
          mutate(ecoff = coalesce(ecoff_mic, ecoff_disk))
      } else if ("ecoff_mic" %in% colnames(ast)) {
        ast <- ast %>% rename(ecoff=ecoff_mic)
      } else if ("ecoff_disk" %in% colnames(ast)) {
        ast <- ast %>% rename(ecoff=ecoff_disk)
      }
    }
  }
  return(ast)
}

#' Import and Process AST Data from an EBI or NCBI antibiogram File
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset in either EBI or NCBI antibiogram format, processes the data, and optionally interprets the results based on MIC or disk diffusion data. It assumes that the input file is a tab-delimited text file (e.g., TSV) or CSV (which may be commpressed) and parses relevant columns (antibiotic names, species names, MIC or disk data) into suitable classes using the AMR package. It optionally can use the AMR package to interpret susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human breakpoints and/or ECOFF). If expected columns are not found warnings will be given, and interpretation may not be possible.
#' @param input A string representing a dataframe, or a path to an input file, containing the AST data in EBI or NCBI antibiogram format. These files can be downloaded from the EBI AMR browser, e.g. https://www.ebi.ac.uk/amr/data/?view=experiments or NCBI browser e.g. https://www.ncbi.nlm.nih.gov/pathogens/ast#Pseudomonas%20aeruginosa.
#' @param format A string indicating the format of the data, either "ebi" (default) or "ncbi". This determines whether the data is passed on to the `import_ebi_ast` or `import_ncbi_ast`()` function to process.
#' @param interpret_eucast A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against ECOFF human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'.
#' @param interpret_clsi A logical value (default is FALSE). If `TRUE`, the function will interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'.
#' @param interpret_ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype).
#' @param species (optional) Name of the species to use for phenotype interpretation. By default, the organism field in the input file will be assumed to specify the species for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single species name via this parameter.
#' @param ab (optional) Name of the antibiotic to use for phenotype interpretation. By default, the antibiotic field in the input file will be assumed to specify the antibiotic for each sample, but if this is missing or you want to override it in the interpretation step, you may provide a single antibiotic name via this parameter.
#' @param source (optional) A single value to record as the source of these data points, e.g. "EBI_browser". By default, the publications field (for EBI data) or BioProject field (for NCBI data) will be used to indicate the source for each row in the input file, but if this is missing or you want to override it with a single value for all samples, you may provide a source name via this parameter.
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of coalesce if_else mutate relocate rename
#' @return A data frame with the processed AST data, including additional columns:
#' - `id`: The biosample identifier.
#' - `spp_pheno`: The species phenotype, formatted using the `as.mo` function.
#' - `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function.
#' - `mic`: The minimum inhibitory concentration (MIC) value, formatted using the `as.mic` function.
#' - `disk`: The disk diffusion measurement (in mm), formatted using the `as.disk` function.
#' - `method`: The AST platform recorded in the input file as the source of the measurement.
#' - `guideline`: The AST standard recorded in the input file as being used for the AST assay.
#' - `pheno_eucast`: The phenotype newly interpreted against EUCAST human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `pheno_clsi`: The phenotype newly interpreted against CLSI human breakpoint standards (as S/I/R), based on the MIC or disk diffusion data.
#' - `ecoff`: The phenotype newly interpreted against the ECOFF (as S/R), based on the MIC or disk diffusion data.
#' - `pheno_provided`: The original phenotype interpretation provided in the input file.
#' - `source`: The source of each data point (from the publications or bioproject field in the input file, or replaced with a single value passed in as the 'source' parameter).
#' @export
#' @examples
#' # small example E. coli AST data from NCBI
#' ecoli_ast_raw
#'
#' # import without re-interpreting resistance
#' pheno <- import_ast(ecoli_ast_raw, format="ncbi")
#' head(pheno)
#'
#' # import and re-interpret resistance (S/I/R) and WT/NWT (vs ECOFF) using AMR package
#' pheno <- import_ast(ecoli_ast_raw, format="ncbi", interpret_eucast = TRUE, interpret_ecoff = TRUE)
#' head(pheno)
import_ast <- function(input, format = "ebi", interpret_eucast = FALSE,
                       interpret_clsi = FALSE, interpret_ecoff = FALSE,
                       species = NULL, ab = NULL, source=NULL) {
  if (format == "ebi") {
    cat("Reading in as EBI AST format\n")
    ast <- import_ebi_ast(input, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, interpret_ecoff = interpret_ecoff, species = species, ab = ab)
  }
  if (format == "ncbi") {
    cat("Reading in as NCBI AST format\n")
    ast <- import_ncbi_ast(input, interpret_eucast = interpret_eucast, interpret_clsi = interpret_clsi, interpret_ecoff = interpret_ecoff, species = species, ab = ab)
  }

  if (!is.null(source)) { ast <- ast %>% mutate(source=source)}
  ast <- ast %>% relocate(any_of(c("id", "drug_agent", "mic", "disk", "pheno_eucast", "pheno_clsi", "ecoff", "guideline", "method", "source", "pheno_provided", "spp_pheno")))

  return(ast)
}
