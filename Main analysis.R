################################################################################
# COMPREHENSIVE MENDELIAN RANDOMIZATION (MR) ANALYSIS PIPELINE
################################################################################
#
# Description:
#   This script performs comprehensive two-sample Mendelian Randomization (MR)
#   analysis, including:
#   - Univariable MR (UVMR) with multiple methods
#   - Two-step mediation analysis using MVMR
#   - Complete sensitivity analyses
#
# Key Features:
#   - Uses chr:pos as SNP identifier (rsID only for LD clumping)
#   - Multiple MR methods: IVW, MR-Egger, Weighted Median, Weighted Mode
#   - MR-PRESSO for outlier detection and correction
#   - Sensitivity analyses: heterogeneity, pleiotropy, instrument strength
#   - Allows mediators without rsID
#
# Instrument Variable Selection (STRICT):
#   - P-value threshold: 5e-8 (genome-wide significance)
#   - LD clumping: r2 < 0.001, window = 10000 kb
#   - Minimum SNPs: 3
#
# Reference:
#   Based on standards from Nature Human Behaviour 2024
#   (Ye et al. Mental well-being and healthy aging)
#   https://github.com/yechaojie/mental_aging
#
################################################################################

# ==== REQUIRED PACKAGES ====
# Installation command (run if needed):
# install.packages(c("data.table","readr","tibble","TwoSampleMR","ieugwasr","MVMR","MRPRESSO","MRlap"))
# Note: MVMR (version 0.4) and MRlap (version 0.0.3.0) may need special installation
# devtools::install_github("WSpiller/MVMR")  # If not on CRAN
# devtools::install_github("n-mounier/MRlap") # If not on CRAN

library(data.table)           # Fast data manipulation
library(readr)                # Read TSV/CSV files
library(tibble)               # Data frames with better defaults
library(TwoSampleMR)          # Main UVMR analysis framework (version 0.5.7)
library(ieugwasr)             # LD clumping utilities

# MVMR package (version 0.4) - REQUIRED for multivariable MR
# Reference: Ye et al. Nat Hum Behav 2024 used this package
suppressWarnings(try(library(MVMR), silent = TRUE))
MVMR_AVAILABLE <- "MVMR" %in% loadedNamespaces()

# MR-PRESSO (version 1.0) - Outlier detection
suppressWarnings(try(library(MRPRESSO), silent = TRUE))
MRPRESSO_AVAILABLE <- "MRPRESSO" %in% loadedNamespaces()

# MRlap (version 0.0.3.0) - Sample overlap correction
suppressWarnings(try(library(MRlap), silent = TRUE))
MRLAP_AVAILABLE <- "MRlap" %in% loadedNamespaces()

# Report package availability
cat("\n[PACKAGE CHECK]\n")
cat("  TwoSampleMR: ✓ Loaded\n")
cat("  MVMR:        ", ifelse(MVMR_AVAILABLE, "✓ Loaded", "✗ NOT loaded - MVMR analysis will be skipped!"), "\n", sep="")
cat("  MRPRESSO:    ", ifelse(MRPRESSO_AVAILABLE, "✓ Loaded", "✗ NOT loaded - MR-PRESSO will be skipped"), "\n", sep="")
cat("  MRlap:       ", ifelse(MRLAP_AVAILABLE, "✓ Loaded", "✗ NOT loaded - MRlap will be skipped"), "\n", sep="")

if(!MVMR_AVAILABLE){
  message("\n[CRITICAL WARNING] MVMR package not loaded!")
  message("  MVMR analysis requires the MVMR package (version 0.4)")
  message("  Install with: devtools::install_github('WSpiller/MVMR')")
  message("  Or: install.packages('MVMR')")
}

cat("\n")

################################################################################
# STEP 1: CONFIGURATION - Paths and Parameters
################################################################################
# 
# Modify BASE_DIR, PLINK_BIN, and EUR_BFILE according to your local setup
# 
################################################################################

# Base directory for the project
BASE_DIR <- "D:/Projects_data&code/MR_pipeline_demo"

# Exposure GWAS data directories
DIR_EXPO01 <- file.path(BASE_DIR, "Standardized Circulating human plasma proteome_Data")
DIR_EXPO02 <- file.path(BASE_DIR, "Standardized Circulating metabolic biomarkers_Data")
DIR_EXPO03 <- file.path(BASE_DIR, "Circulating inflammatory proteins_Data")

# Mediator GWAS data directory
DIR_MEDI   <- file.path(BASE_DIR, "Cerebrospinal fluid metabolomics_Data")

# Covariates/Confounders GWAS data directory (SES indicators)
# For controlling confounding in MVMR analyses
DIR_COVAR  <- file.path(BASE_DIR, "Covariates_SES")
# Expected covariates:
#   - Education: Years of schooling 
#   - Income: Household income 
#   - Occupation: Occupational attainment 

# Outcome GWAS data directory
DIR_OUT    <- file.path(BASE_DIR, "Outcomes")

# Output directory for results
DIR_RES_TRIAL <- file.path(BASE_DIR, "results_trial")
if(!dir.exists(DIR_RES_TRIAL)) dir.create(DIR_RES_TRIAL, recursive = TRUE)

# PLINK executable for local LD clumping
PLINK_BIN <- "D:/Projects_data&code/MR_pipeline_demo/local_clump/plink/plink.exe"

# 1000 Genomes EUR reference panel for LD clumping
# Should point to the file prefix (without .bed/.bim/.fam extension)
EUR_BFILE <- "D:/Projects_data&code/MR_pipeline_demo/local_clump/EUR/1000G.EUR.QC"

# Verify PLINK executable exists
stopifnot(file.exists(PLINK_BIN))

# Auto-detect local 1000G EUR reference panel; fallback to remote clumping if not found
EUR_DIR <- dirname(EUR_BFILE)
USE_LOCAL_CLUMP <- TRUE

# Helper function to check if all PLINK binary files (.bed, .bim, .fam) exist
has_all <- function(prefix) all(file.exists(paste0(prefix, c(".bed",".bim",".fam"))))

# Check and auto-detect reference panel
if(!has_all(EUR_BFILE)){
  message("[INFO] Specified EUR_BFILE incomplete. Attempting auto-detection...")
  cand <- sub("\\.bed$", "", list.files(EUR_DIR, "\\.bed$", full.names = TRUE))
  ok   <- cand[vapply(cand, has_all, logical(1))]
  if(length(ok) >= 1){ 
    EUR_BFILE <- ok[1]
    message("[INFO] EUR_BFILE selected: ", EUR_BFILE)
} else {
    message("[WARNING] No complete .bed/.bim/.fam set found. Will use remote clumping.")
    USE_LOCAL_CLUMP <- FALSE
  }
} else {
  message("[OK] Local reference panel ready: ", EUR_BFILE)
}

# ---- Instrument Variable Selection Thresholds (STRICT) ----
# These are stringent thresholds following genome-wide standards
P_THRESH <- 5e-8      # P-value threshold (genome-wide significance)
CLUMP_R2 <- 0.001     # LD r² threshold for clumping (very strict)
CLUMP_KB <- 10000     # LD window size in kb

# ---- MRlap Configuration ----
# MRlap automatically handles:
# - Sample size extraction from GWAS files
# - LD score regression for genetic parameters
# - Estimation of sample overlap effects
# - Correction for winner's curse and weak instruments
#
# No manual configuration needed!
# MRlap will process full GWAS summary statistics automatically

################################################################################
# STEP 2: FILE DISCOVERY - Scan for GWAS Data Files
################################################################################
#
# Functions to locate and report GWAS summary statistics files
# Supports .tsv, .tsv.gz, and .gz.tsv formats
#
################################################################################

#' List GWAS Files in Directory
#'
#' Supports multiple GWAS file formats: .tsv, .tsv.gz, .txt.gz, .txt
#'
#' @param dir_path Path to directory containing GWAS files
#' @return Character vector of full file paths
list_gwas_files <- function(dir_path){
  if(!dir.exists(dir_path)) return(character())
  list.files(dir_path, pattern="(\\.tsv$)|(\\.tsv\\.gz$)|(\\.gz\\.tsv$)|(\\.txt\\.gz$)|(\\.txt$)|(\\.gz\\.txt$)",
             full.names=TRUE, ignore.case=TRUE, recursive=FALSE)
}

#' Report Files in Directory
#'
#' @param tag Label for the directory type (e.g., "EXPO01", "MEDI")
#' @param path Path to the directory
#' @return Invisibly returns vector of file paths
report_dir <- function(tag, path){
  if(!dir.exists(path)){ 
    message("[MISSING DIRECTORY] ", tag, ": ", path)
    return(invisible(character()))
  }
  files <- list_gwas_files(path)
  message(sprintf("[DIRECTORY SCAN] %s | Files found = %d", tag, length(files)))
  if(length(files)) print(head(basename(files), 5))
  invisible(files)
}

# Scan all exposure, mediator, and covariate directories
invisible(report_dir("EXPO01", DIR_EXPO01))
invisible(report_dir("EXPO02", DIR_EXPO02))
invisible(report_dir("EXPO03", DIR_EXPO03))
invisible(report_dir("MEDI"  , DIR_MEDI   ))
invisible(report_dir("COVAR" , DIR_COVAR  ))  # Covariates (SES indicators)

#' Get Outcome GWAS Files
#'
#' Retrieves outcome GWAS files, optionally filtering for specific outcomes
#' Supports .tsv, .tsv.gz, .txt, .txt.gz formats
#'
#' @param dir_out Path to outcomes directory
#' @return Character vector of outcome file paths
get_outcomes <- function(dir_out){
  if(!dir.exists(dir_out)) stop("Outcome directory does not exist: ", dir_out)
  
  # Specific outcome files to prioritize (if present)
  wanted <- c("MDD37.tsv","LBD38.tsv","FD37.tsv","CP37.tsv","AD37.tsv",
              "UD38.tsv","VD38.tsv","D38.tsv","MD37.tsv","DD38.tsv")
  
  # Find all GWAS files in outcomes directory (TSV and TXT formats)
  cand <- list.files(dir_out, pattern="(\\.tsv$)|(\\.tsv\\.gz$)|(\\.gz\\.tsv$)|(\\.txt$)|(\\.txt\\.gz$)|(\\.gz\\.txt$)",
                     full.names=TRUE, ignore.case=TRUE, recursive=TRUE)
  
  # Filter for wanted files if they exist, otherwise use all files
  keep <- basename(cand) %in% wanted | basename(cand) %in% paste0(wanted, ".gz")
  out  <- if(any(keep)) cand[keep] else cand
  unique(out)
}

# Collect all exposure, mediator, and outcome files
all_expos_files <- unique(c(list_gwas_files(DIR_EXPO01),
                            list_gwas_files(DIR_EXPO02),
                            list_gwas_files(DIR_EXPO03)))
all_medi_files  <- unique(list_gwas_files(DIR_MEDI))
all_out_files   <- unique(get_outcomes(DIR_OUT))

# Covariates files (for MVMR confounder adjustment)
all_covar_files <- if(dir.exists(DIR_COVAR)) {
  unique(list_gwas_files(DIR_COVAR))
} else {
  character(0)
}

# Verify that we have data files for required categories
stopifnot(length(all_expos_files)>0, length(all_medi_files)>0, length(all_out_files)>0)

# Report covariates
if(length(all_covar_files) > 0){
  message("[INFO] Covariates found: ", length(all_covar_files), " files")
  message("  These will be used for confounder adjustment in MVMR")
} else {
  message("[WARNING] No covariates found in DIR_COVAR")
  message("  MVMR will be performed WITHOUT covariate adjustment")
  message("  For proper analysis, add SES covariates (Education, Income, Occupation)")
}

################################################################################
# STEP 3: DATA STANDARDIZATION - Column Mapping and Data Harmonization
################################################################################
#
# This section handles flexible column naming across different GWAS formats
# Supports various naming conventions and derives missing columns when possible
# 
# Key features:
#   - Automatic column name detection and mapping
#   - Derives beta from OR/HR when needed
#   - Derives p-value from Z-score when needed
#   - Derives SE from beta and p-value/Z-score when needed
#   - Allows missing rsID (uses chr:pos instead)
#
################################################################################

# Column name mapping dictionary
# Maps standardized column names to common alternatives found in GWAS files
# Updated to support covariate file formats (Education, Income, Occupation)
COLMAP <- list(
  SNP            = c("SNP","rsid","rsID","variant","marker","rs_number","ID","variant_id","rs_id","MarkerName"),
  effect_allele  = c("effect_allele","EA","A1","a1","alt","Allele1","EFFECT_ALLELE"),
  other_allele   = c("other_allele","OA","A2","a2","ref","Allele2","OTHER_ALLELE"),
  beta           = c("beta","BETA","b","Beta","logOR","lnOR","effect","EFFECT","BETA_COEFF"),
  se             = c("se","SE","sebeta","stderr","standard_error","se_beta","SE_EFFECT","SE"),
  pval           = c("p","P","pval","pvalue","PVAL","p_value","P_VALUE","Pval","LOG.P."),
  eaf            = c("eaf","eaf1","freq","FRQ","maf","ALT_AF","AF","effect_allele_frequency","EAF"),
  chr            = c("chr","CHR","chrom","chromosome","Chromosome"),
  pos            = c("pos","BP","bp","position","base_pair_location","POS"),
  samplesize     = c("n","N","n_total","samplesize","total_sample_size","Neff","Neff_total"),
  or             = c("OR","odds_ratio","oddsratio"),
  hr             = c("HR","hazard_ratio","hazardratio"),
  z              = c("z","Z","zscore","Zscore","z_score","ZSCORE")
)

#' Match First Available Column Name
#'
#' Finds the first matching column name from a list of alternatives (case-insensitive)
#'
#' @param nms Character vector of actual column names in data
#' @param want Character vector of possible names to match
#' @return First matching column name, or NA if none found
match_first <- function(nms, want){
  idx <- match(tolower(want), tolower(nms))
  if(!all(is.na(idx))){ m <- idx[!is.na(idx)]; if(length(m)>0) return(nms[m[1]]) }
  NA_character_
}

#' Convert to Numeric
#'
#' Safely converts a vector to numeric, suppressing warnings
#'
#' @param x Vector to convert
#' @return Numeric vector
to_num <- function(x){ if(is.numeric(x)) x else suppressWarnings(as.numeric(x)) }

#' Clean File Name
#'
#' Removes common GWAS file extensions to get clean phenotype name
#' Supports: .tsv, .tsv.gz, .txt, .txt.gz, .gz
#'
#' @param filename File name with extension
#' @return Clean phenotype name without extension
clean_filename <- function(filename){
  # Remove multiple possible extensions
  filename <- sub("\\.tsv\\.gz$", "", filename, ignore.case = TRUE)
  filename <- sub("\\.txt\\.gz$", "", filename, ignore.case = TRUE)
  filename <- sub("\\.gz\\.tsv$", "", filename, ignore.case = TRUE)
  filename <- sub("\\.gz\\.txt$", "", filename, ignore.case = TRUE)
  filename <- sub("\\.tsv$", "", filename, ignore.case = TRUE)
  filename <- sub("\\.txt$", "", filename, ignore.case = TRUE)
  filename <- sub("\\.gz$", "", filename, ignore.case = TRUE)
  filename
}

#' Standardize Column Names and Derive Missing Columns
#'
#' Main function to standardize GWAS data columns
#' - Maps various column names to standard names
#' - Derives beta from OR/HR if beta is missing
#' - Derives p-value from Z-score if needed
#' - Derives SE from beta and p-value/Z-score if needed
#' - Creates chr:pos identifier for SNP matching
#'
#' @param dt data.table containing GWAS summary statistics
#' @return Standardized data.table with required columns
standardize_cols <- function(dt){
  nm <- names(dt)
  found <- lapply(COLMAP, function(v) match_first(nm, v))
  names(found) <- names(COLMAP)
  
  # Rename columns to standard names (if they exist)
  rename_map <- found[!is.na(found)]
  data.table::setnames(dt, old = unlist(rename_map), new = names(rename_map), skip_absent = TRUE)
  
  # ---- Derive beta from OR or HR if beta is missing ----
  if(!"beta" %in% names(dt)){
    if("or" %in% names(dt)){
      dt[, beta := log(to_num(or))]  # log(OR) = beta
    } else if("hr" %in% names(dt)){
      dt[, beta := log(to_num(hr))]  # log(HR) = beta
    }
  }
  
  # ---- Derive p-value from Z-score if p-value is missing ----
  if(!"pval" %in% names(dt) && "z" %in% names(dt)){
    znum <- to_num(dt$z)
    dt[, pval := 2*pnorm(abs(znum), lower.tail=FALSE)]  # Two-tailed p-value
  }
  
  # ---- Derive standard error (SE) if missing ----
  if(!"se" %in% names(dt) && "beta" %in% names(dt)){
    if("z" %in% names(dt)){
      # SE = beta / Z
      znum <- to_num(dt$z)
      dt[, se := to_num(beta)/znum]
    } else if("pval" %in% names(dt)){
      # SE = |beta| / |Z|, where Z = qnorm(p/2)
      pnum <- pmin(pmax(to_num(dt$pval), .Machine$double.xmin), 1)
      znum <- qnorm(pnum/2, lower.tail = FALSE)
      dt[, se := abs(to_num(beta))/abs(znum)]
    }
  }
  
  # ---- Check for required columns ----
  # Note: SNP column can be missing if chr and pos are available
  need <- c("effect_allele","other_allele","beta","se","pval")
  if(!all(need %in% names(dt))){
    miss <- setdiff(need, names(dt))
    stop("Required columns missing: ", paste(miss, collapse=", "))
  }
  
  # ---- Convert to numeric and standardize ----
  num_cols <- intersect(c("beta","se","pval","eaf","pos","samplesize","chr","z"), names(dt))
  for(cc in num_cols) dt[[cc]] <- to_num(dt[[cc]])
  
  # Handle LOG.P. column if present (convert to regular p-value)
  if("pval" %in% names(dt) && anyNA(dt$pval) && "LOG.P." %in% names(dt)){
    dt[is.na(pval) & !is.na(`LOG.P.`), pval := 10^(-to_num(`LOG.P.`))]
  }
  
  # Standardize allele names to uppercase
  if("effect_allele" %in% names(dt)) dt[, effect_allele := toupper(effect_allele)]
  if("other_allele"  %in% names(dt)) dt[, other_allele  := toupper(other_allele)]
  
  # ---- Create chr:pos identifier if chr and pos are available ----
  # This provides a robust SNP identifier even when rsID is missing
  if(all(c("chr","pos") %in% names(dt))){
    dt[, chr := as.integer(chr)]
    dt[, pos := as.integer(pos)]
    dt[, SNP_cp := paste0(chr, ":", pos)]  # SNP_cp = "chromosome:position"
  }
  
  # Remove duplicates
  unique(dt)
}

#' Read and Standardize GWAS Summary Statistics
#'
#' Reads a GWAS file (TSV or TSV.GZ) and standardizes column names
#' Tries multiple reading strategies to handle different file formats
#'
#' @param file Path to GWAS summary statistics file
#' @return Standardized data.table
read_gwas <- function(file){
  # Try fast reading with fread first
  dt <- try(fread(file, sep="\t", na.strings=c("NA","NaN",""), showProgress=FALSE,
                  nThread = max(1, parallel::detectCores()-1), fill=TRUE), silent = TRUE)
  
  # Fallback to read_tsv if fread fails or produces single column
  if(inherits(dt,"try-error") || (is.data.frame(dt) && ncol(dt) <= 1)){
    dt2 <- try(read_tsv(file, comment="#",
                        col_types = cols(.default = col_character()),
                        progress = FALSE, guess_max = 100000), silent = TRUE)
    if(inherits(dt2,"try-error") || ncol(dt2)==0){
      dt2 <- read_tsv(file, col_types = cols(.default = col_character()),
                      progress = FALSE, guess_max = 100000)
    }
    dt <- as.data.table(dt2)
  }
  
  # Standardize column names and derive missing columns
  standardize_cols(dt)
}

#' Use chr:pos as SNP Identifier
#'
#' Replaces SNP column with chr:pos identifier (if available)
#' Preserves original rsID in SNP_rsid column
#' This approach is more robust than rsID for matching across datasets
#'
#' @param dt data.table with GWAS data
#' @return data.table with SNP set to chr:pos
use_chrpos_id <- function(dt){
  dt <- copy(dt)
  
  if("SNP_cp" %in% names(dt)){
    # Save original rsID if it exists
    if("SNP" %in% names(dt)){
      setnames(dt, "SNP", "SNP_rsid", skip_absent = TRUE)
    }
    # Use chr:pos as SNP identifier
    setnames(dt, "SNP_cp", "SNP", skip_absent = TRUE)
  } else if(!"SNP" %in% names(dt)){
    stop("Cannot create chr:pos identifier (missing chr/pos), and no SNP column exists.")
  }
  
  dt
}

################################################################################
# STEP 4: INSTRUMENT VARIABLE SELECTION - LD Clumping
################################################################################
#
# Selects independent genetic instruments using LD clumping
# Requires rsID for clumping (exposure GWAS only)
# Uses local PLINK or remote API for LD calculation
#
################################################################################

#' Select Independent Genetic Instruments via LD Clumping
#'
#' Performs LD clumping to select independent SNPs as instrumental variables
#' - Filters SNPs by p-value threshold
#' - Removes SNPs in LD (r² threshold)
#' - Tries local clumping first, falls back to remote if needed
#'
#' @param exp_dt data.table with exposure GWAS data
#' @param p_thresh P-value threshold for significance (default: 5e-8)
#' @param r2 LD r² threshold for clumping (default: 0.001)
#' @param kb LD window size in kb (default: 10000)
#' @return data.table with independent SNPs (clumped)
select_instruments <- function(exp_dt, p_thresh = P_THRESH, r2=CLUMP_R2, kb=CLUMP_KB){
  
  # Check if exposure has rsID (required for LD clumping)
  if(!any(grepl("^rs\\d+$", exp_dt$SNP))){
    message("  [WARNING] Exposure lacks rsID (SNP not 'rs...'). Cannot perform LD clumping. Returning empty.")
    return(data.table())
  }
  
  # Filter for genome-wide significant SNPs
  sig <- exp_dt[pval < p_thresh, .(rsid = SNP, pval)]
  if(nrow(sig) == 0) return(data.table())
  
  # Try local LD clumping first (faster and more reliable)
  if (isTRUE(USE_LOCAL_CLUMP)) {
    res_local <- try(ieugwasr::ld_clump_local(
      tibble::as_tibble(transform(sig, id="exposure_gwas")),
      bfile = EUR_BFILE,
      plink_bin = PLINK_BIN,
      clump_kb = kb,
      clump_r2 = r2,
      clump_p = p_thresh
    ), silent = TRUE)
    
    if(!inherits(res_local, "try-error")){
      res_local <- as.data.table(res_local)
      if(nrow(res_local)>0) return(res_local)
    }
    message("  [INFO] Local clumping returned no results. Trying remote clumping...")
  } else {
    message("  [INFO] Local clumping disabled. Trying remote clumping...")
  }
  
  # Fallback to remote LD clumping via API
  res_remote <- try(ieugwasr::ld_clump(
    d = tibble::as_tibble(transform(sig, id="exposure_gwas")),
    clump_kb = kb,
    clump_r2 = r2,
    clump_p = p_thresh,
    pop = "EUR"
  ), silent = TRUE)
  
  if(!inherits(res_remote, "try-error")) return(as.data.table(res_remote))
  
  # Both methods failed
  message("  [WARNING] Remote clumping also failed. Returning empty.")
  data.table()
}

################################################################################
# STEP 5: DATA FORMATTING FOR TWO-SAMPLE MR
################################################################################
#
# Formats data for TwoSampleMR package and harmonizes exposure-outcome data
# Only passes optional columns if they exist in the data
#
################################################################################

#' Format Exposure Data for TwoSampleMR
#'
#' Converts standardized exposure data to TwoSampleMR format
#' Dynamically includes optional columns (eaf, chr, pos, samplesize) if present
#'
#' @param exp_dt data.table with exposure GWAS data
#' @param exposure_label Label for the exposure
#' @return Formatted data frame ready for MR analysis
to_exposure_format <- function(exp_dt, exposure_label){
  exp_df <- as.data.frame(exp_dt)
  
  # Required arguments
  args <- list(
    dat = exp_df,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval"
  )
  
  # Add optional columns if they exist
  if ("eaf" %in% names(exp_df))        args$eaf_col        <- "eaf"
  if ("chr" %in% names(exp_df))        args$chr_col        <- "chr"
  if ("pos" %in% names(exp_df))        args$pos_col        <- "pos"
  if ("samplesize" %in% names(exp_df)) args$samplesize_col <- "samplesize"
  
  # Format and label
  res <- do.call(TwoSampleMR::format_data, args)
  res$exposure <- exposure_label
  res
}

#' Format Outcome Data for TwoSampleMR
#'
#' Converts standardized outcome data to TwoSampleMR format
#' Dynamically includes optional columns (eaf, chr, pos, samplesize) if present
#'
#' @param out_dt data.table with outcome GWAS data
#' @param outcome_label Label for the outcome
#' @return Formatted data frame ready for MR analysis
to_outcome_format <- function(out_dt, outcome_label){
  out_df <- as.data.frame(out_dt)
  
  # Required arguments
  args <- list(
    dat = out_df,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval"
  )
  
  # Add optional columns if they exist
  if ("eaf" %in% names(out_df))        args$eaf_col        <- "eaf"
  if ("chr" %in% names(out_df))        args$chr_col        <- "chr"
  if ("pos" %in% names(out_df))        args$pos_col        <- "pos"
  if ("samplesize" %in% names(out_df)) args$samplesize_col <- "samplesize"
  
  # Format and label
  res <- do.call(TwoSampleMR::format_data, args)
  res$outcome <- outcome_label
  res
}

#' Harmonize Exposure and Outcome Data
#'
#' Aligns SNP effects to the same allele across exposure and outcome
#' Uses TwoSampleMR::harmonise_data with action=2 (default harmonization)
#'
#' @param exp_fmt Formatted exposure data
#' @param out_fmt Formatted outcome data
#' @return Harmonized data frame
harmonise_xy <- function(exp_fmt, out_fmt){
  TwoSampleMR::harmonise_data(exp_fmt, out_fmt, action = 2)
}

################################################################################
# STEP 6: COMPREHENSIVE UNIVARIABLE MR (UVMR) ANALYSIS
################################################################################
#
# Performs complete MR analysis with multiple methods and sensitivity analyses
# Includes:
#   - Multiple MR methods (IVW, MR-Egger, Weighted Median, Weighted Mode)
#   - Heterogeneity testing (Cochran's Q)
#   - Pleiotropy testing (MR-Egger intercept)
#   - MR-PRESSO outlier detection (if available)
#   - Instrument strength assessment (F-statistic)
#
################################################################################

#' Calculate Mean F-Statistic
#'
#' Calculates the mean F-statistic for instrument strength
#' F = (beta / SE)^2
#'
#' @param hdat Harmonized data frame with beta.exposure and se.exposure columns
#' @return Mean F-statistic across all SNPs
calc_mean_F <- function(hdat){
  if(!all(c("beta.exposure","se.exposure") %in% names(hdat))) return(NA_real_)
  mean((hdat$beta.exposure^2)/(hdat$se.exposure^2), na.rm=TRUE)
}

#' Run Comprehensive UVMR Analysis
#'
#' Performs two-sample MR with multiple methods and complete sensitivity analyses
#' This is the main UVMR analysis function
#'
#' @param hdat Harmonized data frame (exposure + outcome)
#' @param exposure_name Name of exposure
#' @param outcome_name Name of outcome
#' @return data.table with results from all methods and sensitivity tests (NULL if <3 SNPs)
run_uvmr_comprehensive <- function(hdat, exposure_name, outcome_name){
  # Minimum 3 SNPs required for MR
  if(nrow(hdat) < 3) return(NULL)
  
  # ---- Define MR methods to use ----
  methods_to_use <- c(
    "mr_ivw",                # Inverse Variance Weighted (main method)
    "mr_egger_regression",   # MR-Egger (detects pleiotropy)
    "mr_weighted_median",    # Weighted Median (robust to outliers)
    "mr_weighted_mode",      # Weighted Mode
    "mr_simple_mode"         # Simple Mode
  )
  
  # ---- Main MR Analysis ----
  main <- tryCatch(
    TwoSampleMR::mr(hdat, method_list = methods_to_use), 
    error = function(e) NULL
  )
  
  # ---- Heterogeneity Test (Cochran's Q) ----
  # Tests whether SNP estimates are more variable than expected by chance
  het <- tryCatch(
    TwoSampleMR::mr_heterogeneity(hdat, method_list = c("mr_ivw", "mr_egger_regression")), 
    error = function(e) NULL
  )
  
  # ---- Pleiotropy Test (MR-Egger Intercept) ----
  # Tests for directional pleiotropy
  # Significant intercept suggests horizontal pleiotropy
  ple <- tryCatch(
    TwoSampleMR::mr_pleiotropy_test(hdat), 
    error = function(e) NULL
  )
  
  # ---- Leave-One-Out Analysis ----
  # Tests if results driven by single SNP
  # (Prepared but not included in final output for brevity)
  loo <- tryCatch(
    TwoSampleMR::mr_leaveoneout(hdat, method = mr_ivw), 
    error = function(e) NULL
  )
  
  # ---- MR-PRESSO: Outlier Detection and Correction ----
  # Detects and corrects for horizontal pleiotropy outliers
  # Requires at least 4 SNPs and MRPRESSO package
  presso_res <- NULL
  presso_global_p <- NA
  presso_outliers <- NA
  presso_corrected_b <- NA
  presso_corrected_se <- NA
  presso_corrected_p <- NA
  
  if(MRPRESSO_AVAILABLE && nrow(hdat) >= 4){
    # Prepare input data for MR-PRESSO
    presso_input <- tryCatch({
      data.frame(
        SNP = hdat$SNP,
        beta.exposure = hdat$beta.exposure,
        beta.outcome = hdat$beta.outcome,
        se.exposure = hdat$se.exposure,
        se.outcome = hdat$se.outcome
      )
    }, error = function(e) NULL)
    
    if(!is.null(presso_input)){
      # Run MR-PRESSO
      # - OUTLIERtest: identify outlier SNPs
      # - DISTORTIONtest: test if removing outliers changes estimate
      # - NbDistribution: number of bootstrap iterations
      presso_res <- tryCatch(
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE,
          data = presso_input,
          NbDistribution = 1000,      # Bootstrap iterations (may be slow)
          SignifThreshold = 0.05
        ),
        error = function(e) NULL
      )
      
      # Extract MR-PRESSO results
      if(!is.null(presso_res)){
        # Global test p-value (overall outlier detection)
        presso_global_p <- presso_res$`MR-PRESSO results`$`Global Test`$Pvalue
        
        # Identified outlier SNPs
        if(!is.null(presso_res$`MR-PRESSO results`$`Outlier Test`)){
          presso_outliers <- paste(presso_res$`MR-PRESSO results`$`Outlier Test`$SNP, collapse = ";")
        }
        
        # Outlier-corrected estimates
        if(!is.null(presso_res$`Main MR results`)){
          presso_corrected_b <- presso_res$`Main MR results`$`Causal Estimate`[presso_res$`Main MR results`$`MR Analysis` == "Outlier-corrected"]
          presso_corrected_se <- presso_res$`Main MR results`$Sd[presso_res$`Main MR results`$`MR Analysis` == "Outlier-corrected"]
          presso_corrected_p <- presso_res$`Main MR results`$`P-value`[presso_res$`Main MR results`$`MR Analysis` == "Outlier-corrected"]
        }
      }
    }
  }
  
  # ---- Calculate Instrument Strength (F-statistic) ----
  F_mean <- calc_mean_F(hdat)
  
  # ---- Organize Results ----
  # Create one row per MR method with all sensitivity statistics
  results_list <- list()
  
  if(!is.null(main) && nrow(main) > 0){
    # Loop through each MR method
    for(i in 1:nrow(main)){
      method_name <- main$method[i]
      
      # ---- Extract heterogeneity statistics for this method ----
      Q_val <- NA
      Q_pval <- NA
      if(!is.null(het)){
        # Match heterogeneity test to method (IVW or Egger)
        if(grepl("Inverse variance weighted", method_name, ignore.case = TRUE)){
          het_row <- het[grepl("Inverse variance weighted", het$method, ignore.case = TRUE), ]
          if(nrow(het_row) > 0){
            Q_val <- het_row$Q[1]
            Q_pval <- het_row$Q_pval[1]
          }
        } else if(grepl("Egger", method_name, ignore.case = TRUE)){
          het_row <- het[grepl("Egger", het$method, ignore.case = TRUE), ]
          if(nrow(het_row) > 0){
            Q_val <- het_row$Q[1]
            Q_pval <- het_row$Q_pval[1]
          }
        }
      }
      
      # ---- Extract Egger intercept (pleiotropy test) ----
      egger_intercept <- NA
      egger_intercept_se <- NA
      egger_intercept_p <- NA
      if(!is.null(ple) && nrow(ple) > 0){
        egger_intercept <- ple$egger_intercept[1]
        egger_intercept_se <- ple$se[1]
        egger_intercept_p <- ple$pval[1]
      }
      
      # ---- Compile results for this method ----
      result_row <- data.table(
        exposure = exposure_name,
        outcome = outcome_name,
        method = method_name,
        nSNP = main$nsnp[i],                 # Number of SNPs
        b = main$b[i],                       # Beta (causal estimate)
        se = main$se[i],                     # Standard error
        pval = main$pval[i],                 # P-value
        Q = Q_val,                           # Cochran's Q statistic
        Q_pval = Q_pval,                     # Heterogeneity p-value
        egger_intercept = egger_intercept,   # MR-Egger intercept
        egger_intercept_se = egger_intercept_se,
        egger_intercept_pval = egger_intercept_p,
        F_statistic = F_mean,                # Mean F-statistic
        presso_global_pval = presso_global_p,  # MR-PRESSO global test
        presso_outliers = presso_outliers,     # Outlier SNPs
        presso_corrected_b = presso_corrected_b,  # Outlier-corrected beta
        presso_corrected_se = presso_corrected_se,
        presso_corrected_pval = presso_corrected_p
      )
      
      results_list[[length(results_list) + 1]] <- result_row
    }
  }
  
  # ---- Return combined results ----
  if(length(results_list) > 0){
    return(rbindlist(results_list))
  } else {
    return(NULL)
  }
}

################################################################################
# STEP 7: MAIN UVMR WORKFLOW - Analyze All Exposure-Outcome Pairs
################################################################################
#
# This function orchestrates the complete UVMR analysis across all combinations
# of exposures and outcomes
#
# Workflow for each exposure-outcome pair:
#   1. Read exposure GWAS
#   2. Select instruments (LD clumping at p<5e-8)
#   3. Read outcome GWAS
#   4. Harmonize exposure-outcome data
#   5. Run comprehensive MR analysis (multiple methods + sensitivity)
#   6. Save results
#
################################################################################

MIN_SNPS <- 3  # Minimum number of SNPs required for MR

#' Run Comprehensive UVMR Analysis Across All Pairs
#'
#' Main function to perform UVMR for all exposure-outcome combinations
#' Saves complete results to CSV file
#'
#' @param expos_files Vector of exposure GWAS file paths
#' @param outcomes_files Vector of outcome GWAS file paths
#' @param out_csv Output CSV file path
#' @return data.table with all UVMR results (invisibly)
run_comprehensive_uvmr <- function(
    expos_files = all_expos_files,
    outcomes_files = all_out_files,
    out_csv = file.path(DIR_RES_TRIAL, "uvmr_comprehensive_results.csv")
){
  message("\n[COMPREHENSIVE UVMR ANALYSIS] Exposures=", length(expos_files), "; Outcomes=", length(outcomes_files),
          "; p<", P_THRESH, "; r2=", CLUMP_R2, "; kb=", CLUMP_KB)
  
  all_results <- list()
  total_pairs <- 0
  successful_pairs <- 0
  
  for(exp_file in expos_files){
    exp_name <- basename(exp_file)
    exp_name <- clean_filename(exp_name)
    message("\n[Exposure ", match(exp_file, expos_files), "/", length(expos_files), "] ", exp_name)
    
    exp_dt <- try(read_gwas(exp_file), silent = TRUE)
    if(inherits(exp_dt, "try-error")) { 
      message("  暴露读取失败，跳过。")
      next 
    }
    
    # 选择工具变量
    ivs <- select_instruments(exp_dt, P_THRESH, CLUMP_R2, CLUMP_KB)
    if(nrow(ivs) < MIN_SNPS){ 
      message("  IV 数=", nrow(ivs), " < ", MIN_SNPS, "，跳过该暴露。")
      next 
    }
    
    message("  选择了 ", nrow(ivs), " 个独立工具变量")
    
    exp_iv <- exp_dt[SNP %in% ivs$rsid]
    # 用 chr:pos 作为 SNP 键对齐
    exp_iv <- use_chrpos_id(exp_iv)
    
    exp_fmt <- to_exposure_format(exp_iv, exp_name)
    
    for(out_file in outcomes_files){
      total_pairs <- total_pairs + 1
      out_name <- basename(out_file)
      out_name <- clean_filename(out_name)
      
      out_dt <- try(read_gwas(out_file), silent = TRUE)
      if(inherits(out_dt, "try-error")) { 
        message("    • ", out_name, "：读取失败，跳过。")
        next 
      }
      out_dt <- use_chrpos_id(out_dt)
      
      out_fmt <- to_outcome_format(out_dt, out_name)
      
      # Harmonize data
      hdat <- harmonise_xy(exp_fmt, out_fmt)
      hdat <- hdat[hdat$remove==FALSE & is.finite(hdat$beta.exposure) & is.finite(hdat$beta.outcome), ]
      
      message("    • ", out_name, ": Harmonized SNPs = ", nrow(hdat))
      
      if(nrow(hdat) >= MIN_SNPS){
        # Run comprehensive MR analysis
        res <- run_uvmr_comprehensive(hdat, exp_name, out_name)
        
        if(!is.null(res) && nrow(res) > 0){
          # Add extra information
          res[, `:=`(
            p_threshold = P_THRESH,
            clump_r2 = CLUMP_R2,
            clump_kb = CLUMP_KB,
            total_ivs = nrow(ivs),
            harmonised_snps = nrow(hdat)
          )]
          
          # Run MRlap analysis for sample overlap correction
          # MRlap uses FULL GWAS files (not just harmonized IVs)
          if(MRLAP_AVAILABLE){
            mrlap_res <- run_mrlap_analysis(exp_file, out_file, exp_name, out_name)
            
            if(!is.null(mrlap_res) && nrow(mrlap_res) > 0){
              # Merge MRlap results with main results
              # Add MRlap columns to all rows (same for all methods)
              res[, `:=`(
                mrlap_status = mrlap_res$mrlap_status[1],
                mrlap_observed_effect = mrlap_res$mrlap_observed_effect[1],
                mrlap_observed_se = mrlap_res$mrlap_observed_se[1],
                mrlap_observed_pval = mrlap_res$mrlap_observed_pval[1],
                mrlap_corrected_effect = mrlap_res$mrlap_corrected_effect[1],
                mrlap_corrected_se = mrlap_res$mrlap_corrected_se[1],
                mrlap_corrected_pval = mrlap_res$mrlap_corrected_pval[1],
                mrlap_n_ivs = mrlap_res$mrlap_n_ivs[1],
                mrlap_diff_pval = mrlap_res$mrlap_diff_pval[1],
                mrlap_h2_exp = mrlap_res$mrlap_h2_exp[1],
                mrlap_h2_out = mrlap_res$mrlap_h2_out[1],
                mrlap_rg = mrlap_res$mrlap_rg[1]
              )]
              
              if(mrlap_res$mrlap_status[1] == "Success"){
                message("      ✓ MRlap completed (IVs used: ", mrlap_res$mrlap_n_ivs[1], ")")
                
                # Report if correction is substantial
                if(!is.na(mrlap_res$mrlap_diff_pval[1]) && mrlap_res$mrlap_diff_pval[1] < 0.05){
                  message("        ⚠️ Significant difference: Observed β=", 
                          round(mrlap_res$mrlap_observed_effect[1], 3),
                          " vs Corrected β=", 
                          round(mrlap_res$mrlap_corrected_effect[1], 3))
                }
              }
            }
          }
          
          all_results[[length(all_results) + 1]] <- res
          successful_pairs <- successful_pairs + 1
          message("      ✓ Analysis completed")
        }
      } else {
        message("      Insufficient SNPs after harmonization, skipping")
      }
    }
  }
  
  # 合并所有结果
  if(length(all_results) > 0){
    final_results <- rbindlist(all_results, fill = TRUE)
    fwrite(final_results, out_csv, row.names = FALSE)
    message("\n[完成] UVMR分析完成！")
    message("  总配对数: ", total_pairs)
    message("  成功分析: ", successful_pairs)
    message("  结果保存至: ", out_csv)
    message("  总行数: ", nrow(final_results))
    return(invisible(final_results))
  } else {
    message("\n[警告] 未生成任何有效结果")
    return(invisible(NULL))
  }
}

################################################################################
# STEP 8: MULTIVARIABLE MR (MVMR) ANALYSIS WITH COVARIATE ADJUSTMENT
################################################################################
#
# Evaluates the independent causal effects of multiple exposures on an outcome
# after adjusting for confounders (e.g., SES indicators)
#
# KEY DIFFERENCE from UVMR:
#   - UVMR: Uses IVs from exposure only
#   - MVMR: Uses COMBINED IVs from exposure OR covariates (union of significant SNPs)
#
# Workflow:
#   1. Select IVs from exposure GWAS (P<5e-8)
#   2. Select IVs from covariate GWAS (P<5e-8)
#   3. Combine IVs (union) and LD clump
#   4. Harmonize exposures + covariates + outcome
#   5. Run MVMR to get confounder-adjusted effects
#
# Reference: 
#   - Burgess S, Thompson SG. Multivariable Mendelian randomization. Am J Epidemiol. 2015.
#   - Ye CJ, et al. Mental well-being and healthy aging. Nat Hum Behav. 2024.
#     "genetic instruments were the combination of SNPs, which had genome-wide 
#      significance in either the GWAS of each exposure or the GWAS of each mediator"
#
################################################################################

MIN_SNPS_MVMR <- 3

#' Select Combined Instrument Variables for MVMR
#'
#' Selects IVs from BOTH exposure and covariate GWAS (union of significant SNPs)
#' This follows the standard MVMR approach where IVs can come from any variable
#'
#' @param gwas_list List of GWAS data.tables (exposures + covariates)
#' @param p_thresh P-value threshold (default: 5e-8)
#' @param r2 LD r² threshold for clumping (default: 0.001)
#' @param kb LD window size (default: 10000)
#' @return Vector of selected SNP IDs (rsID), or empty if insufficient
select_combined_ivs_mvmr <- function(gwas_list, p_thresh = P_THRESH, r2 = CLUMP_R2, kb = CLUMP_KB){
  
  # Collect all significant SNPs from any GWAS
  all_sig_snps <- list()
  
  for(i in seq_along(gwas_list)){
    gwas_dt <- gwas_list[[i]]
    
    # Check if has rsID
    if(!any(grepl("^rs\\d+$", gwas_dt$SNP))){
      next  # Skip if no rsID (can't clump)
    }
    
    # Get significant SNPs from this GWAS
    sig <- gwas_dt[pval < p_thresh, .(rsid = SNP, pval)]
    if(nrow(sig) > 0){
      all_sig_snps[[i]] <- sig
    }
  }
  
  if(length(all_sig_snps) == 0){
    return(character(0))
  }
  
  # Combine all significant SNPs (union)
  combined_snps <- rbindlist(all_sig_snps)
  
  # Keep the minimum p-value for each SNP across all GWAS
  combined_snps <- combined_snps[, .(pval = min(pval)), by = rsid]
  
  if(nrow(combined_snps) == 0){
    return(character(0))
  }
  
  # LD clumping on combined set
  clumped <- try(ieugwasr::ld_clump_local(
    tibble::as_tibble(transform(combined_snps, id="combined_gwas")),
    bfile = EUR_BFILE,
    plink_bin = PLINK_BIN,
    clump_kb = kb,
    clump_r2 = r2,
    clump_p = p_thresh
  ), silent = TRUE)
  
  if(inherits(clumped, "try-error")){
    # Fallback to remote clumping
    clumped <- try(ieugwasr::ld_clump(
      d = tibble::as_tibble(transform(combined_snps, id="combined_gwas")),
      clump_kb = kb,
      clump_r2 = r2,
      clump_p = p_thresh,
      pop = "EUR"
    ), silent = TRUE)
  }
  
  if(inherits(clumped, "try-error") || nrow(clumped) == 0){
    return(character(0))
  }
  
  return(clumped$rsid)
}

#' Prepare Data for MVMR Analysis with Covariates
#'
#' Prepares MVMR data using combined IVs from exposures and covariates
#'
#' @param expo_list List of exposure data.tables (full GWAS data, not just IVs)
#' @param expo_names Vector of exposure names
#' @param covar_list List of covariate data.tables (full GWAS data)
#' @param covar_names Vector of covariate names
#' @param out_dt Outcome data
#' @param out_name Outcome name
#' @return List with harmonized data ready for MVMR (or NULL if insufficient SNPs)
prepare_mvmr_data <- function(expo_list, expo_names, covar_list, covar_names, out_dt, out_name){
  
  # Step 1: Select combined IVs from exposures AND covariates
  all_gwas <- c(expo_list, covar_list)
  combined_ivs <- select_combined_ivs_mvmr(all_gwas)
  
  if(length(combined_ivs) < MIN_SNPS_MVMR){
    return(NULL)
  }
  
  message("    Selected ", length(combined_ivs), " combined IVs (from exposures + covariates)")
  
  # Step 2: Extract these IVs from each exposure and covariate
  # Use chr:pos for matching
  out_dt <- use_chrpos_id(out_dt)
  out_fmt <- to_outcome_format(out_dt, out_name)
  
  harmonised_list <- list()
  var_names <- c(expo_names, covar_names)
  all_var_list <- c(expo_list, covar_list)
  
  for(i in seq_along(all_var_list)){
    var_dt <- all_var_list[[i]]
    
    # Extract SNPs that are in combined IV set
    var_iv <- var_dt[SNP %in% combined_ivs]
    
    if(nrow(var_iv) == 0) next
    
    var_iv <- use_chrpos_id(var_iv)
    var_fmt <- to_exposure_format(var_iv, var_names[i])
    h <- harmonise_xy(var_fmt, out_fmt)
    h <- h[h$remove==FALSE & is.finite(h$beta.exposure) & is.finite(h$beta.outcome), ]
    
    if(nrow(h) > 0){
      harmonised_list[[i]] <- h
    }
  }
  
  if(length(harmonised_list) < 2){
    return(NULL)  # Need at least 2 variables for MVMR
  }
  
  # Step 3: Find common SNPs across all variables (intersection)
  common_snps <- Reduce(intersect, lapply(harmonised_list, function(h) h$SNP))
  
  if(length(common_snps) < MIN_SNPS_MVMR){
    return(NULL)
  }
  
  # Step 4: Subset to common SNPs and align
  harmonised_list <- lapply(harmonised_list, function(h) h[h$SNP %in% common_snps, ])
  harmonised_list <- lapply(harmonised_list, function(h) h[order(h$SNP), ])
  
  # Step 5: Create matrices for MVMR
  bx_matrix <- do.call(cbind, lapply(harmonised_list, function(h) h$beta.exposure))
  bxse_matrix <- do.call(cbind, lapply(harmonised_list, function(h) h$se.exposure))
  by <- harmonised_list[[1]]$beta.outcome
  byse <- harmonised_list[[1]]$se.outcome
  
  # Identify which columns are exposures vs covariates
  n_expo <- length(expo_names)
  expo_indices <- 1:n_expo
  covar_indices <- if(length(covar_names) > 0) (n_expo + 1):(n_expo + length(covar_names)) else integer(0)
  
  list(
    bx = bx_matrix,
    bxse = bxse_matrix,
    by = by,
    byse = byse,
    snps = common_snps,
    n_snps = length(common_snps),
    var_names = var_names[sapply(harmonised_list, function(x) !is.null(x))],
    expo_indices = expo_indices[expo_indices <= ncol(bx_matrix)],
    covar_indices = covar_indices[covar_indices <= ncol(bx_matrix)]
  )
}

#' Run MVMR Analysis with Covariate Adjustment and Sensitivity Analyses
#'
#' Performs multivariable MR using MVMR package (version 0.4)
#' Includes sensitivity analyses: MVMR-IVW, MVMR-Median, MVMR-Egger, MVMR-Lasso
#' Returns results for EXPOSURES only (covariates used for adjustment)
#'
#' Reference: Ye et al. Nat Hum Behav 2024
#'   "we performed the MVMR-median, multivariable MR-Egger (MVMR-Egger), 
#'    and multivariable MR-Lasso (MVMR-Lasso) methods"
#'
#' @param mvmr_data Prepared MVMR data from prepare_mvmr_data()
#' @param out_name Outcome name
#' @return data.table with MVMR results (multiple rows per exposure: one per method)
run_mvmr_analysis <- function(mvmr_data, out_name){
  
  if(!MVMR_AVAILABLE){
    message("      [SKIP] MVMR package not available")
    return(NULL)
  }
  
  var_names <- mvmr_data$var_names
  n_vars <- length(var_names)
  n_expo <- length(mvmr_data$expo_indices)
  n_covar <- length(mvmr_data$covar_indices)
  
  # ==== Prepare data for MVMR package ====
  # MVMR package expects specific format
  mvmr_input <- tryCatch(
    MVMR::format_mvmr(
      BXGs = mvmr_data$bx,     # Matrix: SNPs x Exposures
      BYG = mvmr_data$by,      # Vector: SNP-outcome associations
      seBXGs = mvmr_data$bxse, # Matrix: SNPs x Exposures  
      seBYG = mvmr_data$byse,  # Vector: SE of SNP-outcome
      RSID = mvmr_data$snps    # SNP IDs
    ),
    error = function(e) {
      message("      [ERROR] MVMR format_mvmr failed: ", e$message)
      NULL
    }
  )
  
  if(is.null(mvmr_input)) return(NULL)
  
  # ==== Calculate instrument strength for MVMR ====
  # Conditional F-statistics
  sres <- tryCatch(
    MVMR::strength_mvmr(mvmr_input, gencov=0),
    error = function(e) NULL
  )
  
  # ==== Run MVMR Methods ====
  
  # 1. MVMR-IVW (main method)
  mvmr_ivw <- tryCatch(
    MVMR::ivw_mvmr(mvmr_input),
    error = function(e) NULL
  )
  
  # 2. MVMR-Median (sensitivity)
  mvmr_median <- tryCatch(
    MVMR::median_mvmr(mvmr_input),
    error = function(e) NULL
  )
  
  # 3. MVMR-Egger (sensitivity - detects pleiotropy)
  mvmr_egger <- tryCatch(
    MVMR::egger_mvmr(mvmr_input),
    error = function(e) NULL
  )
  
  # 4. MVMR-Lasso (sensitivity - variable selection)
  mvmr_lasso <- tryCatch(
    MVMR::lasso_mvmr(mvmr_input),
    error = function(e) NULL
  )
  
  # ==== Organize results from all methods ====
  results_list <- list()
  
  # Helper function to extract and organize results
  add_mvmr_results <- function(mvmr_res, method_name){
    if(is.null(mvmr_res)) return(NULL)
    
    # Extract estimates
    if(method_name == "MVMR-Egger"){
      # Egger has intercept as first element
      beta_vals <- mvmr_res$Estimate[-1]  # Remove intercept
      se_vals <- mvmr_res$StdError[-1]
      pval_vals <- mvmr_res$Pvalue[-1]
      egger_intercept <- mvmr_res$Estimate[1]
      egger_intercept_pval <- mvmr_res$Pvalue[1]
    } else {
      beta_vals <- mvmr_res$Estimate
      se_vals <- mvmr_res$StdError
      pval_vals <- mvmr_res$Pvalue
      egger_intercept <- NA
      egger_intercept_pval <- NA
    }
    
    if(length(beta_vals) != n_vars) return(NULL)
    
    # Create results for all variables
    res <- data.table(
      outcome = out_name,
      variable = var_names,
      variable_type = c(rep("Exposure", n_expo), rep("Covariate", n_covar)),
      method = method_name,
      beta_mvmr = beta_vals,
      se_mvmr = se_vals,
      pval_mvmr = pval_vals,
      mvmr_egger_intercept = egger_intercept,
      mvmr_egger_intercept_pval = egger_intercept_pval
    )
    
    # Add conditional F-statistics if available
    if(!is.null(sres) && "exposure" %in% names(sres)){
      for(i in 1:nrow(res)){
        if(i <= nrow(sres)){
          res$conditional_F[i] <- sres$F[i]
        }
      }
    }
    
    # Only return exposure results (covariates used for adjustment)
    res_expo <- res[variable_type == "Exposure"]
    
    if(nrow(res_expo) > 0){
      covar_names_used <- res[variable_type == "Covariate", variable]
      res_expo[, adjusted_for := paste(covar_names_used, collapse = ";")]
      res_expo[, n_covariates := length(covar_names_used)]
    }
    
    res_expo
  }
  
  # Collect results from all methods
  if(!is.null(mvmr_ivw)){
    results_list[[length(results_list) + 1]] <- add_mvmr_results(mvmr_ivw, "MVMR-IVW")
  }
  if(!is.null(mvmr_median)){
    results_list[[length(results_list) + 1]] <- add_mvmr_results(mvmr_median, "MVMR-Median")
  }
  if(!is.null(mvmr_egger)){
    results_list[[length(results_list) + 1]] <- add_mvmr_results(mvmr_egger, "MVMR-Egger")
  }
  if(!is.null(mvmr_lasso)){
    results_list[[length(results_list) + 1]] <- add_mvmr_results(mvmr_lasso, "MVMR-Lasso")
  }
  
  # Combine all results
  if(length(results_list) > 0){
    final_res <- rbindlist(results_list, fill = TRUE)
    final_res[, n_snps := mvmr_data$n_snps]
    final_res[, n_variables := n_vars]
    return(final_res)
  } else {
    return(NULL)
  }
}

#' Run Comprehensive MVMR Analysis with Covariate Adjustment
#'
#' Main function to perform MVMR across exposure groups and outcomes
#' Adjusts for SES covariates (Education, Income, Occupation) if available
#'
#' KEY: Uses COMBINED IVs from exposures OR covariates (following MVMR standard)
#'
#' @param expos_files Vector of exposure file paths
#' @param covar_files Vector of covariate file paths (SES indicators)
#' @param outcomes_files Vector of outcome file paths
#' @param out_csv Output CSV file path
#' @param group_by How to group exposures: "category" or "all"
#' @param max_exposures_per_group Maximum number of exposures in one MVMR model
#' @return data.table with MVMR results (invisibly)
run_comprehensive_mvmr <- function(
    expos_files = all_expos_files,
    covar_files = all_covar_files,
    outcomes_files = all_out_files,
    out_csv = file.path(DIR_RES_TRIAL, "mvmr_comprehensive_results.csv"),
    group_by = "category",
    max_exposures_per_group = 10
){
  message("\n[COMPREHENSIVE MVMR ANALYSIS WITH COVARIATE ADJUSTMENT]")
  message("  Exposures: ", length(expos_files))
  message("  Covariates: ", length(covar_files))
  message("  Outcomes: ", length(outcomes_files))
  message("  Grouping strategy: ", group_by)
  message("  Max exposures per group: ", max_exposures_per_group)
  
  if(length(covar_files) == 0){
    message("\n  [WARNING] No covariates provided!")
    message("  MVMR will be performed WITHOUT confounder adjustment.")
    message("  For proper analysis, provide SES covariates (Education, Income, Occupation)")
  }
  
  all_results <- list()
  total_analyses <- 0
  successful_analyses <- 0
  
  # Group exposures by category
  if(group_by == "category"){
    # Create exposure groups based on directory
    expo_groups <- list(
      Proteome = list_gwas_files(DIR_EXPO01),
      Metabolic = list_gwas_files(DIR_EXPO02),
      Inflammatory = list_gwas_files(DIR_EXPO03)
    )
  } else if(group_by == "all"){
    # Treat all exposures as one group
    expo_groups <- list(All = expos_files)
  } else {
    message("[WARNING] Unknown grouping strategy. Using 'category'.")
    expo_groups <- list(
      Proteome = list_gwas_files(DIR_EXPO01),
      Metabolic = list_gwas_files(DIR_EXPO02),
      Inflammatory = list_gwas_files(DIR_EXPO03)
    )
  }
  
  # Remove empty groups
  expo_groups <- expo_groups[sapply(expo_groups, length) > 0]
  
  # ==== Read covariates (SES indicators) ====
  covar_list <- list()
  covar_names <- character()
  
  if(length(covar_files) > 0){
    message("\n[Reading Covariates for Confounder Adjustment]")
    for(cov_file in covar_files){
      cov_name <- basename(cov_file)
      cov_name <- clean_filename(cov_name)
      
      cov_dt <- try(read_gwas(cov_file), silent = TRUE)
      if(inherits(cov_dt, "try-error")){
        message("  [WARNING] Failed to read covariate: ", cov_name)
        next
      }
      
      covar_list[[length(covar_list) + 1]] <- cov_dt  # Full GWAS data
      covar_names <- c(covar_names, cov_name)
      message("  ✓ Loaded: ", cov_name)
    }
    message("  Total covariates loaded: ", length(covar_list))
  }
  
  # ==== Analyze each exposure group ====
  for(group_name in names(expo_groups)){
    group_files <- expo_groups[[group_name]]
    
    # Limit number of exposures per group
    if(length(group_files) > max_exposures_per_group){
      message("\n[INFO] Group '", group_name, "' has ", length(group_files), 
              " exposures. Using first ", max_exposures_per_group)
      group_files <- head(group_files, max_exposures_per_group)
    }
    
    message("\n[Exposure Group: ", group_name, "] N=", length(group_files))
    
    # Read and prepare all exposures in this group
    expo_list <- list()
    expo_names <- character()
    
    for(exp_file in group_files){
      exp_name <- basename(exp_file)
      exp_name <- clean_filename(exp_name)
      
      exp_dt <- try(read_gwas(exp_file), silent = TRUE)
      if(inherits(exp_dt, "try-error")) next
      
      # Store FULL GWAS data (not just selected IVs)
      # IV selection will be done jointly with covariates
      expo_list[[length(expo_list) + 1]] <- exp_dt
      expo_names <- c(expo_names, exp_name)
    }
    
    if(length(expo_list) < 1){
      message("  [SKIP] No valid exposures in group")
      next
    }
    
    message("  Valid exposures in group: ", length(expo_list))
    
    # Run MVMR for each outcome
    for(out_file in outcomes_files){
      total_analyses <- total_analyses + 1
      out_name <- basename(out_file)
      out_name <- clean_filename(out_name)
      
      out_dt <- try(read_gwas(out_file), silent = TRUE)
      if(inherits(out_dt, "try-error")) next
      
      # Prepare MVMR data WITH covariates
      # This will select combined IVs from exposures + covariates
      mvmr_data <- prepare_mvmr_data(expo_list, expo_names, covar_list, covar_names, out_dt, out_name)
      
      if(is.null(mvmr_data)){
        message("    • ", out_name, ": Insufficient common SNPs")
        next
      }
      
      message("    • ", out_name, ": Common SNPs = ", mvmr_data$n_snps)
      
      # Run MVMR (returns covariate-adjusted results for exposures)
      mvmr_res <- run_mvmr_analysis(mvmr_data, out_name)
      
      if(!is.null(mvmr_res) && nrow(mvmr_res) > 0){
        # Add group information
        mvmr_res[, exposure_group := group_name]
        mvmr_res[, `:=`(
          p_threshold = P_THRESH,
          clump_r2 = CLUMP_R2,
          clump_kb = CLUMP_KB
        )]
        
        all_results[[length(all_results) + 1]] <- mvmr_res
        successful_analyses <- successful_analyses + 1
        message("      ✓ MVMR completed (covariate-adjusted)")
      }
    }
  }
  
  # Combine and save results
  if(length(all_results) > 0){
    final_results <- rbindlist(all_results, fill = TRUE)
    fwrite(final_results, out_csv, row.names = FALSE)
    message("\n[完成] MVMR分析完成！")
    message("  Total analyses attempted: ", total_analyses)
    message("  Successful analyses: ", successful_analyses)
    message("  Results saved to: ", out_csv)
    message("  Total rows: ", nrow(final_results))
    return(invisible(final_results))
  } else {
    message("\n[WARNING] No valid MVMR results generated")
    return(invisible(NULL))
  }
}

################################################################################
# STEP 9: TWO-STEP MEDIATION ANALYSIS (Using MVMR)
################################################################################
#
# - Step1: exp -> med (using exposure IVs; chr:pos as key)
# - Step2: MVMR: [med, exp] -> out; joint SNPs from exposure IVs (mediator can lack rsID)
#
################################################################################
MIN_SNPS_STEP1 <- 3
MIN_SNPS_MVMR  <- 3

# >>> 修复：确保将 harmonise 后的数据转为 data.table 再用 .() 语法
align_for_mvmr <- function(h_exp_out, h_med_out){
  hx <- as.data.table(h_exp_out)
  hm <- as.data.table(h_med_out)
  
  a1 <- hx[, .(SNP,
               eao1 = effect_allele.outcome, oao1 = other_allele.outcome,
               bx_exp = beta.exposure, sx_exp = se.exposure,
               by = beta.outcome, sy = se.outcome)]
  a2 <- hm[, .(SNP,
               eao2 = effect_allele.outcome, oao2 = other_allele.outcome,
               bx_med = beta.exposure, sx_med = se.exposure,
               by2 = beta.outcome, sy2 = se.outcome)]
  M <- merge(a1, a2, by = "SNP")
  if(nrow(M) == 0) return(M[0])
  
  need_flip <- rep(NA_integer_, nrow(M))
  idx_same <- !is.na(M$eao1) & !is.na(M$eao2) & (M$eao1 == M$eao2)
  idx_swap <- !is.na(M$eao1) & !is.na(M$eao2) & (M$eao1 == M$oao2)
  need_flip[idx_same] <- 1L
  need_flip[idx_swap] <- -1L
  idx_unknown <- is.na(need_flip)
  if(any(idx_unknown)){
    s <- sign(M$by[idx_unknown] / M$by2[idx_unknown])
    s[is.na(s) | is.infinite(s)] <- 1
    need_flip[idx_unknown] <- s
  }
  M[, bx_med := bx_med * need_flip]
  out <- M[, .(SNP, bx_med, sx_med, bx_exp, sx_exp, by, sy)]
  out <- out[is.finite(bx_med) & is.finite(bx_exp) & is.finite(by) &
               is.finite(sx_med) & is.finite(sx_exp) & is.finite(sy)]
  out <- unique(out, by = "SNP")
  out
}

delta_ci <- function(b1,s1,b2,s2,bT,sT,alpha=0.05){
  varp <- ( (b2/bT)^2 * s1^2 ) + ( (b1/bT)^2 * s2^2 ) + ( (b1*b2/(bT^2))^2 * sT^2 )
  se_p <- sqrt(varp); z <- qnorm(1 - alpha/2)
  prop <- (b1*b2)/bT
  c(prop=max(0,prop), l=max(0,prop - z*se_p), u=max(0,prop + z*se_p))
}

run_comprehensive_mediation <- function(
    expos_files = all_expos_files,
    medi_files  = all_medi_files,
    outcomes_files = all_out_files,
    out_csv = file.path(DIR_RES_TRIAL, "mediation_comprehensive_results.csv")
){
  message("\n[COMPREHENSIVE MEDIATION ANALYSIS] 暴露=", length(expos_files),
          "；中介=", length(medi_files), "；结局=", length(outcomes_files))
  
  all_results <- list()
  total_triplets <- 0
  successful_triplets <- 0
  
  for(exp_file in expos_files){
    exp_name <- basename(exp_file)
    exp_name <- clean_filename(exp_name)
    message("\n[Exposure ", match(exp_file, expos_files), "/", length(expos_files), "] ", exp_name)
    
    exp_dt <- try(read_gwas(exp_file), silent = TRUE)
    if(inherits(exp_dt, "try-error")) { 
      message("  暴露读取失败，跳过。")
      next 
    }
    
    # 暴露 clump（严格；需要 rsID）
    ivs_exp <- select_instruments(exp_dt, P_THRESH, CLUMP_R2, CLUMP_KB)
    if(nrow(ivs_exp) < MIN_SNPS_STEP1){ 
      message("  IV 数=", nrow(ivs_exp), " < ", MIN_SNPS_STEP1, "，跳过该暴露。")
      next 
    }
    
    message("  选择了 ", nrow(ivs_exp), " 个独立工具变量")
    
    # 暴露 IV（chr:pos 作为 SNP）
    exp_iv <- exp_dt[SNP %in% ivs_exp$rsid]
    exp_iv <- use_chrpos_id(exp_iv)
    exp_fmt <- to_exposure_format(exp_iv, exp_name)
    
    for(med_file in medi_files){
      med_name <- basename(med_file)
      med_name <- clean_filename(med_name)
      message("  [Mediator] ", med_name)
      
      med_dt <- try(read_gwas(med_file), silent = TRUE)
      if(inherits(med_dt, "try-error")) { 
        message("    中介读取失败，跳过。")
        next 
      }
      med_dt <- use_chrpos_id(med_dt)  # 中介无 rsID 也能参与
      
      # Step1: exp -> med（用暴露 IV，chr:pos 键）
      med_fmt_as_out <- to_outcome_format(med_dt, med_name)
      h1 <- harmonise_xy(exp_fmt, med_fmt_as_out)
      h1 <- h1[h1$remove==FALSE & is.finite(h1$beta.exposure) & is.finite(h1$beta.outcome), ]
      
      if(nrow(h1) < MIN_SNPS_STEP1) { 
        message("    协调后SNP不足（Step1: Exp->Med）")
        next 
      }
      
      m1 <- try(TwoSampleMR::mr(h1, method_list = c("mr_ivw_mre")), silent = TRUE)
      if(inherits(m1, "try-error") || nrow(m1)==0) {
        message("    Step1 MR失败")
        next
      }
      
      beta1 <- m1$b[1]
      se1 <- m1$se[1]
      pval1 <- m1$pval[1]
      
      # ==== Reverse MR: Test for bidirectionality (Med -> Exp) ====
      # Required to validate mediation pathway direction
      message("    Testing reverse causation (Med → Exp)...")
      reverse_res <- run_reverse_mr(med_dt, exp_dt, med_name, exp_name)
      
      # Check bidirectionality
      is_bidirectional <- reverse_res$bidirectional == "Yes_Bidirectional"
      
      if(is_bidirectional){
        message("    [WARNING] Bidirectional relationship detected!")
        message("      Forward (Exp→Med): β=", round(beta1,3), ", P=", format(pval1, digits=3))
        message("      Reverse (Med→Exp): β=", round(reverse_res$reverse_beta,3), 
                ", P=", format(reverse_res$reverse_pval, digits=3))
        message("    Mediation validity questionable. Proceeding with caution...")
      } else {
        message("    ✓ Unidirectional (Exp→Med only)")
      }
      
      for(out_file in outcomes_files){
        total_triplets <- total_triplets + 1
        out_name <- basename(out_file)
        out_name <- clean_filename(out_name)
        
        out_dt <- try(read_gwas(out_file), silent = TRUE)
        if(inherits(out_dt, "try-error")) { 
          message("      • 结局读取失败：", out_name)
          next 
        }
        out_dt <- use_chrpos_id(out_dt)
        out_fmt <- to_outcome_format(out_dt, out_name)
        
        # 暴露-结局（仍用暴露 IV；chr:pos）
        h_exp_out <- harmonise_xy(exp_fmt, out_fmt)
        h_exp_out <- h_exp_out[h_exp_out$remove==FALSE & is.finite(h_exp_out$beta.exposure) & is.finite(h_exp_out$beta.outcome), ]
        if(nrow(h_exp_out) < MIN_SNPS_MVMR) {
          message("      • ", out_name, "：协调后SNP不足（Exp->Out）")
          next
        }
        
        # 中介-结局：联合 SNP = 暴露 IV 的 chr:pos
        med_subset <- med_dt[SNP %in% unique(h_exp_out$SNP)]
        if(nrow(med_subset) < MIN_SNPS_MVMR) {
          message("      • ", out_name, "：中介SNP不足")
          next
        }
        
        med_fmt2 <- to_exposure_format(med_subset, med_name)
        h_med_out <- harmonise_xy(med_fmt2, out_fmt)
        h_med_out <- h_med_out[h_med_out$remove==FALSE & is.finite(h_med_out$beta.exposure) & is.finite(h_med_out$beta.outcome), ]
        if(nrow(h_med_out) < MIN_SNPS_MVMR) {
          message("      • ", out_name, "：协调后SNP不足（Med->Out）")
          next
        }
        
        # 对齐到同一等位基因方向并做 MVMR
        M <- align_for_mvmr(h_exp_out, h_med_out)
        if(nrow(M) < MIN_SNPS_MVMR) {
          message("      • ", out_name, "：MVMR对齐后SNP不足")
          next
        }
        
        # MVMR分析
        mvinput <- tryCatch(
          MendelianRandomization::mr_mvinput(
          bx   = as.matrix(M[, .(bx_med, bx_exp)]),
          bxse = as.matrix(M[, .(sx_med, sx_exp)]),
          by   = M$by,
          byse = M$sy,
          exposure = c(med_name, exp_name),
          outcome  = out_name
          ),
          error = function(e) NULL
        )
        
        if(is.null(mvinput)) {
          message("      • ", out_name, "：MVMR输入构建失败")
          next
        }
        
        m2 <- try(
          MendelianRandomization::mr_mvivw(mvinput, model="default", correl=FALSE, distribution="normal"),
          silent = TRUE
        )
        
        if(inherits(m2, "try-error")) {
          message("      • ", out_name, "：MVMR分析失败")
          next
        }
        
        est_tbl <- as.data.frame(m2@Estimate)
        se_tbl <- as.data.frame(m2@StdError)
        pval_tbl <- as.data.frame(m2@Pvalue)
        
        if(nrow(est_tbl) < 2) {
          message("      • ", out_name, "：MVMR结果不完整")
          next
        }
        
        beta2_med <- est_tbl[1,1]  # mediator 的直接效应
        se2_med <- se_tbl[1,1]
        pval2_med <- pval_tbl[1,1]
        
        beta2_exp <- est_tbl[2,1]  # exposure 的直接效应
        se2_exp <- se_tbl[2,1]
        pval2_exp <- pval_tbl[2,1]
        
        # Total effect：exp -> out
        mtot <- try(TwoSampleMR::mr(h_exp_out, method_list = c("mr_ivw_mre")), silent = TRUE)
        if(inherits(mtot, "try-error") || nrow(mtot)==0) {
          message("      • ", out_name, "：总效应计算失败")
          next
        }
        
        betaT <- mtot$b[1]
        seT <- mtot$se[1]
        pvalT <- mtot$pval[1]
        
        # 计算中介效应
        ci <- delta_ci(beta1, se1, beta2_med, se2_med, betaT, seT)
        
        res_row <- data.table(
          exposure = exp_name,
          mediator = med_name,
          outcome = out_name,
          # Step1: Exp -> Med (Forward)
          beta_exp_med = beta1,
          se_exp_med = se1,
          pval_exp_med = pval1,
          n_snps_exp_med = nrow(h1),
          # Reverse MR: Med -> Exp (Bidirectionality test)
          reverse_beta_med_exp = reverse_res$reverse_beta,
          reverse_pval_med_exp = reverse_res$reverse_pval,
          bidirectional = reverse_res$bidirectional,
          # MVMR: Med效应（在Exp存在下）
          beta_med_out_direct = beta2_med,
          se_med_out_direct = se2_med,
          pval_med_out_direct = pval2_med,
          # MVMR: Exp效应（在Med存在下）
          beta_exp_out_direct = beta2_exp,
          se_exp_out_direct = se2_exp,
          pval_exp_out_direct = pval2_exp,
          n_snps_mvmr = nrow(M),
          # Total effect: Exp -> Out
          beta_exp_out_total = betaT,
          se_exp_out_total = seT,
          pval_exp_out_total = pvalT,
          n_snps_total = nrow(h_exp_out),
          # Mediation effect
          mediation_proportion = ci["prop"],
          mediation_prop_lci = ci["l"],
          mediation_prop_uci = ci["u"],
          # Direction check
          exp_med_direction = sign(beta1),
          med_out_direction = sign(beta2_med),
          direction_concordant = (sign(beta1) == sign(beta2_med)),
          # Parameters
          p_threshold = P_THRESH,
          clump_r2 = CLUMP_R2,
          clump_kb = CLUMP_KB
        )
        
        all_results[[length(all_results) + 1]] <- res_row
        successful_triplets <- successful_triplets + 1
        message("      • ", out_name, "：✓ 分析完成")
      }
    }
  }
  
  # 合并所有结果
  if(length(all_results) > 0){
    final_results <- rbindlist(all_results, fill = TRUE)
    fwrite(final_results, out_csv, row.names = FALSE)
    message("\n[完成] 中介MR分析完成！")
    message("  总三元组数: ", total_triplets)
    message("  成功分析: ", successful_triplets)
    message("  结果保存至: ", out_csv)
    message("  总行数: ", nrow(final_results))
    return(invisible(final_results))
  } else {
    message("\n[警告] 未生成任何有效中介分析结果")
    return(invisible(NULL))
  }
}

################################################################################
# STEP 10: POST-PROCESSING AND VALIDATION FUNCTIONS
################################################################################
#
# Functions for result validation and correction:
#   - FDR correction (Benjamini-Hochberg)
#   - Method concordance checking
#   - Reverse MR for bidirectionality
#   - MRlap for sample overlap correction
#
################################################################################

#' Run MRlap Analysis for Sample Overlap Correction
#'
#' MRlap corrects for bias from sample overlap between exposure and outcome GWAS
#' Also robust to winner's curse and weak instruments
#'
#' Reference: 
#'   - Ye et al. Nat Hum Behav 2024
#'   - Mounier N, Kutalik Z. Genet Epidemiol 2023
#'   - GitHub: https://github.com/n-mounier/MRlap
#'
#' MRlap requires full GWAS summary statistics files (not just harmonized IVs)
#'
#' @param exp_file Full path to exposure GWAS file
#' @param out_file Full path to outcome GWAS file
#' @param exposure_name Exposure name for labeling
#' @param outcome_name Outcome name for labeling
#' @return data.table with MRlap results (or NULL if MRlap not available)
run_mrlap_analysis <- function(exp_file, out_file, exposure_name, outcome_name){
  
  if(!MRLAP_AVAILABLE){
    return(NULL)
  }
  
  # MRlap needs FULL summary statistics files
  # It will:
  # 1. Read the files
  # 2. Perform LD score regression for genetic parameters
  # 3. Estimate genetic architecture
  # 4. Correct for sample overlap, winner's curse, weak instruments
  
  message("      Running MRlap (full GWAS correction)...")
  
  # Run MRlap with file paths
  # MRlap::MRlap() takes file paths and does everything internally
  mrlap_res <- tryCatch({
    MRlap::MRlap(
      Exposure = exp_file,           # Full path to exposure GWAS file
      Outcome = out_file,             # Full path to outcome GWAS file
      exposure_name = exposure_name,  # Label
      outcome_name = outcome_name,    # Label
      ld = "1000G_phase3_eur",        # LD reference (default)
      hm3 = NULL,                     # HapMap3 SNPs (default: auto-download)
      save_logfiles = FALSE           # Don't save LDSC log files
    )
  }, error = function(e){
    message("      [MRlap ERROR] ", e$message)
    NULL
  })
  
  if(is.null(mrlap_res)){
    return(data.table(
      exposure = exposure_name,
      outcome = outcome_name,
      mrlap_status = "Failed",
      mrlap_observed_effect = NA,
      mrlap_observed_se = NA,
      mrlap_observed_pval = NA,
      mrlap_corrected_effect = NA,
      mrlap_corrected_se = NA,
      mrlap_corrected_pval = NA,
      mrlap_diff_pval = NA
    ))
  }
  
  # Extract MRlap results
  # MRlap returns a list with:
  # - MRcorrection: observed_effect, corrected_effect, p-values, IVs used
  # - LDSC: heritability, genetic correlation estimates
  # - GeneticArchitecture: polygenicity, per-SNP heritability
  
  mr_correction <- mrlap_res$MRcorrection
  ldsc_res <- mrlap_res$LDSC
  gen_arch <- mrlap_res$GeneticArchitecture
  
  data.table(
    exposure = exposure_name,
    outcome = outcome_name,
    mrlap_status = "Success",
    # MR estimates
    mrlap_observed_effect = mr_correction$observed_effect,
    mrlap_observed_se = mr_correction$observed_effect_se,
    mrlap_observed_pval = mr_correction$observed_effect_p,
    mrlap_corrected_effect = mr_correction$corrected_effect,
    mrlap_corrected_se = mr_correction$corrected_effect_se,
    mrlap_corrected_pval = mr_correction$corrected_effect_p,
    mrlap_n_ivs = mr_correction$m_IVs,
    mrlap_diff_test = mr_correction$test_difference,
    mrlap_diff_pval = mr_correction$p_difference,
    # LDSC estimates
    mrlap_h2_exp = ldsc_res$h2_exp,
    mrlap_h2_out = ldsc_res$h2_out,
    mrlap_gcov = ldsc_res$gcov,
    mrlap_rg = ldsc_res$rg,
    # Genetic architecture
    mrlap_polygenicity = gen_arch$polygenicity,
    mrlap_persnp_h2 = gen_arch$perSNP_heritability
  )
}

################################################################################

#' Apply FDR Correction
#'
#' Calculates Benjamini-Hochberg false discovery rate q-values
#' Required when testing multiple hypotheses
#'
#' Reference: Ye et al. Nat Hum Behav 2024
#'   "false discovery rate q values were calculated by the Benjamini–Hochberg method"
#'
#' @param results_dt Results data.table with pval column
#' @param pval_col Name of p-value column (default: "pval")
#' @return Results with added q_value column
apply_fdr_correction <- function(results_dt, pval_col = "pval"){
  if(nrow(results_dt) == 0) return(results_dt)
  
  if(!pval_col %in% names(results_dt)){
    warning("P-value column not found: ", pval_col)
    return(results_dt)
  }
  
  # Calculate FDR q-values using Benjamini-Hochberg method
  results_dt[, q_value := p.adjust(get(pval_col), method = "BH")]
  
  return(results_dt)
}

#' Check Method Concordance
#'
#' Validates IVW results by checking agreement with sensitivity analyses
#'
#' Reference: Ye et al. Nat Hum Behav 2024
#'   "IVW estimates were considered causal associations only if they had the 
#'    same direction and statistical significance as at least one sensitivity analysis"
#'
#' @param results_dt Results with multiple methods (must have method, b, pval columns)
#' @param p_threshold P-value threshold for significance (default: 0.05)
#' @return Results with added concordance_validated column
check_method_concordance <- function(results_dt, p_threshold = 0.05){
  if(nrow(results_dt) == 0) return(results_dt)
  
  # For each exposure-outcome pair, check if IVW agrees with other methods
  results_dt[, concordance_validated := FALSE]
  
  # Get unique exposure-outcome pairs
  pairs <- unique(results_dt[, .(exposure, outcome)])
  
  for(i in 1:nrow(pairs)){
    exp_i <- pairs$exposure[i]
    out_i <- pairs$outcome[i]
    
    # Get all methods for this pair
    pair_results <- results_dt[exposure == exp_i & outcome == out_i]
    
    # Find IVW result
    ivw_result <- pair_results[grepl("Inverse variance weighted|IVW", method, ignore.case = TRUE)]
    
    if(nrow(ivw_result) == 0) next
    
    ivw_beta <- ivw_result$b[1]
    ivw_pval <- ivw_result$pval[1]
    ivw_sig <- ivw_pval < p_threshold
    ivw_direction <- sign(ivw_beta)
    
    # Check other methods
    other_results <- pair_results[!grepl("Inverse variance weighted|IVW", method, ignore.case = TRUE)]
    
    if(nrow(other_results) > 0){
      # Check for concordance: same direction AND significant
      concordant <- other_results[
        sign(b) == ivw_direction & pval < p_threshold
      ]
      
      # IVW validated if ≥1 other method agrees
      is_validated <- nrow(concordant) >= 1
    } else {
      is_validated <- FALSE
    }
    
    # Mark IVW result
    results_dt[exposure == exp_i & outcome == out_i & 
               grepl("Inverse variance weighted|IVW", method, ignore.case = TRUE),
               concordance_validated := is_validated]
  }
  
  return(results_dt)
}

#' Run Reverse MR
#'
#' Tests for reverse causation (Outcome → Exposure)
#' Critical for validating mediation pathways
#'
#' Reference: Ye et al. Nat Hum Behav 2024
#'   "Reverse MR between the mediator and the well-being spectrum was conducted 
#'    to determine if there was bi-directionality"
#'
#' @param outcome_dt Outcome/Mediator GWAS data (as exposure in reverse)
#' @param exposure_dt Original exposure GWAS data (as outcome in reverse)
#' @param outcome_label Label for what was outcome (now exposure in reverse)
#' @param exposure_label Label for original exposure (now outcome in reverse)
#' @return Reverse MR results (IVW only for speed)
run_reverse_mr <- function(outcome_dt, exposure_dt, outcome_label, exposure_label){
  
  # Select instruments from outcome (now treating as exposure)
  ivs_reverse <- select_instruments(outcome_dt, P_THRESH, CLUMP_R2, CLUMP_KB)
  
  if(nrow(ivs_reverse) < MIN_SNPS){
    return(data.table(
      reverse_direction = paste0(outcome_label, " → ", exposure_label),
      reverse_n_ivs = nrow(ivs_reverse),
      reverse_beta = NA,
      reverse_pval = NA,
      bidirectional = "Insufficient_IVs"
    ))
  }
  
  # Prepare reverse exposure data
  rev_exp_iv <- outcome_dt[SNP %in% ivs_reverse$rsid]
  rev_exp_iv <- use_chrpos_id(rev_exp_iv)
  rev_exp_fmt <- to_exposure_format(rev_exp_iv, outcome_label)
  
  # Prepare reverse outcome data
  exposure_dt_cp <- use_chrpos_id(exposure_dt)
  rev_out_fmt <- to_outcome_format(exposure_dt_cp, exposure_label)
  
  # Harmonize
  h_reverse <- harmonise_xy(rev_exp_fmt, rev_out_fmt)
  h_reverse <- h_reverse[h_reverse$remove==FALSE & 
                         is.finite(h_reverse$beta.exposure) & 
                         is.finite(h_reverse$beta.outcome), ]
  
  if(nrow(h_reverse) < MIN_SNPS){
    return(data.table(
      reverse_direction = paste0(outcome_label, " → ", exposure_label),
      reverse_n_ivs = nrow(ivs_reverse),
      reverse_n_snps = nrow(h_reverse),
      reverse_beta = NA,
      reverse_pval = NA,
      bidirectional = "Insufficient_SNPs"
    ))
  }
  
  # Run reverse MR (IVW only)
  mr_reverse <- tryCatch(
    TwoSampleMR::mr(h_reverse, method_list = c("mr_ivw")),
    error = function(e) NULL
  )
  
  if(is.null(mr_reverse) || nrow(mr_reverse) == 0){
    return(data.table(
      reverse_direction = paste0(outcome_label, " → ", exposure_label),
      reverse_n_ivs = nrow(ivs_reverse),
      reverse_n_snps = nrow(h_reverse),
      reverse_beta = NA,
      reverse_pval = NA,
      bidirectional = "MR_Failed"
    ))
  }
  
  # Return reverse MR results
  data.table(
    reverse_direction = paste0(outcome_label, " → ", exposure_label),
    reverse_n_ivs = nrow(ivs_reverse),
    reverse_n_snps = nrow(h_reverse),
    reverse_beta = mr_reverse$b[1],
    reverse_se = mr_reverse$se[1],
    reverse_pval = mr_reverse$pval[1],
    bidirectional = ifelse(mr_reverse$pval[1] < 0.05, "Yes_Bidirectional", "No_Unidirectional")
  )
}

################################################################################
# MAIN EXECUTION - Complete MR Analysis Pipeline
################################################################################

message("\n", rep("=", 80))
message("COMPREHENSIVE MENDELIAN RANDOMIZATION ANALYSIS PIPELINE")
message(rep("=", 80), "\n")

# ---- Step 1: Univariable MR (UVMR) ----
message(">>> STEP 1: Univariable Mendelian Randomization (UVMR)")
message("    Analyzing individual exposure-outcome associations")
uvmr_results <- run_comprehensive_uvmr()

# Post-process UVMR results
if(!is.null(uvmr_results) && nrow(uvmr_results) > 0){
  message("\n[POST-PROCESSING UVMR RESULTS]")
  
  # Apply method concordance check
  message("  - Checking method concordance...")
  uvmr_results <- check_method_concordance(uvmr_results, p_threshold = 0.05)
  n_validated <- sum(uvmr_results$concordance_validated, na.rm = TRUE)
  message("    IVW results validated by ≥1 sensitivity method: ", n_validated)
  
  # Apply FDR correction for multiple testing
  message("  - Applying FDR correction (Benjamini-Hochberg)...")
  uvmr_results <- apply_fdr_correction(uvmr_results, pval_col = "pval")
  n_fdr_sig <- sum(uvmr_results$q_value < 0.05, na.rm = TRUE)
  message("    Results with q-value < 0.05: ", n_fdr_sig)
  
  # Re-save with validation columns
  fwrite(uvmr_results, file.path(DIR_RES_TRIAL, "uvmr_comprehensive_results.csv"), row.names = FALSE)
}

# ---- Step 2: Multivariable MR (MVMR) ----
message("\n>>> STEP 2: Multivariable Mendelian Randomization (MVMR)")
message("    Analyzing joint independent effects with covariate adjustment")
mvmr_results <- run_comprehensive_mvmr()

# Post-process MVMR results
if(!is.null(mvmr_results) && nrow(mvmr_results) > 0){
  message("\n[POST-PROCESSING MVMR RESULTS]")
  
  # Apply FDR correction
  message("  - Applying FDR correction...")
  mvmr_results <- apply_fdr_correction(mvmr_results, pval_col = "pval_mvmr")
  n_fdr_sig <- sum(mvmr_results$q_value < 0.05, na.rm = TRUE)
  message("    Results with q-value < 0.05: ", n_fdr_sig)
  
  # Re-save with FDR column
  fwrite(mvmr_results, file.path(DIR_RES_TRIAL, "mvmr_comprehensive_results.csv"), row.names = FALSE)
}

# ---- Step 3: Mediation Analysis ----
message("\n>>> STEP 3: Two-Step Mediation Analysis")
message("    Analyzing mediation pathways with reverse MR checking")
mediation_results <- run_comprehensive_mediation()

# Post-process Mediation results
if(!is.null(mediation_results) && nrow(mediation_results) > 0){
  message("\n[POST-PROCESSING MEDIATION RESULTS]")
  
  # Apply FDR correction
  message("  - Applying FDR correction...")
  mediation_results <- apply_fdr_correction(mediation_results, pval_col = "pval_exp_med")
  
  # Re-save with FDR column
  fwrite(mediation_results, file.path(DIR_RES_TRIAL, "mediation_comprehensive_results.csv"), row.names = FALSE)
}

# ---- Summary ----
message("\n", rep("=", 80))
message("ANALYSIS COMPLETE WITH VALIDATION")
message(rep("=", 80), "\n")

message("Output files saved in: ", DIR_RES_TRIAL, "\n")

if(!is.null(uvmr_results)){
  message("✓ UVMR Results: uvmr_comprehensive_results.csv")
  message("  - Total rows: ", nrow(uvmr_results))
  message("  - IVW results validated by sensitivity: ", sum(uvmr_results$concordance_validated, na.rm=TRUE))
  message("  - Results significant after FDR (q<0.05): ", sum(uvmr_results$q_value < 0.05, na.rm=TRUE))
  if(MRLAP_AVAILABLE && "mrlap_status" %in% names(uvmr_results)){
    n_mrlap_success <- sum(uvmr_results$mrlap_status == "Success", na.rm=TRUE)
    n_mrlap_diff <- sum(uvmr_results$mrlap_diff_pval < 0.05, na.rm=TRUE)
    message("  - MRlap analyses completed: ", n_mrlap_success)
    if(n_mrlap_diff > 0){
      message("    (", n_mrlap_diff, " with significant correction - sample overlap/bias detected)")
    }
  }
  message("  - Unique exposures: ", uniqueN(uvmr_results$exposure))
  message("  - Unique outcomes: ", uniqueN(uvmr_results$outcome))
}

if(!is.null(mvmr_results)){
  message("\n✓ MVMR Results: mvmr_comprehensive_results.csv")
  message("  - Total rows: ", nrow(mvmr_results))
  message("  - Results significant after FDR (q<0.05): ", sum(mvmr_results$q_value < 0.05, na.rm=TRUE))
  message("  - Methods used: ", paste(unique(mvmr_results$method), collapse=", "))
  message("  - Covariates adjusted: ", ifelse(length(all_covar_files)>0, "Yes (SES)", "No"))
}

if(!is.null(mediation_results)){
  message("\n✓ Mediation Results: mediation_comprehensive_results.csv")
  message("  - Total rows: ", nrow(mediation_results))
  message("  - Pathways with reverse MR check: ", sum(!is.na(mediation_results$bidirectional), na.rm=TRUE))
  message("  - Unique exposures: ", uniqueN(mediation_results$exposure))
  message("  - Unique mediators: ", uniqueN(mediation_results$mediator))
  message("  - Unique outcomes: ", uniqueN(mediation_results$outcome))
}

message("\n", rep("=", 80))
message("VALIDATION NOTES:")
message("  - FDR q-values calculated using Benjamini-Hochberg method")
message("  - UVMR: IVW validated if ≥1 sensitivity method agrees (concordance_validated column)")
message("  - MVMR: Used MVMR package with IVW, Median, Egger, Lasso methods")
if(length(all_covar_files) > 0){
  message("  - MVMR: Adjusted for SES covariates (check 'adjusted_for' column)")
} else {
  message("  - WARNING: MVMR performed WITHOUT covariate adjustment!")
}
message(rep("=", 80))

message("\nNext steps:")
message("  1. Review results in results_trial/ folder")
message("  2. Check q_value column for FDR-corrected significance")
message("  3. For UVMR: Use only concordance_validated=TRUE results")
message("  4. Filter significant results using Results_Filter_Helper.R")
message("  5. Visualize results with your preferred tools")
message(rep("=", 80), "\n")

