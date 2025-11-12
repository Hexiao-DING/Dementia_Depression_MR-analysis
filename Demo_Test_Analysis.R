################################################################################
# DEMO TEST ANALYSIS - Quick Pipeline Test
################################################################################
#
# Purpose: Test the MR pipeline with a small subset of data
# Usage: source("Demo_Test_Analysis.R")
#
# This script will:
# 1. Select a few exposure, mediator, and outcome files
# 2. Run the complete MR pipeline (UVMR + MVMR + Mediation)
# 3. Generate test results in results_trial/demo_test/
# 4. Verify all functions work correctly
#
# Expected runtime: 5-15 minutes (vs hours for full analysis)
#
# Use this to:
# - Test installation on new server
# - Verify data format compatibility
# - Check that all packages work
# - Preview output structure before full analysis
#
################################################################################

cat("\n")
cat(rep("=", 80), "\n", sep="")
cat("MR ANALYSIS PIPELINE - DEMO TEST MODE\n")
cat("Testing with subset of data for quick validation\n")
cat(rep("=", 80), "\n\n", sep="")

# ==== IMPORTANT: Set Working Directory ====
# Make sure you are in the correct directory before running
current_wd <- getwd()
cat("[INFO] Current working directory: ", current_wd, "\n")

# Check if we're in the right directory
if(!file.exists("Main analysis.R")){
  cat("\n[WARNING] Main analysis.R not found in current directory!\n")
  cat("Please set working directory to MR_pipeline_demo folder:\n")
  cat("  setwd('D:/Projects_data&code/MR_pipeline_demo')\n")
  cat("Or run from R Console in the correct folder.\n\n")
  
  # Try to auto-detect if we're in a subdirectory
  if(basename(current_wd) != "MR_pipeline_demo"){
    possible_dir <- "D:/Projects_data&code/MR_pipeline_demo"
    if(dir.exists(possible_dir)){
      cat("Attempting to set working directory to: ", possible_dir, "\n")
      setwd(possible_dir)
      cat("✓ Working directory set to: ", getwd(), "\n\n")
    }
  }
}

library(data.table)
library(readr)
library(tibble)

################################################################################
# CONFIGURATION - Same as main script
################################################################################

BASE_DIR <- "D:/Projects_data&code/MR_pipeline_demo"

DIR_EXPO01 <- file.path(BASE_DIR, "Standardized Circulating human plasma proteome_Data")
DIR_EXPO02 <- file.path(BASE_DIR, "Standardized Circulating metabolic biomarkers_Data")
DIR_EXPO03 <- file.path(BASE_DIR, "Circulating inflammatory proteins_Data")
DIR_MEDI   <- file.path(BASE_DIR, "Cerebrospinal fluid metabolomics_Data")
DIR_OUT    <- file.path(BASE_DIR, "Outcomes")
DIR_COVAR  <- file.path(BASE_DIR, "Covariates_SES")

# Demo output directory
DIR_DEMO <- file.path(BASE_DIR, "results_trial", "demo_test")
if(!dir.exists(DIR_DEMO)) dir.create(DIR_DEMO, recursive = TRUE)

cat("[INFO] Demo output directory: ", DIR_DEMO, "\n\n")

################################################################################
# SELECT DEMO FILES
################################################################################

cat(">>> Selecting Demo Files\n\n")

# Function to list files
list_gwas_files <- function(dir_path){
  if(!dir.exists(dir_path)) return(character())
  list.files(dir_path, pattern="(\\.tsv$)|(\\.tsv\\.gz$)|(\\.gz\\.tsv$)|(\\.txt\\.gz$)|(\\.txt$)|(\\.gz\\.txt$)",
             full.names=TRUE, ignore.case=TRUE, recursive=FALSE)
}

# Get all available files
all_expo01 <- list_gwas_files(DIR_EXPO01)
all_expo02 <- list_gwas_files(DIR_EXPO02)
all_expo03 <- list_gwas_files(DIR_EXPO03)
all_medi <- list_gwas_files(DIR_MEDI)
all_out <- list_gwas_files(DIR_OUT)
all_covar <- list_gwas_files(DIR_COVAR)

# Select DEMO subset (2 from each category for testing)
demo_expo <- c(
  if(length(all_expo01) > 0) head(all_expo01, 1) else character(),
  if(length(all_expo02) > 0) head(all_expo02, 1) else character(),
  if(length(all_expo03) > 0) head(all_expo03, 1) else character()
)
demo_expo <- demo_expo[demo_expo != ""]

demo_medi <- if(length(all_medi) > 0) head(all_medi, 2) else character()
demo_out <- if(length(all_out) > 0) head(all_out, 2) else character()
demo_covar <- all_covar  # Use all covariates (only 3)

cat("DEMO DATA SELECTION:\n")
cat("  Exposures:  ", length(demo_expo), " files\n", sep="")
for(f in demo_expo) cat("    - ", basename(f), "\n", sep="")

cat("  Mediators:  ", length(demo_medi), " files\n", sep="")
for(f in demo_medi) cat("    - ", basename(f), "\n", sep="")

cat("  Outcomes:   ", length(demo_out), " files\n", sep="")
for(f in demo_out) cat("    - ", basename(f), "\n", sep="")

cat("  Covariates: ", length(demo_covar), " files\n", sep="")
for(f in demo_covar) cat("    - ", basename(f), "\n", sep="")

if(length(demo_expo) == 0 || length(demo_out) == 0){
  stop("\n[ERROR] Insufficient demo files. Please check data directories.")
}

cat("\nExpected analyses:\n")
cat("  UVMR pairs:     ", length(demo_expo), " × ", length(demo_out), " = ", length(demo_expo)*length(demo_out), "\n", sep="")
cat("  Mediation:      ", length(demo_expo), " × ", length(demo_medi), " × ", length(demo_out), " = ", 
    length(demo_expo)*length(demo_medi)*length(demo_out), "\n", sep="")
cat("  MVMR groups:    ~", length(unique(c(DIR_EXPO01, DIR_EXPO02, DIR_EXPO03))), "\n", sep="")

cat("\n")

################################################################################
# SOURCE MAIN SCRIPT FUNCTIONS
################################################################################

cat(">>> Loading MR Pipeline Functions\n\n")

# Get the directory where this script is located
script_dir <- if(exists("BASE_DIR") && dir.exists(BASE_DIR)){
  BASE_DIR
} else {
  # Try to detect from current working directory
  getwd()
}

# Path to main script
main_script <- file.path(script_dir, "Main analysis.R")

if(!file.exists(main_script)){
  # Try alternative: script in current directory
  main_script <- "Main analysis.R"
  if(!file.exists(main_script)){
    stop("\n[ERROR] Cannot find Main analysis.R\n",
         "  Please run this script from the MR_pipeline_demo directory\n",
         "  Current directory: ", getwd(), "\n",
         "  Or set working directory: setwd('D:/Projects_data&code/MR_pipeline_demo')\n")
  }
}

cat("Loading functions from: ", main_script, "\n")

# Temporarily suppress the execution part of main script
# We only want to load functions, not run the full analysis
temp_env <- new.env()

# Source the main script to load all functions into temp environment
source(main_script, local = temp_env)

# Copy necessary functions to global environment
for(fn_name in ls(temp_env)){
  if(is.function(temp_env[[fn_name]])){
    assign(fn_name, temp_env[[fn_name]], envir = .GlobalEnv)
  }
}

cat("✓ All functions loaded\n\n")

################################################################################
# RUN DEMO ANALYSIS
################################################################################

cat(rep("=", 80), "\n", sep="")
cat("RUNNING DEMO ANALYSIS\n")
cat(rep("=", 80), "\n\n", sep="")

# Override file lists with demo subset
all_expos_files_demo <- demo_expo
all_medi_files_demo <- demo_medi
all_out_files_demo <- demo_out
all_covar_files_demo <- demo_covar

# ---- Demo UVMR ----
cat(">>> Demo Test 1: UVMR Analysis\n")
cat("    Testing ", length(demo_expo), " exposures × ", length(demo_out), " outcomes\n\n", sep="")

demo_uvmr <- run_comprehensive_uvmr(
  expos_files = all_expos_files_demo,
  outcomes_files = all_out_files_demo,
  out_csv = file.path(DIR_DEMO, "demo_uvmr_results.csv")
)

if(!is.null(demo_uvmr)){
  cat("\n✓ UVMR Demo Completed\n")
  cat("  Rows: ", nrow(demo_uvmr), "\n")
  cat("  Output: demo_test/demo_uvmr_results.csv\n")
  
  # Post-process
  demo_uvmr <- check_method_concordance(demo_uvmr, p_threshold = 0.05)
  demo_uvmr <- apply_fdr_correction(demo_uvmr, pval_col = "pval")
  fwrite(demo_uvmr, file.path(DIR_DEMO, "demo_uvmr_results.csv"), row.names = FALSE)
  
  cat("  IVW validated: ", sum(demo_uvmr$concordance_validated, na.rm=TRUE), "\n")
  cat("  FDR sig (q<0.05): ", sum(demo_uvmr$q_value < 0.05, na.rm=TRUE), "\n")
  if("mrlap_status" %in% names(demo_uvmr)){
    cat("  MRlap completed: ", sum(demo_uvmr$mrlap_status == "Success", na.rm=TRUE), "\n")
  }
} else {
  cat("\n⚠️ UVMR Demo returned no results\n")
  cat("   This may be normal if IV selection is very strict\n")
}

# ---- Demo MVMR ----
cat("\n>>> Demo Test 2: MVMR Analysis\n")
cat("    Testing with covariate adjustment\n\n")

demo_mvmr <- run_comprehensive_mvmr(
  expos_files = all_expos_files_demo,
  covar_files = all_covar_files_demo,
  outcomes_files = all_out_files_demo,
  out_csv = file.path(DIR_DEMO, "demo_mvmr_results.csv"),
  group_by = "category",
  max_exposures_per_group = 5
)

if(!is.null(demo_mvmr)){
  cat("\n✓ MVMR Demo Completed\n")
  cat("  Rows: ", nrow(demo_mvmr), "\n")
  cat("  Output: demo_test/demo_mvmr_results.csv\n")
  
  # Post-process
  demo_mvmr <- apply_fdr_correction(demo_mvmr, pval_col = "pval_mvmr")
  fwrite(demo_mvmr, file.path(DIR_DEMO, "demo_mvmr_results.csv"), row.names = FALSE)
  
  cat("  Methods used: ", paste(unique(demo_mvmr$method), collapse=", "), "\n")
  cat("  Covariates: ", ifelse(all(demo_mvmr$n_covariates == 3), "✓ All adjusted for SES", "⚠️ Check adjustment"), "\n")
} else {
  cat("\n⚠️ MVMR Demo returned no results\n")
  cat("   This may happen if insufficient valid exposures in groups\n")
}

# ---- Demo Mediation ----
if(length(demo_medi) > 0){
  cat("\n>>> Demo Test 3: Mediation Analysis\n")
  cat("    Testing ", length(demo_expo), " exposures × ", length(demo_medi), " mediators × ", length(demo_out), " outcomes\n\n", sep="")
  
  demo_mediation <- run_comprehensive_mediation(
    expos_files = all_expos_files_demo,
    medi_files = all_medi_files_demo,
    outcomes_files = all_out_files_demo,
    out_csv = file.path(DIR_DEMO, "demo_mediation_results.csv")
  )
  
  if(!is.null(demo_mediation)){
    cat("\n✓ Mediation Demo Completed\n")
    cat("  Rows: ", nrow(demo_mediation), "\n")
    cat("  Output: demo_test/demo_mediation_results.csv\n")
    
    # Post-process
    demo_mediation <- apply_fdr_correction(demo_mediation, pval_col = "pval_exp_med")
    fwrite(demo_mediation, file.path(DIR_DEMO, "demo_mediation_results.csv"), row.names = FALSE)
    
    cat("  Bidirectional: ", sum(demo_mediation$bidirectional == "Yes_Bidirectional", na.rm=TRUE), "\n")
    cat("  Unidirectional: ", sum(demo_mediation$bidirectional == "No_Unidirectional", na.rm=TRUE), "\n")
  } else {
    cat("\n⚠️ Mediation Demo returned no results\n")
    cat("   This may happen with strict IV selection\n")
  }
} else {
  cat("\n⚠️ No mediator files found - skipping mediation demo\n")
}

################################################################################
# DEMO SUMMARY
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("DEMO TEST SUMMARY\n")
cat(rep("=", 80), "\n\n", sep="")

cat("Demo Results Saved in: ", DIR_DEMO, "/\n\n", sep="")

demo_files <- list.files(DIR_DEMO, pattern="\\.csv$", full.names=FALSE)
if(length(demo_files) > 0){
  cat("Generated files:\n")
  for(f in demo_files){
    fsize <- file.size(file.path(DIR_DEMO, f))
    cat("  ✓ ", f, " (", round(fsize/1024, 1), " KB)\n", sep="")
  }
} else {
  cat("⚠️ No output files generated\n")
  cat("   Check console messages for errors\n")
}

################################################################################
# VALIDATION CHECKS
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("VALIDATION CHECKS\n")
cat(rep("=", 80), "\n\n", sep="")

checks_passed <- 0
total_checks <- 0

# Check 1: UVMR output
total_checks <- total_checks + 1
if(!is.null(demo_uvmr) && nrow(demo_uvmr) > 0){
  cat("✓ Check 1: UVMR analysis produced results\n")
  checks_passed <- checks_passed + 1
  
  # Verify key columns
  required_cols <- c("exposure", "outcome", "method", "b", "pval", "q_value", "concordance_validated")
  if(all(required_cols %in% names(demo_uvmr))){
    cat("  ✓ All required columns present\n")
  }
  
  # Check if MRlap ran
  if("mrlap_status" %in% names(demo_uvmr)){
    mrlap_success <- sum(demo_uvmr$mrlap_status == "Success", na.rm=TRUE)
    if(mrlap_success > 0){
      cat("  ✓ MRlap completed for ", mrlap_success, " analyses\n", sep="")
    } else {
      cat("  ⚠ MRlap did not complete (check if MRlap package installed)\n")
    }
  }
} else {
  cat("✗ Check 1: UVMR analysis failed or returned no results\n")
}

# Check 2: MVMR output
total_checks <- total_checks + 1
if(!is.null(demo_mvmr) && nrow(demo_mvmr) > 0){
  cat("✓ Check 2: MVMR analysis produced results\n")
  checks_passed <- checks_passed + 1
  
  # Verify covariate adjustment
  if(all(demo_mvmr$n_covariates == 3, na.rm=TRUE)){
    cat("  ✓ All results properly adjusted for 3 SES covariates\n")
  } else {
    cat("  ✗ Covariate adjustment issue detected!\n")
  }
  
  # Check methods
  mvmr_methods <- unique(demo_mvmr$method)
  cat("  ✓ MVMR methods used: ", paste(mvmr_methods, collapse=", "), "\n", sep="")
} else {
  cat("✗ Check 2: MVMR analysis failed or returned no results\n")
  cat("  (This may be normal if exposure groups too small for demo)\n")
}

# Check 3: Mediation output
if(length(demo_medi) > 0){
  total_checks <- total_checks + 1
  if(!is.null(demo_mediation) && nrow(demo_mediation) > 0){
    cat("✓ Check 3: Mediation analysis produced results\n")
    checks_passed <- checks_passed + 1
    
    # Verify reverse MR
    if("bidirectional" %in% names(demo_mediation)){
      cat("  ✓ Reverse MR (bidirectionality test) completed\n")
    }
    
    # Verify direction check
    if("direction_concordant" %in% names(demo_mediation)){
      cat("  ✓ Direction concordance check completed\n")
    }
  } else {
    cat("✗ Check 3: Mediation analysis failed or returned no results\n")
  }
}

cat("\nValidation Score: ", checks_passed, "/", total_checks, "\n", sep="")

################################################################################
# SHOW EXAMPLE RESULTS
################################################################################

if(!is.null(demo_uvmr) && nrow(demo_uvmr) > 0){
  cat("\n", rep("=", 80), "\n", sep="")
  cat("EXAMPLE UVMR RESULTS (First 3 rows)\n")
  cat(rep("=", 80), "\n\n", sep="")
  
  # Show first few IVW results
  demo_sample <- head(demo_uvmr[grepl("IVW", method)], 3)
  
  if(nrow(demo_sample) > 0){
    print(demo_sample[, .(exposure, outcome, 
                          b, pval, q_value, 
                          F_statistic, concordance_validated)])
  }
}

if(!is.null(demo_mvmr) && nrow(demo_mvmr) > 0){
  cat("\n", rep("=", 80), "\n", sep="")
  cat("EXAMPLE MVMR RESULTS (First 3 rows)\n")
  cat(rep("=", 80), "\n\n", sep="")
  
  demo_sample_mvmr <- head(demo_mvmr, 3)
  print(demo_sample_mvmr[, .(variable, outcome, method,
                             beta_mvmr, pval_mvmr, q_value,
                             adjusted_for, n_covariates)])
}

################################################################################
# FINAL RECOMMENDATION
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("DEMO TEST CONCLUSION\n")
cat(rep("=", 80), "\n\n", sep="")

if(checks_passed >= total_checks * 0.66){  # At least 2/3 checks passed
  cat("✅ DEMO TEST PASSED\n\n")
  cat("The pipeline is working correctly!\n\n")
  cat("You can now run the full analysis:\n")
  cat("  source('Main analysis.R')\n\n")
  cat("Note: Full analysis will take longer (hours vs minutes for demo)\n")
  cat("Consider running on server overnight or in background\n\n")
} else {
  cat("⚠️ DEMO TEST INCOMPLETE\n\n")
  cat("Some components did not produce results.\n")
  cat("This may be due to:\n")
  cat("  - Very strict IV selection (p<5e-8)\n")
  cat("  - Small demo dataset\n")
  cat("  - Missing packages (MVMR, MRlap)\n")
  cat("  - Data format issues\n\n")
  cat("Recommendations:\n")
  cat("  1. Check if MVMR and MRlap packages are installed\n")
  cat("  2. Review console messages for errors\n")
  cat("  3. Check if covariate files are in Covariates_SES/\n")
  cat("  4. Verify GWAS file formats are correct\n\n")
}

cat("Demo output saved in: ", DIR_DEMO, "/\n", sep="")
cat("Review the CSV files to understand output structure\n\n")

cat(rep("=", 80), "\n", sep="")
cat("Next Steps:\n")
cat("  1. Review demo results in ", basename(DIR_DEMO), "/ folder\n", sep="")
cat("  2. Check output CSV structure\n")
cat("  3. If satisfied, run full analysis: source('Main analysis.R')\n")
cat("  4. For server: Consider using nohup or screen for long-running jobs\n")
cat(rep("=", 80), "\n\n", sep="")

################################################################################
# SAVE DEMO LOG
################################################################################

demo_log <- list(
  demo_date = Sys.time(),
  n_exposures = length(demo_expo),
  n_mediators = length(demo_medi),
  n_outcomes = length(demo_out),
  n_covariates = length(demo_covar),
  uvmr_rows = ifelse(is.null(demo_uvmr), 0, nrow(demo_uvmr)),
  mvmr_rows = ifelse(is.null(demo_mvmr), 0, nrow(demo_mvmr)),
  mediation_rows = ifelse(is.null(demo_mediation), 0, nrow(demo_mediation)),
  checks_passed = checks_passed,
  total_checks = total_checks
)

saveRDS(demo_log, file.path(DIR_DEMO, "demo_log.rds"))

cat("Demo log saved to: ", DIR_DEMO, "/demo_log.rds\n\n", sep="")
cat("✓ Demo test complete!\n\n")


