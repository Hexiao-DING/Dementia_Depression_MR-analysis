################################################################################
# MR RESULTS FILTERING AND SUMMARIZATION HELPER
################################################################################
#
# Description:
#   This script provides functions to filter and summarize MR analysis results
#   NO VISUALIZATION - only data filtering and summarization
#
# Usage:
#   source("Results_Filter_Helper.R")
#   generate_summary_report()
#
# Author: MR Analysis Pipeline
# Date: 2025-10-26
#
################################################################################

library(data.table)

# ---- Configuration ----
BASE_DIR <- "D:/Projects_data&code/MR_pipeline_demo"
RES_DIR <- file.path(BASE_DIR, "results_trial")

################################################################################
# UVMR RESULTS FILTERING
################################################################################

#' Filter Significant UVMR Results
#'
#' Filters UVMR results based on statistical significance and quality criteria
#'
#' @param results_file Path to UVMR results CSV file
#' @param p_threshold P-value threshold (default: 0.05)
#' @param method_pattern Method to filter for (default: "Inverse variance weighted")
#' @param min_F Minimum F-statistic for instrument strength (default: 10)
#' @param het_p_min Minimum heterogeneity p-value to exclude heterogeneity (default: NULL, no filter)
#' @return Filtered data.table
filter_significant_uvmr <- function(
    results_file = file.path(RES_DIR, "uvmr_comprehensive_results.csv"),
    p_threshold = 0.05,
    method_pattern = "Inverse variance weighted",
    min_F = 10,
    het_p_min = NULL
){
  if(!file.exists(results_file)){
    stop("Results file does not exist: ", results_file)
  }
  
  cat("\n[INFO] Reading UVMR results...\n")
  results <- fread(results_file)
  cat("[INFO] Total rows:", nrow(results), "\n")
  
  # Filter by method
  if(!is.null(method_pattern)){
    results <- results[grepl(method_pattern, method, ignore.case = TRUE)]
    cat("[INFO] After filtering for method '", method_pattern, "': ", nrow(results), " rows\n", sep="")
  }
  
  # Filter by p-value
  results <- results[pval < p_threshold]
  cat("[INFO] After P <", p_threshold, ":", nrow(results), "rows\n")
  
  # Filter by F-statistic
  if(!is.null(min_F)){
    results <- results[!is.na(F_statistic) & F_statistic >= min_F]
    cat("[INFO] After F >=", min_F, ":", nrow(results), "rows\n")
  }
  
  # Filter by heterogeneity
  if(!is.null(het_p_min)){
    results <- results[is.na(Q_pval) | Q_pval >= het_p_min]
    cat("[INFO] After excluding significant heterogeneity (Q_pval >=", het_p_min, "):", nrow(results), "rows\n")
  }
  
  # Add OR and 95%CI for binary outcomes
  results[, OR := exp(b)]
  results[, OR_lci := exp(b - 1.96*se)]
  results[, OR_uci := exp(b + 1.96*se)]
  results[, beta_lci := b - 1.96*se]
  results[, beta_uci := b + 1.96*se]
  
  cat("\n[INFO] Final filtered results:", nrow(results), "rows\n")
  return(results)
}

#' Summarize UVMR Results
#'
#' Generates summary statistics for UVMR results
#'
#' @param results_dt Filtered UVMR results data.table
summarize_uvmr <- function(results_dt){
  cat("\n========== UVMR RESULTS SUMMARY ==========\n")
  cat("Total analyses:", nrow(results_dt), "\n")
  cat("Unique exposures:", uniqueN(results_dt$exposure), "\n")
  cat("Unique outcomes:", uniqueN(results_dt$outcome), "\n")
  
  cat("\nMethods used:\n")
  print(table(results_dt$method))
  
  cat("\nEffect direction:\n")
  results_dt[, direction := ifelse(b > 0, "Positive (risk-increasing)", "Negative (protective)")]
  print(table(results_dt$direction))
  
  cat("\nF-statistic distribution:\n")
  print(summary(results_dt$F_statistic))
  
  cat("\nP-value distribution:\n")
  print(summary(results_dt$pval))
  
  invisible(results_dt)
}

################################################################################
# MVMR RESULTS FILTERING
################################################################################

#' Filter Significant MVMR Results
#'
#' Filters MVMR results for significant independent effects
#'
#' @param results_file Path to MVMR results CSV file
#' @param p_threshold P-value threshold (default: 0.05)
#' @return Filtered data.table
filter_significant_mvmr <- function(
    results_file = file.path(RES_DIR, "mvmr_comprehensive_results.csv"),
    p_threshold = 0.05
){
  if(!file.exists(results_file)){
    stop("Results file does not exist: ", results_file)
  }
  
  cat("\n[INFO] Reading MVMR results...\n")
  results <- fread(results_file)
  cat("[INFO] Total rows:", nrow(results), "\n")
  
  # Filter by p-value
  results <- results[pval_mvmr < p_threshold]
  cat("[INFO] After P <", p_threshold, ":", nrow(results), "rows\n")
  
  # Add OR and 95%CI
  results[, OR_mvmr := exp(beta_mvmr)]
  results[, OR_mvmr_lci := exp(beta_lci)]
  results[, OR_mvmr_uci := exp(beta_uci)]
  
  cat("\n[INFO] Final filtered results:", nrow(results), "rows\n")
  return(results)
}

#' Summarize MVMR Results
#'
#' Generates summary statistics for MVMR results
#'
#' @param results_dt Filtered MVMR results data.table
summarize_mvmr <- function(results_dt){
  cat("\n========== MVMR RESULTS SUMMARY ==========\n")
  cat("Total analyses:", nrow(results_dt), "\n")
  cat("Unique exposure groups:", uniqueN(results_dt$exposure_group), "\n")
  cat("Unique exposures:", uniqueN(results_dt$exposure), "\n")
  cat("Unique outcomes:", uniqueN(results_dt$outcome), "\n")
  
  cat("\nExposure groups:\n")
  print(table(results_dt$exposure_group))
  
  cat("\nEffect direction:\n")
  results_dt[, direction := ifelse(beta_mvmr > 0, "Positive (risk-increasing)", "Negative (protective)")]
  print(table(results_dt$direction))
  
  cat("\nP-value distribution:\n")
  print(summary(results_dt$pval_mvmr))
  
  cat("\nNumber of exposures per model:\n")
  print(table(results_dt$n_exposures))
  
  invisible(results_dt)
}

################################################################################
# MEDIATION RESULTS FILTERING
################################################################################

#' Filter Significant Mediation Results
#'
#' Filters mediation analysis results for significant pathways
#'
#' @param results_file Path to mediation results CSV file
#' @param p_threshold P-value threshold (default: 0.05)
#' @param min_mediation_prop Minimum mediation proportion (default: 0)
#' @return Filtered data.table
filter_significant_mediation <- function(
    results_file = file.path(RES_DIR, "mediation_comprehensive_results.csv"),
    p_threshold = 0.05,
    min_mediation_prop = 0
){
  if(!file.exists(results_file)){
    stop("Results file does not exist: ", results_file)
  }
  
  cat("\n[INFO] Reading mediation results...\n")
  results <- fread(results_file)
  cat("[INFO] Total rows:", nrow(results), "\n")
  
  # Filter for all three pathways being significant
  results <- results[
    pval_exp_med < p_threshold &           # Exposure -> Mediator
    pval_med_out_direct < p_threshold &    # Mediator -> Outcome (adjusted)
    pval_exp_out_total < p_threshold       # Exposure -> Outcome (total)
  ]
  cat("[INFO] After filtering for all pathways P <", p_threshold, ":", nrow(results), "rows\n")
  
  # Filter by mediation proportion
  if(min_mediation_prop > 0){
    results <- results[mediation_proportion >= min_mediation_prop]
    cat("[INFO] After mediation proportion >=", min_mediation_prop, ":", nrow(results), "rows\n")
  }
  
  # Classify mediation strength
  results[, mediation_strength := fcase(
    mediation_proportion < 0.2, "Weak (<20%)",
    mediation_proportion < 0.5, "Moderate (20-50%)",
    mediation_proportion < 0.8, "Strong (50-80%)",
    default = "Very Strong (>80%)"
  )]
  
  cat("\n[INFO] Final filtered results:", nrow(results), "rows\n")
  return(results)
}

#' Summarize Mediation Results
#'
#' Generates summary statistics for mediation results
#'
#' @param results_dt Filtered mediation results data.table
summarize_mediation <- function(results_dt){
  cat("\n========== MEDIATION RESULTS SUMMARY ==========\n")
  cat("Total analyses:", nrow(results_dt), "\n")
  cat("Unique exposures:", uniqueN(results_dt$exposure), "\n")
  cat("Unique mediators:", uniqueN(results_dt$mediator), "\n")
  cat("Unique outcomes:", uniqueN(results_dt$outcome), "\n")
  
  cat("\nMediation strength distribution:\n")
  print(table(results_dt$mediation_strength))
  
  cat("\nMediation proportion statistics:\n")
  print(summary(results_dt$mediation_proportion))
  
  invisible(results_dt)
}

################################################################################
# EXPORT FUNCTIONS
################################################################################

#' Export Filtered Results
#'
#' Exports filtered results to CSV with organized columns
#'
#' @param results_dt Filtered results data.table
#' @param output_file Output file path
#' @param type Type of results: "uvmr", "mvmr", or "mediation"
export_results <- function(results_dt, output_file, type="uvmr"){
  dt <- copy(results_dt)
  
  if(type == "uvmr"){
    # Organize key columns first for UVMR
    key_cols <- c("exposure", "outcome", "method", "b", "se", "pval", 
                  "OR", "OR_lci", "OR_uci",
                  "nSNP", "F_statistic", "Q", "Q_pval",
                  "egger_intercept_pval", "presso_global_pval")
    other_cols <- setdiff(names(dt), key_cols)
    setcolorder(dt, c(key_cols[key_cols %in% names(dt)], other_cols))
    
  } else if(type == "mvmr"){
    # Organize key columns first for MVMR
    key_cols <- c("exposure_group", "exposure", "outcome", 
                  "n_exposures", "n_snps",
                  "beta_mvmr", "se_mvmr", "pval_mvmr",
                  "beta_lci", "beta_uci",
                  "OR_mvmr", "OR_mvmr_lci", "OR_mvmr_uci")
    other_cols <- setdiff(names(dt), key_cols)
    setcolorder(dt, c(key_cols[key_cols %in% names(dt)], other_cols))
    
  } else if(type == "mediation"){
    # Organize key columns first for mediation
    key_cols <- c("exposure", "mediator", "outcome",
                  "beta_exp_med", "pval_exp_med",
                  "beta_med_out_direct", "pval_med_out_direct",
                  "beta_exp_out_total", "pval_exp_out_total",
                  "mediation_proportion", "mediation_prop_lci", "mediation_prop_uci",
                  "mediation_strength")
    other_cols <- setdiff(names(dt), key_cols)
    setcolorder(dt, c(key_cols[key_cols %in% names(dt)], other_cols))
  }
  
  fwrite(dt, output_file, row.names = FALSE)
  cat("[INFO] Results exported to:", output_file, "\n")
}

################################################################################
# ONE-CLICK SUMMARY REPORT
################################################################################

#' Generate Complete Summary Report
#'
#' One-click function to filter and export UVMR, MVMR, and mediation results
#'
#' @param uvmr_file Path to UVMR results file
#' @param mvmr_file Path to MVMR results file
#' @param mediation_file Path to mediation results file
#' @param output_dir Output directory for filtered results
generate_summary_report <- function(
    uvmr_file = file.path(RES_DIR, "uvmr_comprehensive_results.csv"),
    mvmr_file = file.path(RES_DIR, "mvmr_comprehensive_results.csv"),
    mediation_file = file.path(RES_DIR, "mediation_comprehensive_results.csv"),
    output_dir = file.path(RES_DIR, "filtered_results")
){
  # Create output directory
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("\n", rep("=", 60), "\n", sep="")
  cat("GENERATING MR ANALYSIS SUMMARY REPORT\n")
  cat(rep("=", 60), "\n\n", sep="")
  
  # ========== UVMR ANALYSIS ==========
  if(file.exists(uvmr_file)){
    cat(">>> Analyzing UVMR Results...\n")
    
    # Filter significant results
    uvmr_sig <- filter_significant_uvmr(
      results_file = uvmr_file,
      p_threshold = 0.05,
      min_F = 10,
      het_p_min = 0.05  # Exclude significant heterogeneity
    )
    
    # Summarize
    summarize_uvmr(uvmr_sig)
    
    # Export
    if(nrow(uvmr_sig) > 0){
      export_results(
        uvmr_sig,
        file.path(output_dir, "uvmr_significant_filtered.csv"),
        type = "uvmr"
      )
    } else {
      cat("[WARNING] No significant UVMR results after filtering.\n")
    }
  } else {
    cat("[WARNING] UVMR results file not found.\n")
  }
  
  # ========== MVMR ANALYSIS ==========
  if(file.exists(mvmr_file)){
    cat("\n>>> Analyzing MVMR Results...\n")
    
    # Filter significant results
    mvmr_sig <- filter_significant_mvmr(
      results_file = mvmr_file,
      p_threshold = 0.05
    )
    
    # Summarize
    summarize_mvmr(mvmr_sig)
    
    # Export
    if(nrow(mvmr_sig) > 0){
      export_results(
        mvmr_sig,
        file.path(output_dir, "mvmr_significant_filtered.csv"),
        type = "mvmr"
      )
    } else {
      cat("[WARNING] No significant MVMR results after filtering.\n")
    }
  } else {
    cat("[WARNING] MVMR results file not found.\n")
  }
  
  # ========== MEDIATION ANALYSIS ==========
  if(file.exists(mediation_file)){
    cat("\n>>> Analyzing Mediation Results...\n")
    
    # Filter significant results
    med_sig <- filter_significant_mediation(
      results_file = mediation_file,
      p_threshold = 0.05,
      min_mediation_prop = 0.1  # At least 10% mediated
    )
    
    # Summarize
    summarize_mediation(med_sig)
    
    # Export
    if(nrow(med_sig) > 0){
      export_results(
        med_sig,
        file.path(output_dir, "mediation_significant_filtered.csv"),
        type = "mediation"
      )
    } else {
      cat("[WARNING] No significant mediation results after filtering.\n")
    }
  } else {
    cat("[WARNING] Mediation results file not found.\n")
  }
  
  cat("\n", rep("=", 60), "\n", sep="")
  cat("SUMMARY REPORT COMPLETE\n")
  cat("Output directory:", output_dir, "\n")
  cat(rep("=", 60), "\n\n", sep="")
}

################################################################################
# USAGE INSTRUCTIONS
################################################################################

cat("\n==========================================================\n")
cat("MR Results Filtering Helper Loaded\n")
cat("==========================================================\n\n")

cat("Available Functions:\n")
cat("  - filter_significant_uvmr()      : Filter significant UVMR results\n")
cat("  - filter_significant_mvmr()      : Filter significant MVMR results\n")
cat("  - filter_significant_mediation() : Filter significant mediation results\n")
cat("  - summarize_uvmr()               : Summarize UVMR results\n")
cat("  - summarize_mvmr()               : Summarize MVMR results\n")
cat("  - summarize_mediation()          : Summarize mediation results\n")
cat("  - export_results()               : Export filtered results to CSV\n")
cat("  - generate_summary_report()      : One-click complete report (all 3 analyses)\n\n")

cat("Quick Start:\n")
cat("  source(\"Results_Filter_Helper.R\")\n")
cat("  generate_summary_report()\n\n")

cat("==========================================================\n\n")

