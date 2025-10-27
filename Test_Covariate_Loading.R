################################################################################
# COVARIATE LOADING TEST SCRIPT
################################################################################
#
# Purpose: Test if your covariate files can be correctly loaded and standardized
# Run this BEFORE running the main analysis to verify everything is set up correctly
#
# Usage:
#   source("Test_Covariate_Loading.R")
#
################################################################################

library(data.table)
library(readr)

# Configuration
BASE_DIR <- "D:/Projects_data&code/MR_pipeline_demo"
DIR_COVAR <- file.path(BASE_DIR, "Covariates_SES")

cat("\n", rep("=", 80), "\n", sep="")
cat("TESTING COVARIATE FILE LOADING\n")
cat(rep("=", 80), "\n\n", sep="")

# ---- Check directory exists ----
cat("Step 1: Checking covariate directory...\n")
if(!dir.exists(DIR_COVAR)){
  cat("  [ERROR] Directory does not exist: ", DIR_COVAR, "\n")
  cat("  Please create this directory and place your covariate files there.\n")
  cat("\n  PowerShell command:\n")
  cat("  New-Item -Path \"", DIR_COVAR, "\" -ItemType Directory -Force\n\n", sep="")
  stop("Directory not found")
} else {
  cat("  ✓ Directory exists: ", DIR_COVAR, "\n\n")
}

# ---- List files ----
cat("Step 2: Scanning for files...\n")
files <- list.files(DIR_COVAR, pattern="(\\.tsv$)|(\\.tsv\\.gz$)|(\\.txt$)|(\\.txt\\.gz$)",
                    full.names=TRUE, ignore.case=TRUE, recursive=FALSE)

if(length(files) == 0){
  cat("  [ERROR] No covariate files found!\n")
  cat("  Expected to find .tsv.gz or .txt.gz files\n")
  cat("  Please place your covariate files in: ", DIR_COVAR, "\n\n")
  stop("No files found")
}

cat("  Found ", length(files), " file(s):\n", sep="")
for(f in files){
  cat("    • ", basename(f), "\n", sep="")
}
cat("\n")

# ---- Test reading each file ----
cat("Step 3: Testing file reading and column mapping...\n\n")

success_count <- 0

for(file_path in files){
  file_name <- basename(file_path)
  cat("  Testing: ", file_name, "\n", sep="")
  
  # Try reading
  df <- tryCatch({
    fread(file_path, nrows=10)  # Read first 10 rows only for testing
  }, error = function(e){
    # Fallback to read_tsv
    tryCatch({
      as.data.table(read_tsv(file_path, n_max=10, show_col_types=FALSE))
    }, error = function(e2){
      NULL
    })
  })
  
  if(is.null(df)){
    cat("    [ERROR] Failed to read file\n")
    next
  }
  
  cat("    ✓ File read successfully\n")
  cat("    Columns found: ", ncol(df), "\n")
  cat("    Column names: ", paste(names(df)[1:min(5, ncol(df))], collapse=", "), "...\n")
  
  # Check for required columns
  has_snp <- any(tolower(names(df)) %in% c("snp","rsid","rs_id","markername","variant"))
  has_ea <- any(tolower(names(df)) %in% c("effect_allele","ea","a1"))
  has_oa <- any(tolower(names(df)) %in% c("other_allele","oa","a2"))
  has_beta <- any(tolower(names(df)) %in% c("beta","b"))
  has_se <- any(tolower(names(df)) %in% c("se","standard_error"))
  has_p <- any(tolower(names(df)) %in% c("p","pval","p_value"))
  has_chr <- any(tolower(names(df)) %in% c("chr","chromosome"))
  has_pos <- any(tolower(names(df)) %in% c("pos","bp","position","base_pair_location"))
  
  cat("    Column detection:\n")
  cat("      SNP identifier: ", ifelse(has_snp, "✓", "✗"), "\n", sep="")
  cat("      Effect allele:  ", ifelse(has_ea, "✓", "✗"), "\n", sep="")
  cat("      Other allele:   ", ifelse(has_oa, "✓", "✗"), "\n", sep="")
  cat("      Beta:           ", ifelse(has_beta, "✓", "✗"), "\n", sep="")
  cat("      SE:             ", ifelse(has_se, "✓", "✗"), "\n", sep="")
  cat("      P-value:        ", ifelse(has_p, "✓", "✗"), "\n", sep="")
  cat("      Chr:            ", ifelse(has_chr, "✓", "✗"), "\n", sep="")
  cat("      Pos:            ", ifelse(has_pos, "✓", "✗"), "\n", sep="")
  
  # Check if all required columns present
  required_ok <- has_snp && has_ea && has_oa && has_beta && has_se && has_p
  
  if(required_ok){
    cat("    ✓ ALL REQUIRED COLUMNS PRESENT\n")
    success_count <- success_count + 1
  } else {
    cat("    [WARNING] Some required columns missing\n")
  }
  
  cat("\n")
}

# ---- Summary ----
cat(rep("=", 80), "\n", sep="")
cat("SUMMARY\n")
cat(rep("=", 80), "\n\n", sep="")

cat("Files scanned: ", length(files), "\n")
cat("Files with complete columns: ", success_count, "\n\n")

if(success_count == 3){
  cat("✓✓✓ EXCELLENT! All 3 covariate files are ready!\n\n")
  cat("Expected covariates:\n")
  cat("  1. Education (years of schooling)\n")
  cat("  2. Income (household income)\n")
  cat("  3. Occupation (occupational attainment)\n\n")
  cat("You can now run the main analysis:\n")
  cat("  source(\"MR_Debug_Script.R\")\n\n")
} else if(success_count > 0){
  cat("⚠️  PARTIAL SUCCESS: ", success_count, " file(s) OK, but ", 3-success_count, " file(s) have issues\n", sep="")
  cat("Review the errors above and fix the problematic files.\n\n")
} else {
  cat("❌ ERROR: No files could be loaded successfully\n")
  cat("Please check:\n")
  cat("  1. File formats are correct (.tsv.gz or .txt.gz)\n")
  cat("  2. Files are not corrupted\n")
  cat("  3. Files have required columns\n\n")
}

cat(rep("=", 80), "\n\n", sep="")

# ---- Test actual standardization ----
if(success_count > 0){
  cat("Testing column standardization with actual code...\n\n")
  
  # Source just the necessary functions from main script
  source("MR_Debug_Script.R", local=TRUE)
  
  for(file_path in files){
    file_name <- basename(file_path)
    cat("  Processing: ", file_name, "\n", sep="")
    
    result <- tryCatch({
      df <- read_gwas(file_path)
      cat("    ✓ Successfully standardized\n")
      cat("    Rows: ", nrow(df), "\n", sep="")
      cat("    Standardized columns: ", paste(names(df)[1:min(8, ncol(df))], collapse=", "), "...\n")
      
      # Check if chr:pos created
      if("SNP_cp" %in% names(df)){
        cat("    ✓ chr:pos identifier created\n")
      }
      
      # Count significant SNPs
      if("pval" %in% names(df)){
        n_sig <- sum(df$pval < 5e-8, na.rm=TRUE)
        cat("    Significant SNPs (P<5e-8): ", n_sig, "\n", sep="")
      }
      
      TRUE
    }, error = function(e){
      cat("    [ERROR] ", e$message, "\n", sep="")
      FALSE
    })
    
    cat("\n")
  }
}

cat(rep("=", 80), "\n")
cat("TEST COMPLETE\n")
cat(rep("=", 80), "\n\n")

