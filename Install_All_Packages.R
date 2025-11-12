################################################################################
# COMPLETE R PACKAGE INSTALLATION SCRIPT
# For MR Analysis Pipeline v2.5
################################################################################
#
# Purpose: Install ALL required packages for the MR analysis pipeline
# Usage: source("Install_All_Packages.R")
# 
# This script will:
# 1. Check and install CRAN packages
# 2. Check and install GitHub packages
# 3. Verify all installations
# 4. Report any installation failures
#
# Recommended: Run this on your server before running the main analysis
#
################################################################################

cat("\n")
cat(rep("=", 80), "\n", sep="")
cat("MR ANALYSIS PIPELINE - PACKAGE INSTALLATION\n")
cat("Version 2.5 - Complete Installation\n")
cat(rep("=", 80), "\n\n", sep="")

# Function to install package if not already installed
install_if_missing <- function(pkg, repo = "CRAN"){
  if(!require(pkg, character.only = TRUE, quietly = TRUE)){
    cat("[INSTALLING] ", pkg, " from ", repo, "...\n", sep="")
    if(repo == "CRAN"){
      install.packages(pkg, dependencies = TRUE)
    }
    # Verify installation
    if(require(pkg, character.only = TRUE, quietly = TRUE)){
      cat("  ✓ ", pkg, " installed successfully\n", sep="")
      return(TRUE)
    } else {
      cat("  ✗ ", pkg, " installation FAILED\n", sep="")
      return(FALSE)
    }
  } else {
    cat("[OK] ", pkg, " already installed\n", sep="")
    return(TRUE)
  }
}

# Function to install GitHub package
install_github_if_missing <- function(pkg_name, github_repo){
  if(!require(pkg_name, character.only = TRUE, quietly = TRUE)){
    cat("[INSTALLING] ", pkg_name, " from GitHub (", github_repo, ")...\n", sep="")
    
    # Check if devtools is available
    if(!require("devtools", quietly = TRUE)){
      install.packages("devtools")
    }
    
    tryCatch({
      devtools::install_github(github_repo, upgrade = "never")
      
      # Verify installation
      if(require(pkg_name, character.only = TRUE, quietly = TRUE)){
        cat("  ✓ ", pkg_name, " installed successfully from GitHub\n", sep="")
        return(TRUE)
      } else {
        cat("  ✗ ", pkg_name, " installation FAILED\n", sep="")
        return(FALSE)
      }
    }, error = function(e){
      cat("  ✗ Error installing ", pkg_name, ": ", e$message, "\n", sep="")
      return(FALSE)
    })
  } else {
    cat("[OK] ", pkg_name, " already installed\n", sep="")
    return(TRUE)
  }
}

################################################################################
# STEP 1: ESSENTIAL CRAN PACKAGES
################################################################################

cat("\n>>> STEP 1: Installing Essential CRAN Packages\n\n")

essential_packages <- c(
  "data.table",      # Fast data manipulation
  "readr",           # Read TSV/CSV files
  "tibble",          # Modern data frames
  "devtools"         # For GitHub installations
)

success_count <- 0
for(pkg in essential_packages){
  if(install_if_missing(pkg, "CRAN")){
    success_count <- success_count + 1
  }
}

cat("\nEssential packages: ", success_count, "/", length(essential_packages), " installed\n", sep="")

################################################################################
# STEP 2: MR ANALYSIS PACKAGES (CRAN)
################################################################################

cat("\n>>> STEP 2: Installing MR Analysis Packages from CRAN\n\n")

mr_cran_packages <- c(
  "TwoSampleMR",     # Main UVMR framework (v0.5.7)
  "ieugwasr",        # LD clumping utilities
  "MVMR",            # Multivariable MR (v0.4) - CRITICAL!
  "MRPRESSO"         # Outlier detection (v1.0)
)

mr_success <- 0
mr_failed <- character()

for(pkg in mr_cran_packages){
  if(install_if_missing(pkg, "CRAN")){
    mr_success <- mr_success + 1
  } else {
    mr_failed <- c(mr_failed, pkg)
  }
}

cat("\nMR packages (CRAN): ", mr_success, "/", length(mr_cran_packages), " installed\n", sep="")

# Special handling for packages that might not be on CRAN
if(length(mr_failed) > 0){
  cat("\n[INFO] Some packages failed from CRAN. Will try GitHub installation...\n")
}

################################################################################
# STEP 3: GITHUB PACKAGES
################################################################################

cat("\n>>> STEP 3: Installing Packages from GitHub\n\n")

github_packages <- list(
  list(name = "MRlap", repo = "n-mounier/MRlap"),           # Sample overlap correction
  list(name = "MVMR", repo = "WSpiller/MVMR"),              # If not on CRAN
  list(name = "MRPRESSO", repo = "rondolab/MR-PRESSO")      # If not on CRAN
)

github_success <- 0

for(pkg_info in github_packages){
  pkg_name <- pkg_info$name
  
  # Skip if already installed successfully from CRAN
  if(require(pkg_name, character.only = TRUE, quietly = TRUE)){
    cat("[OK] ", pkg_name, " already installed (skipping GitHub)\n", sep="")
    github_success <- github_success + 1
    next
  }
  
  # Try GitHub installation
  if(install_github_if_missing(pkg_name, pkg_info$repo)){
    github_success <- github_success + 1
  }
}

cat("\nGitHub packages: ", github_success, "/", length(github_packages), " available\n", sep="")

################################################################################
# STEP 4: VERIFY ALL CRITICAL PACKAGES
################################################################################

cat("\n>>> STEP 4: Verifying Critical Package Installation\n\n")

critical_packages <- list(
  list(name = "data.table", version_required = "1.0.0"),
  list(name = "TwoSampleMR", version_required = "0.5.0"),
  list(name = "MVMR", version_required = "0.3"),
  list(name = "ieugwasr", version_required = "0.1.0")
)

optional_packages <- list(
  list(name = "MRPRESSO", version_required = "1.0"),
  list(name = "MRlap", version_required = "0.0.3")
)

# Check critical packages
cat("CRITICAL PACKAGES (MUST HAVE):\n")
all_critical_ok <- TRUE

for(pkg_info in critical_packages){
  pkg_name <- pkg_info$name
  if(require(pkg_name, character.only = TRUE, quietly = TRUE)){
    version <- as.character(packageVersion(pkg_name))
    cat("  ✓ ", pkg_name, " (v", version, ")\n", sep="")
  } else {
    cat("  ✗ ", pkg_name, " NOT INSTALLED - CRITICAL!\n", sep="")
    all_critical_ok <- FALSE
  }
}

# Check optional packages
cat("\nOPTIONAL PACKAGES (Highly Recommended):\n")

for(pkg_info in optional_packages){
  pkg_name <- pkg_info$name
  if(require(pkg_name, character.only = TRUE, quietly = TRUE)){
    version <- as.character(packageVersion(pkg_name))
    cat("  ✓ ", pkg_name, " (v", version, ")\n", sep="")
  } else {
    cat("  ⚠  ", pkg_name, " NOT INSTALLED - Analysis will proceed without it\n", sep="")
  }
}

################################################################################
# STEP 5: LOAD AND TEST ALL PACKAGES
################################################################################

cat("\n>>> STEP 5: Loading and Testing All Packages\n\n")

test_packages <- c("data.table", "readr", "tibble", "TwoSampleMR", 
                   "MVMR", "ieugwasr", "MRPRESSO", "MRlap")

loaded_count <- 0
for(pkg in test_packages){
  result <- suppressWarnings(try(library(pkg, character.only = TRUE), silent = TRUE))
  if(!inherits(result, "try-error")){
    cat("  ✓ ", pkg, " loaded successfully\n", sep="")
    loaded_count <- loaded_count + 1
  } else {
    cat("  ✗ ", pkg, " failed to load\n", sep="")
  }
}

cat("\nPackages loaded: ", loaded_count, "/", length(test_packages), "\n", sep="")

################################################################################
# STEP 6: INSTALLATION SUMMARY
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("INSTALLATION SUMMARY\n")
cat(rep("=", 80), "\n\n", sep="")

if(all_critical_ok){
  cat("✅ ALL CRITICAL PACKAGES INSTALLED\n")
  cat("   You can proceed with the MR analysis!\n\n")
} else {
  cat("❌ SOME CRITICAL PACKAGES MISSING\n")
  cat("   Please install missing packages before running analysis\n\n")
}

# Check specific critical packages
mvmr_ok <- require("MVMR", quietly = TRUE)
mrlap_ok <- require("MRlap", quietly = TRUE)

cat("Package Status:\n")
cat("  TwoSampleMR: ", ifelse(require("TwoSampleMR", quietly = TRUE), "✓ OK", "✗ MISSING"), "\n", sep="")
cat("  MVMR:        ", ifelse(mvmr_ok, "✓ OK", "✗ MISSING - CRITICAL!"), "\n", sep="")
cat("  MRPRESSO:    ", ifelse(require("MRPRESSO", quietly = TRUE), "✓ OK", "⚠ Missing (optional)"), "\n", sep="")
cat("  MRlap:       ", ifelse(mrlap_ok, "✓ OK", "⚠ Missing (optional but recommended)"), "\n", sep="")

################################################################################
# STEP 7: TROUBLESHOOTING TIPS
################################################################################

if(!mvmr_ok || !mrlap_ok){
  cat("\n", rep("=", 80), "\n", sep="")
  cat("TROUBLESHOOTING TIPS\n")
  cat(rep("=", 80), "\n\n", sep="")
  
  if(!mvmr_ok){
    cat("MVMR Installation Failed:\n")
    cat("  Try: install.packages('MVMR')\n")
    cat("  Or:  devtools::install_github('WSpiller/MVMR')\n\n")
  }
  
  if(!mrlap_ok){
    cat("MRlap Installation Failed:\n")
    cat("  Try: devtools::install_github('n-mounier/MRlap')\n")
    cat("  Note: MRlap is only on GitHub, not CRAN\n\n")
  }
  
  cat("Common Issues:\n")
  cat("  1. devtools not installed: install.packages('devtools')\n")
  cat("  2. GitHub rate limit: Wait an hour or use GitHub token\n")
  cat("  3. Dependencies missing: Install dependencies manually\n")
  cat("  4. R version too old: Update R to ≥4.0.0\n\n")
}

################################################################################
# STEP 8: NEXT STEPS
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("NEXT STEPS\n")
cat(rep("=", 80), "\n\n", sep="")

if(all_critical_ok && mvmr_ok){
  cat("✅ Environment ready for MR analysis!\n\n")
  cat("Next:\n")
  cat("  1. Verify covariate files are in Covariates_SES/ folder\n")
  cat("  2. (Optional) Test with: source('Demo_Test_Analysis.R')\n")
  cat("  3. Run full analysis: source('Main analysis.R')\n\n")
  
  if(!mrlap_ok){
    cat("⚠️  MRlap not installed:\n")
    cat("   Analysis will run without MRlap (sample overlap correction)\n")
    cat("   To install: devtools::install_github('n-mounier/MRlap')\n\n")
  }
} else {
  cat("❌ Please install missing critical packages first\n")
  cat("   Re-run this script after installation\n\n")
}

cat(rep("=", 80), "\n\n", sep="")

################################################################################
# SAVE INSTALLATION REPORT
################################################################################

# Create installation report
report <- data.frame(
  Package = c("data.table", "readr", "tibble", "TwoSampleMR", "MVMR", 
              "ieugwasr", "MRPRESSO", "MRlap"),
  Installed = c(
    require("data.table", quietly = TRUE),
    require("readr", quietly = TRUE),
    require("tibble", quietly = TRUE),
    require("TwoSampleMR", quietly = TRUE),
    require("MVMR", quietly = TRUE),
    require("ieugwasr", quietly = TRUE),
    require("MRPRESSO", quietly = TRUE),
    require("MRlap", quietly = TRUE)
  ),
  Version = sapply(c("data.table", "readr", "tibble", "TwoSampleMR", "MVMR", 
                     "ieugwasr", "MRPRESSO", "MRlap"), function(pkg){
    if(require(pkg, character.only = TRUE, quietly = TRUE)){
      as.character(packageVersion(pkg))
    } else {
      "Not installed"
    }
  }),
  Critical = c("Yes", "Yes", "Yes", "Yes", "YES - MVMR!", "Yes", "Optional", "Optional")
)

write.csv(report, "package_installation_report.csv", row.names = FALSE)
cat("Installation report saved to: package_installation_report.csv\n\n")

# Print report
print(report)

cat("\n")
cat(rep("=", 80), "\n", sep="")
cat("INSTALLATION COMPLETE\n")
cat(rep("=", 80), "\n\n", sep="")

# Return installation status
invisible(all(report$Installed[report$Critical == "Yes" | report$Critical == "YES - MVMR!"]))


