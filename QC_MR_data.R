library(tidyverse)
library(MungeSumstats)
library(future.apply)  # For parallel processing
library(fs)

plan(multisession, workers = 5)  # Reserve 2 cores for system

input_dir <- "D:/Projects_data&code/Circulating human plasma proteome_Data"
output_dir <- "D:/Projects_data&code/Standardized Circulating human plasma proteome_Data"
dir.create(output_dir, showWarnings = FALSE)

gwas_files <- fs::dir_ls(input_dir, regexp = "\\.tsv(\\.gz)?$")

process_one_file <- function(file_path) {
  tryCatch({
    fname <- fs::path_file(file_path)
    output_path <- file.path(output_dir, fname)
    
    # Process file
    MungeSumstats::format_sumstats(
      path = file_path,
      ref_genome = "GRCh37",
      nThread = 22,
      dbSNP = 144,
      return_data = FALSE,
      check_dups = TRUE,
      save_path = output_path
    )
    
    # Remove original file after successful processing
    file.remove(file_path)
    message("✅ Processed and removed: ", fname)
    
    # Return success message
    return(paste("Success:", fname))
    
  }, error = function(e) {
    message("❌ Error processing ", fs::path_file(file_path), ": ", e$message)
    return(paste("Error:", fs::path_file(file_path)))
  })
}

# Multi-threaded parallel processing of all files
results <- future_lapply(gwas_files, process_one_file, future.seed = TRUE)

# Print processing summary
success_count <- sum(grepl("^Success", unlist(results)))
error_count <- length(gwas_files) - success_count
message("\n===== Processing Summary =====")
message("Total files: ", length(gwas_files))
message("Successfully processed: ", success_count)
message("Failed: ", error_count)

# Explicitly close parallel worker threads
plan(sequential)
gc()


