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
# PERFORMANCE OPTIMIZATION: GLOBAL CACHE SYSTEMS
################################################################################
#
# These caches significantly improve performance by avoiding redundant computations:
# - HARMONIZE_CACHE: Stores harmonization results (SNP matching is expensive)
# - GWAS_READ_CACHE: Stores standardized GWAS data to avoid re-reading/re-standardizing
# - CLUMP_CACHE: Stores LD clumping results (PLINK calls are expensive)
#
################################################################################

# Harmonization cache (stores harmonized exposure-outcome pairs)
HARMONIZE_CACHE <- new.env(hash = TRUE, parent = emptyenv())
HARMONIZE_CACHE_HITS <- 0L
HARMONIZE_CACHE_MISSES <- 0L

# GWAS read cache (stores standardized GWAS data)
GWAS_READ_CACHE <- new.env(hash = TRUE, parent = emptyenv())
GWAS_READ_CACHE_HITS <- 0L
GWAS_READ_CACHE_MISSES <- 0L

# Clumping cache (stores LD clumping results)
CLUMP_CACHE <- new.env(hash = TRUE, parent = emptyenv())
CLUMP_CACHE_HITS <- 0L
CLUMP_CACHE_MISSES <- 0L

# MVMR combined IV cache (stores combined IV sets per exposure/covariate group)
MVMR_IVS_CACHE <- new.env(hash = TRUE, parent = emptyenv())
MVMR_IVS_CACHE_HITS <- 0L
MVMR_IVS_CACHE_MISSES <- 0L

# Reverse MR IV cache (stores clumped IVs for mediator/outcome in reverse MR)
REVERSE_IVS_CACHE <- new.env(hash = TRUE, parent = emptyenv())
REVERSE_IVS_CACHE_HITS <- 0L
REVERSE_IVS_CACHE_MISSES <- 0L

# Mediator format cache (stores TwoSampleMR formatted mediators as outcomes)
MEDIATOR_FMT_CACHE <- new.env(hash = TRUE, parent = emptyenv())
MEDIATOR_FMT_CACHE_HITS <- 0L
MEDIATOR_FMT_CACHE_MISSES <- 0L

# Optional disk cache stats (for read_gwas persistent cache)
DISK_CACHE_HITS <- 0L
DISK_CACHE_MISSES <- 0L

# Optional disk cache stats (for preload_exposures persistent cache)
EXPOSURE_CACHE_HITS <- 0L
EXPOSURE_CACHE_MISSES <- 0L

# Optional disk cache stats (for reverse MR IVs persistent cache)
REVERSE_IVS_DISK_CACHE_HITS <- 0L
REVERSE_IVS_DISK_CACHE_MISSES <- 0L

# Helper: Generate cache key from data
make_cache_key <- function(...){
  # Use digest for fast hashing (falls back to base if not available)
  if(requireNamespace("digest", quietly = TRUE)){
    digest::digest(list(...), algo = "xxhash64")
  } else {
    paste0(sapply(list(...), function(x) paste(head(x, 100), collapse = ":")), collapse = "_")
  }
}

# Build a cache key from file paths + mtimes and clumping params
make_file_cache_key <- function(paths, p_thresh, r2, kb){
  if(length(paths) == 0) return(NA_character_)
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = FALSE))
  paths <- sort(paths)
  mtimes <- vapply(paths, function(p) as.numeric(file.info(p)$mtime), 0)
  make_cache_key(paths, mtimes, p_thresh, r2, kb)
}

# Cache statistics reporter
report_cache_stats <- function(){
  total_harmonize <- HARMONIZE_CACHE_HITS + HARMONIZE_CACHE_MISSES
  total_gwas <- GWAS_READ_CACHE_HITS + GWAS_READ_CACHE_MISSES
  total_clump <- CLUMP_CACHE_HITS + CLUMP_CACHE_MISSES
  total_mvmr_ivs <- MVMR_IVS_CACHE_HITS + MVMR_IVS_CACHE_MISSES
  total_reverse_ivs <- REVERSE_IVS_CACHE_HITS + REVERSE_IVS_CACHE_MISSES
  total_med_fmt <- MEDIATOR_FMT_CACHE_HITS + MEDIATOR_FMT_CACHE_MISSES
  total_disk <- DISK_CACHE_HITS + DISK_CACHE_MISSES
  total_expo_disk <- EXPOSURE_CACHE_HITS + EXPOSURE_CACHE_MISSES
  total_reverse_disk <- REVERSE_IVS_DISK_CACHE_HITS + REVERSE_IVS_DISK_CACHE_MISSES

  if(total_harmonize > 0){
    hit_rate_harmonize <- round(100 * HARMONIZE_CACHE_HITS / total_harmonize, 1)
    message("[CACHE STATS - Harmonization] Hits: ", HARMONIZE_CACHE_HITS,
            " | Misses: ", HARMONIZE_CACHE_MISSES, " | Hit rate: ", hit_rate_harmonize, "%")
  }

  if(total_gwas > 0){
    hit_rate_gwas <- round(100 * GWAS_READ_CACHE_HITS / total_gwas, 1)
    message("[CACHE STATS - GWAS Read] Hits: ", GWAS_READ_CACHE_HITS,
            " | Misses: ", GWAS_READ_CACHE_MISSES, " | Hit rate: ", hit_rate_gwas, "%")
  }

  if(total_clump > 0){
    hit_rate_clump <- round(100 * CLUMP_CACHE_HITS / total_clump, 1)
    message("[CACHE STATS - LD Clumping] Hits: ", CLUMP_CACHE_HITS,
            " | Misses: ", CLUMP_CACHE_MISSES, " | Hit rate: ", hit_rate_clump, "%")
  }

  if(total_mvmr_ivs > 0){
    hit_rate_mvmr_ivs <- round(100 * MVMR_IVS_CACHE_HITS / total_mvmr_ivs, 1)
    message("[CACHE STATS - MVMR Combined IVs] Hits: ", MVMR_IVS_CACHE_HITS,
            " | Misses: ", MVMR_IVS_CACHE_MISSES, " | Hit rate: ", hit_rate_mvmr_ivs, "%")
  }

  if(total_reverse_ivs > 0){
    hit_rate_reverse_ivs <- round(100 * REVERSE_IVS_CACHE_HITS / total_reverse_ivs, 1)
    message("[CACHE STATS - Reverse MR IVs] Hits: ", REVERSE_IVS_CACHE_HITS,
            " | Misses: ", REVERSE_IVS_CACHE_MISSES, " | Hit rate: ", hit_rate_reverse_ivs, "%")
  }

  if(total_med_fmt > 0){
    hit_rate_med_fmt <- round(100 * MEDIATOR_FMT_CACHE_HITS / total_med_fmt, 1)
    message("[CACHE STATS - Mediator Format] Hits: ", MEDIATOR_FMT_CACHE_HITS,
            " | Misses: ", MEDIATOR_FMT_CACHE_MISSES, " | Hit rate: ", hit_rate_med_fmt, "%")
  }

  if(total_disk > 0){
    hit_rate_disk <- round(100 * DISK_CACHE_HITS / total_disk, 1)
    message("[CACHE STATS - Disk GWAS Cache] Hits: ", DISK_CACHE_HITS,
            " | Misses: ", DISK_CACHE_MISSES, " | Hit rate: ", hit_rate_disk, "%")
  }

  if(total_expo_disk > 0){
    hit_rate_expo_disk <- round(100 * EXPOSURE_CACHE_HITS / total_expo_disk, 1)
    message("[CACHE STATS - Disk Exposure Cache] Hits: ", EXPOSURE_CACHE_HITS,
            " | Misses: ", EXPOSURE_CACHE_MISSES, " | Hit rate: ", hit_rate_expo_disk, "%")
  }

  if(total_reverse_disk > 0){
    hit_rate_reverse_disk <- round(100 * REVERSE_IVS_DISK_CACHE_HITS / total_reverse_disk, 1)
    message("[CACHE STATS - Disk Reverse IVs] Hits: ", REVERSE_IVS_DISK_CACHE_HITS,
            " | Misses: ", REVERSE_IVS_DISK_CACHE_MISSES, " | Hit rate: ", hit_rate_reverse_disk, "%")
  }
}

################################################################################
# PARALLEL EXECUTION SETTINGS
################################################################################
#
# Controls worker count and backend selection for parallelized loops below.
# Defaults:
#   - Enable parallel execution when >1 worker is available.
#   - Use multicore on Unix-alike platforms; multisession (PSOCK) on Windows.
#   - Falls back to sequential if future.apply is missing on Windows.
#
################################################################################

safe_detect_cores <- function(){
  cores <- NA_integer_
  cores <- try(parallel::detectCores(logical = TRUE), silent = TRUE)
  if(inherits(cores, "try-error") || is.na(cores) || cores < 1) return(1L)
  as.integer(cores)
}

env_int <- function(var, default = NA_integer_){
  val <- Sys.getenv(var, unset = NA_character_)
  if(is.na(val) || val == "") return(default)
  suppressWarnings(as.integer(val))
}

env_bool <- function(var, default = FALSE){
  val <- Sys.getenv(var, unset = NA_character_)
  if(is.na(val) || val == "") return(default)
  tolower(val) %in% c("1","true","yes","y","on")
}

# User-configurable knobs (can be overridden before sourcing the script)
ENABLE_PARALLEL <- getOption("mr_parallel_enable", env_bool("MR_PARALLEL_ENABLE", TRUE))
PARALLEL_WORKERS <- env_int(
  "MR_PARALLEL_WORKERS",
  getOption("mr_parallel_workers", safe_detect_cores())
)
PARALLEL_WORKERS <- max(1L, min(safe_detect_cores(), PARALLEL_WORKERS))
# 优化：智能配置data.table线程数，避免过度订阅
# 当使用多个R worker时，限制DT线程；否则可以使用更多线程
default_dt_threads <- getOption(
  "mr_dt_threads",
  if(ENABLE_PARALLEL && PARALLEL_WORKERS > 1L) 1L else min(4L, safe_detect_cores())
)
DT_THREADS <- env_int("MR_DT_THREADS", default_dt_threads)
DT_THREADS <- max(1L, min(safe_detect_cores(), DT_THREADS))
default_future_strategy <- if(.Platform$OS.type != "windows") "multicore" else "multisession"
options(mr_future_strategy = getOption("mr_future_strategy", default_future_strategy))
# Optional: keep/restore previous future plan (on for interactive, off for speed)
RESET_FUTURE_PLAN <- getOption("mr_reset_future_plan", env_bool("MR_RESET_FUTURE_PLAN", FALSE))
# Optional: run MR leave-one-out (off by default; expensive and not used in output)
RUN_LOO <- getOption("mr_run_loo", env_bool("MR_RUN_LOO", FALSE))
# Optional: run MRlap sample-overlap correction (on by default)
RUN_MRLAP <- getOption("mr_run_mrlap", env_bool("MR_RUN_MRLAP", TRUE))
# Optional: skip reverse MR (bidirectionality check) to save time
SKIP_REVERSE_MR <- getOption("mr_skip_reverse_mr", env_bool("MR_SKIP_REVERSE_MR", FALSE))

# data.table threading (decoupled from PARALLEL_WORKERS to avoid oversubscription)
data.table::setDTthreads(DT_THREADS)

# Option: preload exposures once and reuse across UVMR/MVMR/Mediation to save time
PRELOAD_EXPOSURES_ONCE <- getOption(
  "mr_preload_exposures_once",
  env_bool("MR_PRELOAD_EXPOSURES_ONCE", TRUE)
)

# Option: skip MVMR analysis (useful when covariates have insufficient SNPs)
SKIP_MVMR <- getOption(
  "mr_skip_mvmr",
  env_bool("MR_SKIP_MVMR", FALSE)
)

# Helper: parallel map with safe fallbacks
.FUTURE_PLAN_STATE <- new.env(parent = emptyenv())
.FUTURE_PLAN_STATE$active <- FALSE
.FUTURE_PLAN_STATE$workers <- NA_integer_
.FUTURE_PLAN_STATE$strategy <- NA_character_

parallel_map <- function(X, FUN, cores = PARALLEL_WORKERS, enable = ENABLE_PARALLEL, ...){
  if(!enable || length(X) <= 1L || cores <= 1L){
    return(lapply(X, FUN, ...))
  }
  
  cores <- max(1L, min(cores, length(X)))
  
  # Prefer future.apply for cross-platform support
  if(requireNamespace("future.apply", quietly = TRUE)){
    # Default strategy can be overridden by option/env; multicore on Unix is faster.
    strategy <- getOption("mr_future_strategy", default_future_strategy)
    if(isTRUE(RESET_FUTURE_PLAN)){
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(strategy, workers = cores)
    } else {
      if(!isTRUE(.FUTURE_PLAN_STATE$active) ||
         .FUTURE_PLAN_STATE$workers != cores ||
         .FUTURE_PLAN_STATE$strategy != strategy){
        future::plan(strategy, workers = cores)
        .FUTURE_PLAN_STATE$active <- TRUE
        .FUTURE_PLAN_STATE$workers <- cores
        .FUTURE_PLAN_STATE$strategy <- strategy
      }
    }
    return(future.apply::future_lapply(X, FUN, future.seed = TRUE, ...))
  }
  
  # Unix fallback: mclapply (fork-based; not available on Windows)
  if(.Platform$OS.type != "windows"){
    return(parallel::mclapply(X, FUN, mc.cores = cores, ...))
  }
  
  message("[INFO] future.apply not installed; running sequentially (Windows PSOCK fallback).")
  lapply(X, FUN, ...)
}

# Preload exposures (read + clump) in parallel to avoid serial PLINK bottleneck
# keep_exp_dt controls whether to retain full GWAS (needed for mediation reverse MR)
preload_exposures <- function(exp_files, keep_exp_dt = FALSE, report_step_prefix = ""){
  total_files <- length(exp_files)
  message("[DEBUG] preload_exposures 开始: ", total_files, " 个文件, keep_exp_dt=", keep_exp_dt)

  # Report initial progress if prefix provided
  if(nzchar(report_step_prefix)){
    report_progress(paste0(report_step_prefix, "_clump"), 0, total_files)
  }

  # 激进优化：增加批量大小以最大化PLINK调用效率
  # 更大的批量可以显著减少进程启动开销
  batch_size <- env_int("MR_CLUMP_BATCH_SIZE", getOption("mr_clump_batch_size", 100))
  batch_size <- max(1L, batch_size)

  # 智能限制批量大小，避免过大导致卡顿
  # 当批量超过200时，可能导致PLINK卡住或内存不足
  if(batch_size > 200){
    message("[WARNING] MR_CLUMP_BATCH_SIZE=", batch_size, " 过大，自动限制为200")
    batch_size <- 200L
  }

  # Optional disk cache for preloaded exposures (persistent across runs)
  use_disk_cache <- isTRUE(EXPOSURE_CACHE_ENABLED) &&
    requireNamespace("digest", quietly = TRUE) &&
    (!isTRUE(keep_exp_dt) || isTRUE(EXPOSURE_CACHE_FULL))

  cache_path <- NA_character_
  if(use_disk_cache && total_files > 0){
    exp_files_norm <- sort(unique(normalizePath(exp_files, winslash = "/", mustWork = FALSE)))
    mtimes <- vapply(exp_files_norm, function(p) as.numeric(file.info(p)$mtime), 0)
    cache_key <- digest::digest(
      list(
        "preload_exposures",
        EXPOSURE_CACHE_VERSION,
        exp_files_norm,
        mtimes,
        P_THRESH,
        CLUMP_R2,
        CLUMP_KB,
        MIN_SNPS,
        batch_size,
        keep_exp_dt,
        USE_LOCAL_CLUMP,
        EUR_BFILE,
        PLINK_BIN
      ),
      algo = "xxhash64"
    )
    cache_path <- file.path(EXPOSURE_CACHE_DIR, paste0(cache_key, ".rds"))
    lock_path <- paste0(cache_path, ".lock")

    # 尝试读取缓存（带重试机制，避免并发冲突）
    if(file.exists(cache_path)){
      for(retry in 1:3){
        cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
        if(!is.null(cached)){
          EXPOSURE_CACHE_HITS <<- EXPOSURE_CACHE_HITS + 1L
          if(nzchar(report_step_prefix)){
            report_progress(paste0(report_step_prefix, "_clump"), total_files, total_files)
          }
          return(cached)
        }
        if(retry < 3) Sys.sleep(0.5)  # 等待其他进程完成写入
      }
    }

    # 检查是否有其他进程正在生成缓存
    if(file.exists(lock_path)){
      lock_age <- difftime(Sys.time(), file.info(lock_path)$mtime, units = "secs")
      if(lock_age < 60){  # 1分钟内的锁文件认为有效
        message("[INFO] 检测到其他进程正在生成缓存，等待5秒...")
        for(wait in 1:10){  # 最多等待5秒
          Sys.sleep(0.5)
          if(file.exists(cache_path)){
            cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
            if(!is.null(cached)){
              EXPOSURE_CACHE_HITS <<- EXPOSURE_CACHE_HITS + 1L
              if(nzchar(report_step_prefix)){
                report_progress(paste0(report_step_prefix, "_clump"), total_files, total_files)
              }
              return(cached)
            }
          }
        }
        message("[WARNING] 等待超时，继续独立生成缓存")
      } else {
        # 锁文件过期，删除
        message("[WARNING] 发现过期锁文件 (", round(lock_age), "秒)，删除并继续")
        unlink(lock_path)
      }
    }

    EXPOSURE_CACHE_MISSES <<- EXPOSURE_CACHE_MISSES + 1L
  }

  reader <- function(i){
    exp_file <- exp_files[i]
    exp_name <- clean_filename(basename(exp_file))
    exp_dt <- try(read_gwas(exp_file), silent = TRUE)
    if(inherits(exp_dt, "try-error")){
      message("  [SKIP] Failed to read exposure: ", exp_name)
      return(NULL)
    }
    list(
      exp_file = exp_file,
      exp_name = exp_name,
      exp_dt = exp_dt
    )
  }

  # Read all exposure files (in parallel)
  # 智能限制并行度：在多进程环境下减少每个进程的worker数量
  # 避免过度并行导致资源竞争和死锁
  effective_workers <- PARALLEL_WORKERS
  if(PARALLEL_WORKERS > 8 && total_files < 50){
    # 如果文件数量少，减少worker避免过度并行
    effective_workers <- min(8L, max(1L, ceiling(total_files / 4)))
    message("[INFO] 自动调整并行度: ", PARALLEL_WORKERS, " -> ", effective_workers, " (文件数=", total_files, ")")
  }
  read_list <- parallel_map(seq_along(exp_files), reader, cores = effective_workers)

  processed <- 0L
  res_list <- list()

  # Count and report failures early
  for(i in seq_along(read_list)){
    if(is.null(read_list[[i]])){
      processed <- processed + 1L
      if(nzchar(report_step_prefix)){
        report_progress(paste0(report_step_prefix, "_clump"), processed, total_files)
      }
    }
  }

  valid_idx <- which(!sapply(read_list, is.null))
  if(length(valid_idx) == 0) return(list())

  chunks <- split(valid_idx, ceiling(seq_along(valid_idx) / batch_size))
  message("[DEBUG] 分成 ", length(chunks), " 个批次，每批最多 ", batch_size, " 个文件")

  for(chunk_idx in seq_along(chunks)){
    chunk <- chunks[[chunk_idx]]
    message("[DEBUG] 处理批次 ", chunk_idx, "/", length(chunks), " (", length(chunk), " 个文件)")

    exp_dt_list <- lapply(chunk, function(i) read_list[[i]]$exp_dt)
    exp_labels <- vapply(chunk, function(i) read_list[[i]]$exp_file, "", USE.NAMES = FALSE)

    message("[DEBUG] 开始 batch_clump_exposures...")
    clump_list <- batch_clump_exposures(exp_dt_list, exp_labels, P_THRESH, CLUMP_R2, CLUMP_KB)
    message("[DEBUG] batch_clump_exposures 完成")

    for(idx in chunk){
      info <- read_list[[idx]]
      ivs <- clump_list[[info$exp_file]]
      ivs_n <- if(is.null(ivs)) 0 else nrow(ivs)

      if(ivs_n < MIN_SNPS){
        message("  [SKIP] ", info$exp_name, " IVs=", ivs_n, " < ", MIN_SNPS)
        processed <- processed + 1L
        if(nzchar(report_step_prefix)){
          report_progress(paste0(report_step_prefix, "_clump"), processed, total_files)
        }
        next
      }

      exp_iv <- info$exp_dt[SNP %in% ivs$rsid]
      exp_iv <- use_chrpos_id(exp_iv)
      exp_fmt <- to_exposure_format(exp_iv, info$exp_name)

      res_list[[length(res_list) + 1]] <- list(
        exp_file = info$exp_file,
        exp_name = info$exp_name,
        exp_fmt = exp_fmt,
        total_ivs = ivs_n,
        exp_dt = if(keep_exp_dt) info$exp_dt else NULL
      )

      if(!keep_exp_dt){
        read_list[[idx]]$exp_dt <- NULL
      }

      processed <- processed + 1L
      if(nzchar(report_step_prefix)){
        report_progress(paste0(report_step_prefix, "_clump"), processed, total_files)
      }
    }
  }

  if(length(res_list) == 0) return(list())
  names(res_list) <- vapply(res_list, function(x) x$exp_file, "", USE.NAMES = FALSE)

  if(use_disk_cache && !is.na(cache_path)){
    if(!file.exists(cache_path)){
      # 创建锁文件，防止其他进程同时写入
      lock_path <- paste0(cache_path, ".lock")
      writeLines(as.character(Sys.getpid()), lock_path)

      tmp_path <- paste0(cache_path, ".", Sys.getpid(), ".tmp")
      tryCatch({
        saveRDS(res_list, tmp_path, compress = EXPOSURE_CACHE_COMPRESS)
        if(!file.exists(cache_path)){
          file.rename(tmp_path, cache_path)
          message("[INFO] 缓存已保存: ", basename(cache_path))
        }
      }, error = function(e) {
        message("[WARNING] 缓存写入失败: ", e$message)
      })
      if(file.exists(tmp_path)) unlink(tmp_path)

      # 删除锁文件
      if(file.exists(lock_path)) unlink(lock_path)
    }
  }
  res_list
}

# Preload outcomes once to avoid repeated disk I/O and parsing across exposures
preload_outcomes <- function(out_files, report_step_prefix = ""){
  total_files <- length(out_files)

  # Report initial progress if prefix provided
  if(nzchar(report_step_prefix)){
    report_progress(paste0(report_step_prefix, "_load"), 0, total_files)
  }

  loader <- function(i){
    out_file <- out_files[i]
    out_name <- clean_filename(basename(out_file))
    out_dt <- try(read_gwas(out_file), silent = TRUE)

    # Report progress after each file (works in both sequential and parallel)
    if(nzchar(report_step_prefix)){
      report_progress(paste0(report_step_prefix, "_load"), i, total_files)
    }

    if(inherits(out_dt, "try-error")) return(NULL)
    out_dt <- use_chrpos_id(out_dt)
    out_fmt <- to_outcome_format(out_dt, out_name)

    list(
      out_file = out_file,
      out_name = out_name,
      out_dt = out_dt,
      out_fmt = out_fmt
    )
  }

  # Process by index instead of file to enable progress tracking
  res_list <- parallel_map(seq_along(out_files), loader, cores = PARALLEL_WORKERS)

  res_list <- res_list[!sapply(res_list, is.null)]
  if(length(res_list) == 0) return(list())
  names(res_list) <- vapply(res_list, function(x) x$out_file, "", USE.NAMES = FALSE)
  res_list
}

################################################################################
# STEP 1: CONFIGURATION - Paths and Parameters
################################################################################
# 
# Modify BASE_DIR, PLINK_BIN, and EUR_BFILE according to your local setup
# 
################################################################################

# Skip auto-run and strict file checks when requested (used by Python runner)
SKIP_AUTORUN <- isTRUE(getOption("mr_skip_autorun", FALSE)) ||
  identical(tolower(Sys.getenv("MR_SKIP_AUTORUN", "false")), "true")

# Progress logging throttle (reduce I/O overhead for large batch runs)
.PROGRESS_STATE <- new.env(parent = emptyenv())
PROGRESS_MAX_UPDATES <- env_int("MR_PROGRESS_MAX_UPDATES", 200L)
PROGRESS_MIN_INTERVAL_MS <- env_int("MR_PROGRESS_INTERVAL_MS", 250L)

# Report progress back to Python runner (if PY_PROGRESS_FILE is set)
report_progress <- function(step, done = 0, total = 0, msg = ""){
  pf <- Sys.getenv("PY_PROGRESS_FILE", unset = "")
  if(pf == "") return(invisible())
  
  done_num <- suppressWarnings(as.integer(done))
  total_num <- suppressWarnings(as.integer(total))
  now <- Sys.time()
  
  # Always report start/end; throttle intermediate updates to reduce I/O
  force <- FALSE
  if(is.na(done_num) || is.na(total_num)){
    force <- TRUE
  } else if(done_num <= 0L || (total_num > 0L && done_num >= total_num)){
    force <- TRUE
  }
  
  if(!force){
    max_updates <- if(is.na(PROGRESS_MAX_UPDATES) || PROGRESS_MAX_UPDATES < 1L) 200L else PROGRESS_MAX_UPDATES
    min_delta <- if(!is.na(total_num) && total_num > 0L) max(1L, floor(total_num / max_updates)) else 1L
    min_interval <- if(is.na(PROGRESS_MIN_INTERVAL_MS) || PROGRESS_MIN_INTERVAL_MS < 0L) 0 else PROGRESS_MIN_INTERVAL_MS / 1000
    state <- .PROGRESS_STATE[[step]]
    if(!is.null(state)){
      last_done <- state$done
      last_time <- state$time
      if(!is.na(last_done) && !is.na(done_num)){
        if((done_num - last_done) < min_delta &&
           as.numeric(difftime(now, last_time, units = "secs")) < min_interval){
          return(invisible())
        }
      } else if(as.numeric(difftime(now, last_time, units = "secs")) < min_interval){
        return(invisible())
      }
    }
  }
  
  .PROGRESS_STATE[[step]] <- list(done = done_num, time = now)
  ts <- format(now, "%Y-%m-%d %H:%M:%S")
  line <- paste("PROGRESS", step, done, total, msg, ts, sep = "\t")
  con <- file(pf, open = "a")
  on.exit(close(con), add = TRUE)
  writeLines(line, con = con, sep = "\n", useBytes = TRUE)
  invisible()
}

# Report errors back to Python runner (if PY_ERROR_LOG is set)
# Universal error logging function for all types of errors
report_error <- function(error_type, details = "", exposure = "", outcome = "", location = ""){
  ef <- Sys.getenv("PY_ERROR_LOG", unset = "")
  if(ef == "") return(invisible())

  tryCatch({
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    # Format: ERROR | timestamp | error_type | exposure | outcome | location | details
    line <- paste("ERROR", ts, error_type, exposure, outcome, location, details, sep = "\t")
    con <- file(ef, open = "a")
    on.exit(close(con), add = TRUE)
    writeLines(line, con = con, sep = "\n", useBytes = TRUE)
  }, error = function(e){
    # Silently fail if logging fails - don't want logging errors to break the analysis
    invisible()
  })
  invisible()
}

# Report general log messages (DEBUG, INFO, WARNING) to run.log
report_log <- function(level = "INFO", message_text = ""){
  lf <- Sys.getenv("PY_RUN_LOG", unset = "")
  if(lf == "") return(invisible())

  tryCatch({
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    # Format: LEVEL | timestamp | message
    line <- paste(level, ts, message_text, sep = "\t")
    con <- file(lf, open = "a")
    on.exit(close(con), add = TRUE)
    writeLines(line, con = con, sep = "\n", useBytes = TRUE)
  }, error = function(e){
    # Silently fail if logging fails - don't want logging errors to break the analysis
    invisible()
  })
  invisible()
}

# Append a single line to a task-level log (e.g., uvmr.log / mvmr.log)
append_task_log <- function(log_file, message_text){
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste(ts, message_text, sep = "\t")

  # Simple cross-process lock using a temporary directory
  lock_path <- paste0(log_file, ".lock")
  start <- Sys.time()
  acquired <- FALSE
  lock_timeout <- env_int("MR_LOG_LOCK_TIMEOUT", 5L)
  if(is.na(lock_timeout) || lock_timeout < 0L) lock_timeout <- 5L
  while(!acquired && as.numeric(difftime(Sys.time(), start, units = "secs")) < lock_timeout){
    acquired <- dir.create(lock_path, showWarnings = FALSE)
    if(!acquired) Sys.sleep(runif(1, 0.01, 0.05))
  }
  if(!acquired){
    return(invisible())
  }
  on.exit(unlink(lock_path, recursive = TRUE, force = TRUE), add = TRUE)

  tryCatch({
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    con <- file(log_file, open = "a")
    on.exit(close(con), add = TRUE)
    writeLines(line, con = con, sep = "\n", useBytes = TRUE)
    flush(con)
  }, error = function(e) invisible())
  invisible()
}

# Clean up temporary folders created under an output directory
cleanup_temp_dirs <- function(output_dir){
  tmp_base <- file.path(output_dir, "tmp")
  if(!dir.exists(tmp_base)) return(invisible())

  # Remove Rtmp* subdirectories
  tmp_dirs <- list.dirs(tmp_base, recursive = FALSE, full.names = TRUE)
  for(d in tmp_dirs){
    if(grepl("Rtmp", basename(d), fixed = TRUE)){
      unlink(d, recursive = TRUE, force = TRUE)
    }
  }

  # Remove tmp/ itself if empty
  remaining <- list.files(tmp_base, all.files = TRUE, no.. = TRUE)
  if(length(remaining) == 0){
    unlink(tmp_base, recursive = TRUE, force = TRUE)
  }
  invisible()
}

# Base directory for the project
BASE_DIR <- "/home/Dementia_Depression_MR-analysis/data/MR_pipeline_demo"
base_override <- Sys.getenv("BASE_DIR", "")
if(nzchar(base_override)) BASE_DIR <- base_override

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

# Optional disk cache for GWAS reads (persistent across runs)
DISK_CACHE_ENABLED <- getOption("mr_disk_cache", env_bool("MR_DISK_CACHE", FALSE))
DISK_CACHE_DIR <- getOption(
  "mr_disk_cache_dir",
  Sys.getenv("MR_DISK_CACHE_DIR", file.path(DIR_RES_TRIAL, "cache"))
)
DISK_CACHE_COMPRESS <- getOption("mr_disk_cache_compress", env_bool("MR_DISK_CACHE_COMPRESS", FALSE))
if(isTRUE(DISK_CACHE_ENABLED)){
  dir.create(DISK_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Optional disk cache for preloaded exposures (persistent across runs)
EXPOSURE_CACHE_ENABLED <- getOption("mr_exposure_cache", env_bool("MR_EXPOSURE_CACHE", FALSE))
EXPOSURE_CACHE_DIR <- getOption(
  "mr_exposure_cache_dir",
  Sys.getenv("MR_EXPOSURE_CACHE_DIR", file.path(DIR_RES_TRIAL, "cache", "exposures"))
)
EXPOSURE_CACHE_COMPRESS <- getOption("mr_exposure_cache_compress", env_bool("MR_EXPOSURE_CACHE_COMPRESS", FALSE))
# By default avoid caching full GWAS tables (keep_exp_dt=TRUE) to limit disk usage
EXPOSURE_CACHE_FULL <- getOption("mr_exposure_cache_full", env_bool("MR_EXPOSURE_CACHE_FULL", FALSE))
EXPOSURE_CACHE_VERSION <- getOption("mr_exposure_cache_version", "v1")
if(isTRUE(EXPOSURE_CACHE_ENABLED)){
  dir.create(EXPOSURE_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Optional disk cache for reverse MR IVs (persistent across runs/batches)
REVERSE_IVS_DISK_CACHE_ENABLED <- getOption("mr_reverse_ivs_disk_cache", env_bool("MR_REVERSE_IVS_DISK_CACHE", FALSE))
REVERSE_IVS_DISK_CACHE_DIR <- getOption(
  "mr_reverse_ivs_disk_cache_dir",
  Sys.getenv("MR_REVERSE_IVS_CACHE_DIR", file.path(DIR_RES_TRIAL, "cache", "reverse_ivs"))
)
REVERSE_IVS_DISK_CACHE_COMPRESS <- getOption("mr_reverse_ivs_disk_cache_compress", env_bool("MR_REVERSE_IVS_DISK_COMPRESS", FALSE))
REVERSE_IVS_DISK_CACHE_VERSION <- getOption("mr_reverse_ivs_disk_cache_version", "v1")
if(isTRUE(REVERSE_IVS_DISK_CACHE_ENABLED)){
  dir.create(REVERSE_IVS_DISK_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Optional cache for mediator outcome formatting (fast reuse across batches)
MEDIATOR_FMT_CACHE_ENABLED <- getOption("mr_mediator_format_cache", env_bool("MR_MEDIATOR_FORMAT_CACHE", TRUE))

# PLINK executable for local LD clumping
PLINK_BIN <- "/home/Dementia_Depression_MR-analysis/data/MR_pipeline_demo/local_clump/plink/plink"
plink_override <- Sys.getenv("PLINK_BIN", "")
if(nzchar(plink_override)) PLINK_BIN <- plink_override

# 1000 Genomes EUR reference panel for LD clumping
# Should point to the file prefix (without .bed/.bim/.fam extension)
EUR_BFILE <- "/home/Dementia_Depression_MR-analysis/data/MR_pipeline_demo/local_clump/EUR/1000G.EUR.QC"
eur_override <- Sys.getenv("EUR_BFILE", "")
if(nzchar(eur_override)) EUR_BFILE <- eur_override

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
  if(!dir.exists(dir_out)){
    if(SKIP_AUTORUN){
      message("[INFO] Outcome directory missing: ", dir_out, " (MR_SKIP_AUTORUN=TRUE). Returning empty set.")
      return(character())
    }
    stop("Outcome directory does not exist: ", dir_out)
  }
  
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

# Verify that we have data files for required categories (skip when running programmatically)
if(!SKIP_AUTORUN){
  stopifnot(length(all_expos_files)>0, length(all_medi_files)>0, length(all_out_files)>0)
} else {
  if(length(all_expos_files)==0 || length(all_medi_files)==0 || length(all_out_files)==0){
    message("[INFO] MR_SKIP_AUTORUN=TRUE; data presence checks bypassed.")
  }
}

# Report covariates
if(length(all_covar_files) > 0){
  message("[INFO] Covariates found: ", length(all_covar_files), " files")
  report_log("INFO", paste("Covariates found:", length(all_covar_files), "files"))
  message("  These will be used for confounder adjustment in MVMR")
} else {
  message("[WARNING] No covariates found in DIR_COVAR")
  report_log("WARNING", "No covariates found in DIR_COVAR")
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
  # 激进优化：向量化批量转换，避免循环
  num_cols <- intersect(c("beta","se","pval","eaf","pos","samplesize","chr","z"), names(dt))
  if(length(num_cols) > 0){
    dt[, (num_cols) := lapply(.SD, to_num), .SDcols = num_cols]
  }
  
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
  dt <- unique(dt)

  # ---- PERFORMANCE OPTIMIZATION: Create SNP index ----
  # Setting a key dramatically speeds up SNP lookups (from O(n) to O(log n))
  # This is crucial when filtering SNPs (e.g., exp_dt[SNP %in% ivs$rsid])
  if("SNP" %in% names(dt)){
    setkey(dt, SNP)
  }

  dt
}

#' Read and Standardize GWAS Summary Statistics
#'
#' Reads a GWAS file (TSV or TSV.GZ) and standardizes column names
#' Tries multiple reading strategies to handle different file formats
#'
#' PERFORMANCE: Results are cached to avoid re-reading/re-standardizing
#' Cache key is based on file path + modification time
#'
#' @param file Path to GWAS summary statistics file
#' @return Standardized data.table
read_gwas <- function(file){
  # ---- PERFORMANCE OPTIMIZATION: Check cache first ----
  # Cache key: file path + modification time (detects file changes)
  file_mtime <- file.info(file)$mtime
  cache_key <- paste0(file, "_", as.numeric(file_mtime))

  if(!is.null(GWAS_READ_CACHE[[cache_key]])){
    GWAS_READ_CACHE_HITS <<- GWAS_READ_CACHE_HITS + 1L
    return(copy(GWAS_READ_CACHE[[cache_key]]))  # Return a copy to avoid reference issues
  }

  GWAS_READ_CACHE_MISSES <<- GWAS_READ_CACHE_MISSES + 1L

  # ---- Optional disk cache (persistent across runs) ----
  disk_path <- NA_character_
  if(isTRUE(DISK_CACHE_ENABLED) && requireNamespace("digest", quietly = TRUE)){
    disk_key <- digest::digest(
      list(normalizePath(file, winslash = "/", mustWork = FALSE), as.numeric(file_mtime)),
      algo = "xxhash64"
    )
    disk_path <- file.path(DISK_CACHE_DIR, paste0(disk_key, ".rds"))
    if(file.exists(disk_path)){
      dt_disk <- tryCatch(readRDS(disk_path), error = function(e) NULL)
      if(!is.null(dt_disk)){
        DISK_CACHE_HITS <<- DISK_CACHE_HITS + 1L
        GWAS_READ_CACHE[[cache_key]] <- dt_disk
        return(copy(dt_disk))
      }
    }
    DISK_CACHE_MISSES <<- DISK_CACHE_MISSES + 1L
  }

  # ---- Read and standardize (cache miss) ----
  # 激进优化：最大化fread性能
  # 1. 使用所有可用线程（除非在多worker模式）
  # 2. 禁用verbose以减少I/O
  # 3. 使用data.table的原生操作避免拷贝
  read_threads <- if(ENABLE_PARALLEL && PARALLEL_WORKERS > 1L) 1L else DT_THREADS
  dt <- try(fread(file, sep="\t", na.strings=c("NA","NaN",""), showProgress=FALSE,
                  nThread = read_threads, fill=TRUE, verbose=FALSE,
                  data.table=TRUE), silent = TRUE)

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
  dt <- standardize_cols(dt)

  # Store in cache
  GWAS_READ_CACHE[[cache_key]] <- dt

  # Persist to disk cache for reuse across batches
  if(!is.na(disk_path)){
    if(!file.exists(disk_path)){
      tmp_path <- paste0(disk_path, ".", Sys.getpid(), ".tmp")
      tryCatch({
        saveRDS(dt, tmp_path, compress = DISK_CACHE_COMPRESS)
        if(!file.exists(disk_path)){
          file.rename(tmp_path, disk_path)
        }
      }, error = function(e) NULL)
      if(file.exists(tmp_path)) unlink(tmp_path)
    }
  }

  dt
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

#' Batch LD Clumping for Multiple Exposures (PERFORMANCE OPTIMIZATION)
#'
#' Performs LD clumping for multiple exposures in a single PLINK call
#' This dramatically reduces overhead by:
#' - Running PLINK only once instead of N times
#' - Loading LD reference panel once instead of N times
#' - Amortizing process startup costs
#'
#' @param exp_dt_list List of exposure data.tables
#' @param exp_labels Character vector of exposure labels (same length as exp_dt_list)
#' @param p_thresh P-value threshold for significance (default: 5e-8)
#' @param r2 LD r² threshold for clumping (default: 0.001)
#' @param kb LD window size in kb (default: 10000)
#' @return Named list of clumped SNP data.tables (one per exposure)
batch_clump_exposures <- function(exp_dt_list, exp_labels, p_thresh = P_THRESH, r2 = CLUMP_R2, kb = CLUMP_KB){

  if(length(exp_dt_list) != length(exp_labels)){
    stop("exp_dt_list and exp_labels must have the same length")
  }

  # ---- Step 1: Collect significant SNPs from all exposures ----
  all_sig_snps <- list()
  valid_indices <- c()

  for(i in seq_along(exp_dt_list)){
    exp_dt <- exp_dt_list[[i]]

    # Check if exposure has rsID
    if(!any(grepl("^rs\\d+$", exp_dt$SNP))){
      next  # Skip exposures without rsID
    }

    # Get significant SNPs
    sig <- exp_dt[pval < p_thresh, .(rsid = SNP, pval, exp_id = i)]

    if(nrow(sig) > 0){
      all_sig_snps[[length(all_sig_snps) + 1]] <- sig
      valid_indices <- c(valid_indices, i)
    }
  }

  if(length(all_sig_snps) == 0){
    message("  [Batch Clumping] No valid exposures with significant SNPs")
    return(setNames(vector("list", length(exp_dt_list)), exp_labels))
  }

  # ---- Step 2: Combine all significant SNPs ----
  combined_snps <- rbindlist(all_sig_snps)

  # Keep minimum p-value for each SNP across all exposures
  combined_snps <- combined_snps[, .(pval = min(pval), exp_id = first(exp_id)), by = rsid]

  message("  [Batch Clumping] Combined ", nrow(combined_snps), " significant SNPs from ",
          length(valid_indices), " exposures")

  # ---- Step 3: Single PLINK clumping call ----
  # Check cache first
  cache_key <- make_cache_key(combined_snps$rsid, p_thresh, r2, kb)

  if(!is.null(CLUMP_CACHE[[cache_key]])){
    CLUMP_CACHE_HITS <<- CLUMP_CACHE_HITS + 1L
    clumped <- CLUMP_CACHE[[cache_key]]
    message("  [Batch Clumping] Cache hit! Skipping PLINK call")
  } else {
    CLUMP_CACHE_MISSES <<- CLUMP_CACHE_MISSES + 1L

    # Try local clumping
    if(isTRUE(USE_LOCAL_CLUMP)){
      message("  [Batch Clumping] 开始 PLINK clumping (", nrow(combined_snps), " SNPs)...")
      clumped <- try(ieugwasr::ld_clump_local(
        tibble::as_tibble(transform(combined_snps, id = "batch_clump")),
        bfile = EUR_BFILE,
        plink_bin = PLINK_BIN,
        clump_kb = kb,
        clump_r2 = r2,
        clump_p = p_thresh
      ), silent = TRUE)
      message("  [Batch Clumping] PLINK clumping 完成")

      if(inherits(clumped, "try-error")){
        message("  [Batch Clumping] Local clumping failed, trying remote...")
        clumped <- try(ieugwasr::ld_clump(
          d = tibble::as_tibble(transform(combined_snps, id = "batch_clump")),
          clump_kb = kb,
          clump_r2 = r2,
          clump_p = p_thresh,
          pop = "EUR"
        ), silent = TRUE)
      }
    } else {
      clumped <- try(ieugwasr::ld_clump(
        d = tibble::as_tibble(transform(combined_snps, id = "batch_clump")),
        clump_kb = kb,
        clump_r2 = r2,
        clump_p = p_thresh,
        pop = "EUR"
      ), silent = TRUE)
    }

    if(inherits(clumped, "try-error")){
      message("  [Batch Clumping] Both local and remote clumping failed")
      return(setNames(vector("list", length(exp_dt_list)), exp_labels))
    }

    clumped <- as.data.table(clumped)

    # Store in cache
    CLUMP_CACHE[[cache_key]] <- clumped
  }

  message("  [Batch Clumping] Clumped to ", nrow(clumped), " independent SNPs")

  # ---- Step 4: Assign clumped SNPs back to each exposure ----
  results <- setNames(vector("list", length(exp_dt_list)), exp_labels)

  for(i in seq_along(exp_dt_list)){
    exp_dt <- exp_dt_list[[i]]
    exp_snps <- clumped[rsid %in% exp_dt$SNP]
    results[[i]] <- exp_snps
  }

  results
}

#' Select Independent Genetic Instruments via LD Clumping
#'
#' Performs LD clumping to select independent SNPs as instrumental variables
#' - Filters SNPs by p-value threshold
#' - Removes SNPs in LD (r² threshold)
#' - Tries local clumping first, falls back to remote if needed
#'
#' NOTE: For batch processing of multiple exposures, use batch_clump_exposures() instead
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
  sig <- sig[!is.na(rsid) & nzchar(rsid)]
  if(nrow(sig) == 0) return(data.table())

  # Cache clumping results for repeated calls (e.g., reverse MR)
  sig_rsids <- sort(unique(sig$rsid))
  cache_key <- make_cache_key("select_instruments", sig_rsids, p_thresh, r2, kb,
                              USE_LOCAL_CLUMP, EUR_BFILE, PLINK_BIN)
  if(!is.null(CLUMP_CACHE[[cache_key]])){
    CLUMP_CACHE_HITS <<- CLUMP_CACHE_HITS + 1L
    return(copy(CLUMP_CACHE[[cache_key]]))
  }
  CLUMP_CACHE_MISSES <<- CLUMP_CACHE_MISSES + 1L
  
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
      if(nrow(res_local)>0){
        CLUMP_CACHE[[cache_key]] <- res_local
        return(res_local)
      }
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
  
  if(!inherits(res_remote, "try-error")){
    res_remote <- as.data.table(res_remote)
    if(nrow(res_remote)>0){
      CLUMP_CACHE[[cache_key]] <- res_remote
    } else {
      # Cache empty result to avoid repeated remote calls
      CLUMP_CACHE[[cache_key]] <- res_remote
    }
    return(res_remote)
  }
  
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

# Cached mediator formatting helper (reduces repeated format_data calls)
get_mediator_out_fmt <- function(med_dt, med_name, med_file = NULL){
  if(!isTRUE(MEDIATOR_FMT_CACHE_ENABLED) || is.null(med_file) || !nzchar(med_file)){
    return(to_outcome_format(med_dt, med_name))
  }
  med_path <- normalizePath(med_file, winslash = "/", mustWork = FALSE)
  med_mtime <- as.numeric(file.info(med_path)$mtime)
  cache_key <- make_cache_key("med_out_fmt", med_path, med_mtime, med_name)
  if(!is.null(MEDIATOR_FMT_CACHE[[cache_key]])){
    MEDIATOR_FMT_CACHE_HITS <<- MEDIATOR_FMT_CACHE_HITS + 1L
    return(MEDIATOR_FMT_CACHE[[cache_key]])
  }
  MEDIATOR_FMT_CACHE_MISSES <<- MEDIATOR_FMT_CACHE_MISSES + 1L
  res <- to_outcome_format(med_dt, med_name)
  MEDIATOR_FMT_CACHE[[cache_key]] <- res
  res
}

#' Harmonize Exposure and Outcome Data
#'
#' Aligns SNP effects to the same allele across exposure and outcome
#' Uses TwoSampleMR::harmonise_data with action=2 (default harmonization)
#'
#' PERFORMANCE: Results are cached based on SNP sets to avoid redundant harmonization
#' Especially beneficial for mediation analysis where same pairs are reused
#'
#' @param exp_fmt Formatted exposure data
#' @param out_fmt Formatted outcome data
#' @return Harmonized data frame
harmonise_xy <- function(exp_fmt, out_fmt){
  # ---- PERFORMANCE OPTIMIZATION: Check cache first ----
  # 激进优化：简化缓存键生成，避免不必要的排序和去重
  # 使用更快的哈希方法
  cache_key <- make_cache_key(
    exp_fmt$exposure[1],
    out_fmt$outcome[1],
    length(exp_fmt$SNP),
    length(out_fmt$SNP)
  )

  if(!is.null(HARMONIZE_CACHE[[cache_key]])){
    HARMONIZE_CACHE_HITS <<- HARMONIZE_CACHE_HITS + 1L
    return(HARMONIZE_CACHE[[cache_key]])
  }

  HARMONIZE_CACHE_MISSES <<- HARMONIZE_CACHE_MISSES + 1L

  # ---- Perform harmonization (cache miss) ----
  result <- TwoSampleMR::harmonise_data(exp_fmt, out_fmt, action = 2)

  # Store in cache
  HARMONIZE_CACHE[[cache_key]] <- result

  result
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
  loo <- NULL
  if(isTRUE(RUN_LOO)){
    loo <- tryCatch(
      TwoSampleMR::mr_leaveoneout(hdat, method = mr_ivw), 
      error = function(e) NULL
    )
  }
  
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
#' @param preloaded_outs Optional preloaded outcomes list from preload_outcomes()
#' @param preloaded_expos Optional preloaded exposures list from preload_exposures()
#' @param preloaded_expos Optional preloaded exposures list from preload_exposures()
#' @return data.table with all UVMR results (invisibly)
run_comprehensive_uvmr <- function(
    expos_files = all_expos_files,
    outcomes_files = all_out_files,
    out_csv = file.path(DIR_RES_TRIAL, "uvmr_comprehensive_results.csv"),
    preloaded_outs = NULL,
    preloaded_expos = NULL
){
  message("\n[COMPREHENSIVE UVMR ANALYSIS] Exposures=", length(expos_files), "; Outcomes=", length(outcomes_files),
          "; p<", P_THRESH, "; r2=", CLUMP_R2, "; kb=", CLUMP_KB)

  preloaded <- NULL
  if(!is.null(preloaded_expos)){
    preloaded <- preloaded_expos
    if(length(expos_files) > 0){
      preloaded <- preloaded[intersect(expos_files, names(preloaded_expos))]
    }
    if(length(preloaded) > 0){
      message("[INFO] Using preloaded exposures (", length(preloaded), ")")
      report_progress("uvmr_clump", length(preloaded), length(expos_files))
    } else {
      preloaded <- NULL
    }
  }
  if(is.null(preloaded)){
    message("[INFO] Preloading exposures and running PLINK clumping in parallel...")
    preloaded <- preload_exposures(expos_files, keep_exp_dt = FALSE, report_step_prefix = "uvmr")
  }
  if(length(preloaded) == 0){
    message("[ERROR] No exposures with valid IVs after clumping.")
    report_progress("uvmr_pairs", 0, 1)
    return(invisible(NULL))
  }

  if(is.null(preloaded_outs)){
    message("[INFO] Preloading outcomes once to avoid repeated reads...")
    preloaded_outs <- preload_outcomes(outcomes_files, report_step_prefix = "uvmr")
  } else if(length(outcomes_files) > 0){
    loaded_n <- sum(outcomes_files %in% names(preloaded_outs))
    report_progress("uvmr_load", loaded_n, length(outcomes_files))
  }
  if(length(preloaded_outs) == 0){
    message("[ERROR] No outcomes could be loaded.")
    report_progress("uvmr_pairs", 0, 1)
    return(invisible(NULL))
  }

  # Use only exposures/outcomes that successfully preloaded for progress accounting
  valid_expos_files <- Filter(function(f) f %in% names(preloaded), expos_files)
  valid_out_files <- Filter(function(f) f %in% names(preloaded_outs), outcomes_files)
  # Calculate actual pairs that will be processed
  total_pairs_planned <- length(valid_expos_files) * length(valid_out_files)
  if(total_pairs_planned == 0) total_pairs_planned <- 1  # Avoid division by zero

  all_results <- list()
  total_pairs <- 0
  successful_pairs <- 0
  pairs_done <- 0
  # Report initial progress with correct total
  report_progress("uvmr_pairs", pairs_done, total_pairs_planned)
  
  for(exp_file in valid_expos_files){
    pre <- preloaded[[exp_file]]
    exp_name <- pre$exp_name
    exp_fmt <- pre$exp_fmt
    total_ivs_exp <- pre$total_ivs
    
    message("\n[Exposure ", match(exp_file, expos_files), "/", length(expos_files), "] ", exp_name)
    message("  选择了 ", total_ivs_exp, " 个独立工具变量")

    # Per-outcome worker function to enable parallel processing
    process_outcome <- function(out_file){
      if(!(out_file %in% names(preloaded_outs))){
        message("    • ", basename(out_file), "：未预加载，跳过。")
        return(NULL)
      }

      out_pre <- preloaded_outs[[out_file]]
      out_name <- out_pre$out_name
      out_fmt <- out_pre$out_fmt

      # Harmonize data
      hdat <- harmonise_xy(exp_fmt, out_fmt)
      hdat <- hdat[hdat$remove==FALSE & is.finite(hdat$beta.exposure) & is.finite(hdat$beta.outcome), ]

      message("    • ", out_name, ": Harmonized SNPs = ", nrow(hdat))

      if(nrow(hdat) < MIN_SNPS){
        message("      Insufficient SNPs after harmonization, skipping")
        return(NULL)
      }

      # Run comprehensive MR analysis
      res <- run_uvmr_comprehensive(hdat, exp_name, out_name)

      if(is.null(res) || nrow(res) == 0){
        return(NULL)
      }

      # Add extra information
      res[, `:=`(
        p_threshold = P_THRESH,
        clump_r2 = CLUMP_R2,
        clump_kb = CLUMP_KB,
        total_ivs = total_ivs_exp,
        harmonised_snps = nrow(hdat)
      )]

    # Run MRlap analysis for sample overlap correction
    # MRlap uses FULL GWAS files (not just harmonized IVs)
    if(MRLAP_AVAILABLE && isTRUE(RUN_MRLAP)){
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

      message("      ✓ Analysis completed")
      res
    }

    # Only process valid outcome files that were successfully preloaded
    out_files_run <- valid_out_files
    total_pairs <- total_pairs + length(out_files_run)
    outcome_results <- parallel_map(out_files_run, process_outcome, cores = PARALLEL_WORKERS)

    # Update progress for each outcome processed
    for(or in outcome_results){
      if(!is.null(or) && nrow(or) > 0){
        all_results[[length(all_results) + 1]] <- or
        successful_pairs <- successful_pairs + 1
      }
      pairs_done <- pairs_done + 1
      report_progress("uvmr_pairs", pairs_done, total_pairs_planned)
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
    report_progress("uvmr_complete", total_pairs, total_pairs_planned)
    return(invisible(final_results))
  } else {
    message("\n[警告] 未生成任何有效结果")
    report_progress("uvmr_complete", pairs_done, total_pairs_planned)
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
prepare_mvmr_data <- function(expo_list, expo_names, covar_list, covar_names, out_dt, out_name,
                              combined_ivs = NULL, out_fmt = NULL){
  
  # Step 1: Select combined IVs from exposures AND covariates (once per group if provided)
  computed_here <- is.null(combined_ivs)
  if(computed_here){
    all_gwas <- c(expo_list, covar_list)
    combined_ivs <- select_combined_ivs_mvmr(all_gwas)
  }
  
  if(length(combined_ivs) < MIN_SNPS_MVMR){
    return(NULL)
  }
  
  if(computed_here){
    message("    Selected ", length(combined_ivs), " combined IVs (from exposures + covariates)")
  }
  
  # Step 2: Extract these IVs from each exposure and covariate
  # Use chr:pos for matching
  if(is.null(out_fmt)){
    if(is.null(out_dt)) stop("prepare_mvmr_data requires out_dt or out_fmt")
    out_dt <- use_chrpos_id(out_dt)
    out_fmt <- to_outcome_format(out_dt, out_name)
  }
  
  harmonised_list <- list()
  kept_var_names <- character()
  kept_var_types <- character()
  var_names <- c(expo_names, covar_names)
  all_var_list <- c(expo_list, covar_list)
  
  for(i in seq_along(all_var_list)){
    var_dt <- all_var_list[[i]]
    var_name <- var_names[i]
    var_type <- if(i <= length(expo_names)) "Exposure" else "Covariate"
    
    # Extract SNPs that are in combined IV set
    var_iv <- var_dt[SNP %in% combined_ivs]
    
    if(nrow(var_iv) == 0) next
    
    var_iv <- use_chrpos_id(var_iv)
    var_fmt <- to_exposure_format(var_iv, var_name)
    h <- harmonise_xy(var_fmt, out_fmt)
    h <- h[h$remove==FALSE & is.finite(h$beta.exposure) & is.finite(h$beta.outcome), ]
    
    if(nrow(h) > 0){
      harmonised_list[[length(harmonised_list) + 1]] <- h
      kept_var_names <- c(kept_var_names, var_name)
      kept_var_types <- c(kept_var_types, var_type)
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
  
  list(
    bx = bx_matrix,
    bxse = bxse_matrix,
    by = by,
    byse = byse,
    snps = common_snps,
    n_snps = length(common_snps),
    var_names = kept_var_names,
    var_types = kept_var_types,
    expo_indices = which(kept_var_types == "Exposure"),
    covar_indices = which(kept_var_types == "Covariate")
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
  var_types <- mvmr_data$var_types
  if(is.null(var_types) || length(var_types) != n_vars){
    n_expo <- length(mvmr_data$expo_indices)
    n_covar <- length(mvmr_data$covar_indices)
    var_types <- c(rep("Exposure", n_expo), rep("Covariate", n_covar))
  }
  n_expo <- sum(var_types == "Exposure")
  n_covar <- sum(var_types == "Covariate")
  
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
    if(is.atomic(mvmr_res)){
      message("      [WARNING] ", method_name, " returned atomic vector; skipping.")
      return(NULL)
    }
    if(!is.list(mvmr_res) || is.null(mvmr_res$Estimate) || is.null(mvmr_res$StdError) || is.null(mvmr_res$Pvalue)){
      message("      [WARNING] ", method_name, " result structure unexpected; skipping.")
      return(NULL)
    }
    
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
      variable_type = var_types,
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
    } else if(MRLAP_AVAILABLE && !isTRUE(RUN_MRLAP)){
      message("      • MRlap skipped (MR_RUN_MRLAP=FALSE)")
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
#' @param preloaded_outs Optional preloaded outcomes list from preload_outcomes()
#' @return data.table with MVMR results (invisibly)
run_comprehensive_mvmr <- function(
    expos_files = all_expos_files,
    covar_files = all_covar_files,
    outcomes_files = all_out_files,
    out_csv = file.path(DIR_RES_TRIAL, "mvmr_comprehensive_results.csv"),
    group_by = "category",
    max_exposures_per_group = 10,
    preloaded_outs = NULL,
    preloaded_expos = NULL
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
  mvmr_done <- 0

  # Initial progress report - will update as we process
  report_progress("mvmr_start", 0, 1)  # Start with placeholder
  
  # Group exposures by category
  if(group_by == "category"){
    if(length(expos_files) > 0){
      expo_groups <- list(User = expos_files)
    } else {
      # Create exposure groups based on directory
      expo_groups <- list(
        Proteome = list_gwas_files(DIR_EXPO01),
        Metabolic = list_gwas_files(DIR_EXPO02),
        Inflammatory = list_gwas_files(DIR_EXPO03)
      )
    }
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
  covar_files_used <- character()
  
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
      covar_files_used <- c(covar_files_used, cov_file)
      message("  ✓ Loaded: ", cov_name)
    }
    message("  Total covariates loaded: ", length(covar_list))
  }

  if(is.null(preloaded_outs)){
    message("\n[Preloading Outcomes for MVMR]")
    preloaded_outs <- preload_outcomes(outcomes_files, report_step_prefix = "mvmr")
  } else if(length(outcomes_files) > 0){
    loaded_n <- sum(outcomes_files %in% names(preloaded_outs))
    report_progress("mvmr_load", loaded_n, length(outcomes_files))
  }
  if(length(preloaded_outs) == 0){
    message("  [SKIP] No outcomes could be loaded.")
    report_progress("mvmr_pairs", 0, 1)
    return(invisible(NULL))
  }
  out_files_run <- intersect(outcomes_files, names(preloaded_outs))
  total_outcomes <- length(out_files_run)
  total_analyses_planned <- total_outcomes * max(1, length(expo_groups))
  report_progress("mvmr_pairs", mvmr_done, total_analyses_planned)
  
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
    expo_files_used <- character()
    
    for(exp_file in group_files){
      exp_name <- basename(exp_file)
      exp_name <- clean_filename(exp_name)

      exp_dt <- NULL
      if(!is.null(preloaded_expos) &&
         exp_file %in% names(preloaded_expos) &&
         !is.null(preloaded_expos[[exp_file]]$exp_dt)){
        exp_dt <- preloaded_expos[[exp_file]]$exp_dt
      } else {
        exp_dt <- try(read_gwas(exp_file), silent = TRUE)
      }
      if(inherits(exp_dt, "try-error")) next
      
      # Store FULL GWAS data (not just selected IVs)
      # IV selection will be done jointly with covariates
      expo_list[[length(expo_list) + 1]] <- exp_dt
      expo_names <- c(expo_names, exp_name)
      expo_files_used <- c(expo_files_used, exp_file)
    }
    
    if(length(expo_list) < 1){
      message("  [SKIP] No valid exposures in group")
      next
    }
    
    message("  Valid exposures in group: ", length(expo_list))

    # Pre-compute combined IVs once per exposure group (shared across outcomes)
    cache_key <- make_file_cache_key(c(expo_files_used, covar_files_used), P_THRESH, CLUMP_R2, CLUMP_KB)
    if(!is.na(cache_key) && !is.null(MVMR_IVS_CACHE[[cache_key]])){
      MVMR_IVS_CACHE_HITS <<- MVMR_IVS_CACHE_HITS + 1L
      combined_ivs <- MVMR_IVS_CACHE[[cache_key]]
    } else {
      MVMR_IVS_CACHE_MISSES <<- MVMR_IVS_CACHE_MISSES + 1L
      combined_ivs <- select_combined_ivs_mvmr(c(expo_list, covar_list))
      if(!is.na(cache_key)){
        MVMR_IVS_CACHE[[cache_key]] <- combined_ivs
      }
    }
    if(length(combined_ivs) < MIN_SNPS_MVMR){
      message("  [SKIP] Combined IVs < ", MIN_SNPS_MVMR, " for group '", group_name, "'")
      # Advance progress for this group's outcomes to keep totals consistent
      mvmr_done <- mvmr_done + length(out_files_run)
      report_progress("mvmr_pairs", mvmr_done, total_analyses_planned)
      next
    }
    message("  Selected ", length(combined_ivs), " combined IVs (from exposures + covariates)")

    # Worker for each outcome (enables parallel execution)
    process_outcome <- function(out_file){
      if(!(out_file %in% names(preloaded_outs))){
        message("    • ", basename(out_file), "：未预加载，跳过。")
        return(NULL)
      }

      out_pre <- preloaded_outs[[out_file]]
      out_name <- out_pre$out_name
      out_dt <- out_pre$out_dt

      # Prepare MVMR data WITH covariates
      # This will select combined IVs from exposures + covariates
      mvmr_data <- prepare_mvmr_data(
        expo_list, expo_names, covar_list, covar_names,
        out_dt, out_name,
        combined_ivs = combined_ivs,
        out_fmt = out_pre$out_fmt
      )
      
      if(is.null(mvmr_data)){
        message("    • ", out_name, ": Insufficient common SNPs")
        return(NULL)
      }
      
      message("    • ", out_name, ": Common SNPs = ", mvmr_data$n_snps)
      
      # Run MVMR (returns covariate-adjusted results for exposures)
      mvmr_res <- run_mvmr_analysis(mvmr_data, out_name)
      
      if(is.null(mvmr_res) || nrow(mvmr_res) == 0){
        return(NULL)
      }
      
      # Add group information
      mvmr_res[, exposure_group := group_name]
      mvmr_res[, `:=`(
        p_threshold = P_THRESH,
        clump_r2 = CLUMP_R2,
        clump_kb = CLUMP_KB
      )]
      
      message("      ✓ MVMR completed (covariate-adjusted)")
      mvmr_res
    }
    
    # Note: out_files_run is already defined above
    total_analyses <- total_analyses + length(out_files_run)
    mvmr_list <- parallel_map(out_files_run, process_outcome, cores = PARALLEL_WORKERS)

    for(res in mvmr_list){
      if(!is.null(res) && nrow(res) > 0){
        all_results[[length(all_results) + 1]] <- res
        successful_analyses <- successful_analyses + 1
      }
      mvmr_done <- mvmr_done + 1
      report_progress("mvmr_pairs", mvmr_done, total_analyses_planned)
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
    report_progress("mvmr_complete", successful_analyses, total_analyses)
    return(invisible(final_results))
  } else {
    message("\n[WARNING] No valid MVMR results generated")
    report_progress("mvmr_complete", 0, 1)
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
    out_csv = file.path(DIR_RES_TRIAL, "mediation_comprehensive_results.csv"),
    preloaded_outs = NULL,
    preloaded_expos = NULL
){
  message("\n[COMPREHENSIVE MEDIATION ANALYSIS] 暴露=", length(expos_files),
          "；中介=", length(medi_files), "；结局=", length(outcomes_files))

  preloaded <- NULL
  if(!is.null(preloaded_expos)){
    preloaded <- preloaded_expos
    if(length(expos_files) > 0){
      preloaded <- preloaded[intersect(expos_files, names(preloaded_expos))]
    }
    if(length(preloaded) > 0){
      # Ensure exp_dt is available for reverse MR
      missing_dt <- vapply(preloaded, function(x) is.null(x$exp_dt), logical(1))
      if(any(missing_dt)){
        message("[INFO] Preloaded exposures missing full GWAS; reloading for mediation...")
        preloaded <- NULL
      } else {
        message("[INFO] Using preloaded exposures (", length(preloaded), ")")
        report_progress("mediation_clump", length(preloaded), length(expos_files))
      }
    } else {
      preloaded <- NULL
    }
  }
  if(is.null(preloaded)){
    message("[INFO] Preloading exposures (with full GWAS) and running PLINK clumping in parallel...")
    preloaded <- preload_exposures(expos_files, keep_exp_dt = TRUE, report_step_prefix = "mediation")
  }
  if(length(preloaded) == 0){
    message("[ERROR] No exposures with valid IVs after clumping.")
    return(invisible(NULL))
  }

  if(is.null(preloaded_outs)){
    message("[INFO] Preloading outcomes once to avoid repeated reads...")
    preloaded_outs <- preload_outcomes(outcomes_files, report_step_prefix = "mediation")
  } else if(length(outcomes_files) > 0){
    loaded_n <- sum(outcomes_files %in% names(preloaded_outs))
    report_progress("mediation_load", loaded_n, length(outcomes_files))
  }
  if(length(preloaded_outs) == 0){
    message("[ERROR] No outcomes could be loaded.")
    return(invisible(NULL))
  }

  # Calculate planned triplets based on valid preloaded files
  valid_expos_files <- Filter(function(f) f %in% names(preloaded), expos_files)
  valid_out_files <- Filter(function(f) f %in% names(preloaded_outs), outcomes_files)
  total_triplets_planned <- length(valid_expos_files) * length(medi_files) * length(valid_out_files)
  if(total_triplets_planned == 0) total_triplets_planned <- 1  # Avoid division by zero

  all_results <- list()
  total_triplets <- 0
  processed_triplets <- 0
  successful_triplets <- 0
  report_progress("mediation_start", 0, total_triplets_planned)
  
  for(exp_file in expos_files){
    if(!(exp_file %in% names(preloaded))){
      message("\n[Exposure] Skipping (no IVs): ", basename(exp_file))
      next
    }
    pre <- preloaded[[exp_file]]
    exp_name <- pre$exp_name
    exp_dt <- pre$exp_dt
    exp_fmt <- pre$exp_fmt
    ivs_exp_n <- pre$total_ivs
    exp_out_fmt <- NULL
    if(!is.null(exp_dt) && nrow(exp_dt) > 0){
      exp_out_fmt <- to_outcome_format(use_chrpos_id(exp_dt), exp_name)
    }
    
    message("\n[Exposure ", match(exp_file, expos_files), "/", length(expos_files), "] ", exp_name)
    message("  选择了 ", ivs_exp_n, " 个独立工具变量")

    # Precompute Exp -> Out harmonization and total effects (shared across mediators)
    compute_exp_out <- function(out_file){
      out_pre <- preloaded_outs[[out_file]]
      if(is.null(out_pre)) return(NULL)

      out_name <- out_pre$out_name
      out_fmt <- out_pre$out_fmt

      h_exp_out <- harmonise_xy(exp_fmt, out_fmt)
      h_exp_out <- h_exp_out[h_exp_out$remove==FALSE & is.finite(h_exp_out$beta.exposure) & is.finite(h_exp_out$beta.outcome), ]

      if(nrow(h_exp_out) < MIN_SNPS_MVMR){
        return(NULL)
      }

      mtot <- try(TwoSampleMR::mr(h_exp_out, method_list = c("mr_ivw_mre")), silent = TRUE)
      if(inherits(mtot, "try-error") || nrow(mtot)==0) {
        return(NULL)
      }

      list(
        out_file = out_file,
        out_name = out_name,
        out_fmt = out_fmt,
        h_exp_out = h_exp_out,
        exp_snps = unique(h_exp_out$SNP),
        betaT = mtot$b[1],
        seT = mtot$se[1],
        pvalT = mtot$pval[1]
      )
    }

    exp_out_list <- parallel_map(valid_out_files, compute_exp_out, cores = PARALLEL_WORKERS)
    exp_out_cache <- setNames(vector("list", length(valid_out_files)), valid_out_files)
    for(item in exp_out_list){
      if(!is.null(item)){
        exp_out_cache[[item$out_file]] <- item
      }
    }
    
    for(med_idx in seq_along(medi_files)){
      med_file <- medi_files[[med_idx]]
      med_name <- basename(med_file)
      med_name <- clean_filename(med_name)
      if(length(medi_files) > 0){
        # Report the current mediator while keeping the fraction based on completed mediators.
        report_progress("mediation_mediators", med_idx - 1L, length(medi_files), med_name)
      }
      message("  [Mediator] ", med_name)
      
      med_dt <- try(read_gwas(med_file), silent = TRUE)
      if(inherits(med_dt, "try-error")) { 
        message("    中介读取失败，跳过。")
        processed_triplets <- processed_triplets + length(valid_out_files)
        report_progress("mediation_triplets", processed_triplets, total_triplets_planned, med_name)
        report_progress("mediation_mediators", med_idx, length(medi_files), med_name)
        next 
      }
      med_dt <- use_chrpos_id(med_dt)  # 中介无 rsID 也能参与
      
      # Step1: exp -> med（用暴露 IV，chr:pos 键）
      med_fmt_as_out <- get_mediator_out_fmt(med_dt, med_name, med_file)
      h1 <- harmonise_xy(exp_fmt, med_fmt_as_out)
      h1 <- h1[h1$remove==FALSE & is.finite(h1$beta.exposure) & is.finite(h1$beta.outcome), ]
      
      if(nrow(h1) < MIN_SNPS_STEP1) { 
        message("    协调后SNP不足（Step1: Exp->Med）")
        processed_triplets <- processed_triplets + length(valid_out_files)
        report_progress("mediation_triplets", processed_triplets, total_triplets_planned, med_name)
        report_progress("mediation_mediators", med_idx, length(medi_files), med_name)
        next 
      }
      
      m1 <- try(TwoSampleMR::mr(h1, method_list = c("mr_ivw_mre")), silent = TRUE)
      if(inherits(m1, "try-error") || nrow(m1)==0) {
        message("    Step1 MR失败")
        processed_triplets <- processed_triplets + length(valid_out_files)
        report_progress("mediation_triplets", processed_triplets, total_triplets_planned, med_name)
        report_progress("mediation_mediators", med_idx, length(medi_files), med_name)
        next
      }
      
      beta1 <- m1$b[1]
      se1 <- m1$se[1]
      pval1 <- m1$pval[1]
      
      # ==== Reverse MR: Test for bidirectionality (Med -> Exp) ====
      # Required to validate mediation pathway direction
      if(isTRUE(SKIP_REVERSE_MR)){
        message("    Reverse MR skipped (MR_SKIP_REVERSE_MR=TRUE)")
        reverse_res <- data.table(
          reverse_direction = paste0(med_name, " → ", exp_name),
          reverse_n_ivs = NA_integer_,
          reverse_n_snps = NA_integer_,
          reverse_beta = NA_real_,
          reverse_pval = NA_real_,
          bidirectional = "Skipped"
        )
      } else {
        message("    Testing reverse causation (Med → Exp)...")
        reverse_res <- run_reverse_mr(med_dt, exp_dt, med_name, exp_name,
                                      outcome_file = med_file, exposure_out_fmt = exp_out_fmt)
      }
      
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

      # Worker per outcome to enable parallel execution
      process_mediation_outcome <- function(out_file){
        cache <- exp_out_cache[[out_file]]
        if(is.null(cache)){
          message("      • ", basename(out_file), "：Exp->Out预计算失败，跳过。")
          return(NULL)
        }

        out_name <- cache$out_name
        out_fmt <- cache$out_fmt
        h_exp_out <- cache$h_exp_out
        betaT <- cache$betaT
        seT <- cache$seT
        pvalT <- cache$pvalT

        # 中介-结局：联合 SNP = 暴露 IV 的 chr:pos
        med_subset <- med_dt[SNP %in% cache$exp_snps]
        if(nrow(med_subset) < MIN_SNPS_MVMR) {
          message("      • ", out_name, "：中介SNP不足")
          return(NULL)
        }

        med_fmt2 <- to_exposure_format(med_subset, med_name)
        h_med_out <- harmonise_xy(med_fmt2, out_fmt)
        h_med_out <- h_med_out[h_med_out$remove==FALSE & is.finite(h_med_out$beta.exposure) & is.finite(h_med_out$beta.outcome), ]
        if(nrow(h_med_out) < MIN_SNPS_MVMR) {
          message("      • ", out_name, "：协调后SNP不足（Med->Out）")
          return(NULL)
        }

        # 对齐到同一等位基因方向并做 MVMR
        M <- align_for_mvmr(h_exp_out, h_med_out)
        if(nrow(M) < MIN_SNPS_MVMR) {
          message("      • ", out_name, "：MVMR对齐后SNP不足")
          return(NULL)
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
          return(NULL)
        }

        m2 <- try(
          MendelianRandomization::mr_mvivw(mvinput, model="default", correl=FALSE, distribution="normal"),
          silent = TRUE
        )

        if(inherits(m2, "try-error")) {
          message("      • ", out_name, "：MVMR分析失败")
          return(NULL)
        }

        est_tbl <- as.data.frame(m2@Estimate)
        se_tbl <- as.data.frame(m2@StdError)
        pval_tbl <- as.data.frame(m2@Pvalue)

        if(nrow(est_tbl) < 2) {
          message("      • ", out_name, "：MVMR结果不完整")
          return(NULL)
        }

        beta2_med <- est_tbl[1,1]  # mediator 的直接效应
        se2_med <- se_tbl[1,1]
        pval2_med <- pval_tbl[1,1]

        beta2_exp <- est_tbl[2,1]  # exposure 的直接效应
        se2_exp <- se_tbl[2,1]
        pval2_exp <- pval_tbl[2,1]

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

        message("      • ", out_name, "：✓ 分析完成")
        res_row
      }

      # Only process valid outcome files
      out_files_run <- valid_out_files
      total_triplets <- total_triplets + length(out_files_run)
      mediation_list <- parallel_map(out_files_run, process_mediation_outcome, cores = PARALLEL_WORKERS)

      for(res_row in mediation_list){
        if(!is.null(res_row) && nrow(res_row) > 0){
          all_results[[length(all_results) + 1]] <- res_row
          successful_triplets <- successful_triplets + 1
        }
      }
      processed_triplets <- processed_triplets + length(out_files_run)
      report_progress("mediation_triplets", processed_triplets, total_triplets_planned, med_name)
      report_progress("mediation_mediators", med_idx, length(medi_files), med_name)
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
    report_progress("mediation_complete", processed_triplets, total_triplets_planned,
                    paste0("successful=", successful_triplets))
    return(invisible(final_results))
  } else {
    message("\n[警告] 未生成任何有效中介分析结果")
    report_progress("mediation_complete", processed_triplets, total_triplets_planned,
                    paste0("successful=", successful_triplets))
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
  report_log("INFO", paste("Running MRlap for", exposure_name, "->", outcome_name))

  # Store error message for diagnostics
  error_msg <- NA_character_

  # Temporary files for standardized GWAS (cleanup on exit)
  # MRlap expects files named as {exposure_name}.sumstats and {outcome_name}.sumstats
  # So we create temp files with appropriate names in a temp directory
  mrlap_temp_dir <- tempfile(pattern = "mrlap_dir_")
  dir.create(mrlap_temp_dir, recursive = TRUE)

  exp_temp <- file.path(mrlap_temp_dir, paste0(exposure_name, ".sumstats"))
  out_temp <- file.path(mrlap_temp_dir, paste0(outcome_name, ".sumstats"))
  exp_temp_gz <- paste0(exp_temp, ".gz")
  out_temp_gz <- paste0(out_temp, ".gz")

  on.exit({
    if(file.exists(exp_temp)) unlink(exp_temp)
    if(file.exists(out_temp)) unlink(out_temp)
    if(file.exists(exp_temp_gz)) unlink(exp_temp_gz)
    if(file.exists(out_temp_gz)) unlink(out_temp_gz)
    if(dir.exists(mrlap_temp_dir)) unlink(mrlap_temp_dir, recursive = TRUE)
  }, add = TRUE)

  # Preprocess exposure and outcome files for MRlap compatibility
  # MRlap expects specific column names (uppercase: POS, CHR, SNP, etc.)
  # read_gwas() already standardizes columns to lowercase (pos, chr, etc.)
  # but MRlap may require uppercase versions

  # Debug: Log file paths
  report_log("DEBUG", paste("MRlap preprocessing | Exposure file:", basename(exp_file)))
  report_log("DEBUG", paste("MRlap preprocessing | Outcome file:", basename(out_file)))
  report_log("DEBUG", paste("MRlap temp dir:", mrlap_temp_dir))

  tryCatch({
    # Read and standardize exposure file
    exp_dt <- read_gwas(exp_file)

    # Ensure required columns for MRlap (uppercase versions)
    # read_gwas() produces: pos, chr, beta, se, pval, effect_allele, other_allele, samplesize
    # MRlap expects: POS, CHR, BETA, SE, P (or pval), A1, A2, SNP, N

    # Rename lowercase standardized columns to MRlap-expected uppercase
    col_map <- c(
      "pos" = "POS",
      "chr" = "CHR",
      "beta" = "BETA",
      "se" = "SE",
      "pval" = "P",
      "effect_allele" = "A1",
      "other_allele" = "A2",
      "samplesize" = "N"
    )

    for(old_col in names(col_map)){
      new_col <- col_map[[old_col]]
      if(old_col %in% names(exp_dt)){
        setnames(exp_dt, old_col, new_col)
      }
    }

    # Add sample size N if missing (required by MRlap for LDSC)
    # Use fixed N=3301 for proteome/metabolome/inflammatory protein exposures
    if(!"N" %in% names(exp_dt)){
      message("      [MRlap] Exposure missing N column, adding N=3301...")
      report_log("DEBUG", paste("MRlap: Exposure", exposure_name, "missing N column, adding N=3301"))
      exp_dt[, N := 3301]
    }

    # Check for required MRlap columns and provide helpful error message
    # N (sample size) is required by MRlap for LDSC calculations
    required_cols <- c("POS", "CHR", "BETA", "SE", "P", "A1", "A2", "SNP", "N")
    missing_exp <- setdiff(required_cols, names(exp_dt))
    if(length(missing_exp) > 0){
      error_msg <<- paste0("Exposure missing required columns for MRlap: ", paste(missing_exp, collapse=", "),
                          ". Available columns: ", paste(names(exp_dt), collapse=", "))
      message("      [MRlap SKIP] ", error_msg)
      report_log("ERROR", paste("MRlap preprocessing:", error_msg))
      report_error(
        error_type = "MRlap_missing_columns",
        details = error_msg,
        exposure = exposure_name,
        outcome = outcome_name,
        location = "run_mrlap_analysis"
      )
      return(NULL)
    }

    # Write standardized exposure to temp file
    fwrite(exp_dt, exp_temp, sep = "\t", quote = FALSE)
    # Also create gzipped copy for tools that expect .sumstats.gz
    if("compress" %in% names(formals(data.table::fwrite))){
      fwrite(exp_dt, exp_temp_gz, sep = "\t", quote = FALSE, compress = "gzip")
    } else if(requireNamespace("R.utils", quietly = TRUE)){
      R.utils::gzip(exp_temp, destname = exp_temp_gz, overwrite = TRUE)
    }

    # Read and standardize outcome file
    out_dt <- read_gwas(out_file)

    # Apply same column mapping for outcome
    for(old_col in names(col_map)){
      new_col <- col_map[[old_col]]
      if(old_col %in% names(out_dt)){
        setnames(out_dt, old_col, new_col)
      }
    }

    # Check for required MRlap columns in outcome
    missing_out <- setdiff(required_cols, names(out_dt))
    if(length(missing_out) > 0){
      error_msg <<- paste0("Outcome missing required columns for MRlap: ", paste(missing_out, collapse=", "),
                          ". Available columns: ", paste(names(out_dt), collapse=", "))
      message("      [MRlap SKIP] ", error_msg)
      report_log("ERROR", paste("MRlap preprocessing:", error_msg))
      report_error(
        error_type = "MRlap_missing_columns",
        details = error_msg,
        exposure = exposure_name,
        outcome = outcome_name,
        location = "run_mrlap_analysis"
      )
      return(NULL)
    }

    # Write standardized outcome to temp file
    fwrite(out_dt, out_temp, sep = "\t", quote = FALSE)
    # Also create gzipped copy for tools that expect .sumstats.gz
    if("compress" %in% names(formals(data.table::fwrite))){
      fwrite(out_dt, out_temp_gz, sep = "\t", quote = FALSE, compress = "gzip")
    } else if(requireNamespace("R.utils", quietly = TRUE)){
      R.utils::gzip(out_temp, destname = out_temp_gz, overwrite = TRUE)
    }

  }, error = function(e){
    error_msg <<- paste("File preprocessing failed:", e$message)
    message("      [MRlap PREPROCESSING ERROR] ", error_msg)
    report_log("ERROR", paste("MRlap preprocessing failed:", e$message))
    report_error(
      error_type = "MRlap_preprocessing_failed",
      details = error_msg,
      exposure = exposure_name,
      outcome = outcome_name,
      location = "run_mrlap_analysis"
    )
    return(NULL)
  })

  # Run MRlap with preprocessed files
  # Pass explicit file paths to the standardized .sumstats files we just wrote
  # Note: Suppress harmless package loading warnings (e.g., GenomicSEM namespace conflicts)

  # Use gz files if they exist, otherwise fall back to uncompressed
  exp_input <- if(file.exists(exp_temp_gz)) exp_temp_gz else exp_temp
  out_input <- if(file.exists(out_temp_gz)) out_temp_gz else out_temp

  # Change working directory to the temp dir so MRlap finds the files even if it uses relative paths
  orig_wd <- getwd()
  setwd(mrlap_temp_dir)
  on.exit(setwd(orig_wd), add = TRUE)

  # Debug: Verify files exist before calling MRlap
  report_log("DEBUG", paste("MRlap: Checking file existence"))
  report_log("DEBUG", paste("  Exposure temp file exists:", file.exists(exp_temp)))
  report_log("DEBUG", paste("  Outcome temp file exists:", file.exists(out_temp)))
  report_log("DEBUG", paste("  Exposure temp gz exists:", file.exists(exp_temp_gz)))
  report_log("DEBUG", paste("  Outcome temp gz exists:", file.exists(out_temp_gz)))
  report_log("DEBUG", paste("  Exposure input path:", exp_input))
  report_log("DEBUG", paste("  Outcome input path:", out_input))

  mrlap_res <- tryCatch({
    suppressWarnings({
      MRlap::MRlap(
        exposure = exp_input,           # Path to exposure .sumstats/.gz file
        outcome = out_input,            # Path to outcome .sumstats/.gz file
        exposure_name = exposure_name,  # Labels in the output
        outcome_name = outcome_name,    # Labels in the output
        ld            = "/home/Dementia_Depression_MR-analysis/data/MR_pipeline_demo/ldsc/eur_w_ld_chr",
        hm3           = "/home/Dementia_Depression_MR-analysis/data/MR_pipeline_demo/ldsc/w_hm3.snplist",
        save_logfiles = FALSE           # Don't save LDSC log files
      )
    })
  }, error = function(e){
    error_msg <<- e$message
    message("      [MRlap ERROR] ", e$message)
    message("      [MRlap ERROR] Exposure: ", basename(exp_file))
    message("      [MRlap ERROR] Outcome: ", basename(out_file))

    report_log("ERROR", paste("MRlap execution failed:", e$message, "| Exposure:", basename(exp_file), "| Outcome:", basename(out_file)))

    # Log error using universal error logging function
    report_error(
      error_type = "MRlap_execution_failed",
      details = e$message,
      exposure = exposure_name,
      outcome = outcome_name,
      location = "run_mrlap_analysis"
    )

    NULL
  })

  if(is.null(mrlap_res)){
    return(data.table(
      exposure = exposure_name,
      outcome = outcome_name,
      mrlap_status = "Failed",
      mrlap_error = error_msg,
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

  message("      [MRlap SUCCESS] Correction completed")
  report_log("INFO", paste("MRlap SUCCESS:", exposure_name, "->", outcome_name, "| IVs:", mrlap_res$MRcorrection$m_IVs))

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
#' @param outcome_file Optional file path used for caching reverse IVs
#' @return Reverse MR results (IVW only for speed)
run_reverse_mr <- function(outcome_dt, exposure_dt, outcome_label, exposure_label,
                           outcome_file = NULL, exposure_out_fmt = NULL){
  # Prefer rsID source for clumping if available (SNP_rsid or SNP)
  outcome_dt_sel <- copy(outcome_dt)
  if("SNP_rsid" %in% names(outcome_dt_sel) && any(grepl("^rs\\d+$", outcome_dt_sel$SNP_rsid))){
    outcome_dt_sel[, SNP := SNP_rsid]
  }

  # Select instruments from outcome (now treating as exposure)
  cache_key <- NA_character_
  disk_path <- NA_character_
  lock_path <- NA_character_
  ivs_reverse <- NULL

  if(!is.null(outcome_file) && nzchar(outcome_file)){
    cache_key <- make_file_cache_key(outcome_file, P_THRESH, CLUMP_R2, CLUMP_KB)
    if(!is.na(cache_key) && !is.null(REVERSE_IVS_CACHE[[cache_key]])){
      REVERSE_IVS_CACHE_HITS <<- REVERSE_IVS_CACHE_HITS + 1L
      ivs_reverse <- copy(REVERSE_IVS_CACHE[[cache_key]])
    } else {
      REVERSE_IVS_CACHE_MISSES <<- REVERSE_IVS_CACHE_MISSES + 1L
    }
  }

  # Try disk cache for reverse IVs (persistent across batches)
  if(is.null(ivs_reverse) &&
     isTRUE(REVERSE_IVS_DISK_CACHE_ENABLED) &&
     !is.null(outcome_file) && nzchar(outcome_file) &&
     requireNamespace("digest", quietly = TRUE)){
    out_norm <- normalizePath(outcome_file, winslash = "/", mustWork = FALSE)
    out_mtime <- as.numeric(file.info(out_norm)$mtime)
    disk_key <- digest::digest(
      list("reverse_ivs", REVERSE_IVS_DISK_CACHE_VERSION, out_norm, out_mtime, P_THRESH, CLUMP_R2, CLUMP_KB),
      algo = "xxhash64"
    )
    disk_path <- file.path(REVERSE_IVS_DISK_CACHE_DIR, paste0(disk_key, ".rds"))
    lock_path <- paste0(disk_path, ".lock")

    if(file.exists(disk_path)){
      ivs_disk <- tryCatch(readRDS(disk_path), error = function(e) NULL)
      if(!is.null(ivs_disk)){
        REVERSE_IVS_DISK_CACHE_HITS <<- REVERSE_IVS_DISK_CACHE_HITS + 1L
        ivs_reverse <- ivs_disk
        if(!is.na(cache_key)){
          REVERSE_IVS_CACHE[[cache_key]] <- ivs_reverse
        }
      }
    }

    if(is.null(ivs_reverse)){
      REVERSE_IVS_DISK_CACHE_MISSES <<- REVERSE_IVS_DISK_CACHE_MISSES + 1L
      # If another process is building the same cache, wait briefly
      if(file.exists(lock_path)){
        lock_age <- difftime(Sys.time(), file.info(lock_path)$mtime, units = "secs")
        if(lock_age < 60){
          for(wait in 1:10){
            Sys.sleep(0.5)
            if(file.exists(disk_path)){
              ivs_disk <- tryCatch(readRDS(disk_path), error = function(e) NULL)
              if(!is.null(ivs_disk)){
                REVERSE_IVS_DISK_CACHE_HITS <<- REVERSE_IVS_DISK_CACHE_HITS + 1L
                ivs_reverse <- ivs_disk
                if(!is.na(cache_key)){
                  REVERSE_IVS_CACHE[[cache_key]] <- ivs_reverse
                }
                break
              }
            }
          }
        } else {
          unlink(lock_path)
        }
      }
    }
  }

  # Compute if not cached
  if(is.null(ivs_reverse)){
    ivs_reverse <- select_instruments(outcome_dt_sel, P_THRESH, CLUMP_R2, CLUMP_KB)
    if(!is.na(cache_key)){
      REVERSE_IVS_CACHE[[cache_key]] <- ivs_reverse
    }

    # Persist to disk cache (best effort)
    if(!is.na(disk_path) && isTRUE(REVERSE_IVS_DISK_CACHE_ENABLED)){
      if(!file.exists(disk_path)){
        # Create lock file
        writeLines(as.character(Sys.getpid()), lock_path)
        tmp_path <- paste0(disk_path, ".", Sys.getpid(), ".tmp")
        tryCatch({
          saveRDS(ivs_reverse, tmp_path, compress = REVERSE_IVS_DISK_CACHE_COMPRESS)
          if(!file.exists(disk_path)){
            file.rename(tmp_path, disk_path)
          }
        }, error = function(e) NULL)
        if(file.exists(tmp_path)) unlink(tmp_path)
        if(file.exists(lock_path)) unlink(lock_path)
      }
    }
  }
  
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
  # Use original outcome_dt so chr:pos can be restored (if available)
  if("SNP_rsid" %in% names(outcome_dt) && any(grepl("^rs\\d+$", outcome_dt$SNP_rsid))){
    rev_exp_iv <- outcome_dt[SNP_rsid %in% ivs_reverse$rsid]
  } else {
    rev_exp_iv <- outcome_dt[SNP %in% ivs_reverse$rsid]
  }
  rev_exp_iv <- use_chrpos_id(rev_exp_iv)
  rev_exp_fmt <- to_exposure_format(rev_exp_iv, outcome_label)
  
  # Prepare reverse outcome data
  if(is.null(exposure_out_fmt)){
    exposure_dt_cp <- use_chrpos_id(exposure_dt)
    exposure_out_fmt <- to_outcome_format(exposure_dt_cp, exposure_label)
  }
  rev_out_fmt <- exposure_out_fmt
  
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
# PROGRAMMATIC ENTRYPOINT (for Python/parallel runs)
################################################################################
#
# Exposes a single-batch interface that mirrors the full pipeline but writes
# outputs into a caller-specified directory. Auto-run at the bottom of this
# file can be disabled by setting option mr_skip_autorun=TRUE or environment
# variable MR_SKIP_AUTORUN=true.
#
run_main_analysis <- function(
  exposure_path,
  mediator_path = NULL,
  outcome_path,
  covariate_paths = NULL,
  output_dir = file.path(DIR_RES_TRIAL, "batch_run"),
  ...
){
  report_progress("run_start", 0, 1)
  if(missing(exposure_path) || missing(outcome_path)){
    stop("run_main_analysis requires exposure_path and outcome_path")
  }
  
  normalize_files <- function(x){
    if(is.null(x) || length(x) == 0) return(character())
    x <- unlist(x, use.names = FALSE)
    x <- x[nzchar(x)]
    if(length(x) == 0) return(character())
    unique(normalizePath(x, winslash = "/", mustWork = FALSE))
  }
  
  exp_files   <- normalize_files(exposure_path)
  med_files   <- normalize_files(mediator_path)
  out_files   <- normalize_files(outcome_path)
  covar_files <- normalize_files(covariate_paths)
  
  missing_req <- c(exp_files[!file.exists(exp_files)], out_files[!file.exists(out_files)])
  if(length(missing_req)){
    stop("Required GWAS file(s) not found: ", paste(missing_req, collapse = ", "))
  }
  
  if(length(med_files)){
    missing_med <- med_files[!file.exists(med_files)]
    if(length(missing_med)){
      warning("Mediator files not found and will be skipped: ", paste(missing_med, collapse = ", "))
      med_files <- med_files[file.exists(med_files)]
    }
  }
  
  if(length(covar_files)){
    missing_cov <- covar_files[!file.exists(covar_files)]
    if(length(missing_cov)){
      warning("Covariate files not found and will be skipped: ", paste(missing_cov, collapse = ", "))
      covar_files <- covar_files[file.exists(covar_files)]
    }
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # Logs for batches live one level up from the batch directory
  log_base <- normalizePath(dirname(output_dir), winslash = "/", mustWork = FALSE)
  
  uvmr_out <- file.path(output_dir, "uvmr_results.csv")
  mvmr_out <- file.path(output_dir, "mvmr_results.csv")
  med_out  <- file.path(output_dir, "mediation_results.csv")

  # Preload outcomes once and reuse across UVMR/MVMR/Mediation
  preloaded_outs <- NULL
  if(length(out_files) > 0){
    preloaded_outs <- preload_outcomes(out_files, report_step_prefix = "")
  }

  # Preload exposures once (with full GWAS) if mediation is requested
  preloaded_expos <- NULL
  if(PRELOAD_EXPOSURES_ONCE && length(exp_files) > 0 && length(med_files) > 0){
    message("[INFO] Preloading exposures once for UVMR/MVMR/Mediation...")
    preloaded_expos <- preload_exposures(exp_files, keep_exp_dt = TRUE, report_step_prefix = "")
  }
  
  uvmr_results <- run_comprehensive_uvmr(
    expos_files = exp_files,
    outcomes_files = out_files,
    out_csv = uvmr_out,
    preloaded_outs = preloaded_outs,
    preloaded_expos = preloaded_expos
  )
  
  if(!is.null(uvmr_results) && nrow(uvmr_results) > 0){
    uvmr_results <- check_method_concordance(uvmr_results, p_threshold = 0.05)
    uvmr_results <- apply_fdr_correction(uvmr_results, pval_col = "pval")
    fwrite(uvmr_results, uvmr_out, row.names = FALSE)
    uvmr_out_abs <- normalizePath(uvmr_out, winslash = "/", mustWork = FALSE)
    append_task_log(file.path(log_base, "uvmr.log"),
                    paste("uvmr_results.csv", uvmr_out_abs, "rows", nrow(uvmr_results)))
  }
  
  mvmr_results <- NULL
  message("[DEBUG] Covariate files count: ", length(covar_files))
  report_log("DEBUG", paste("Covariate files count:", length(covar_files)))

  # 检查是否应该跳过 MVMR
  if(SKIP_MVMR){
    message("[INFO] MVMR analysis skipped (MR_SKIP_MVMR=TRUE)")
    report_log("INFO", "MVMR analysis skipped (MR_SKIP_MVMR=TRUE)")
  } else if(length(covar_files) > 0){
    message("[DEBUG] Covariate files: ", paste(basename(covar_files), collapse=", "))
    report_log("DEBUG", paste("Covariate files:", paste(basename(covar_files), collapse=", ")))
    mvmr_results <- run_comprehensive_mvmr(
      expos_files = exp_files,
      covar_files = covar_files,
      outcomes_files = out_files,
      out_csv = mvmr_out,
      preloaded_outs = preloaded_outs,
      preloaded_expos = preloaded_expos
    )

    if(!is.null(mvmr_results) && nrow(mvmr_results) > 0){
      mvmr_results <- apply_fdr_correction(mvmr_results, pval_col = "pval_mvmr")
      fwrite(mvmr_results, mvmr_out, row.names = FALSE)
      message("[INFO] MVMR results saved to: ", mvmr_out)
      report_log("INFO", paste("MVMR results saved to:", mvmr_out))
      mvmr_out_abs <- normalizePath(mvmr_out, winslash = "/", mustWork = FALSE)
      append_task_log(file.path(log_base, "mvmr.log"),
                      paste("mvmr_results.csv", mvmr_out_abs, "rows", nrow(mvmr_results)))
    } else {
      message("[WARNING] MVMR returned NULL or empty results")
      report_log("WARNING", "MVMR returned NULL or empty results")
    }
  } else {
    message("[INFO] No covariate files supplied; skipping MVMR.")
    report_log("INFO", "No covariate files supplied; skipping MVMR.")
    message("[DEBUG] covariate_paths parameter: ", if(is.null(covariate_paths)) "NULL" else paste(length(covariate_paths), "paths provided"))
    report_log("DEBUG", paste("covariate_paths parameter:", if(is.null(covariate_paths)) "NULL" else paste(length(covariate_paths), "paths provided")))
  }
  
  mediation_results <- NULL
  if(length(med_files) > 0){
    mediation_results <- run_comprehensive_mediation(
      expos_files = exp_files,
      medi_files = med_files,
      outcomes_files = out_files,
      out_csv = med_out,
      preloaded_outs = preloaded_outs,
      preloaded_expos = preloaded_expos
    )
    
    if(!is.null(mediation_results) && nrow(mediation_results) > 0){
      mediation_results <- apply_fdr_correction(mediation_results, pval_col = "pval_exp_med")
      fwrite(mediation_results, med_out, row.names = FALSE)
    }
  } else {
    message("[INFO] No mediator files supplied; skipping mediation.")
    report_log("INFO", "No mediator files supplied; skipping mediation.")
  }
  
  report_progress("run_end", 1, 1)
  # Clean temp folders (e.g., Rtmp* created under output_dir/tmp)
  cleanup_temp_dirs(output_dir)
  
  invisible(list(
    uvmr = uvmr_results,
    mvmr = mvmr_results,
    mediation = mediation_results
  ))
}

################################################################################
# MAIN EXECUTION - Complete MR Analysis Pipeline
################################################################################

if(!SKIP_AUTORUN){
  message("\n", rep("=", 80))
  message("COMPREHENSIVE MENDELIAN RANDOMIZATION ANALYSIS PIPELINE")
  message(rep("=", 80), "\n")

  # Optional: preload exposures once for UVMR + Mediation
  preloaded_expos_main <- NULL
  if(PRELOAD_EXPOSURES_ONCE && length(all_medi_files) > 0){
    message("[INFO] Preloading exposures once for UVMR/MVMR/Mediation...")
    preloaded_expos_main <- preload_exposures(all_expos_files, keep_exp_dt = TRUE, report_step_prefix = "")
  }
  
  # ---- Step 1: Univariable MR (UVMR) ----
  message(">>> STEP 1: Univariable Mendelian Randomization (UVMR)")
  message("    Analyzing individual exposure-outcome associations")
  uvmr_results <- run_comprehensive_uvmr(preloaded_expos = preloaded_expos_main)
  
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
  mvmr_results <- run_comprehensive_mvmr(preloaded_expos = preloaded_expos_main)
  
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
  mediation_results <- run_comprehensive_mediation(preloaded_expos = preloaded_expos_main)
  
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

  # ---- PERFORMANCE OPTIMIZATION: Report cache statistics ----
  message("\n[PERFORMANCE OPTIMIZATION REPORT]")
  report_cache_stats()
  message(rep("=", 80))

  message("\nNext steps:")
  message("  1. Review results in results_trial/ folder")
  message("  2. Check q_value column for FDR-corrected significance")
  message("  3. For UVMR: Use only concordance_validated=TRUE results")
  message("  4. Filter significant results using Results_Filter_Helper.R")
  message("  5. Visualize results with your preferred tools")
  message(rep("=", 80), "\n")
} else {
  message("[INFO] MR_SKIP_AUTORUN set - skipping automatic pipeline execution.")
}

