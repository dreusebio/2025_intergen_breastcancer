#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 08: Module and Sample Diagnostics
#
# Pipeline: comethyl WGBS network analysis
#
# PURPOSE
#   - Load one or more module objects from Script 07
#   - Optionally align module eigengenes to sample metadata
#   - Generate:
#       1) module-module structure diagnostics from module eigengenes
#       2) sample-sample structure diagnostics based on module eigengenes
#       3) sample-by-module eigengene heatmap
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --modules_v1   : path to v1 Modules.rds from Script 07
#
# OPTIONAL INPUTS
#   --modules_v2            : optional path to v2 Modules.rds
#   --modules_v3            : optional path to v3 Modules.rds
#   --sample_info            : optional metadata xlsx with sample IDs as rownames
#   --module_cor            : correlation for module structure:
#                             bicor or pearson [default = pearson]
#   --sample_dendro_distance: distance for sample dendrogram:
#                             euclidean, pearson, or bicor [default = euclidean]
#   --max_p_outliers        : maxPOutliers for bicor calls [default = 0.1]
#   --enforce_min_module_size : TRUE/FALSE. If TRUE, Script 08 reads the
#                              min_module_size recorded by Script 07 and
#                              reassigns smaller non-grey modules to grey
#                              [default = TRUE]
#   --min_module_size_check  : optional manual override; otherwise auto-read
#                              from Script 07 run_parameters.txt
#
# OUTPUTS
#   project_root/comethyl_output/08_module_and_sample_diagnostics/<cpg_label>/<region_label>/<variant>/
#       Module_ME_Dendrogram.pdf
#       Module_Correlation_Heatmap.pdf
#       Module_Correlation_Stats.tsv
#       Sample_ME_Dendrogram.pdf
#       Sample_Correlation_Heatmap.pdf
#       Sample_ME_Heatmap.pdf
#       Modules_minSizeQC.rds
#       Modules_for_downstream.rds
#       Module_Size_QC_before.tsv
#       Module_Size_QC_after.tsv
#       Module_Size_QC_summary.txt
#       ME_Region_Consistency_QC.txt
#       run_parameters.txt
#
# NOTES
#   - The post-detection minimum module-size threshold is read from Script 07
#     run_parameters.txt. If Script 07 used min_module_size: 10, Script 08 uses 10.
#   - Small non-grey modules are reassigned to grey before downstream use.
#   - Orphan ME columns are removed. This handles both ME-prefixed columns
#     such as MEdarkmagenta and unprefixed columns such as darkmagenta.
#   - Grey module columns are removed before plotting, whether named grey or MEgrey.
#   - All plot calls are wrapped in tryCatch() so one failure does not halt the script.
#   - Sample-level plots use transpose = FALSE when computing dendrograms and
#     correlations from the ME matrix because MEs are stored as samples x modules.
# ================================================================
message("Starting Script 08")

suppressPackageStartupMessages({
  library(optparse)
  library(comethyl)
  library(WGCNA)
  library(AnnotationHub)
  library(openxlsx)
  library(dplyr)
})

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
write_log_lines <- function(lines, file) {
  writeLines(as.character(lines), con = file)
}

safe_close_pdf <- function() {
  cur_dev <- grDevices::dev.cur()
  if (!is.null(cur_dev) && cur_dev > 1) {
    try(grDevices::dev.off(), silent = TRUE)
  }
}

validate_modules_object <- function(x, label) {
  if (!is.list(x)) stop(label, " must be a list-like module object.")
  if (is.null(x$MEs)) stop(label, " is missing $MEs.")
  if (!(is.matrix(x$MEs) || is.data.frame(x$MEs))) stop(label, "$MEs must be matrix-like.")

  MEs <- as.matrix(x$MEs)
  if (!is.numeric(MEs)) stop(label, "$MEs must be numeric.")
  if (nrow(MEs) < 2) stop(label, "$MEs must have at least 2 samples.")
  if (ncol(MEs) < 1) stop(label, "$MEs must have at least 1 module.")
  if (is.null(rownames(MEs))) stop(label, "$MEs must have sample IDs in rownames.")
  if (is.null(colnames(MEs))) stop(label, "$MEs must have module names in colnames.")

  if (anyDuplicated(rownames(MEs))) {
    dup_ids <- unique(rownames(MEs)[duplicated(rownames(MEs))])
    stop(label, "$MEs has duplicated sample IDs. Example: ",
         paste(head(dup_ids, 10), collapse = ", "))
  }

  if (anyDuplicated(colnames(MEs))) {
    dup_mods <- unique(colnames(MEs)[duplicated(colnames(MEs))])
    stop(label, "$MEs has duplicated module names. Example: ",
         paste(head(dup_mods, 10), collapse = ", "))
  }

  x$MEs <- MEs
  x
}

parse_run_parameters <- function(file) {
  if (is.null(file) || is.na(file) || !file.exists(file)) return(list())
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  out <- list()
  for (ln in lines) {
    # Supports both "key: value" and "key\tvalue" formats.
    if (grepl("\t", ln)) {
      parts <- strsplit(ln, "\t", fixed = FALSE)[[1]]
      key <- trimws(parts[1])
      val <- trimws(paste(parts[-1], collapse = "\t"))
    } else if (grepl(":", ln, fixed = TRUE)) {
      key <- trimws(sub(":.*$", "", ln))
      val <- trimws(sub("^[^:]*:", "", ln))
    } else {
      next
    }
    if (nzchar(key)) out[[key]] <- val
  }
  out
}

infer_script07_run_parameters_file <- function(modules_file, explicit_file = NULL) {
  if (!is.null(explicit_file) && nzchar(explicit_file)) {
    if (!file.exists(explicit_file)) stop("script07_run_parameters_file not found: ", explicit_file)
    return(explicit_file)
  }

  variant_dir <- dirname(modules_file)
  step_dir <- dirname(variant_dir)

  candidates <- c(
    file.path(variant_dir, "run_parameters.txt"),
    file.path(step_dir, "run_parameters.txt")
  )
  candidates[file.exists(candidates)][1]
}

get_script07_min_module_size <- function(modules_file,
                                         explicit_params_file = NULL,
                                         fallback = NA_integer_) {
  params_file <- infer_script07_run_parameters_file(modules_file, explicit_params_file)

  if (is.na(params_file) || is.null(params_file) || !nzchar(params_file)) {
    return(list(min_module_size = fallback, params_file = NA_character_))
  }

  params <- parse_run_parameters(params_file)
  val <- suppressWarnings(as.integer(params$min_module_size))

  if (is.na(val)) {
    return(list(min_module_size = fallback, params_file = params_file))
  }

  list(min_module_size = val, params_file = params_file)
}

get_module_column <- function(modules) {
  if (is.null(modules$regions) || !is.data.frame(modules$regions)) {
    stop("Modules object is missing $regions as a data.frame; cannot enforce module-size QC.")
  }

  candidate_cols <- c("module", "Module", "color", "colour", "module_color", "moduleColor")
  hit <- candidate_cols[candidate_cols %in% colnames(modules$regions)][1]

  if (is.na(hit)) {
    stop("modules$regions does not contain a recognizable module/color column. Expected one of: ",
         paste(candidate_cols, collapse = ", "))
  }

  hit
}

write_module_size_table <- function(counts, file) {
  df <- data.frame(
    module = names(counts),
    n_regions = as.integer(counts),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$n_regions, df$module), , drop = FALSE]
  write.table(df, file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(df)
}

me_col_to_module <- function(me_cols) {
  # Handles both common formats:
  #   WGCNA-style: MEblue, MEgrey
  #   comethyl-style in this pipeline: blue, grey
  sub("^ME", "", as.character(me_cols))
}

synchronize_me_columns_with_regions <- function(modules,
                                                variant_name,
                                                variant_out_dir,
                                                modules_file = NA_character_) {
  module_col <- get_module_column(modules)

  valid_modules <- unique(as.character(modules$regions[[module_col]]))
  valid_modules <- valid_modules[!is.na(valid_modules) & valid_modules != ""]

  me_cols <- colnames(modules$MEs)
  me_modules <- me_col_to_module(me_cols)

  keep_cols <- me_modules %in% valid_modules
  orphan_me_cols <- me_cols[!keep_cols]
  orphan_me_modules <- me_modules[!keep_cols]

  missing_me_modules <- setdiff(valid_modules, me_modules[keep_cols])

  if (length(orphan_me_cols) > 0) {
    message("[", variant_name, "] Dropping orphan ME columns not present in modules$regions$", module_col, ": ",
            paste(orphan_me_cols, collapse = ", "))
    modules$MEs <- as.matrix(modules$MEs[, keep_cols, drop = FALSE])
  } else {
    message("[", variant_name, "] No orphan ME columns found; MEs and region modules are consistent.")
  }

  writeLines(
    c(
      paste("variant_name:", variant_name),
      paste("modules_file_input:", modules_file),
      paste("module_column_used:", module_col),
      paste("n_region_modules:", length(valid_modules)),
      paste("n_ME_columns_before_sync:", length(me_cols)),
      paste("n_ME_columns_after_sync:", ncol(modules$MEs)),
      paste("orphan_ME_columns_dropped:",
            ifelse(length(orphan_me_cols) == 0, "NONE", paste(orphan_me_cols, collapse = ", "))),
      paste("orphan_ME_modules_dropped:",
            ifelse(length(orphan_me_modules) == 0, "NONE", paste(orphan_me_modules, collapse = ", "))),
      paste("region_modules_without_ME_column:",
            ifelse(length(missing_me_modules) == 0, "NONE", paste(missing_me_modules, collapse = ", "))),
      paste("date:", as.character(Sys.time()))
    ),
    con = file.path(variant_out_dir, "ME_Region_Consistency_QC.txt")
  )

  list(
    modules = modules,
    module_col = module_col,
    orphan_me_cols = orphan_me_cols,
    orphan_me_modules = orphan_me_modules,
    missing_me_modules = missing_me_modules
  )
}

enforce_modules_min_size <- function(modules,
                                     min_module_size,
                                     variant_name,
                                     variant_out_dir,
                                     modules_file,
                                     overwrite_modules_rds = FALSE) {
  if (is.na(min_module_size) || !is.finite(min_module_size)) {
    stop("Cannot enforce minimum module size because min_module_size is NA. Check Script 07 run_parameters.txt.")
  }
  min_module_size <- as.integer(min_module_size)
  if (min_module_size < 2) stop("min_module_size must be >= 2 for module-size QC.")

  module_col <- get_module_column(modules)
  module_vec <- as.character(modules$regions[[module_col]])
  module_vec[is.na(module_vec) | module_vec == ""] <- "grey"

  counts_before <- sort(table(module_vec), decreasing = TRUE)
  write_module_size_table(counts_before,
                          file.path(variant_out_dir, "Module_Size_QC_before.tsv"))

  non_grey <- names(counts_before)[!tolower(names(counts_before)) %in% "grey"]
  small_modules <- non_grey[as.integer(counts_before[non_grey]) < min_module_size]

  if (length(small_modules) > 0) {
    message("[", variant_name, "] Reassigning ", length(small_modules),
            " module(s) smaller than Script 07 min_module_size=", min_module_size,
            " to grey: ", paste(small_modules, collapse = ", "))

    module_vec[module_vec %in% small_modules] <- "grey"
    modules$regions[[module_col]] <- module_vec

    # If a separate color vector exists and has the same length as regions, keep it consistent.
    if (!is.null(modules$colors) && length(modules$colors) == length(module_vec)) {
      colors_vec <- as.character(modules$colors)
      colors_vec[colors_vec %in% small_modules] <- "grey"
      modules$colors <- colors_vec
    }

    # If another common region-level color column exists, keep it consistent too.
    extra_cols <- intersect(c("color", "colour", "module_color", "moduleColor"), colnames(modules$regions))
    for (ec in extra_cols) {
      ec_vec <- as.character(modules$regions[[ec]])
      ec_vec[ec_vec %in% small_modules] <- "grey"
      modules$regions[[ec]] <- ec_vec
    }
  } else {
    message("[", variant_name, "] No non-grey modules below Script 07 min_module_size=", min_module_size)
  }

  counts_after <- sort(table(as.character(modules$regions[[module_col]])), decreasing = TRUE)
  write_module_size_table(counts_after,
                          file.path(variant_out_dir, "Module_Size_QC_after.tsv"))

  # Final consistency guard: remove any ME columns whose module no longer exists
  # in modules$regions after small modules have been reassigned. This catches both
  # ME-prefixed columns (e.g., MEdarkmagenta) and unprefixed columns (e.g., darkmagenta).
  sync_res <- synchronize_me_columns_with_regions(
    modules = modules,
    variant_name = variant_name,
    variant_out_dir = variant_out_dir,
    modules_file = modules_file
  )
  modules <- sync_res$modules

  writeLines(
    c(
      paste("script07_min_module_size:", min_module_size),
      paste("module_column_used:", module_col),
      paste("small_modules_reassigned_to_grey:",
            ifelse(length(small_modules) == 0, "NONE", paste(small_modules, collapse = ", "))),
      paste("n_small_modules_reassigned:", length(small_modules)),
      paste("orphan_ME_columns_dropped:",
            ifelse(length(sync_res$orphan_me_cols) == 0, "NONE", paste(sync_res$orphan_me_cols, collapse = ", "))),
      paste("n_orphan_ME_columns_dropped:", length(sync_res$orphan_me_cols)),
      paste("modules_file_input:", modules_file),
      paste("overwrite_modules_rds:", overwrite_modules_rds),
      paste("date:", as.character(Sys.time()))
    ),
    con = file.path(variant_out_dir, "Module_Size_QC_summary.txt")
  )

  corrected_rds <- file.path(variant_out_dir, "Modules_minSizeQC.rds")
  downstream_rds <- file.path(variant_out_dir, "Modules_for_downstream.rds")

  saveRDS(modules, corrected_rds)
  saveRDS(modules, downstream_rds)

  if (isTRUE(overwrite_modules_rds)) {
    backup_file <- paste0(modules_file, ".pre_minSizeQC_backup")
    if (!file.exists(backup_file)) file.copy(modules_file, backup_file)
    saveRDS(modules, modules_file)
  }

  list(
    modules = modules,
    min_module_size = min_module_size,
    module_col = module_col,
    small_modules = small_modules,
    corrected_rds = corrected_rds,
    downstream_rds = downstream_rds
  )
}

# ------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--modules_v1", type = "character",
              help = "Path to v1 Modules.rds from Script 07"),

  make_option("--modules_v2", type = "character", default = NULL,
              help = "Optional path to v2 Modules.rds from Script 07"),

  make_option("--modules_v3", type = "character", default = NULL,
              help = "Optional path to v3 Modules.rds from Script 07"),

  make_option("--sample_info", type = "character", default = NULL,
              help = "Optional metadata xlsx with sample IDs as rownames for sample alignment/logging"),

  make_option("--module_cor", type = "character", default = "pearson",
              help = "Correlation for module structure: bicor or pearson [default = pearson]"),

  make_option("--sample_dendro_distance", type = "character", default = "euclidean",
              help = "Distance for sample dendrogram: euclidean, pearson, or bicor [default = euclidean]"),

  make_option("--max_p_outliers", type = "double", default = 0.1,
              help = "maxPOutliers for bicor calls [default = 0.1]"),

  make_option("--enforce_min_module_size", type = "logical", default = TRUE,
              help = "After loading Script 07 Modules.rds, reassign non-grey modules smaller than the min_module_size recorded in Script 07 run_parameters.txt to grey [default = TRUE]"),

  make_option("--min_module_size_check", type = "integer", default = NA,
              help = "Optional override for post-detection minimum module-size QC. If not provided, Script 08 reads min_module_size from Script 07 run_parameters.txt [default = auto]"),

  make_option("--script07_run_parameters_file", type = "character", default = NULL,
              help = "Optional explicit path to Script 07 run_parameters.txt. Usually not needed because Script 08 searches next to Modules.rds and one folder above."),

  make_option("--overwrite_modules_rds", type = "logical", default = FALSE,
              help = "If TRUE, overwrite the input Modules.rds after saving a .pre_minSizeQC_backup. Recommended: FALSE [default = FALSE]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$modules_v1))   stop("--modules_v1 is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$modules_v1))  stop("modules_v1 not found: ", opt$modules_v1)

if (!is.null(opt$modules_v2) && !file.exists(opt$modules_v2)) {
  stop("modules_v2 not found: ", opt$modules_v2)
}
if (!is.null(opt$modules_v3) && !file.exists(opt$modules_v3)) {
  stop("modules_v3 not found: ", opt$modules_v3)
}
if (!is.null(opt$sample_info) && !file.exists(opt$sample_info)) {
  stop("sample_info not found: ", opt$sample_info)
}

module_cor <- tolower(opt$module_cor)
if (!module_cor %in% c("pearson", "bicor")) {
  stop("--module_cor must be 'pearson' or 'bicor'")
}

sample_dendro_distance <- tolower(opt$sample_dendro_distance)
if (!sample_dendro_distance %in% c("euclidean", "pearson", "bicor")) {
  stop("--sample_dendro_distance must be one of: euclidean, pearson, bicor")
}

if (!is.numeric(opt$max_p_outliers) || opt$max_p_outliers < 0 || opt$max_p_outliers > 1) {
  stop("--max_p_outliers must be between 0 and 1")
}

if (!is.na(opt$min_module_size_check) && opt$min_module_size_check < 2) {
  stop("--min_module_size_check must be >= 2 when provided")
}

if (!is.null(opt$script07_run_parameters_file) && !file.exists(opt$script07_run_parameters_file)) {
  stop("script07_run_parameters_file not found: ", opt$script07_run_parameters_file)
}

# ------------------------------------------------------------
# Configure cache and threads
# ------------------------------------------------------------
AnnotationHub::setAnnotationHubOption(
  "CACHE",
  value = file.path(opt$project_root, ".cache")
)

WGCNA::enableWGCNAThreads()

# ------------------------------------------------------------
# Optional sample metadata
# ------------------------------------------------------------
colData_num <- NULL

if (!is.null(opt$sample_info)) {
  colData <- openxlsx::read.xlsx(opt$sample_info, rowNames = TRUE)
  message("Loaded sample_info: ", opt$sample_info)
  message("colData dimensions: ", nrow(colData), " samples x ", ncol(colData), " traits")

  colData_num <- colData %>%
    dplyr::mutate(dplyr::across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
    dplyr::select(where(~ is.numeric(.) && is.atomic(.)))

  if (nrow(colData_num) == 0) {
    stop("sample_info was loaded but has zero rows after processing.")
  }

  message("colData_num dimensions: ", nrow(colData_num), " samples x ", ncol(colData_num), " numeric traits")
}

# ------------------------------------------------------------
# Derive lineage from modules_v1
# Expected input:
#   .../07_module_detection/<cpg_label>/<region_label>/v1_all_pcs/Modules.rds
# ------------------------------------------------------------
v1_variant_dir <- dirname(opt$modules_v1)
v1_region_dir  <- dirname(v1_variant_dir)
region_label   <- basename(v1_region_dir)
cpg_label      <- basename(dirname(v1_region_dir))

pipeline_root <- file.path(opt$project_root, "comethyl_output")
step_dir      <- file.path(pipeline_root, "08_module_and_sample_diagnostics")
out_dir       <- file.path(step_dir, cpg_label, region_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Build variant input list
# ------------------------------------------------------------
variant_inputs <- list(v1_all_pcs = opt$modules_v1)

if (!is.null(opt$modules_v2)) {
  variant_inputs[[basename(dirname(opt$modules_v2))]] <- opt$modules_v2
}
if (!is.null(opt$modules_v3)) {
  variant_inputs[[basename(dirname(opt$modules_v3))]] <- opt$modules_v3
}

# ------------------------------------------------------------
# Run per variant
# ------------------------------------------------------------
for (variant_name in names(variant_inputs)) {
  message("\n==============================")
  message("Running module/sample diagnostics for variant: ", variant_name)
  message("==============================\n")

  modules_file <- variant_inputs[[variant_name]]

  modules <- validate_modules_object(
    readRDS(modules_file),
    paste0(variant_name, " modules object")
  )

  variant_out_dir <- file.path(out_dir, variant_name)
  dir.create(variant_out_dir, recursive = TRUE, showWarnings = FALSE)

  script07_min_info <- get_script07_min_module_size(
    modules_file,
    explicit_params_file = opt$script07_run_parameters_file,
    fallback = opt$min_module_size_check
  )

  min_module_size_used_for_qc <- script07_min_info$min_module_size
  script07_params_file_used <- script07_min_info$params_file

  if (isTRUE(opt$enforce_min_module_size)) {
    if (is.na(min_module_size_used_for_qc)) {
      stop("[", variant_name, "] Could not determine min_module_size for QC. ",
           "Provide --min_module_size_check or ensure Script 07 run_parameters.txt exists next to Modules.rds.")
    }

    qc_res <- enforce_modules_min_size(
      modules = modules,
      min_module_size = min_module_size_used_for_qc,
      variant_name = variant_name,
      variant_out_dir = variant_out_dir,
      modules_file = modules_file,
      overwrite_modules_rds = opt$overwrite_modules_rds
    )

    modules <- validate_modules_object(qc_res$modules, paste0(variant_name, " modules object after min-size QC"))
  } else {
    message("[", variant_name, "] Module-size QC disabled by --enforce_min_module_size false")

    # Even when minimum-size reassignment is disabled, still enforce internal
    # consistency so downstream scripts do not analyze ME columns that are not
    # represented in modules$regions.
    sync_res <- synchronize_me_columns_with_regions(
      modules = modules,
      variant_name = variant_name,
      variant_out_dir = variant_out_dir,
      modules_file = modules_file
    )
    modules <- validate_modules_object(sync_res$modules, paste0(variant_name, " modules object after ME-region sync"))
    saveRDS(modules, file.path(variant_out_dir, "Modules_for_downstream.rds"))
  }

  MEs <- as.data.frame(modules$MEs)

  # ----------------------------------------------------------
  # Remove grey module before all plotting
  # ----------------------------------------------------------
  grey_col <- colnames(MEs)[tolower(me_col_to_module(colnames(MEs))) == "grey"]
  if (length(grey_col) > 0) {
    MEs <- MEs[, !colnames(MEs) %in% grey_col, drop = FALSE]
    message("[", variant_name, "] Grey ME column removed: ", paste(grey_col, collapse = ", "),
            " (", ncol(MEs), " real MEs remaining)")
  } else {
    message("[", variant_name, "] No grey ME column found (", ncol(MEs), " real MEs remaining)")
  }

  if (ncol(MEs) < 2) {
    message("[", variant_name, "] Only ", ncol(MEs), " real ME(s) after grey removal — skipping all plots")
    writeLines(
      paste("Skipped: only", ncol(MEs), "real ME(s) after removing grey."),
      con = file.path(variant_out_dir, paste0(variant_name, "_plots_SKIPPED.txt"))
    )
    next
  }

  # ----------------------------------------------------------
  # Align samples to metadata if provided
  # ----------------------------------------------------------
  if (!is.null(colData_num)) {
    common_samples <- intersect(rownames(MEs), rownames(colData_num))

    if (length(common_samples) == 0) {
      stop("[", variant_name, "] No overlapping samples between module eigengenes and sample_info metadata.")
    }

    MEs_use <- MEs[common_samples, , drop = FALSE]
    colData_use <- colData_num[common_samples, , drop = FALSE]

    message("[", variant_name, "] Sample alignment to metadata: ",
            length(common_samples), " overlapping samples retained from ",
            nrow(MEs), " ME samples and ", nrow(colData_num), " metadata samples")
  } else {
    MEs_use <- MEs
    colData_use <- NULL
    common_samples <- rownames(MEs_use)
    message("[", variant_name, "] No sample_info provided; using all ", nrow(MEs_use), " ME samples")
  }

  if (nrow(MEs_use) < 2) {
    message("[", variant_name, "] Fewer than 2 samples after alignment — skipping plots")
    writeLines(
      "Skipped: fewer than 2 samples after metadata alignment.",
      con = file.path(variant_out_dir, paste0(variant_name, "_plots_SKIPPED.txt"))
    )
    next
  }

  message("[", variant_name, "] MEs used: ", nrow(MEs_use), " samples x ", ncol(MEs_use), " modules")

  f_module_me_dendro <- file.path(variant_out_dir, "Module_ME_Dendrogram.pdf")
  f_module_cor_hm    <- file.path(variant_out_dir, "Module_Correlation_Heatmap.pdf")
  f_module_stats     <- file.path(variant_out_dir, "Module_Correlation_Stats.tsv")
  f_sample_me_dendro <- file.path(variant_out_dir, "Sample_ME_Dendrogram.pdf")
  f_sample_cor_hm    <- file.path(variant_out_dir, "Sample_Correlation_Heatmap.pdf")
  f_sample_me_hm     <- file.path(variant_out_dir, "Sample_ME_Heatmap.pdf")

  # ==========================================================
  # MODULE-LEVEL diagnostics
  # ==========================================================
  moduleDendro <- tryCatch({
    getDendro(MEs_use, distance = module_cor)
  }, error = function(e) {
    message("[", variant_name, "] Module dendrogram failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(moduleDendro)) {
    tryCatch({
      plotDendro(moduleDendro, labelSize = 4, nBreaks = 5, file = f_module_me_dendro)
      message("[", variant_name, "] Saved module ME dendrogram")
    }, error = function(e) {
      safe_close_pdf()
      message("[", variant_name, "] Failed saving module dendrogram: ", conditionMessage(e))
    })
  }

  moduleCor <- tryCatch({
    getCor(MEs_use, corType = module_cor, maxPOutliers = opt$max_p_outliers)
  }, error = function(e) {
    message("[", variant_name, "] Module correlation failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(moduleCor) && !is.null(moduleDendro)) {
    tryCatch({
      plotHeatmap(
        moduleCor,
        rowDendro = moduleDendro,
        colDendro = moduleDendro,
        file      = f_module_cor_hm
      )
      message("[", variant_name, "] Saved module correlation heatmap")
    }, error = function(e) {
      safe_close_pdf()
      message("[", variant_name, "] Failed module correlation heatmap: ", conditionMessage(e))
    })
  }

  tryCatch({
    getMEtraitCor(
      MEs_use,
      colData  = MEs_use,
      corType  = module_cor,
      robustY  = TRUE,
      file     = f_module_stats
    )
    message("[", variant_name, "] Saved module correlation stats")
  }, error = function(e) {
    message("[", variant_name, "] Failed module correlation stats: ", conditionMessage(e))
    writeLines(
      paste("Failed:", conditionMessage(e)),
      con = file.path(variant_out_dir, "Module_Correlation_Stats_FAILED.txt")
    )
  })

  # ==========================================================
  # SAMPLE-LEVEL diagnostics
  # IMPORTANT: transpose = FALSE because MEs are samples x modules
  # ==========================================================
  sampleDendro <- tryCatch({
    getDendro(MEs_use, transpose = FALSE, distance = sample_dendro_distance)
  }, error = function(e) {
    message("[", variant_name, "] Sample dendrogram failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(sampleDendro)) {
    tryCatch({
      plotDendro(sampleDendro, labelSize = 3, nBreaks = 5, file = f_sample_me_dendro)
      message("[", variant_name, "] Saved sample ME dendrogram")
    }, error = function(e) {
      safe_close_pdf()
      message("[", variant_name, "] Failed saving sample dendrogram: ", conditionMessage(e))
    })
  }

  sampleCor <- tryCatch({
    getCor(MEs_use, transpose = TRUE, corType = module_cor, maxPOutliers = opt$max_p_outliers)
  }, error = function(e) {
    message("[", variant_name, "] Sample correlation failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(sampleCor) && !is.null(sampleDendro)) {
    tryCatch({
      plotHeatmap(
        sampleCor,
        rowDendro = sampleDendro,
        colDendro = sampleDendro,
        file      = f_sample_cor_hm
      )
      message("[", variant_name, "] Saved sample correlation heatmap")
    }, error = function(e) {
      safe_close_pdf()
      message("[", variant_name, "] Failed sample correlation heatmap: ", conditionMessage(e))
    })
  }

  tryCatch({
    plotHeatmap(
      MEs_use,
      rowDendro       = sampleDendro,
      colDendro       = moduleDendro,
      legend.title    = "Module\nEigengene",
      legend.position = c(0.37, 0.89),
      file            = f_sample_me_hm
    )
    message("[", variant_name, "] Saved sample x module eigengene heatmap")
  }, error = function(e) {
    safe_close_pdf()
    message("[", variant_name, "] Sample x module heatmap failed: ", conditionMessage(e))
    writeLines(
      paste("Failed:", conditionMessage(e)),
      con = file.path(variant_out_dir, "Sample_ME_Heatmap_FAILED.txt")
    )
  })

  # ----------------------------------------------------------
  # Per-variant run parameters
  # ----------------------------------------------------------
  write_log_lines(
    c(
      paste("variant_name:", variant_name),
      paste("modules_file:", variant_inputs[[variant_name]]),
      paste("modules_for_downstream:", file.path(variant_out_dir, "Modules_for_downstream.rds")),
      paste("sample_info:", ifelse(is.null(opt$sample_info), "NULL", opt$sample_info)),
      paste("module_cor:", module_cor),
      paste("sample_dendro_distance:", sample_dendro_distance),
      paste("max_p_outliers:", opt$max_p_outliers),
      paste("enforce_min_module_size:", opt$enforce_min_module_size),
      paste("script07_run_parameters_file_used:", ifelse(is.na(script07_params_file_used), "NA", script07_params_file_used)),
      paste("script07_min_module_size_used_for_qc:", ifelse(is.na(min_module_size_used_for_qc), "NA", min_module_size_used_for_qc)),
      paste("overwrite_modules_rds:", opt$overwrite_modules_rds),
      paste("n_samples_before_metadata_alignment:", nrow(MEs)),
      paste("n_samples_after_metadata_alignment:", nrow(MEs_use)),
      paste("n_traits_numeric:", ifelse(is.null(colData_use), 0, ncol(colData_use))),
      paste("n_modules:", ncol(MEs_use)),
      paste("module_me_dendrogram:", f_module_me_dendro),
      paste("module_correlation_heatmap:", f_module_cor_hm),
      paste("module_correlation_stats:", f_module_stats),
      paste("sample_me_dendrogram:", f_sample_me_dendro),
      paste("sample_correlation_heatmap:", f_sample_cor_hm),
      paste("sample_me_heatmap:", f_sample_me_hm),
      paste("date:", as.character(Sys.time()))
    ),
    file.path(variant_out_dir, "run_parameters.txt")
  )

  message("Finished variant: ", variant_name)
  message("  Outputs in: ", variant_out_dir)
}

# ------------------------------------------------------------
# Top-level run parameters
# ------------------------------------------------------------
write_log_lines(
  c(
    paste("project_root:", opt$project_root),
    paste("modules_v1:", opt$modules_v1),
    paste("modules_v2:", ifelse(is.null(opt$modules_v2), "NULL", opt$modules_v2)),
    paste("modules_v3:", ifelse(is.null(opt$modules_v3), "NULL", opt$modules_v3)),
    paste("sample_info:", ifelse(is.null(opt$sample_info), "NULL", opt$sample_info)),
    paste("module_cor:", module_cor),
    paste("sample_dendro_distance:", sample_dendro_distance),
    paste("max_p_outliers:", opt$max_p_outliers),
    paste("enforce_min_module_size:", opt$enforce_min_module_size),
    paste("min_module_size_check_override:", ifelse(is.na(opt$min_module_size_check), "auto_from_script07", opt$min_module_size_check)),
    paste("script07_run_parameters_file:", ifelse(is.null(opt$script07_run_parameters_file), "auto", opt$script07_run_parameters_file)),
    paste("overwrite_modules_rds:", opt$overwrite_modules_rds),
    paste("cpg_label:", cpg_label),
    paste("region_label:", region_label),
    paste("variants_run:", paste(names(variant_inputs), collapse = ", ")),
    paste("date:", as.character(Sys.time()))
  ),
  file.path(out_dir, "run_parameters.txt")
)

message("Script 08 complete: module and sample diagnostics finished")
message("Outputs saved under:\n  ", out_dir)
