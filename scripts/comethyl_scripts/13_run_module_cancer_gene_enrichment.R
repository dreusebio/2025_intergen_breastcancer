#!/usr/bin/env Rscript

# ============================================================
# Cancer gene overlap/enrichment analysis for comethyl modules
#
# Purpose:
#   Compare genes annotated to comethyl modules against a cancer
#   gene census/database and test whether each module is enriched
#   for cancer-associated genes.
#
# Input:
#   1. Cancer gene database CSV/XLSX/TSV
#   2. comethyl Annotated_Regions.xlsx
#
# Output:
#   Automatically creates:
#   comethyl_output/13_cancer_gene_enrichment/<same subfolders as 12_annotation>/
#
# Example:
#   If annot_file is:
#   comethyl_output/12_annotation/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/Annotated_Regions.xlsx
#
#   Output will be:
#   comethyl_output/13_cancer_gene_enrichment/cov3_75pct/covMin10_methSD0p07/v2_exclude_protected_pcs/
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# qvalue is a Bioconductor package — load gracefully so the script
# still runs if not installed, but warn clearly.
qvalue_available <- requireNamespace("qvalue", quietly = TRUE)
if (!qvalue_available) {
  warning(
    "Package 'qvalue' is not installed. Q-value columns will be NA.\n",
    "Install with: BiocManager::install('qvalue')"
  )
}

# ============================================================
# Arguments
# ============================================================

option_list <- list(

  make_option(
    "--cancer_file",
    type = "character",
    help = "Path to cancer gene census/database file. Supported: csv, tsv, txt, xlsx, xls"
  ),

  make_option(
    "--annot_file",
    type = "character",
    help = "Path to comethyl Annotated_Regions.xlsx"
  ),

  make_option(
    "--out_dir",
    type = "character",
    default = NULL,
    help = "Optional exact output directory. If not provided, script infers output directory from annot_file"
  ),

  make_option(
    "--output_step_name",
    type = "character",
    default = "13_cancer_gene_enrichment",
    help = "Output step folder to create under comethyl_output [default: %default]"
  ),

  make_option(
    "--cancer_fdr_cutoff",
    type = "double",
    default = 0.05,
    help = "FDR cutoff for cancer gene enrichment significance [default: %default]"
  ),

  make_option(
    "--module_col",
    type = "character",
    default = "module",
    help = "Column name in Annotated_Regions.xlsx containing module assignment [default: %default]"
  ),

  make_option(
    "--gene_col",
    type = "character",
    default = "gene_symbol",
    help = "Column name in Annotated_Regions.xlsx containing gene symbols [default: %default]"
  ),

  make_option(
    "--module_trait_file",
    type = "character",
    default = NULL,
    help = "Optional module-trait association file. Supported: csv, tsv, txt, xlsx, xls"
  ),

  make_option(
    "--trait_module_col",
    type = "character",
    default = "Module",
    help = "Column name in module-trait file containing module names [default: %default]"
  ),

  make_option(
    "--trait_col",
    type = "character",
    default = "trait",
    help = "Column name in module-trait file containing trait/outcome names [default: %default]"
  ),

  make_option(
    "--trait_fdr_col",
    type = "character",
    default = "NONE",
    help = "Column name in module-trait file containing adjusted p-value/FDR [default: %default]"
  ),

  make_option(
    "--trait_p_col",
    type = "character",
    default = NULL,
    help = "Optional p-value column name in module-trait file"
  ),

  make_option(
    "--trait_effect_col",
    type = "character",
    default = NULL,
    help = "Optional effect/correlation column name in module-trait file"
  ),

  make_option(
    "--trait_fdr_cutoff",
    type = "double",
    default = 0.05,
    help = "FDR cutoff for defining outcome-associated modules [default: %default]"
  ),

  make_option(
    "--cancer_keyword",
    type = "character",
    default = "breast",
    help = "Cancer keyword to flag in tumour type columns, e.g. breast, ovarian, colorectal [default: %default]"
  ),

  make_option(
    "--top_cancer_types",
    type = "integer",
    default = 15,
    help = "Number of top cancer types to show in bubble plot [default: %default]"
  ),

  make_option(
    "--focus_modules",
    type = "character",
    default = "yellow,brown,tan,salmon,magenta,darkgrey",
    help = "Comma-separated list of modules for the focused Figure 4D plots [default: %default]"
  ),

  make_option(
    "--plot_width",
    type = "double",
    default = 10,
    help = "Plot width in inches [default: %default]"
  ),

  make_option(
    "--plot_height",
    type = "double",
    default = 7,
    help = "Plot height in inches [default: %default]"
  ),

  make_option(
    "--dpi",
    type = "integer",
    default = 600,
    help = "DPI for PNG/TIFF outputs [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

required_args <- c("cancer_file", "annot_file")
missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]

if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

# Parse focus modules from comma-separated string
focus_modules_vec <- str_trim(str_split(opt$focus_modules, ",")[[1]])
focus_modules_vec <- focus_modules_vec[focus_modules_vec != ""]
message("Focus modules: ", paste(focus_modules_vec, collapse = ", "))

# ============================================================
# Helper functions
# ============================================================

read_any_table <- function(path) {

  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }

  ext <- tolower(tools::file_ext(path))

  if (ext == "csv") {
    read_csv(path, show_col_types = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    read_tsv(path, show_col_types = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    read_excel(path)
  } else {
    stop("Unsupported file extension: ", ext)
  }
}

clean_gene_symbol <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    toupper()
}

safe_neglog10 <- function(x) {
  x <- as.numeric(x)
  x <- ifelse(is.na(x), NA_real_, x)
  x <- ifelse(x == 0, .Machine$double.xmin, x)
  -log10(x)
}

infer_out_dir <- function(annot_file, output_step_name) {

  annot_dir <- dirname(normalizePath(annot_file, mustWork = TRUE))

  path_parts <- strsplit(annot_dir, .Platform$file.sep, fixed = TRUE)[[1]]

  comethyl_idx <- which(path_parts == "comethyl_output")

  if (length(comethyl_idx) == 0) {
    stop(
      "Could not infer output directory because 'comethyl_output' was not found in annot_file path.\n",
      "Please provide --out_dir manually."
    )
  }

  annot_idx <- which(path_parts == "12_annotation")

  if (length(annot_idx) == 0) {
    stop(
      "Could not infer output subdirectory because '12_annotation' was not found in annot_file path.\n",
      "Please provide --out_dir manually."
    )
  }

  comethyl_root <- paste(
    path_parts[1:comethyl_idx[1]],
    collapse = .Platform$file.sep
  )

  analysis_subpath <- paste(
    path_parts[(annot_idx[1] + 1):length(path_parts)],
    collapse = .Platform$file.sep
  )

  file.path(comethyl_root, output_step_name, analysis_subpath)
}

make_plot_exports <- function(plot, filename_base, out_dir, width, height, dpi) {

  ggsave(
    filename = file.path(out_dir, paste0(filename_base, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )

  ggsave(
    filename = file.path(out_dir, paste0(filename_base, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )

  ggsave(
    filename = file.path(out_dir, paste0(filename_base, ".tiff")),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw"
  )
}

# ============================================================
# Module color mapping
#
# WGCNA module names are themselves color names, so we can map
# each module directly to its corresponding R color. A small
# set of adjustments are made for modules whose standard names
# are too light or ambiguous for bar fills on a white background:
#   - "white" -> medium grey (otherwise invisible)
#   - "grey"  -> removed before plotting (unassigned regions)
#   - "darkgrey" -> a readable dark grey distinct from "grey"
# ============================================================

module_fill_color <- function(module_names) {
  # For each module name, try to use the name directly as an R color.
  # Fall back to a neutral mid-grey for names that are not valid R colors.
  sapply(module_names, function(m) {
    m_lower <- tolower(m)
    if (m_lower == "white") return("#AAAAAA")      # white is invisible on white bg
    if (m_lower == "lightcyan") return("#A8D8EA")  # slightly deeper for visibility
    if (m_lower == "lightyellow") return("#F7E060") # slightly deeper for visibility
    # Test whether the name is a known R color
    is_valid <- tryCatch({
      col2rgb(m_lower)
      TRUE
    }, error = function(e) FALSE)
    if (is_valid) return(m_lower)
    return("#888888")  # fallback grey for unknown module names
  }, USE.NAMES = FALSE)
}

# ============================================================
# Create output directory
# ============================================================

if (is.null(opt$out_dir)) {
  opt$out_dir <- infer_out_dir(
    annot_file = opt$annot_file,
    output_step_name = opt$output_step_name
  )
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("============================================================")
message("Cancer gene enrichment analysis")
message("============================================================")
message("Cancer file: ", opt$cancer_file)
message("Annotation file: ", opt$annot_file)
message("Output directory: ", opt$out_dir)
message("Cancer keyword: ", opt$cancer_keyword)
message("============================================================")

# ============================================================
# Read cancer census/database
# ============================================================

message("Reading cancer file...")

cancer <- read_any_table(opt$cancer_file)

required_cancer_cols <- c("Gene Symbol")
missing_cancer_cols <- setdiff(required_cancer_cols, names(cancer))

if (length(missing_cancer_cols) > 0) {
  stop(
    "Cancer file is missing required columns: ",
    paste(missing_cancer_cols, collapse = ", "),
    "\nAvailable columns are:\n",
    paste(names(cancer), collapse = ", ")
  )
}

optional_cancer_cols <- c(
  "Name",
  "Tier",
  "Hallmark",
  "Somatic",
  "Germline",
  "Tumour Types(Somatic)",
  "Tumour Types(Germline)",
  "Cancer Syndrome",
  "Role in Cancer",
  "Mutation Types",
  "Synonyms"
)

for (cc in optional_cancer_cols) {
  if (!cc %in% names(cancer)) {
    cancer[[cc]] <- NA_character_
  }
}

cancer_clean <- cancer %>%
  mutate(
    cancer_gene_symbol = clean_gene_symbol(`Gene Symbol`),
    Synonyms = ifelse(is.na(Synonyms), "", Synonyms)
  )

cancer_symbols <- cancer_clean %>%
  mutate(match_symbol = cancer_gene_symbol)

cancer_synonyms <- cancer_clean %>%
  separate_rows(Synonyms, sep = ",") %>%
  mutate(match_symbol = clean_gene_symbol(Synonyms)) %>%
  filter(!is.na(match_symbol), match_symbol != "")

cancer_lookup <- bind_rows(cancer_symbols, cancer_synonyms) %>%
  filter(!is.na(match_symbol), match_symbol != "") %>%
  distinct(match_symbol, cancer_gene_symbol, .keep_all = TRUE)

message("Unique cancer genes in database: ", n_distinct(cancer_clean$cancer_gene_symbol))
message("Cancer lookup symbols including synonyms: ", n_distinct(cancer_lookup$match_symbol))

# ============================================================
# Read comethyl annotation file
# ============================================================

message("Reading comethyl annotation file...")

if (!file.exists(opt$annot_file)) {
  stop("Annotation file does not exist: ", opt$annot_file)
}

annot <- read_excel(opt$annot_file)

message("Available columns in annotation file:")
message(paste(names(annot), collapse = ", "))

if (!opt$module_col %in% names(annot)) {
  stop(
    "module_col not found: ", opt$module_col,
    "\nAvailable columns are:\n",
    paste(names(annot), collapse = ", ")
  )
}

if (!opt$gene_col %in% names(annot)) {
  stop(
    "gene_col not found: ", opt$gene_col,
    "\nAvailable columns are:\n",
    paste(names(annot), collapse = ", ")
  )
}

module_genes <- annot %>%
  select(
    Module = all_of(opt$module_col),
    Gene_raw = all_of(opt$gene_col)
  ) %>%
  filter(!is.na(Module), !is.na(Gene_raw)) %>%
  mutate(
    Module = as.character(Module),
    Gene_raw = as.character(Gene_raw)
  ) %>%
  # ---- Remove grey module (unassigned/background regions) ----
  filter(tolower(Module) != "grey") %>%
  separate_rows(Gene_raw, sep = "[,;| ]+") %>%
  mutate(Gene = clean_gene_symbol(Gene_raw)) %>%
  filter(!is.na(Gene), Gene != "", Gene != "NA") %>%
  distinct(Module, Gene)

if (nrow(module_genes) == 0) {
  stop("No module-gene pairs were extracted. Check --module_col and --gene_col.")
}

n_grey_removed <- sum(tolower(annot[[opt$module_col]]) == "grey", na.rm = TRUE)
message("Grey module regions excluded (unassigned): ", n_grey_removed)
message("Unique modules (grey excluded): ", n_distinct(module_genes$Module))
message("Unique annotated comethyl genes: ", n_distinct(module_genes$Gene))

# ============================================================
# Define universe/background
# ============================================================

universe_genes <- unique(module_genes$Gene)

cancer_genes_in_universe <- cancer_lookup %>%
  filter(match_symbol %in% universe_genes) %>%
  distinct(match_symbol)

n_universe <- length(universe_genes)
n_cancer_universe <- nrow(cancer_genes_in_universe)

message("Cancer genes present in comethyl universe: ", n_cancer_universe)

if (n_cancer_universe == 0) {
  stop(
    "No cancer genes from the cancer database were found in the comethyl gene universe.\n",
    "Check whether gene symbols are in the expected columns."
  )
}

# ============================================================
# Module-level cancer gene overlap details
# ============================================================

module_overlap_details <- module_genes %>%
  left_join(cancer_lookup, by = c("Gene" = "match_symbol")) %>%
  filter(!is.na(cancer_gene_symbol)) %>%
  distinct(
    Module,
    Gene,
    cancer_gene_symbol,
    Name,
    Tier,
    Hallmark,
    Somatic,
    Germline,
    `Tumour Types(Somatic)`,
    `Tumour Types(Germline)`,
    `Cancer Syndrome`,
    `Role in Cancer`,
    `Mutation Types`
  ) %>%
  mutate(
    tumour_types_combined = paste(
      `Tumour Types(Somatic)`,
      `Tumour Types(Germline)`,
      sep = ", "
    ),
    cancer_keyword_flag = str_detect(
      str_to_lower(tumour_types_combined),
      str_to_lower(opt$cancer_keyword)
    )
  )

message("Total module-cancer gene matches: ", nrow(module_overlap_details))

# ============================================================
# Helper: compute q-values from a vector of p-values
#
# With very few tests (e.g. 6 focus modules), Storey's pi0
# estimator can be unstable. safe_qvalue() falls back to BH
# FDR with a warning rather than crashing.
# ============================================================

safe_qvalue <- function(p_vec, label = "") {
  if (!qvalue_available) {
    return(rep(NA_real_, length(p_vec)))
  }
  tryCatch({
    qvalue::qvalue(p_vec)$qvalues
  }, error = function(e) {
    warning(
      "qvalue estimation failed", if (nchar(label) > 0) paste0(" (", label, ")"),
      ": ", conditionMessage(e),
      "\nFalling back to BH FDR for q-value column."
    )
    p.adjust(p_vec, method = "BH")
  })
}

# ============================================================
# Fisher enrichment per module
# ============================================================

message("Running Fisher enrichment per module...")

module_enrichment <- module_genes %>%
  group_by(Module) %>%
  summarise(
    module_gene_count = n_distinct(Gene),
    module_cancer_gene_count = n_distinct(
      Gene[Gene %in% cancer_genes_in_universe$match_symbol]
    ),
    cancer_genes = paste(
      sort(unique(Gene[Gene %in% cancer_genes_in_universe$match_symbol])),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    a = module_cancer_gene_count,
    b = module_gene_count - module_cancer_gene_count,
    c = n_cancer_universe - module_cancer_gene_count,
    d = n_universe - a - b - c,

    fisher_p = fisher.test(
      matrix(c(a, b, c, d), nrow = 2),
      alternative = "greater"
    )$p.value,

    odds_ratio = unname(
      fisher.test(
        matrix(c(a, b, c, d), nrow = 2),
        alternative = "greater"
      )$estimate
    )
  ) %>%
  ungroup() %>%
  mutate(
    Cancer_Enrichment_P = fisher_p,

    # --- Global BH FDR: adjusted across ALL modules ---
    Cancer_Enrichment_FDR_global = p.adjust(Cancer_Enrichment_P, method = "BH"),

    # --- Global q-value: Storey's method across ALL modules ---
    Cancer_Enrichment_Qvalue_global = safe_qvalue(Cancer_Enrichment_P, label = "all modules"),

    expected_cancer_genes = module_gene_count * n_cancer_universe / n_universe,

    # Convenience -log10 columns for plotting (use global FDR as default)
    minus_log10_P       = safe_neglog10(Cancer_Enrichment_P),
    minus_log10_FDR     = safe_neglog10(Cancer_Enrichment_FDR_global),
    minus_log10_Qvalue  = safe_neglog10(Cancer_Enrichment_Qvalue_global),

    # Significance label using global FDR (used in all-module plot 4A)
    cancer_enrichment_result_global = case_when(
      Cancer_Enrichment_FDR_global <= opt$cancer_fdr_cutoff & odds_ratio > 1 ~
        paste0("FDR <= ", opt$cancer_fdr_cutoff),
      Cancer_Enrichment_P < 0.05 & odds_ratio > 1 ~
        "Nominal p < 0.05",
      TRUE ~ "Not significant"
    )
  ) %>%
  select(
    Module,
    module_gene_count,
    module_cancer_gene_count,
    expected_cancer_genes,
    odds_ratio,
    Cancer_Enrichment_P,
    Cancer_Enrichment_FDR_global,
    Cancer_Enrichment_Qvalue_global,
    minus_log10_P,
    minus_log10_FDR,
    minus_log10_Qvalue,
    cancer_enrichment_result_global,
    cancer_genes
  ) %>%
  arrange(Cancer_Enrichment_P)

# ============================================================
# Cancer type summary
# ============================================================

cancer_type_summary <- module_overlap_details %>%
  mutate(
    tumour_types_combined = paste(
      `Tumour Types(Somatic)`,
      `Tumour Types(Germline)`,
      sep = ", "
    )
  ) %>%
  separate_rows(tumour_types_combined, sep = ",") %>%
  mutate(
    tumour_types_combined = str_trim(tumour_types_combined)
  ) %>%
  filter(
    !is.na(tumour_types_combined),
    tumour_types_combined != "",
    tumour_types_combined != "NA"
  ) %>%
  group_by(Module, tumour_types_combined) %>%
  summarise(
    n_genes = n_distinct(Gene),
    genes = paste(sort(unique(Gene)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Module, desc(n_genes))

# ============================================================
# Keyword-specific cancer overlap
# Example: breast cancer
# ============================================================

keyword_specific_overlap <- module_overlap_details %>%
  filter(cancer_keyword_flag) %>%
  arrange(Module, Gene)

keyword_module_summary <- keyword_specific_overlap %>%
  group_by(Module) %>%
  summarise(
    n_keyword_cancer_genes = n_distinct(Gene),
    keyword_cancer_genes = paste(sort(unique(Gene)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_keyword_cancer_genes))

message(
  "Modules with ",
  opt$cancer_keyword,
  " cancer-associated genes: ",
  n_distinct(keyword_module_summary$Module)
)

# ============================================================
# Optional outcome-associated module summary
# ============================================================

outcome_module_cancer_summary <- NULL

# ============================================================
# Save tables
# ============================================================

message("Saving output tables...")

wb <- createWorkbook()

addWorksheet(wb, "Cancer gene overlap details")
writeData(wb, "Cancer gene overlap details", module_overlap_details)

addWorksheet(wb, "Module enrichment")
writeData(wb, "Module enrichment", module_enrichment)

addWorksheet(wb, "Cancer type summary")
writeData(wb, "Cancer type summary", cancer_type_summary)

keyword_sheet_1 <- paste0(opt$cancer_keyword, "_cancer_overlap")
keyword_sheet_2 <- paste0(opt$cancer_keyword, "_module_summary")

addWorksheet(wb, keyword_sheet_1)
writeData(wb, keyword_sheet_1, keyword_specific_overlap)

addWorksheet(wb, keyword_sheet_2)
writeData(wb, keyword_sheet_2, keyword_module_summary)

if (!is.null(outcome_module_cancer_summary)) {
  addWorksheet(wb, "Outcome modules cancer")
  writeData(wb, "Outcome modules cancer", outcome_module_cancer_summary)
}

workbook_file <- file.path(
  opt$out_dir,
  "Comethyl_Module_CancerGene_Overlap_and_Enrichment.xlsx"
)

saveWorkbook(
  wb,
  workbook_file,
  overwrite = TRUE
)

write_csv(
  module_overlap_details,
  file.path(opt$out_dir, "Module_CancerGene_Overlap_Details.csv")
)

write_csv(
  module_enrichment,
  file.path(opt$out_dir, "Module_CancerGene_Enrichment.csv")
)

write_csv(
  cancer_type_summary,
  file.path(opt$out_dir, "Module_CancerType_Summary.csv")
)

write_csv(
  keyword_specific_overlap,
  file.path(opt$out_dir, paste0("Module_", opt$cancer_keyword, "_CancerGene_Overlap.csv"))
)

write_csv(
  keyword_module_summary,
  file.path(opt$out_dir, paste0("Module_", opt$cancer_keyword, "_CancerGene_Summary.csv"))
)

if (!is.null(outcome_module_cancer_summary)) {
  write_csv(
    outcome_module_cancer_summary,
    file.path(opt$out_dir, "Outcome_Associated_Module_CancerGene_Summary.csv")
  )
}

# ============================================================
# Figure 4A: Module cancer-gene enrichment (all modules, colored)
# ============================================================

message("Generating Figure 4A...")

plot_enrichment <- module_enrichment %>%
  filter(module_cancer_gene_count > 0) %>%
  mutate(
    Module   = fct_reorder(Module, minus_log10_FDR),
    bar_fill = module_fill_color(as.character(Module)),
    sig_label = case_when(
      Cancer_Enrichment_FDR_global <= opt$cancer_fdr_cutoff ~
        paste0("FDR=", signif(Cancer_Enrichment_FDR_global, 2)),
      Cancer_Enrichment_P < 0.05 ~
        paste0("p=", signif(Cancer_Enrichment_P, 2)),
      TRUE ~ ""
    )
  )

if (nrow(plot_enrichment) > 0) {

  fill_colors_4A <- setNames(plot_enrichment$bar_fill, as.character(plot_enrichment$Module))

  p1 <- ggplot(
    plot_enrichment,
    aes(
      x = Module,
      y = minus_log10_FDR,
      fill = Module
    )
  ) +
    geom_col(color = "black", linewidth = 0.25) +
    scale_fill_manual(values = fill_colors_4A, guide = "none") +
    coord_flip() +
    geom_hline(
      yintercept = -log10(opt$cancer_fdr_cutoff),
      linetype = "dashed",
      linewidth = 0.4,
      color = "red"
    ) +
    geom_text(
      aes(label = sig_label),
      hjust = -0.1,
      size = 3,
      color = "black"
    ) +
    labs(
      title = "Cancer gene enrichment across comethyl modules",
      subtitle = "Fisher's exact test; BH FDR adjusted across all modules; grey module excluded",
      x = "Comethyl module",
      y = expression(-log[10]("BH FDR, all modules")),
      caption = paste0("Dashed line: FDR = ", opt$cancer_fdr_cutoff)
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title  = element_text(face = "bold"),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11)
    )

  make_plot_exports(
    plot = p1,
    filename_base = "Figure4A_Module_CancerGene_Enrichment",
    out_dir = opt$out_dir,
    width = opt$plot_width,
    height = opt$plot_height,
    dpi = opt$dpi
  )
}

# ============================================================
# Figure 4B: Module x cancer type bubble plot (colored by module)
# ============================================================

message("Generating Figure 4B...")

if (nrow(cancer_type_summary) > 0) {

  top_types <- cancer_type_summary %>%
    group_by(tumour_types_combined) %>%
    summarise(total_genes = sum(n_genes), .groups = "drop") %>%
    arrange(desc(total_genes)) %>%
    slice_head(n = opt$top_cancer_types) %>%
    pull(tumour_types_combined)

  bubble_df <- cancer_type_summary %>%
    filter(tumour_types_combined %in% top_types) %>%
    left_join(module_enrichment, by = "Module") %>%
    mutate(
      Module = fct_reorder(Module, module_cancer_gene_count),
      tumour_types_combined = fct_reorder(tumour_types_combined, n_genes),
      point_color = module_fill_color(as.character(Module))
    )

  fill_colors_4B <- setNames(
    bubble_df$point_color,
    as.character(bubble_df$Module)
  )

  p2 <- ggplot(
    bubble_df,
    aes(
      x = tumour_types_combined,
      y = Module,
      size = n_genes,
      fill = Module
    )
  ) +
    geom_point(shape = 21, color = "black", stroke = 0.3, alpha = 0.85) +
    scale_fill_manual(values = fill_colors_4B, guide = "none") +
    scale_size_continuous(range = c(2, 10), name = "Gene count") +
    labs(
      title = "Cancer types represented among cancer genes in comethyl modules",
      subtitle = paste0(
        "Top ", opt$top_cancer_types, " tumour types by overlap count\n",
        "(grey module excluded)"
      ),
      x = "Cancer type",
      y = "Comethyl module"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title  = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 11)
    )

  make_plot_exports(
    plot = p2,
    filename_base = "Figure4B_Module_CancerType_BubblePlot",
    out_dir = opt$out_dir,
    width  = opt$plot_width + 2,
    height = opt$plot_height,
    dpi    = opt$dpi
  )
}

# ============================================================
# Figure 4C: Keyword-specific cancer gene count by module (colored)
# ============================================================

message("Generating Figure 4C...")

if (nrow(keyword_module_summary) > 0) {

  keyword_plot_df <- keyword_module_summary %>%
    mutate(
      Module     = fct_reorder(Module, n_keyword_cancer_genes),
      bar_fill   = module_fill_color(as.character(Module))
    )

  fill_colors_4C <- setNames(keyword_plot_df$bar_fill, as.character(keyword_plot_df$Module))

  p3 <- ggplot(
    keyword_plot_df,
    aes(
      x = Module,
      y = n_keyword_cancer_genes,
      fill = Module
    )
  ) +
    geom_col(color = "black", linewidth = 0.25) +
    scale_fill_manual(values = fill_colors_4C, guide = "none") +
    coord_flip() +
    labs(
      title = paste0("Modules containing ", opt$cancer_keyword, " cancer-associated genes"),
      subtitle = "Cancer keyword detected in tumour type annotations\n(grey module excluded)",
      x = "Comethyl module",
      y = paste0("Number of ", opt$cancer_keyword, " cancer genes")
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title  = element_text(face = "bold"),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11)
    )

  make_plot_exports(
    plot = p3,
    filename_base = paste0("Figure4C_Module_", opt$cancer_keyword, "_CancerGene_Count"),
    out_dir = opt$out_dir,
    width  = opt$plot_width,
    height = opt$plot_height,
    dpi    = opt$dpi
  )
}

# ============================================================
# Figure 4D: Focused module analysis
#
# Uses the full gene universe for Fisher's test (same p-values
# as all-module analysis), but applies multiple-testing
# correction WITHIN the focus modules only — fewer tests means
# less penalty, which is statistically appropriate when modules
# were pre-selected for biological reasons.
#
# Three panels:
#   4D-i   : Enrichment bar chart colored by -log10(focus BH FDR)
#   4D-ii  : Enrichment bar chart colored by -log10(focus q-value)
#   4D-iii : Bubble plot of top cancer types for focus modules
# ============================================================

message("Generating Figure 4D (focused modules: ", paste(focus_modules_vec, collapse = ", "), ")...")

# Match focus modules case-insensitively
present_modules <- unique(module_enrichment$Module)
focus_matched   <- present_modules[tolower(present_modules) %in% tolower(focus_modules_vec)]

not_found <- focus_modules_vec[!tolower(focus_modules_vec) %in% tolower(present_modules)]
if (length(not_found) > 0) {
  message(
    "Warning: the following focus modules were not found in the data and will be skipped: ",
    paste(not_found, collapse = ", ")
  )
}

if (length(focus_matched) == 0) {
  message("No focus modules matched the data; skipping Figure 4D.")
} else {

  # -------------------------------------------------------
  # Compute focus-adjusted corrections on the focus subset
  # p-values are IDENTICAL to the all-module analysis —
  # only the multiple-testing adjustment changes.
  # -------------------------------------------------------
  focus_enrichment <- module_enrichment %>%
    filter(Module %in% focus_matched) %>%
    mutate(
      # BH FDR adjusted within focus modules only
      Cancer_Enrichment_FDR_focus = p.adjust(Cancer_Enrichment_P, method = "BH"),

      # Q-value within focus modules only
      Cancer_Enrichment_Qvalue_focus = safe_qvalue(
        Cancer_Enrichment_P,
        label = paste0("focus modules (n=", length(focus_matched), ")")
      ),

      minus_log10_FDR_focus    = safe_neglog10(Cancer_Enrichment_FDR_focus),
      minus_log10_Qvalue_focus = safe_neglog10(Cancer_Enrichment_Qvalue_focus),

      cancer_enrichment_result_focus = case_when(
        Cancer_Enrichment_FDR_focus <= opt$cancer_fdr_cutoff & odds_ratio > 1 ~
          paste0("FDR <= ", opt$cancer_fdr_cutoff),
        Cancer_Enrichment_P < 0.05 & odds_ratio > 1 ~
          "Nominal p < 0.05",
        TRUE ~ "Not significant"
      )
    )

  message(
    "Focus module corrections computed (n=", length(focus_matched), " modules):\n",
    "  Raw p-values are unchanged from the all-module analysis.\n",
    "  BH FDR and q-value are now adjusted within focus modules only."
  )

  # -------------------------------------------------------
  # Save focus enrichment table separately
  # -------------------------------------------------------
  write_csv(
    focus_enrichment %>%
      select(
        Module,
        module_gene_count,
        module_cancer_gene_count,
        expected_cancer_genes,
        odds_ratio,
        Cancer_Enrichment_P,
        Cancer_Enrichment_FDR_global,
        Cancer_Enrichment_FDR_focus,
        Cancer_Enrichment_Qvalue_global,
        Cancer_Enrichment_Qvalue_focus,
        minus_log10_P,
        minus_log10_FDR,
        minus_log10_FDR_focus,
        minus_log10_Qvalue,
        minus_log10_Qvalue_focus,
        cancer_enrichment_result_global,
        cancer_enrichment_result_focus,
        cancer_genes
      ) %>%
      arrange(Cancer_Enrichment_P),
    file.path(opt$out_dir, "FocusModules_CancerGene_Enrichment_AllCorrections.csv")
  )

  # ---- 4D-i: BH FDR (focus-adjusted) bar chart ----
  focus_enr_plot <- focus_enrichment %>%
    filter(module_cancer_gene_count > 0) %>%
    mutate(
      Module    = fct_reorder(Module, minus_log10_FDR_focus),
      bar_fill  = module_fill_color(as.character(Module)),
      sig_label = case_when(
        Cancer_Enrichment_FDR_focus <= opt$cancer_fdr_cutoff ~
          paste0("FDR=", signif(Cancer_Enrichment_FDR_focus, 2)),
        Cancer_Enrichment_P < 0.05 ~
          paste0("p=", signif(Cancer_Enrichment_P, 2)),
        TRUE ~ ""
      )
    )

  if (nrow(focus_enr_plot) > 0) {

    fill_4Di <- setNames(focus_enr_plot$bar_fill, as.character(focus_enr_plot$Module))

    p4i <- ggplot(
      focus_enr_plot,
      aes(x = Module, y = minus_log10_FDR_focus, fill = Module)
    ) +
      geom_col(color = "black", linewidth = 0.25) +
      scale_fill_manual(values = fill_4Di, guide = "none") +
      coord_flip() +
      geom_hline(
        yintercept = -log10(opt$cancer_fdr_cutoff),
        linetype = "dashed", linewidth = 0.4, color = "red"
      ) +
      geom_text(aes(label = sig_label), hjust = -0.1, size = 3.5, color = "black") +
      labs(
        title = "Cancer gene enrichment – modules of interest",
        subtitle = paste0(
          "BH FDR adjusted within focus modules only (n=", length(focus_matched), ")\n",
          "Full universe background; grey module excluded"
        ),
        x = "Comethyl module",
        y = expression(-log[10]("BH FDR, focus-adjusted")),
        caption = paste0("Dashed line: FDR = ", opt$cancer_fdr_cutoff)
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title  = element_text(face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11)
      )

    make_plot_exports(
      plot = p4i,
      filename_base = "Figure4D-i_FocusModules_FDR_Enrichment",
      out_dir = opt$out_dir,
      width  = opt$plot_width,
      height = max(3, length(focus_matched) * 0.7 + 2),
      dpi    = opt$dpi
    )
  }

  # ---- 4D-ii: Q-value (focus-adjusted) bar chart ----
  focus_qval_plot <- focus_enrichment %>%
    filter(module_cancer_gene_count > 0) %>%
    mutate(
      Module    = fct_reorder(Module, minus_log10_Qvalue_focus),
      bar_fill  = module_fill_color(as.character(Module)),
      sig_label = case_when(
        !is.na(Cancer_Enrichment_Qvalue_focus) &
          Cancer_Enrichment_Qvalue_focus <= opt$cancer_fdr_cutoff ~
          paste0("q=", signif(Cancer_Enrichment_Qvalue_focus, 2)),
        Cancer_Enrichment_P < 0.05 ~
          paste0("p=", signif(Cancer_Enrichment_P, 2)),
        TRUE ~ ""
      )
    )

  if (nrow(focus_qval_plot) > 0) {

    fill_4Dii <- setNames(focus_qval_plot$bar_fill, as.character(focus_qval_plot$Module))

    qval_subtitle <- if (qvalue_available) {
      paste0(
        "Storey q-value adjusted within focus modules only (n=", length(focus_matched), ")\n",
        "Full universe background; grey module excluded"
      )
    } else {
      "NOTE: qvalue package not available — q-value column reflects BH FDR fallback"
    }

    p4ii <- ggplot(
      focus_qval_plot,
      aes(x = Module, y = minus_log10_Qvalue_focus, fill = Module)
    ) +
      geom_col(color = "black", linewidth = 0.25) +
      scale_fill_manual(values = fill_4Dii, guide = "none") +
      coord_flip() +
      geom_hline(
        yintercept = -log10(opt$cancer_fdr_cutoff),
        linetype = "dashed", linewidth = 0.4, color = "red"
      ) +
      geom_text(aes(label = sig_label), hjust = -0.1, size = 3.5, color = "black") +
      labs(
        title = "Cancer gene enrichment – modules of interest (q-value)",
        subtitle = qval_subtitle,
        x = "Comethyl module",
        y = expression(-log[10]("q-value, focus-adjusted")),
        caption = paste0("Dashed line: q-value = ", opt$cancer_fdr_cutoff)
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title  = element_text(face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11)
      )

    make_plot_exports(
      plot = p4ii,
      filename_base = "Figure4D-ii_FocusModules_Qvalue_Enrichment",
      out_dir = opt$out_dir,
      width  = opt$plot_width,
      height = max(3, length(focus_matched) * 0.7 + 2),
      dpi    = opt$dpi
    )
  }

  # ---- 4D-iii: Bubble plot of top cancer types for focus modules ----
  focus_top_types <- cancer_type_summary %>%
    filter(Module %in% focus_matched) %>%
    group_by(tumour_types_combined) %>%
    summarise(total_genes = sum(n_genes), .groups = "drop") %>%
    arrange(desc(total_genes)) %>%
    slice_head(n = opt$top_cancer_types) %>%
    pull(tumour_types_combined)

  focus_bubble_df <- cancer_type_summary %>%
    filter(
      Module %in% focus_matched,
      tumour_types_combined %in% focus_top_types
    ) %>%
    left_join(
      focus_enrichment %>%
        select(Module, module_cancer_gene_count,
               Cancer_Enrichment_FDR_focus, Cancer_Enrichment_Qvalue_focus),
      by = "Module"
    ) %>%
    mutate(
      Module = fct_reorder(Module, module_cancer_gene_count),
      tumour_types_combined = fct_reorder(tumour_types_combined, n_genes, .fun = sum),
      point_color = module_fill_color(as.character(Module))
    )

  if (nrow(focus_bubble_df) == 0) {
    message("No cancer type data found for focus modules; skipping Figure 4D-iii.")
  } else {

    fill_4Diii <- setNames(
      focus_bubble_df$point_color,
      as.character(focus_bubble_df$Module)
    )

    p4iii <- ggplot(
      focus_bubble_df,
      aes(
        x = tumour_types_combined,
        y = Module,
        size = n_genes,
        fill = Module
      )
    ) +
      geom_point(shape = 21, color = "black", stroke = 0.4, alpha = 0.88) +
      scale_fill_manual(values = fill_4Diii, guide = "none") +
      scale_size_continuous(
        range  = c(2, 12),
        name   = "Gene count",
        breaks = pretty_breaks(n = 4)
      ) +
      labs(
        title = "Cancer types associated with modules of interest",
        subtitle = paste0(
          "Top ", opt$top_cancer_types, " cancer types ranked within focus modules\n",
          "Modules: ", paste(sort(focus_modules_vec), collapse = ", "),
          "; grey module excluded"
        ),
        x = "Cancer type",
        y = "Comethyl module"
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title  = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 12),
        legend.position = "right"
      )

    make_plot_exports(
      plot = p4iii,
      filename_base = "Figure4D-iii_FocusModules_CancerType_BubblePlot",
      out_dir = opt$out_dir,
      width  = opt$plot_width + 2,
      height = opt$plot_height,
      dpi    = opt$dpi
    )

    message(
      "Figure 4D-iii saved: top ", opt$top_cancer_types,
      " cancer types across ", length(focus_matched), " focus modules."
    )
  }
}

# ============================================================
# Done
# ============================================================

message("============================================================")
message("Done.")
message("Results saved to:")
message(opt$out_dir)
message("")
message("Main workbook:")
message(workbook_file)
message("")
message("Main figure files:")
message(file.path(opt$out_dir, "Figure4A_Module_CancerGene_Enrichment.pdf"))
message(file.path(opt$out_dir, "Figure4B_Module_CancerType_BubblePlot.pdf"))
message(file.path(opt$out_dir, paste0("Figure4C_Module_", opt$cancer_keyword, "_CancerGene_Count.pdf")))
message(file.path(opt$out_dir, "Figure4D-i_FocusModules_FDR_Enrichment.pdf"))
message(file.path(opt$out_dir, "Figure4D-ii_FocusModules_Qvalue_Enrichment.pdf"))
message(file.path(opt$out_dir, "Figure4D-iii_FocusModules_CancerType_BubblePlot.pdf"))
message(file.path(opt$out_dir, "FocusModules_CancerGene_Enrichment_AllCorrections.csv"))
message("============================================================")