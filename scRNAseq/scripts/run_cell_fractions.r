#' ============================================================================
#' Cell Type Composition Analysis
#' ============================================================================
#'
#' Manuscript: Multi-omics Reveals How Diet and Exercise Improve Peripheral Neuropathy
#' Authors:    Eid SA, Guo K, Noureldein MH, Jang D-G, Allouch AM, Miller CM, Kiriluk CP,
#   Lentz W, Boutary S, Mule JJ, Pennathur S, Bhatt DK, Feldman EL, Hur J
#'
#' Contact:
#   Dr. Kai Guo        <guokai8@gmail.com>
#   Dr. Junguk Hur     <junguk.hur@med.UND.edu>
#'
#' Description:
#'   Summarizes cell type composition across experimental groups. Calculates cell counts
#   and fractions per group and generates stacked barplots showing absolute and
#   percentage composition.
#'
#' Input:
#'   steph_major.rdata or step_scwann.rdata (Seurat objects)
#'
#' Output:
#'   Composition tables, stacked barplots (absolute and percentage)
#'
#' Related figures:
#'   Fig. S3, Fig. S4
#'
#' ============================================================================
#!/usr/bin/env Rscript
# Overall cell-type composition tables and plots.
# - Loads a Seurat object from an .rdata file
# - Writes counts/fractions and stacked barplots by Group x celltype (and major if available)
# Output directory: ./251208_overall-cell-fractions

#!/usr/bin/env Rscript
# Overall cell-type composition tables and plots.
# - Loads a Seurat object from an .rdata file (defaults: steph_major.rdata, then step_scwann.rdata)
# - Writes counts/fractions and stacked barplots by Group x celltype (and major if available)
# Output directory: ./251208_overall-cell-fractions

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

# Set working directory to the script location (works in Rscript or RStudio)
script_dir <- tryCatch({
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- rstudioapi::getActiveDocumentContext()
    if (!is.null(ctx$path) && nzchar(ctx$path)) dirname(ctx$path) else getwd()
  } else {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", cmd_args[grepl("^--file=", cmd_args)])[1]
    if (!is.na(file_arg) && nzchar(file_arg)) dirname(normalizePath(file_arg)) else getwd()
  }
}, error = function(e) getwd())
setwd(script_dir)

# Config / CLI: pass input .rdata path as first argument; otherwise try defaults.
args <- commandArgs(trailingOnly = TRUE)
defaults <- c("steph_major.rdata", "step_scwann.rdata", "sc.rdata")
rdata_file <- if (length(args) >= 1) args[1] else defaults[1]
if (!file.exists(rdata_file)) {
  alt <- defaults[defaults != rdata_file]
  found <- alt[dir.exists(dirname(alt)) | file.exists(alt)]
  if (length(found)) rdata_file <- found[1]
}
if (!file.exists(rdata_file)) {
  stop("Could not find an input .rdata file. Tried: ", paste(c(args, defaults), collapse = ", "))
}

output_dir <- "251208_overall-cell-fractions"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load Seurat object
after_load <- load(rdata_file)
# Prefer object named sam; else first Seurat object found
got <- intersect(c("sam", after_load), after_load)
if (length(got) && inherits(get(got[1]), "Seurat")) {
  sam <- get(got[1])
} else {
  seurat_names <- Filter(function(x) inherits(get(x), "Seurat"), after_load)
  if (!length(seurat_names)) {
    stop("No Seurat object found in ", rdata_file)
  }
  sam <- get(seurat_names[1])
}

meta <- sam@meta.data
if (!all(c("Group", "celltype") %in% names(meta))) {
  stop("Metadata must contain Group and celltype columns.")
}

# Helper to plot and write composition
write_comp <- function(meta_tbl, xcol, group_col = "Group", label = "celltype") {
  df <- meta_tbl %>%
    count(.data[[group_col]], .data[[xcol]], name = "n") %>%
    group_by(.data[[group_col]]) %>%
    mutate(frac = n / sum(n)) %>%
    ungroup() %>%
    rename(Group = !!group_col, Category = !!xcol)

  csv_fn <- file.path(output_dir, paste0(label, "_counts_and_fractions.csv"))
  write_csv(df, csv_fn)

  p_pct <- ggplot(df, aes(x = Group, y = frac, fill = Category)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Fraction", title = paste0(label, " composition by group"), fill = label) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_cnt <- ggplot(df, aes(x = Group, y = n, fill = Category)) +
    geom_bar(stat = "identity") +
    labs(x = NULL, y = "Cells", title = paste0(label, " counts by group"), fill = label) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(output_dir, paste0(label, "_fraction_by_group.pdf")), p_pct, width = 8, height = 4)
  ggsave(file.path(output_dir, paste0(label, "_counts_by_group.pdf")),   p_cnt, width = 8, height = 4)
}

# Celltype-level composition
write_comp(meta, xcol = "celltype", label = "celltype")

# Major-level composition if present or derivable
# Optional: user-provided mapping from celltype -> major (named character vector)
# Default here mirrors run_step_final.r major groups:
# SC: mySC, SC3, nmSC, ImmSC
# Fib: EndoFib1, EndoFib2, EpiFib
# Mac: Mac
# Perineurial: Perineurial
# SMC: VSMC
# Pericytes: Pericytes
# Endo: Endo
major_map <- c(
  mySC = "SC", SC3 = "SC", nmSC = "SC", ImmSC = "SC",
  EndoFib1 = "Fib", EndoFib2 = "Fib", EpiFib = "Fib",
  Mac = "Mac",
  Perineurial = "Perineurial",
  VSMC = "SMC",
  Pericytes = "Pericytes",
  Endo = "Endo"
)
# If no explicit major column or mapping, derive from celltype by collapsing common patterns
auto_major_if_missing <- TRUE
derive_major <- function(ct) {
  out <- ct
  out[ct %in% c("ImmSC", "mySC", "nmSC")] <- "SC"
  out <- sub("_.*$", "", out)        # collapse suffix after first underscore
  out <- sub("-.*$", "", out)        # collapse suffix after first hyphen
  out
}

major_cols <- c("major", "Major", "major_cell", "MajorCell")
maj_col <- intersect(major_cols, names(meta))
if (length(maj_col)) {
  meta$MajorPlot <- meta[[maj_col[1]]]
  write_comp(meta, xcol = "MajorPlot", label = "major")
} else if (!is.null(major_map)) {
  meta$MajorPlot <- dplyr::coalesce(unname(major_map[meta$celltype]), meta$celltype)
  write_comp(meta, xcol = "MajorPlot", label = "major")
} else if (auto_major_if_missing) {
  meta$MajorPlot <- derive_major(meta$celltype)
  write_comp(meta, xcol = "MajorPlot", label = "major")
} else {
  message("No major column found and auto-derivation disabled; skipping major-level outputs.")
}

# Sample-level (replicate-level) composition if a sample/replicate column exists
sample_cols    <- c("sample", "Sample", "sample_id", "sampleID", "orig.ident", "orig_ident", "orig.ident2")
replicate_cols <- c("replicate", "Replicate", "rep", "Rep", "mouse", "Mouse", "animal", "Animal", "subject", "Subject")
sample_col  <- intersect(sample_cols, names(meta))
rep_col     <- intersect(replicate_cols, names(meta))

use_col <- NULL
label_colname <- NULL
if (length(rep_col)) {
  use_col <- rep_col[1]; label_colname <- use_col
} else if (length(sample_col)) {
  use_col <- sample_col[1]; label_colname <- use_col
}

if (!is.null(use_col)) {
  # order by Group then replicate for nicer plotting
  sample_levels <- meta %>%
    distinct(.data[[use_col]], Group) %>%
    arrange(Group, .data[[use_col]]) %>%
    pull(.data[[use_col]])

  write_comp_sample <- function(meta_tbl, xcol, label) {
    df <- meta_tbl %>%
      count(.data[[use_col]], .data[[xcol]], Group, name = "n") %>%
      group_by(.data[[use_col]]) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup() %>%
      rename(Replicate = !!use_col, Category = !!xcol)
    df$Replicate <- factor(df$Replicate, levels = sample_levels)

    csv_fn <- file.path(output_dir, paste0(label, "_by-replicate_counts_and_fractions.csv"))
    write_csv(df, csv_fn)

    group_df <- df %>%
      distinct(Replicate, Group)
    spans <- tibble::tibble(
      Replicate = sample_levels,
      idx = seq_along(sample_levels)
    ) %>%
      dplyr::left_join(group_df, by = "Replicate") %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(
        xmin = min(idx) - 0.5,
        xmax = max(idx) + 0.5,
        xmid = (min(idx) + max(idx)) / 2,
        .groups = "drop"
      )

    p_pct_core <- ggplot(df, aes(x = Replicate, y = frac, fill = Category)) +
      geom_col(position = "fill") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Replicate", y = "Fraction", title = paste0(label, " composition by replicate"), fill = label) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = "bottom") +
      guides(fill = guide_legend(order = 1, ncol = 1), color = guide_legend(order = 2))

    p_cnt_core <- ggplot(df, aes(x = Replicate, y = n, fill = Category)) +
      geom_col() +
      labs(x = "Replicate", y = "Cells", title = paste0(label, " counts by replicate"), fill = label) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = "bottom") +
      guides(fill = guide_legend(order = 1, ncol = 1), color = guide_legend(order = 2))

    # Thin segmented group bar below x-axis labels with one label per group span
    group_bar <- ggplot(spans) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.25, ymax = 0.75, fill = Group), color = NA) +
      geom_text(aes(x = xmid, y = 0.5, label = Group), size = 2.8, color = "black") +
      scale_x_continuous(expand = c(0, 0), limits = c(0.5, length(sample_levels) + 0.5)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
      guides(fill = "none") +
      theme_void(base_size = 9) +
      theme(plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5))

    p_pct <- patchwork::wrap_plots(p_pct_core, group_bar, ncol = 1, heights = c(0.88, 0.12), guides = "collect")
    p_cnt <- patchwork::wrap_plots(p_cnt_core, group_bar, ncol = 1, heights = c(0.88, 0.12), guides = "collect")

    ggsave(file.path(output_dir, paste0(label, "_fraction_by_replicate.pdf")), p_pct, width = 10, height = 5)
    ggsave(file.path(output_dir, paste0(label, "_counts_by_replicate.pdf")),   p_cnt, width = 10, height = 5)
  }

  write_comp_sample(meta, xcol = "celltype", label = "celltype")
  if ("MajorPlot" %in% names(meta)) {
    write_comp_sample(meta, xcol = "MajorPlot", label = "major")
  }
} else {
  message("No sample/replicate column found (tried sample: ", paste(sample_cols, collapse = ", "),
          "; replicate: ", paste(replicate_cols, collapse = ", "), "). Skipping replicate-level outputs.")
}

message("Done. Used file: ", rdata_file,
        " | Seurat object: ", if (exists("sam")) deparse(substitute(sam)) else "unknown",
        " | Outputs in: ", output_dir)
