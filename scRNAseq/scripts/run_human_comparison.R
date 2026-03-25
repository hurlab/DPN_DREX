#' ============================================================================
#' Human-Mouse Transcriptomic Comparison
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
#'   Compares mouse DREX scRNA-seq DEGs with human sural nerve transcriptomics from
#   diabetic peripheral neuropathy patients. Maps human genes to mouse orthologs
#   via homologene and identifies shared DEGs and enriched pathways.
#'
#' Input:
#'   Mouse DEG CSVs, human sural nerve DEG CSV (HSuralRNASeq_Human),
#   homologene package for ortholog mapping
#'
#' Output:
#'   Venn diagrams, overlapping DEG lists, heatmaps, Excel workbooks
#'
#' Related figures:
#'   Fig. S8, Tables S8-S9
#'
#' ============================================================================
################################################################################
# Initialize the workspace
################################################################################
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  setwd(dirname(script_path))
  base_dir <- dirname(script_path)
}

# Required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, VennDetail, tools, readxl, stringr,
               org.Mm.eg.db, org.Hs.eg.db, AnnotationDbi, homologene,
               pheatmap, tibble, openxlsx)

################################################################################
# Settings
################################################################################
human_deg_path <- "./HSuralRNASeq_Human/Group1vsGroup2_Table S3 - DEGs.csv"

# Folders to process
mouse_deg_dirs <- c(
  "../DEG",
  "../DEG_Major"
)

# Corresponding output folders (one per input)
output_roots <- c(
  "Output_DEG_vsHuman",
  "Output_DEG_Major_vsHuman"
)


################################################################################
# Initialize the workspace
################################################################################
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  setwd(dirname(script_path))
  base_dir <- dirname(script_path)
}

# Required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, VennDetail, tools, readxl, stringr,
               org.Mm.eg.db, org.Hs.eg.db, AnnotationDbi, homologene,
               pheatmap, tibble)

################################################################################
# Settings
################################################################################
human_deg_path <- "./HSuralRNASeq_Human/Group1vsGroup2_Table S3 - DEGs.csv"

# Folders to process
mouse_deg_dirs <- c(
  "../DEG",
  "../DEG_Major"
)

# Corresponding output folders (one per input)
output_roots <- c(
  "Output_DEG_vsHuman",
  "Output_DEG_Major_vsHuman"
)

# Create outputs
invisible(mapply(function(o) dir.create(o, showWarnings = FALSE, recursive = TRUE), output_roots))

################################################################################
# Helpers: safe file base, plotting, and I/O
################################################################################
safe_base <- function(x, limit = 150) {
  # shorten and remove problematic characters for long Windows paths
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  if (nchar(x) > limit) substr(x, 1, limit) else x
}

save_df <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, path, row.names = FALSE)
}

save_plot <- function(p, path, width = 7, height = 6, dpi = 300) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi, limitsize = FALSE)
}

sanitize_sheet <- function(x) {
  x <- gsub('[:\\\\/?*\\[\\]]', "_", x)   # invalid chars
  x <- ifelse(nchar(x) == 0, "Sheet", x)
  substr(x, 1, 31)                        # Excel 31-char limit
}


################################################################################
# Load human DEGs and map Human -> Mouse using homologene
################################################################################
# Human table columns (as given):
# "Symbol  Log2FoldChange  Pvalue  Adjust P value  Entrez ID  Ensembl ID  Annotation"
read_human_deg <- function(human_csv) {
  df <- read.csv(human_csv, check.names = FALSE, stringsAsFactors = FALSE)
  # normalize column names for robustness
  cn <- tolower(gsub("\\s+", "_", colnames(df)))
  colnames(df) <- cn
  # Expect: symbol, log2foldchange, adjust_p_value (or similar)
  sym_col <- "symbol"
  lfc_col <- if ("log2foldchange" %in% cn) "log2foldchange" else stop("Cannot find 'Log2FoldChange' in human file.")
  padj_col <- if ("adjust_p_value" %in% cn) "adjust_p_value" else if ("adjust_pvalue" %in% cn) "adjust_pvalue" else stop("Cannot find 'Adjust P value' in human file.")
  
  df %>%
    dplyr::filter(!is.na(.data[[sym_col]]), .data[[sym_col]] != "") %>%
    dplyr::select(Human_Symbol = all_of(sym_col),
                  Human_log2FC = all_of(lfc_col),
                  Human_padj   = all_of(padj_col))
}

human_deg <- read_human_deg(human_deg_path)

# Build Human -> Mouse symbol mapping (keep one mouse symbol per human)
build_human_to_mouse_map <- function() {
  hom <- homologene::homologeneData
  human <- hom %>%
    dplyr::filter(Taxonomy == 9606) %>%
    dplyr::select(HID, Human_Symbol = Gene.Symbol)
  mouse <- hom %>%
    dplyr::filter(Taxonomy == 10090) %>%
    dplyr::select(HID, Mouse_Entrez = Gene.ID, Mouse_Symbol = Gene.Symbol)
  joined <- dplyr::left_join(human, mouse, by = "HID") %>%
    dplyr::filter(!is.na(Mouse_Symbol))
  # pick the one with the lowest Entrez per human to make it deterministic
  joined %>%
    dplyr::arrange(Human_Symbol, Mouse_Entrez) %>%
    dplyr::group_by(Human_Symbol) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(Human_Symbol, Mouse_Symbol)
}

human2mouse <- build_human_to_mouse_map()

# Human significant set mapped to mouse symbols
human_sig_mapped <- human_deg %>%
  dplyr::inner_join(human2mouse, by = "Human_Symbol") %>%
  dplyr::distinct(Human_Symbol, .keep_all = TRUE)

################################################################################
# Read one mouse DEG file and return significant genes and log2FC
################################################################################
read_mouse_deg_file <- function(file_path) {
  df <- suppressWarnings(read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE))
  
  # First column contains gene symbols even if unnamed
  if (is.na(colnames(df)[1]) || colnames(df)[1] == "" || colnames(df)[1] %in% c("X", "...1")) {
    colnames(df)[1] <- "Gene"
  } else {
    # enforce common name
    colnames(df)[1] <- "Gene"
  }
  
  # normalize column names to locate p_val_adj and avg_log2FC
  cn <- tolower(gsub("\\s+", "_", colnames(df)))
  colnames(df) <- cn
  
  padj_col <- "p_val_adj"
  lfc_col  <- if ("avg_log2fc" %in% cn) "avg_log2fc" else if ("log2foldchange" %in% cn) "log2foldchange" else "avg_log2fc"
  if (!padj_col %in% cn) stop(paste("p_val_adj column not found in", basename(file_path)))
  if (!lfc_col %in% cn)  stop(paste("log2FC column not found in", basename(file_path)))
  
  df %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::filter(.data[[padj_col]] < 0.05) %>%
    dplyr::select(Gene = gene, Mouse_log2FC = all_of(lfc_col), Mouse_padj = all_of(padj_col)) %>%
    dplyr::distinct(Gene, .keep_all = TRUE)
}

################################################################################
# Compare human vs a single mouse file
################################################################################
compare_one_file <- function(mouse_file, out_dir, human_sig_mapped_tbl) {
  message("Processing: ", mouse_file)
  base <- safe_base(tools::file_path_sans_ext(basename(mouse_file)))
  
  # parse "comp_group" and "cell_type" from the base name split by "_"
  parts <- strsplit(base, "_", fixed = TRUE)[[1]]
  comp_group <- if (length(parts) >= 1) parts[1] else NA_character_
  cell_type  <- if (length(parts) >= 2) parts[2] else NA_character_
  
  # read mouse significant genes
  m_df <- read_mouse_deg_file(mouse_file)
  m_genes <- unique(m_df$Gene)
  
  # human mapped to mouse
  h_df <- human_sig_mapped_tbl %>%
    dplyr::filter(Human_padj < 0.05) %>%
    dplyr::select(Mouse_Symbol, Human_Symbol, Human_log2FC, Human_padj) %>%
    dplyr::distinct(Mouse_Symbol, .keep_all = TRUE)
  h_genes <- unique(h_df$Mouse_Symbol)
  
  # overlap and concordance
  both <- dplyr::inner_join(
    m_df %>% dplyr::select(Gene, Mouse_log2FC),
    h_df %>% dplyr::select(Mouse_Symbol, Human_log2FC, Human_Symbol),
    by = c("Gene" = "Mouse_Symbol")
  ) %>%
    dplyr::mutate(
      sign_mouse = sign(Mouse_log2FC),
      sign_human = sign(Human_log2FC),
      concordant = sign_mouse == sign_human
    )
  
  # summary row WITH comp_group and cell_type
  summary_row <- tibble::tibble(
    file            = basename(mouse_file),
    comp_group      = comp_group,
    cell_type       = cell_type,
    n_human_sig     = length(unique(h_df$Human_Symbol)),
    n_human_mapped  = length(h_genes),
    n_mouse_sig     = length(unique(m_df$Gene)),
    n_overlap       = nrow(both),
    jaccard         = ifelse(length(union(m_genes, h_genes)) == 0, NA_real_,
                             length(intersect(m_genes, h_genes)) / length(union(m_genes, h_genes))),
    n_concordant    = sum(both$concordant, na.rm = TRUE),
    n_discordant    = sum(!both$concordant, na.rm = TRUE)
  )
  
  # save overlap tables (unchanged)
  overlap_tbl <- both %>% dplyr::arrange(desc(concordant), desc(abs(Mouse_log2FC)))
  save_df(overlap_tbl, file.path(out_dir, paste0(base, "_Overlap_HumanMapped_vsMouse.csv")))
  save_df(tibble::tibble(MouseGene = sort(h_genes)), file.path(out_dir, paste0(base, "_HumanMapped_MouseGeneList.csv")))
  save_df(tibble::tibble(MouseGene = sort(unique(m_df$Gene))), file.path(out_dir, paste0(base, "_Mouse_SignifGeneList.csv")))
  save_df(tibble::tibble(MouseGene = sort(intersect(m_genes, h_genes))),  file.path(out_dir, paste0(base, "_Overlap_Genes.csv")))
  
  # Venn plot (record once, replay to files) — your latest working approach
  vd <- VennDetail::venndetail(list(Human_to_Mouse = h_genes, Mouse = m_genes))
  plot(vd, main = paste0("Overlap: ", base))
  rp <- recordPlot()
  dev.off()
  png(file.path(out_dir, paste0(base, "_Venn.png")), width = 2000, height = 1600, res = 300)
  replayPlot(rp); dev.off()
  pdf(file.path(out_dir, paste0(base, "_Venn.pdf")), width = 8, height = 6)
  replayPlot(rp); dev.off()
  
  return(summary_row)
}

# Helpers to match selected comparisons and cell types
normalize_key <- function(x) tolower(gsub("[^A-Za-z0-9]+", "", x))

is_selected <- function(comp_group, cell_type) {
  if (is.na(comp_group) || is.na(cell_type)) return(FALSE)
  allowed_comp <- c("HFD vs SD", "DR vs HFD", "EX vs HFD", "DREX vs HFD")
  allowed_cell <- c("majorSC", "mySC", "nmSC", "ImmSC")
  normalize_key(comp_group) %in% normalize_key(allowed_comp) &&
    normalize_key(cell_type) %in% normalize_key(allowed_cell)
}

################################################################################
# Compare for a whole folder and write a master summary (+ selected subset)
################################################################################
compare_folder <- function(input_dir, output_dir, human_sig_mapped_tbl) {
  message("=== Folder: ", input_dir, " ===")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
  files <- files[!basename(files) %>% startsWith(".")]
  if (length(files) == 0) {
    warning("No CSV files found in: ", input_dir)
    return(invisible(NULL))
  }
  
  # ---------- Full run (unchanged behavior) ----------
  per_file_out <- file.path(output_dir, "per_file")
  dir.create(per_file_out, showWarnings = FALSE, recursive = TRUE)
  
  all_summaries <- lapply(files, function(f) {
    tryCatch(
      compare_one_file(f, per_file_out, human_sig_mapped_tbl),
      error = function(e) {
        warning("Failed on ", basename(f), ": ", conditionMessage(e))
        parts <- strsplit(safe_base(tools::file_path_sans_ext(basename(f))), "_", fixed = TRUE)[[1]]
        comp_group <- if (length(parts) >= 1) parts[1] else NA_character_
        cell_type  <- if (length(parts) >= 2) parts[2] else NA_character_
        tibble::tibble(
          file = basename(f),
          comp_group = comp_group,
          cell_type = cell_type,
          n_human_sig = NA_integer_, n_human_mapped = NA_integer_, n_mouse_sig = NA_integer_,
          n_overlap = NA_integer_, jaccard = NA_real_,
          n_concordant = NA_integer_, n_discordant = NA_integer_
        )
      }
    )
  })
  
  master <- dplyr::bind_rows(all_summaries) %>%
    dplyr::arrange(dplyr::desc(n_overlap))
  
  save_df(master, file.path(output_dir, "Master_Summary.csv"))
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Summary")
  openxlsx::writeData(wb, "Summary", master)
  
  cell_types <- sort(unique(master$cell_type))
  cell_types <- cell_types[!is.na(cell_types)]
  for (ct in cell_types) {
    sheet <- sanitize_sheet(ct)
    openxlsx::addWorksheet(wb, sheet)
    df_ct <- master %>% dplyr::filter(cell_type == ct) %>%
      dplyr::arrange(dplyr::desc(n_overlap))
    openxlsx::writeData(wb, sheet, df_ct)
  }
  
  xlsx_path <- file.path(output_dir, "Master_Summary.xlsx")
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  message("Master Excel written: ", xlsx_path)
  
  # ---------- Selected run (new) ----------
  # pick files by parsing comp_group and cell_type from filename
  pick_files <- vapply(files, function(f) {
    base <- safe_base(tools::file_path_sans_ext(basename(f)))
    parts <- strsplit(base, "_", fixed = TRUE)[[1]]
    comp_group <- if (length(parts) >= 1) parts[1] else NA_character_
    cell_type  <- if (length(parts) >= 2) parts[2] else NA_character_
    is_selected(comp_group, cell_type)
  }, logical(1))
  
  sel_files <- files[pick_files]
  
  if (length(sel_files) > 0) {
    per_file_sel_out <- file.path(output_dir, "per_file_select")
    dir.create(per_file_sel_out, showWarnings = FALSE, recursive = TRUE)
    
    sel_summaries <- lapply(sel_files, function(f) {
      tryCatch(
        compare_one_file(f, per_file_sel_out, human_sig_mapped_tbl),
        error = function(e) {
          warning("Failed on ", basename(f), ": ", conditionMessage(e))
          parts <- strsplit(safe_base(tools::file_path_sans_ext(basename(f))), "_", fixed = TRUE)[[1]]
          comp_group <- if (length(parts) >= 1) parts[1] else NA_character_
          cell_type  <- if (length(parts) >= 2) parts[2] else NA_character_
          tibble::tibble(
            file = basename(f),
            comp_group = comp_group,
            cell_type = cell_type,
            n_human_sig = NA_integer_, n_human_mapped = NA_integer_, n_mouse_sig = NA_integer_,
            n_overlap = NA_integer_, jaccard = NA_real_,
            n_concordant = NA_integer_, n_discordant = NA_integer_
          )
        }
      )
    })
    
    master_sel <- dplyr::bind_rows(sel_summaries) %>%
      dplyr::arrange(dplyr::desc(n_overlap))
    
    save_df(master_sel, file.path(output_dir, "Master_Summary_selected.csv"))
    
    wb2 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb2, "Summary_selected")
    openxlsx::writeData(wb2, "Summary_selected", master_sel)
    
    sel_cell_types <- sort(unique(master_sel$cell_type))
    sel_cell_types <- sel_cell_types[!is.na(sel_cell_types)]
    for (ct in sel_cell_types) {
      sheet <- sanitize_sheet(ct)
      openxlsx::addWorksheet(wb2, sheet)
      df_ct <- master_sel %>% dplyr::filter(cell_type == ct) %>%
        dplyr::arrange(dplyr::desc(n_overlap))
      openxlsx::writeData(wb2, sheet, df_ct)
    }
    
    xlsx_path_sel <- file.path(output_dir, "Master_Summary_selected.xlsx")
    openxlsx::saveWorkbook(wb2, xlsx_path_sel, overwrite = TRUE)
    message("Selected Excel written: ", xlsx_path_sel)
  } else {
    message("No files matched the selected comparisons and cell types in: ", input_dir)
  }
  
  invisible(master)
}



################################################################################
# Run for each folder
################################################################################
invisible(mapply(function(in_dir, out_dir) {
  compare_folder(in_dir, out_dir, human_sig_mapped)
}, mouse_deg_dirs, output_roots))

message("All done.")


# 
# # Folders to process
# mouse_deg_dirs <- c(
#   "../DEG",
#   "../DEG_Major"
# )
# 
# # Corresponding output folders (one per input)
# output_roots <- c(
#   "Output_DEG_vsHuman",
#   "Output_DEG_Major_vsHuman"
# )


# 
# ################################################################################
# # DROP-IN: filters and selected-only runner (no dot-prefixed objects)
# ################################################################################
# 
# # Normalize strings for robust matching (remove non-alnum, lowercase)
# normalize_key <- function(x) {
#   tolower(gsub("[^A-Za-z0-9]+", "", x))
# }
# 
# # Allowed sets
# allowed_comp_raw  <- c("HFD vs SD", "DR vs HFD", "EX vs HFD", "DREX vs HFD")
# allowed_comp_keys <- normalize_key(allowed_comp_raw)
# 
# allowed_cell_raw  <- c("SC", "mySC", "nmSC", "ImmSC")
# allowed_cell_keys <- normalize_key(allowed_cell_raw)
# 
# # Return TRUE if comp_group and cell_type are in the allowed sets
# is_selected <- function(comp_group, cell_type) {
#   if (is.na(comp_group) || is.na(cell_type)) return(FALSE)
#   cg_ok <- normalize_key(comp_group) %in% allowed_comp_keys
#   ct_ok <- normalize_key(cell_type)  %in% allowed_cell_keys
#   cg_ok && ct_ok
# }
# 
# # Read a folder, run only on selected files, and write selected-only outputs
# compare_folder_selected <- function(input_dir, output_dir, human_sig_mapped_tbl) {
#   message("=== Folder (selected): ", input_dir, " ===")
#   dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
#   
#   files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
#   files <- files[!basename(files) %>% startsWith(".")]
#   if (length(files) == 0) {
#     warning("No CSV files found in: ", input_dir)
#     return(invisible(NULL))
#   }
#   
#   # Decide which files are selected based on parsed comp_group and cell_type
#   pick_files <- vapply(files, function(f) {
#     base <- safe_base(tools::file_path_sans_ext(basename(f)))
#     parts <- strsplit(base, "_", fixed = TRUE)[[1]]
#     comp_group <- if (length(parts) >= 1) parts[1] else NA_character_
#     cell_type  <- if (length(parts) >= 2) parts[2] else NA_character_
#     is_selected(comp_group, cell_type)
#   }, logical(1))
#   
#   sel_files <- files[pick_files]
#   if (length(sel_files) == 0) {
#     warning("No files matched the selected comparisons and cell types in: ", input_dir)
#     return(invisible(NULL))
#   }
#   
#   per_file_out <- file.path(output_dir, "per_file_select")
#   dir.create(per_file_out, showWarnings = FALSE, recursive = TRUE)
#   
#   all_summaries <- lapply(sel_files, function(f) {
#     tryCatch(
#       compare_one_file(f, per_file_out, human_sig_mapped_tbl),
#       error = function(e) {
#         warning("Failed on ", basename(f), ": ", conditionMessage(e))
#         parts <- strsplit(safe_base(tools::file_path_sans_ext(basename(f))), "_", fixed = TRUE)[[1]]
#         comp_group <- if (length(parts) >= 1) parts[1] else NA_character_
#         cell_type  <- if (length(parts) >= 2) parts[2] else NA_character_
#         tibble::tibble(
#           file = basename(f),
#           comp_group = comp_group,
#           cell_type = cell_type,
#           n_human_sig = NA_integer_, n_human_mapped = NA_integer_, n_mouse_sig = NA_integer_,
#           n_overlap = NA_integer_, jaccard = NA_real_,
#           n_concordant = NA_integer_, n_discordant = NA_integer_
#         )
#       }
#     )
#   })
#   
#   master_sel <- dplyr::bind_rows(all_summaries) %>%
#     dplyr::arrange(dplyr::desc(n_overlap))
#   
#   # Save selected-only master tables
#   save_df(master_sel, file.path(output_dir, "Master_Summary_selected.csv"))
#   
#   # Excel with per-cell-type sheets (selected only)
#   wb <- openxlsx::createWorkbook()
#   openxlsx::addWorksheet(wb, "Summary_selected")
#   openxlsx::writeData(wb, "Summary_selected", master_sel)
#   
#   cell_types <- sort(unique(master_sel$cell_type))
#   cell_types <- cell_types[!is.na(cell_types)]
#   for (ct in cell_types) {
#     sheet <- sanitize_sheet(ct)
#     openxlsx::addWorksheet(wb, sheet)
#     df_ct <- master_sel %>% dplyr::filter(cell_type == ct) %>%
#       dplyr::arrange(dplyr::desc(n_overlap))
#     openxlsx::writeData(wb, sheet, df_ct)
#   }
#   
#   xlsx_path <- file.path(output_dir, "Master_Summary_selected.xlsx")
#   openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
#   message("Selected-only Excel written: ", xlsx_path)
#   
#   invisible(master_sel)
# }
# 
# ################################################################################
# # DROP-IN: run selected-only passes for each folder
# ################################################################################
# invisible(mapply(function(in_dir, out_dir) {
#   compare_folder_selected(in_dir, out_dir, human_sig_mapped)
# }, mouse_deg_dirs, output_roots))
# 
