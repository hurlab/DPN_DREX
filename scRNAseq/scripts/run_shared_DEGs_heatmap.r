#' ============================================================================
#' Shared DEGs Venn Diagrams and Heatmaps
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
#'   Identifies differentially expressed genes shared across pairwise contrasts. Creates
#   Venn diagrams of DEG overlaps and row-z-scored heatmaps of shared genes using
#   the Seurat object expression data.
#'
#' Input:
#'   steph_major.rdata, step_scwann.rdata (Seurat objects with DEG lists)
#'
#' Output:
#'   Venn diagrams, row-z-scored heatmaps of shared DEGs
#'
#' Related figures:
#'   Fig. S4
#'
#' ============================================================================
#####################  Set Current Working Directory ###########################
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(current_path)
base_dir <- dirname(current_path)

#####################  Packages  ###############################################
# Install pacman for package management if not already installed
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

# Required CRAN packages
cran_packages <- c(
  "readxl", "writexl", "dplyr", "ggplot2", "ggVennDiagram", "tibble"
)

# Required Bioconductor packages
bioc_packages <- c("VennDetail")

# Automatically install and load CRAN and Bioconductor packages
p_load(char = cran_packages, install = TRUE)
p_load(char = bioc_packages, install = TRUE, repos = BiocManager::repositories())

#####################  Define Output Directory #################################
output_dir <- file.path(base_dir, "250818_StackedHeatmaps/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#####################  Load Saved Data #########################################
if (file.exists("steph_major.rdata")) load("steph_major.rdata")
if (file.exists("step_scwann.rdata")) load("step_scwann.rdata")


# =========================================================
# Helper 1: Venn + shared-DEG heatmap (+ optional dotplot)
# =========================================================
sc_celltype_shared_venn_heatmap <- function(
    sam,                     # Seurat object
    deg,                     # list of contrasts -> per-celltype DE tables
    celltype    = "ImmSC",   # "ImmSC", "mySC", "nmSC"
    p_cutoff    = 0.01,
    group_order = c("SD","HFD","DR","EX","DREX"),
    output_dir  = ".",
    label       = "ImmSC",   # <-- default label now ImmSC
    contrasts_for_venn = names(deg)[1:4],
    make_dotplot = TRUE
) {
  pacman::p_load(Seurat, VennDetail, pheatmap, ggplot2, grid)
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Ensure "cond" exists (Group_Celltype) for AverageExpression
  if (!"cond" %in% colnames(sam@meta.data)) {
    sam$cond <- paste0(sam$Group, "_", sam$celltype)
  }
  
  # Collect DEG rownames for the requested celltype
  deg_cell <- lapply(deg, function(x) {
    if (!is.null(x[[celltype]])) subset(x[[celltype]], p_val < p_cutoff) else data.frame()
  })
  
  # Venn over selected contrasts
  ven <- VennDetail::venndetail(lapply(deg_cell[contrasts_for_venn], function(x) rownames(x)))
  pdf(file.path(output_dir, paste0(label, "_", celltype, "_shared_deg_venn_p", sub("\\.","", p_cutoff), ".pdf")))
  plot(ven)
  dev.off()
  
  # Shared genes
  shared <- VennDetail::getSet(ven, "Shared")$Detail
  
  # Average expression matrix and subset to this celltype’s group columns
  avg  <- Seurat::AverageExpression(sam, assays = "RNA", group.by = "condition")
  avgd <- as.data.frame(avg$RNA)
  want_cols <- paste0(group_order, "-", celltype)
  have_cols <- intersect(want_cols, colnames(avgd))
  if (length(have_cols) == 0) stop("No columns matching *_", celltype, " in AverageExpression output.")
  
  mat <- avgd[intersect(shared, rownames(avgd)), have_cols, drop = FALSE]
  colnames(mat) <- have_cols
  
  # Heatmap (same styling as your script)
  if (nrow(mat) >= 2) {
    hm <- pheatmap::pheatmap(
      mat, scale = "row", cluster_col = FALSE,
      color = colorRampPalette(c("darkblue", "white", "red"))(128),
      border_color = "white", main = paste0(celltype, " shared DEGs"),
      silent = TRUE
    )
  } else {
    hm <- pheatmap::pheatmap(
      mat, scale = "row", cluster_col = FALSE, cluster_row = FALSE, # disable row clustering
      color = colorRampPalette(c("darkblue", "white", "red"))(128),
      border_color = "white", main = paste0(celltype, " shared DEGs"),
      silent = TRUE
    )
  }
  
  pdf(file.path(output_dir, paste0(label, "_", celltype, "_shared_deg_heatmap.pdf")))
  grid::grid.newpage(); grid::grid.draw(hm$gtable)
  dev.off()
  
  # Optional dotplot for this subtype
  dot_pdf <- NA_character_
  if (make_dotplot && length(shared) > 0) {
    Idents(sam) <- "celltype"
    dp <- Seurat::DotPlot(subset(sam, celltype == !!celltype), features = shared) +
      coord_flip() +
      scale_color_viridis_c() +
      theme_light(base_size = 15) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      labs(x = NULL, y = NULL, title = paste0(celltype, " shared DEGs"))
    dot_pdf <- file.path(output_dir, paste0(label, "_", celltype, "_shared_deg_dotplot.pdf"))
    ggsave(dot_pdf, plot = dp, width = 6.5, height = 8)
  }
  
  invisible(list(
    venn        = ven,
    shared      = shared,
    mat         = mat,
    heatmap     = hm,
    heatmap_pdf = file.path(output_dir, paste0(label, "_", celltype, "_shared_deg_heatmap.pdf")),
    venn_pdf    = file.path(output_dir, paste0(label, "_", celltype, "_shared_deg_venn_p", sub("\\.","", p_cutoff), ".pdf")),
    dot_pdf     = dot_pdf
  ))
}

# =========================================================
# Helper 2: stack three pheatmap objects vertically
# =========================================================
stack_sc_heatmaps <- function(
    hm_list,                   # list(ImmSC_hm, mySC_hm, nmSC_hm)
    labels     = c("ImmSC","mySC","nmSC"),
    output_dir = ".",
    label      = "SC",         # prefix for combined file
    rel_heights = c(1,1,1)
) {
  pacman::p_load(ggplotify, cowplot)
  stopifnot(length(hm_list) == 3)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  gg_list <- lapply(seq_along(hm_list), function(i) {
    g <- ggplotify::as.ggplot(hm_list[[i]]$gtable)
    cowplot::ggdraw() +
      cowplot::draw_label(labels[i], x = 0.02, y = 0.99, hjust = 0, vjust = 1, size = 12) +
      cowplot::draw_plot(g, x = 0, y = 0, width = 1, height = 1)
  })
  
  out_pdf <- file.path(output_dir, paste0(label, "_SC_shared_deg_heatmaps_stacked.pdf"))
  pdf(out_pdf, width = 8.5, height = 9)
  print(cowplot::plot_grid(plotlist = gg_list, ncol = 1, align = "v", rel_heights = rel_heights))
  dev.off()
  invisible(out_pdf)
}


# =========================================================
immSC_out <- sc_celltype_shared_venn_heatmap(
  sam        = sam,
  deg        = deg,
  celltype   = "ImmSC",
  p_cutoff   = 0.01,
  output_dir = output_dir,
  label      = "ImmSC"  
)

mySC_out <- sc_celltype_shared_venn_heatmap(
  sam        = sam,
  deg        = deg,
  celltype   = "mySC",
  p_cutoff   = 0.01,
  output_dir = output_dir,
  label      = "mySC"
)

nmSC_out <- sc_celltype_shared_venn_heatmap(
  sam        = sam,
  deg        = deg,
  celltype   = "nmSC",
  p_cutoff   = 0.01,
  output_dir = output_dir,
  label      = "nmSC"
)

# Combine the three heatmaps (ImmSC, mySC, nmSC)
stack_sc_heatmaps(
  hm_list    = list(immSC_out$heatmap, mySC_out$heatmap, nmSC_out$heatmap),
  labels     = c("ImmSC shared DEGs", "mySC shared DEGs", "nmSC shared DEGs"),
  output_dir = output_dir,
  label      = "SC"
)










#####################  Single combined heatmap builder #########################
build_single_sc_heatmap <- function(
    outs,                       
    sets         = c("ImmSC","mySC","nmSC"),
    group_order  = c("SD","HFD","DR","EX","DREX"),
    output_dir   = ".",
    filename_pdf = "SC_allinone_shared_deg_heatmap.pdf",
    filename_csv = "SC_allinone_shared_deg_matrix.csv",
    within_set_cluster = TRUE
) {
  pacman::p_load(pheatmap, grDevices)
  
  # 1) Collect matrices per set, strip cell type suffixes
  mats <- lapply(sets, function(s) {
    m <- outs[[s]]$mat
    if (is.null(m) || nrow(m) == 0) return(NULL)
    # Remove suffix "-ImmSC" / "-mySC" / "-nmSC" from column names
    clean_names <- sub(paste0("-", s, "$"), "", colnames(m))
    colnames(m) <- clean_names
    # Reorder columns to match group_order
    m <- m[, intersect(group_order, clean_names), drop = FALSE]
    m
  })
  names(mats) <- sets
  mats <- mats[!vapply(mats, is.null, logical(1))]
  
  if (length(mats) == 0) stop("No matrices to combine. Nothing to plot.")
  
  # 2) Cluster rows within each block if requested
  mats_ordered <- lapply(mats, function(m) {
    if (nrow(m) >= 2 && within_set_cluster) {
      o <- tryCatch(hclust(dist(m))$order, error = function(e) seq_len(nrow(m)))
      m[o, , drop = FALSE]
    } else {
      m
    }
  })
  
  # 3) Bind blocks row-wise, now columns match
  comb_mat <- do.call(rbind, mats_ordered)
  
  # 4) Annotation and gaps
  set_sizes <- vapply(mats_ordered, nrow, integer(1))
  gaps_row  <- cumsum(head(set_sizes, -1))
  gaps_col  <- NULL  # no column gaps since now columns are common groups
  
  annot_row <- data.frame(
    Set = factor(rep(names(mats_ordered), times = set_sizes), levels = sets)
  )
  rownames(annot_row) <- rownames(comb_mat)
  
  annot_col <- data.frame(
    Group = factor(colnames(comb_mat), levels = group_order)
  )
  rownames(annot_col) <- colnames(comb_mat)
  
  # 5) Save combined matrix
  write.csv(comb_mat, file.path(output_dir, filename_csv), row.names = TRUE)
  
  # 6) Heatmap with single global scale
  pdf(file.path(output_dir, filename_pdf), width = 10, height = 10)
  hm <- pheatmap::pheatmap(
    comb_mat,
    scale          = "row",
    cluster_row    = FALSE,
    cluster_col    = FALSE,
    gaps_row       = if (length(gaps_row)) gaps_row else NULL,
    annotation_row = annot_row,
    annotation_col = annot_col,
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    color          = colorRampPalette(c("darkblue","white","red"))(128),
    border_color   = "white",
    main           = "Shared DEGs across SC sets (single heatmap)"
  )
  dev.off()
  
  invisible(list(
    combined_matrix = comb_mat,
    gaps_row        = gaps_row,
    annotation_row  = annot_row,
    annotation_col  = annot_col,
    pdf             = file.path(output_dir, filename_pdf),
    csv             = file.path(output_dir, filename_csv)
  ))
}




combined_out <- build_single_sc_heatmap(
  outs = list(ImmSC = immSC_out, mySC = mySC_out, nmSC = nmSC_out),
  sets = c("ImmSC", "mySC", "nmSC"),
  group_order = c("SD", "HFD", "DR", "EX", "DREX"),
  output_dir = output_dir,
  filename_pdf = "SC_allinone_shared_deg_heatmap_v2.pdf",
  filename_csv = "SC_allinone_shared_deg_matrix_v2.csv",
  within_set_cluster = TRUE
)


