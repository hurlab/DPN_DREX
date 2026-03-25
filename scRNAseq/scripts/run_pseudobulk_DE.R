#' ============================================================================
#' Comprehensive Pseudobulk DE with Pathway Validation
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
#'   Builds pseudobulk matrices per cell type and major type from Seurat object, runs
#   DESeq2 differential expression, computes GSVA pathway scores, and performs
#   fGSEA on ranked gene lists for pathway-level validation.
#'
#' Input:
#'   steph_major.rdata (Seurat object), senescence markers CSV, GSVA/fGSEA
#   annotations
#'
#' Output:
#'   DESeq2 results, GSVA pathway scores, fGSEA enrichment, heatmaps
#'
#' Related figures:
#'   Pseudobulk DE used for Fig. S8, Tables S1-S4
#'
#' ============================================================================
#!/usr/bin/env Rscript
# Pseudobulk DE + pathway validation aligned to existing workflow (steph_major.rdata)
# Date: 2026-01-29

suppressPackageStartupMessages({
  pkgs <- c("Seurat", "Matrix", "DESeq2", "dplyr", "readr", "ggplot2", "tibble", "tidyr",
            "limma", "edgeR", "GSVA", "fgsea", "msigdbr", "patchwork", "pheatmap", "scales", "ashr")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  }
  lapply(pkgs, library, character.only = TRUE)
})

# ------------------------- Config -------------------------------------------
rdata_path   <- "steph_major.rdata"
assay_name   <- "RNA"          # use raw counts
sample_col   <- "orig.ident"   # per-sample id
group_col    <- "Group"        # biological group (SD/HFD/DR/EX/DREX)
celltype_col <- "celltype"     # fine cell type
major_col    <- "CellType"     # major cell type (defined in run_step_final.r)
min_cells_pb <- 20              # min cells per sample x grouping to keep (reuse prior code default)
min_genes    <- 10              # filter genes with <10 total counts
marker_csv   <- "../../../Senescence_JH/markgetSets/Senescence-Stephanie-revised.csv"
# contrasts: disease vs control; moderate vs severe not available in metadata (Group only)
contrasts <- list(
  HFD_vs_SD   = c("HFD", "SD"),
  DR_vs_HFD   = c("DR", "HFD"),
  EX_vs_HFD   = c("EX", "HFD"),
  DREX_vs_HFD = c("DREX", "HFD")
)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path("260129_pseudobulk_DE", ts)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- Load data ----------------------------------------
if (!file.exists(rdata_path)) stop("Cannot find ", rdata_path)
load(rdata_path)
if (!exists("sam")) stop("`sam` not found in ", rdata_path)
DefaultAssay(sam) <- assay_name
meta <- sam@meta.data
stopifnot(all(c(sample_col, group_col) %in% colnames(meta)))

# ------------------------- Helper: major labels -----------------------------
if (!major_col %in% colnames(meta)) {
  message("major_col not present; deriving from celltype (SC collapse etc.)")
  meta[[major_col]] <- meta[[celltype_col]]
  meta[[major_col]] <- ifelse(meta[[celltype_col]] %in% c("mySC","SC3","nmSC","ImmSC"), "SC", meta[[major_col]])
  meta[[major_col]] <- ifelse(meta[[celltype_col]] %in% c("EndoFib1","EndoFib2","EpiFib"), "Fib", meta[[major_col]])
}

# ------------------------- Pseudobulk builder -------------------------------
get_counts <- function(obj, assay) {
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = "counts")
  )
}

build_pb <- function(group_by_col) {
  counts <- get_counts(sam, assay_name)
  key <- paste(meta[[sample_col]], meta[[group_by_col]], sep = "|")
  col_sets <- split(seq_len(ncol(counts)), key)
  cell_counts <- vapply(col_sets, length, integer(1))
  keep <- names(col_sets)[cell_counts >= min_cells_pb]
  col_sets <- col_sets[keep]
  if (!length(col_sets)) stop("No pseudobulk samples after filtering for ", group_by_col)

  pb_mat <- do.call(cbind, lapply(col_sets, function(idx) Matrix::rowSums(counts[, idx, drop = FALSE])))
  colnames(pb_mat) <- names(col_sets)
  pb_coldata <- do.call(rbind, lapply(names(col_sets), function(nm) {
    parts <- strsplit(nm, "\\|" )[[1]]
    data.frame(sample = parts[1], level = parts[2], stringsAsFactors = FALSE)
  }))
  pb_coldata$group <- meta[[group_col]][match(pb_coldata$sample, meta[[sample_col]])]
  pb_coldata$cells <- cell_counts[match(paste(pb_coldata$sample, pb_coldata$level, sep="|"), names(cell_counts))]

  # gene filter
  keep_genes <- Matrix::rowSums(pb_mat) >= min_genes
  pb_mat <- pb_mat[keep_genes, , drop = FALSE]
  list(mat = pb_mat, coldata = pb_coldata)
}

# per-sample (all cells together)
counts_all <- get_counts(sam, assay_name)
pb_all_cols <- split(seq_len(ncol(counts_all)), meta[[sample_col]])
pb_all_counts <- vapply(pb_all_cols, length, integer(1))
keep_all <- names(pb_all_cols)[pb_all_counts >= min_cells_pb]
pb_all_mat <- do.call(cbind, lapply(pb_all_cols[keep_all], function(idx) Matrix::rowSums(counts_all[, idx, drop = FALSE])))
colnames(pb_all_mat) <- paste(keep_all, "all", sep = "|")
pb_all_coldata <- data.frame(sample = keep_all,
                             group = meta[[group_col]][match(keep_all, meta[[sample_col]])],
                             cells = pb_all_counts[keep_all], level = "all", stringsAsFactors = FALSE)
keep_genes_all <- Matrix::rowSums(pb_all_mat) >= min_genes
pb_all_mat <- pb_all_mat[keep_genes_all, , drop = FALSE]
pb_all <- list(mat = pb_all_mat, coldata = pb_all_coldata)

# per-celltype and per-major
pb_cell  <- build_pb(group_by_col = celltype_col)
pb_major <- build_pb(group_by_col = major_col)

saveRDS(pb_all,   file.path(out_dir, "pseudobulk_all.rds"))
saveRDS(pb_cell,  file.path(out_dir, "pseudobulk_celltype.rds"))
saveRDS(pb_major, file.path(out_dir, "pseudobulk_major.rds"))

# summary stats
summary_cells <- meta %>% count(.data[[sample_col]]) %>% rename(cells = n)
summary_groups <- meta %>% count(.data[[group_col]])
write_csv(summary_cells, file.path(out_dir, "cells_per_sample.csv"))
write_csv(summary_groups, file.path(out_dir, "cells_per_group.csv"))

# ------------------------- DESeq2 runner ------------------------------------
run_deseq <- function(pb_obj, level_label) {
  mat <- pb_obj$mat; cd <- pb_obj$coldata
  cd$sample_level <- paste(cd$sample, cd$level, sep = "|")
  rownames(cd) <- cd$sample_level
  mat <- mat[, cd$sample_level, drop = FALSE]

  vst_mat <- NULL
  de_list <- list(); res_summ <- list()

  for (cn in names(contrasts)) {
    gpair <- contrasts[[cn]]
    cd_sub <- cd %>% filter(group %in% gpair)
    if (n_distinct(cd_sub$group) < 2) next
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat[, cd_sub$sample_level, drop = FALSE],
                                          colData   = cd_sub,
                                          design    = ~ group)
    dds <- dds[rowSums(counts(dds)) >= min_genes, ]
    if (nrow(dds) < 10) next
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    vst_obj <- DESeq2::vst(dds, blind = TRUE)
    vst_mat <- SummarizedExperiment::assay(vst_obj)
    res <- DESeq2::results(dds, contrast = c("group", gpair[1], gpair[2]))
    res <- DESeq2::lfcShrink(dds, contrast = c("group", gpair[1], gpair[2]), type = "ashr")
    res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene")
    res_df$contrast <- cn; res_df$level <- level_label
    de_list[[cn]] <- res_df

    # MA plot (ggplot)
    ma_df <- res_df %>% mutate(sig = ifelse(!is.na(padj) & padj < 0.05, "sig", "ns"))
    ma_p <- ggplot(ma_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = sig)) +
      geom_point(alpha = 0.5, size = 1) +
      scale_color_manual(values = c(ns = "grey70", sig = "firebrick")) +
      labs(title = sprintf("MA %s %s", level_label, cn), x = "log10 baseMean", y = "log2FC") +
      theme_bw(10)
    ggsave(file.path(out_dir, sprintf("MA_%s_%s.pdf", level_label, cn)), ma_p, width = 6, height = 5)

    # volcano
    volc <- res_df %>% mutate(sig = ifelse(!is.na(padj) & padj < 0.05, "sig", "ns"))
    volc_p <- ggplot(volc, aes(log2FoldChange, -log10(pvalue), color = sig)) +
      geom_point(alpha = 0.6, size = 1.2) +
      scale_color_manual(values = c(ns = "grey70", sig = "firebrick")) +
      labs(title = sprintf("Volcano %s %s", level_label, cn)) + theme_bw(10)
    ggsave(file.path(out_dir, sprintf("volcano_%s_%s.pdf", level_label, cn)), volc_p, width = 5, height = 4)

    res_summ[[cn]] <- data.frame(level = level_label, contrast = cn,
                                 n = nrow(res_df),
                                 sig_padj_0.05 = sum(!is.na(res_df$padj) & res_df$padj < 0.05),
                                 sig_p_0.05    = sum(!is.na(res_df$pvalue) & res_df$pvalue < 0.05),
                                 sig_p_0.01    = sum(!is.na(res_df$pvalue) & res_df$pvalue < 0.01))
  }
  list(de = bind_rows(de_list), vst = vst_mat, summary = bind_rows(res_summ))
}

res_all   <- run_deseq(pb_all,   "all")
res_major <- run_deseq(pb_major, "major")
res_cell  <- run_deseq(pb_cell,  "celltype")

# save DE results
write_csv(res_all$de,   file.path(out_dir, "deseq2_all.csv"))
if (!is.null(res_major$de)) write_csv(res_major$de, file.path(out_dir, "deseq2_major.csv"))
if (!is.null(res_cell$de)) write_csv(res_cell$de,  file.path(out_dir, "deseq2_celltype.csv"))
write_csv(bind_rows(res_all$summary, res_major$summary, res_cell$summary),
          file.path(out_dir, "deseq2_summary.csv"))

saveRDS(res_all,   file.path(out_dir, "deseq2_all.rds"))
saveRDS(res_major, file.path(out_dir, "deseq2_major.rds"))
saveRDS(res_cell,  file.path(out_dir, "deseq2_celltype.rds"))

# PCA of samples (vst from all-level)
if (!is.null(res_all$vst)) {
  pca <- prcomp(t(res_all$vst))
  pca_df <- data.frame(pca$x[,1:2], sample = rownames(pca$x))
  pca_df$Group <- meta[[group_col]][match(sub("|.*", "", pca_df$sample), meta[[sample_col]])]
  p <- ggplot(pca_df, aes(PC1, PC2, color = Group, label = sample)) + geom_point(size=2) + theme_bw(10)
  ggsave(file.path(out_dir, "PCA_pseudobulk_all.pdf"), p, width = 6, height = 5)
}

# Heatmap of top DE genes (all-level first contrast)
if (!is.null(res_all$de) && nrow(res_all$de)) {
  topgenes <- res_all$de %>% arrange(padj) %>% filter(!is.na(padj)) %>% slice_head(n = 50) %>% pull(gene)
  hm_mat <- res_all$vst[topgenes, , drop = FALSE]
  pheatmap::pheatmap(hm_mat, scale = "row", show_rownames = TRUE,
                     filename = file.path(out_dir, "heatmap_top50_all.pdf"))
}

# ------------------------- Marker set scores --------------------------------
marker_sets <- list()
if (file.exists(marker_csv)) {
  markers_raw <- readr::read_csv(marker_csv, col_types = readr::cols())
  if ("Source" %in% colnames(markers_raw)) {
    sources <- unique(markers_raw$Source)
    for (s in sources) {
      genes <- unique(na.omit(as.character(markers_raw$Mouse[markers_raw$Source == s])))
      marker_sets[[s]] <- genes
    }
    marker_sets[["ALL"]] <- unique(na.omit(as.character(markers_raw$Mouse)))
  }
}
if (length(marker_sets) == 0) warning("Marker sets not found; skipping scores.")

if (length(marker_sets) > 0) {
  logcpm <- edgeR::cpm(pb_all$mat, log = TRUE, prior.count = 1)
  scores_mat <- tryCatch(
    {
      GSVA::gsva(expr = as.matrix(logcpm),
                 gset.idx.list = marker_sets,
                 method = "ssgsea",
                 ssgsea.norm = TRUE,
                 parallel.sz = 1)
    },
    error = function(e) {
      message("GSVA failed (", e$message, "); using mean logCPM per set instead.")
      sapply(marker_sets, function(gs) {
        genes <- intersect(gs, rownames(logcpm))
        if (!length(genes)) return(rep(NA_real_, ncol(logcpm)))
        colMeans(logcpm[genes, , drop = FALSE])
      })
    }
  )
  scores <- t(scores_mat) %>% as.data.frame() %>% tibble::rownames_to_column("sample_level")
  scores$sample <- sub("|.*", "", scores$sample_level)
  scores$group  <- meta[[group_col]][match(scores$sample, meta[[sample_col]])]
  write_csv(scores, file.path(out_dir, "marker_scores_by_sample.csv"))

  scores_long <- scores %>% pivot_longer(-c(sample_level, sample, group), names_to = "set", values_to = "score")
  p_box <- ggplot(scores_long, aes(x = group, y = score, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) + facet_wrap(~ set, scales = "free_y") +
    theme_bw(10) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(file.path(out_dir, "marker_scores_boxplots.pdf"), p_box, width = 10, height = 6)

  score_tests <- scores_long %>% group_by(set) %>%
    summarise(p_anova = if (dplyr::n_distinct(group) < 2) NA_real_ else summary(aov(score ~ group))[["Pr(>F)"]][1,1])
  write_csv(score_tests, file.path(out_dir, "marker_scores_anova.csv"))
}

# ------------------------- Pathway enrichment (fgsea Hallmark) --------------
pathways_h <- NULL
try({ pathways_h <- msigdbr::msigdbr(species = "Mus musculus", category = "H") %>%
        split(.$gene_symbol, .$gs_name) }, silent = TRUE)

run_fgsea <- function(de_df, label) {
  if (is.null(pathways_h) || nrow(de_df) == 0) return(NULL)
  if (!"stat" %in% colnames(de_df) || all(is.na(de_df$stat))) {
    de_df$stat <- de_df$log2FoldChange / de_df$lfcSE
  }
  de_df <- de_df %>% filter(!is.na(stat))
  if (nrow(de_df) == 0) return(NULL)
  ranks <- de_df$stat; names(ranks) <- de_df$gene
  fg <- fgsea::fgsea(pathways = pathways_h, stats = ranks, eps = 1e-6)
  fg$contrast <- unique(de_df$contrast)
  fg$level <- label
  fg
}

fg_list <- list()
for (lv in c("all", "major")) {
  de_df <- if (lv == "all") res_all$de else res_major$de
  if (!is.null(de_df) && nrow(de_df)) {
    for (cn in unique(de_df$contrast)) {
      fg <- run_fgsea(de_df %>% filter(contrast == cn), lv)
      if (!is.null(fg)) fg_list[[paste(lv, cn, sep = "_")]] <- fg
    }
  }
}
fg_all <- bind_rows(fg_list)
if (nrow(fg_all) > 0) {
  write_csv(fg_all, file.path(out_dir, "fgsea_results.csv"))
}

# ------------------------- Save RDS -----------------------------------------
saveRDS(list(pb_all = pb_all, pb_cell = pb_cell, pb_major = pb_major), file.path(out_dir, "pseudobulk_objects.rds"))
saveRDS(fg_all, file.path(out_dir, "fgsea_results.rds"))

# ------------------------- Summary print ------------------------------------
summary_text <- list(
  paste0("Output dir: ", out_dir),
  paste0("Samples: ", dplyr::n_distinct(meta[[sample_col]]),
         " | Groups: ", paste(paste0(summary_groups[[group_col]], "=", summary_groups$n), collapse = "; ")),
  paste0("Cells per sample median: ", stats::median(summary_cells$cells)),
  paste0("DEGs (padj<0.05) all-level: ",
         if (!is.null(res_all$summary) && nrow(res_all$summary)) {
           paste(res_all$summary$contrast, res_all$summary$sig_padj_0.05, sep = ":", collapse = "; ")
         } else "none"),
  paste0("DEGs (p<0.05) all-level: ",
         if (!is.null(res_all$summary) && nrow(res_all$summary)) {
           paste(res_all$summary$contrast, res_all$summary$sig_p_0.05, sep = ":", collapse = "; ")
         } else "none"),
  paste0("DEGs (p<0.01) all-level: ",
         if (!is.null(res_all$summary) && nrow(res_all$summary)) {
           paste(res_all$summary$contrast, res_all$summary$sig_p_0.01, sep = ":", collapse = "; ")
         } else "none"),
  paste0("Top fgsea pathways (all): ",
         if (!is.null(fg_all) && nrow(fg_all)) {
           paste(head(fg_all %>% arrange(padj) %>% pull(pathway), 5), collapse = ", ")
         } else "none")
)
summary_text <- vapply(summary_text, as.character, FUN.VALUE = character(1))
writeLines(summary_text)
