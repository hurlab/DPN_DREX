#' ============================================================================
#' Gene Set Variation Analysis (GSVA)
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
#'   Computes GSVA pathway activity scores per cell using KEGG and MSigDB HALLMARK
#   gene sets. Tests for significant pathway differences across experimental groups
#   using Wilcoxon rank-sum tests with Benjamini-Hochberg correction.
#'
#' Input:
#'   sc.rdata (Seurat object), KEGG mouse annotations, MSigDB HALLMARK gene sets
#'
#' Output:
#'   GSVA scores embedded in Seurat object, significant pathway lists
#   (padj < 0.05), violin plots, cell-type-specific heatmaps
#'
#' Related figures:
#'   Fig. S5-S7, Tables S5-S7
#'
#' ============================================================================
# ==============================================================================
# Gene Set Variation Analysis (GSVA) for Pathway Enrichment
# Methods: "Pathway enrichment scores were computed using gene set variation 
#          analysis (GSVA)"
# ==============================================================================

library(scGSVA)
library(Seurat)

# ------------------------------------------------------------------------------
# 1. Load Single-cell Data
# ------------------------------------------------------------------------------
load("sc.rdata")
head(sam@meta.data)

# ------------------------------------------------------------------------------
# 2. Build KEGG Pathway Annotations for Mouse
# Methods: "A custom annotation file was used to map metabolites to pathways"
# Methods: "GSVA scores were calculated for each dataset using the scgsva 
#          function (kcdf = 'Gaussian', min.sz = 2)"
# ------------------------------------------------------------------------------
mmko <- buildAnnot(species="mouse", keytype="SYMBOL", anntype="KEGG", builtin=F)

# ------------------------------------------------------------------------------
# 3. Run GSVA Analysis
# Methods: "GSVA scores were calculated for each dataset"
# Using batch processing for computational efficiency
# cores=20: parallel processing with 20 cores
# ------------------------------------------------------------------------------
res <- scgsva(sam, mmko, batch=4000, cores=20)

# ------------------------------------------------------------------------------
# 4. Identify Significantly Altered Pathways Between Groups
# Methods: "Group comparisons were performed using the Wilcoxon rank-sum 
#          applied to GSVA scores, and p-values were adjusted for multiple 
#          comparisons via the Benjamini-Hochberg correction"
# Methods: "Pathways with an adjusted p-value<0.05 were considered 
#          significantly different between groups"
# ------------------------------------------------------------------------------

# HFD vs SD (Standard Diet) comparison
hfvssd <- sigPathway(subset(res, Group %in% c("HFD", "SD")),
                     group="Group", ref="SD")

# DR (Diet Restriction) vs HFD comparison
drvshfd <- sigPathway(subset(res, Group %in% c("HFD", "DR")),
                      group="Group", ref="HFD")

# EX (Exercise) vs HFD comparison
exvshfd <- sigPathway(subset(res, Group %in% c("EX", "HFD")),
                      group="Group", ref="HFD")

# DREX (Diet Restriction + Exercise) vs HFD comparison
drexvshfd <- sigPathway(subset(res, Group %in% c("DREX", "HFD")),
                        group="Group", ref="HFD")

# DR vs EX comparison
drvsex <- sigPathway(subset(res, Group %in% c("DR", "EX")),
                     group="Group", ref="EX")

# DR vs DREX comparison
drvsdrex <- sigPathway(subset(res, Group %in% c("DR", "DREX")),
                       group="Group", ref="DREX")

# EX vs DREX comparison
exvsdrex <- sigPathway(subset(res, Group %in% c("EX", "DREX")),
                       group="Group", ref="DREX")

# Find pathways specific to HFD group
hfd <- findPathway(res, group="Group", ref="HFD")

# ------------------------------------------------------------------------------
# 5. GSVA Analysis with MSigDB HALLMARK Gene Sets
# Use hallmark pathways for biological process interpretation
# ------------------------------------------------------------------------------
mmmg <- buildMSIGDB(species="mouse",
                    keytype="SYMBOL",
                    anntype="HALLMARK")

# Run GSVA with HALLMARK gene sets
resmg <- scgsva(sam, mmmg, batch=4000, cores=20)

# ------------------------------------------------------------------------------
# 6. Identify Significant Pathways with HALLMARK Gene Sets
# Same comparisons as above, but using HALLMARK pathways
# ------------------------------------------------------------------------------
hfvssdmg <- sigPathway(subset(resmg, Group %in% c("HFD", "SD")),
                       group="Group", ref="SD")
drvshfdmg <- sigPathway(subset(resmg, Group %in% c("HFD", "DR")),
                        group="Group", ref="HFD")
exvshfdmg <- sigPathway(subset(resmg, Group %in% c("EX", "HFD")),
                        group="Group", ref="HFD")
drexvshfdmg <- sigPathway(subset(resmg, Group %in% c("DREX", "HFD")),
                          group="Group", ref="HFD")
drvsexmg <- sigPathway(subset(resmg, Group %in% c("DR", "EX")),
                       group="Group", ref="EX")
drvsdrexmg <- sigPathway(subset(resmg, Group %in% c("DR", "DREX")),
                         group="Group", ref="DREX")
exvsdrexmg <- sigPathway(subset(resmg, Group %in% c("EX", "DREX")),
                         group="Group", ref="DREX")

# ------------------------------------------------------------------------------
# 7. Organize Results and Sort by Significance
# ------------------------------------------------------------------------------

# HALLMARK pathway results
mg <- list(HFDvsSD=hfvssdmg, DRvsHFD=drvshfdmg, EXvsHFD=exvshfdmg, 
           DREXvsHFD=drexvshfdmg, DRvsEX=drvsexmg, DRvsDREX=drvsdrexmg,
           EXvsDREX=exvsdrexmg)
mgs <- lapply(mg, function(x) x[order(x$p), ])

# KEGG pathway results
kg <- list(HFDvsSD=hfvssd, DRvsHFD=drvshfd, EXvsHFD=exvshfd, 
           DREXvsHFD=drexvshfd, DRvsEX=drvsex, DRvsDREX=drvsdrex,
           EXvsDREX=exvsdrex)
kgs <- lapply(kg, function(x) x[order(x$p), ])

# ------------------------------------------------------------------------------
# 8. Cell Type-Specific Pathway Analysis
# Identify significant pathways for each cell type separately
# ------------------------------------------------------------------------------

# HFD vs SD for each cell type
hfvssdc <- lapply(unique(sam$celltype),
                  function(x) sigPathway(
                    subset(res, condition %in% c(paste0(c("HFD_", "SD_"), x))),
                    group="condition"))
names(hfvssdc) <- unique(sam$celltype)

# DR vs HFD for each cell type
drvshfdc <- lapply(unique(sam$celltype),
                   function(x) sigPathway(
                     subset(res, condition %in% c(paste0(c("DR_", "HFD_"), x))),
                     group="condition"))
names(drvshfdc) <- unique(sam$celltype)

# EX vs HFD for each cell type
exvshfdc <- lapply(unique(sam$celltype),
                   function(x) sigPathway(
                     subset(res, condition %in% c(paste0(c("EX_", "HFD_"), x))),
                     group="condition"))
names(exvshfdc) <- unique(sam$celltype)

# DREX vs HFD for each cell type
drexvshfdc <- lapply(unique(sam$celltype),
                     function(x) sigPathway(
                       subset(res, condition %in% c(paste0(c("DREX_", "HFD_"), x))),
                       group="condition"))
names(drexvshfdc) <- unique(sam$celltype)

# ------------------------------------------------------------------------------
# 9. Cell Type-Specific HALLMARK Pathway Analysis
# ------------------------------------------------------------------------------
hfvssdcmg <- lapply(unique(sam$celltype),
                    function(x) sigPathway(
                      subset(resmg, condition %in% c(paste0(c("HFD_", "SD_"), x))),
                      group="condition"))
names(hfvssdcmg) <- unique(sam$celltype)

drvshfdcmg <- lapply(unique(sam$celltype),
                     function(x) sigPathway(
                       subset(resmg, condition %in% c(paste0(c("DR_", "HFD_"), x))),
                       group="condition"))
names(drvshfdcmg) <- unique(sam$celltype)

exvshfdcmg <- lapply(unique(sam$celltype),
                     function(x) sigPathway(
                       subset(resmg, condition %in% c(paste0(c("EX_", "HFD_"), x))),
                       group="condition"))
names(exvshfdcmg) <- unique(sam$celltype)

drexvshfdcmg <- lapply(unique(sam$celltype),
                       function(x) sigPathway(
                         subset(resmg, condition %in% c(paste0(c("DREX_", "HFD_"), x))),
                         group="condition"))
names(drexvshfdcmg) <- unique(sam$celltype)

# ------------------------------------------------------------------------------
# 10. Filter and Sort Cell Type-Specific Results
# Select significant pathways (p.adj < 0.05) and sort by p-value
# ------------------------------------------------------------------------------

# HALLMARK pathways
hfvssdcmgs <- lapply(hfvssdcmg, function(x){
  tmp <- subset(x, p.adj < 0.05)
  tmp <- tmp[order(tmp$p), ]
  return(tmp)
})

drvshfdcmgs <- lapply(drvshfdcmg, function(x){
  tmp <- subset(x, p.adj < 0.05)
  tmp <- tmp[order(tmp$p), ]
  return(tmp)
})

exvshfdcmgs <- lapply(exvshfdcmg, function(x){
  tmp <- subset(x, p.adj < 0.05)
  tmp <- tmp[order(tmp$p), ]
  return(tmp)
})

drexvshfdcmgs <- lapply(drexvshfdcmg, function(x){
  tmp <- subset(x, p.adj < 0.05)
  tmp <- tmp[order(tmp$p), ]
  return(tmp)
})

# KEGG pathways
hfvssdcs <- lapply(hfvssdc, function(x){
  tmp <- subset(x, p.adj < 0.05)
  tmp <- tmp[order(tmp$p), ]
  return(tmp)
})

drvshfdcs <- lapply(drvshfdc, function(x){
  tmp <- subset(x, p.adj < 0.05)
  tmp <- tmp[order(tmp$p), ]
  return(tmp)
})

exvshfdcs <- lapply(exvshfdc, function(x){
  tmp <- subset(x, p.adj < 0.05)
  tmp <- tmp[order(tmp$p), ]
  return(tmp)
})

drexvshfdcs <- lapply(drexvshfdc, function(x){
  tmp <- subset(x, p.adj < 0.05)
  tmp <- tmp[order(tmp$p), ]
  return(tmp)
})

# ------------------------------------------------------------------------------
# 11. Visualization of GSVA Results
# Create violin plots and heatmaps for pathway scores
# ------------------------------------------------------------------------------

# Create output directory
dir.create("GSVA")
setwd("GSVA")

# Violin plots for overall KEGG pathway comparisons (top 12 pathways)
for(i in names(kgs)){
  nam <- unlist(strsplit(i, "vs"))
  vlnPlot(subset(res, Group %in% nam), 
          features=kgs[[i]]$Path[1:12], 
          group_by="Group")
  ggsave(file=paste0(i, "_overall_gsva_top12.pdf"), width=9, height=6)
}

# Violin plots for overall HALLMARK pathway comparisons (top 12 pathways)
for(i in names(mgs)){
  nam <- unlist(strsplit(i, "vs"))
  vlnPlot(subset(resmg, Group %in% nam), 
          features=mgs[[i]]$Path[1:12], 
          group_by="Group")
  ggsave(file=paste0(i, "_overall_gsva_hallmark_top12.pdf"), width=9, height=6)
}

# ------------------------------------------------------------------------------
# 12. Heatmaps for Cell Type-Specific Pathway Scores
# Show top 10 pathways for each cell type
# ------------------------------------------------------------------------------

# KEGG pathways - HFD vs SD
Heatmap(subset(res, Group %in% c("HFD", "SD")),
        features=na.omit(unique(unlist(lapply(hfvssdcs, 
                                              function(x) x$Path[1:10])))),
        show_colnames=F,
        group_by=c("Group", "celltype"),
        cluster_cols=F,
        order_by="celltype",
        fontsize_row=8)
dev.print(pdf, file="HFDvsSD_celltype_heatmap_top10.pdf")

# KEGG pathways - DR vs HFD
Heatmap(subset(res, Group %in% c("DR", "HFD")),
        features=na.omit(unique(unlist(lapply(drvshfdcs, 
                                              function(x) x$Path[1:10])))),
        show_colnames=F,
        group_by=c("Group", "celltype"),
        cluster_cols=F,
        order_by="celltype",
        fontsize_row=8)
dev.print(pdf, file="DRvsHFD_celltype_heatmap_top10.pdf")

# KEGG pathways - EX vs HFD
Heatmap(subset(res, Group %in% c("EX", "HFD")),
        features=na.omit(unique(unlist(lapply(exvshfdcs, 
                                              function(x) x$Path[1:10])))),
        show_colnames=F,
        group_by=c("Group", "celltype"),
        cluster_cols=F,
        order_by="celltype",
        fontsize_row=8)
dev.print(pdf, file="EXvsHFD_celltype_heatmap_top10.pdf")

# KEGG pathways - DREX vs HFD
Heatmap(subset(res, Group %in% c("DREX", "HFD")),
        features=na.omit(unique(unlist(lapply(drexvshfdcs, 
                                              function(x) x$Path[1:10])))),
        show_colnames=F,
        group_by=c("Group", "celltype"),
        cluster_cols=F,
        order_by="celltype",
        fontsize_row=8)
dev.print(pdf, file="DREXvsHFD_celltype_heatmap_top10.pdf")

# HALLMARK pathways - HFD vs SD
Heatmap(subset(resmg, Group %in% c("HFD", "SD")),
        features=na.omit(unique(unlist(lapply(hfvssdcmgs, 
                                              function(x) x$Path[1:10])))),
        show_colnames=F,
        group_by=c("Group", "celltype"),
        cluster_cols=F,
        order_by="celltype",
        fontsize_row=8)
dev.print(pdf, file="HFDvsSD_celltype_heatmap_hallmark_top10.pdf")

# HALLMARK pathways - DR vs HFD
Heatmap(subset(resmg, Group %in% c("DR", "HFD")),
        features=na.omit(unique(unlist(lapply(drvshfdcmgs, 
                                              function(x) x$Path[1:10])))),
        show_colnames=F,
        group_by=c("Group", "celltype"),
        cluster_cols=F,
        order_by="celltype",
        fontsize_row=8)
dev.print(pdf, file="DRvsHFD_celltype_heatmap_hallmark_top10.pdf")

# HALLMARK pathways - EX vs HFD
Heatmap(subset(resmg, Group %in% c("EX", "HFD")),
        features=na.omit(unique(unlist(lapply(exvshfdcmgs, 
                                              function(x) x$Path[1:10])))),
        show_colnames=F,
        group_by=c("Group", "celltype"),
        cluster_cols=F,
        order_by="celltype",
        fontsize_row=8)
dev.print(pdf, file="EXvsHFD_celltype_heatmap_hallmark_top10.pdf")

# HALLMARK pathways - DREX vs HFD
Heatmap(subset(resmg, Group %in% c("DREX", "HFD")),
        features=na.omit(unique(unlist(lapply(drexvshfdcmgs, 
                                              function(x) x$Path[1:10])))),
        show_colnames=F,
        group_by=c("Group", "celltype"),
        cluster_cols=F,
        order_by="celltype",
        fontsize_row=8)
dev.print(pdf, file="DREXvsHFD_celltype_heatmap_hallmark_top10.pdf")
