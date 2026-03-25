#' ============================================================================
#' scRNA-seq Processing Pipeline
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
#'   Complete single-cell RNA-seq processing: Cell Ranger data import, QC filtering
#   (200-6000 genes, <5% mitochondrial), SCTransform normalization, Harmony batch
#   integration, UMAP visualization (resolution 0.3), cell type annotation with
#   CellMarker 2.0/PanglaoDB markers, and differential expression (FindMarkers)
#   for all pairwise comparisons across 5 experimental groups.
#'
#' Input:
#'   Cell Ranger filtered_feature_bc_matrix (h5 files), sample_metadata.csv
#'
#' Output:
#'   sc.rdata (annotated Seurat object), UMAP plots (PDF), DEG CSVs per
#   cell type, cell composition pie charts
#'
#' Related figures:
#'   Fig. 3A-B, Fig. S3, Fig. S4, Tables S1-S4
#'
#' ============================================================================
# ==============================================================================
# scRNA-seq Data Processing and Quality Control
# Processing gene-cell matrices from Cell Ranger output (Seurat v4)
# Reference: Methods - "scRNA-seq data processing and quality control (QC)"
# ==============================================================================

setwd("/extData/NGS/UM_R01/res/new")
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 500000 * 1024^2)

# ------------------------------------------------------------------------------
# 1. Data Import and Seurat Object Creation
# FASTQ files aligned to mm10 using Cell Ranger v7.1.0
# ------------------------------------------------------------------------------
filename <- list.files("..", pattern="filtered_feature_bc_matrix", 
                       recursive=T, full.names=T)
obj <- lapply(filename, function(x) Read10X_h5(x))
names(obj) <- sub('\\-','', sub('\\-\\SC\\-.*','', 
                                sub('.*8826-MN','MN', filename)))

# Create Seurat objects for each sample
# min.cells=3: retain genes expressed in at least 3 cells
# min.features=200: retain cells with at least 200 detected genes (initial filter)
object <- lapply(names(obj), function(x) 
  CreateSeuratObject(counts=obj[[x]], project=x, min.cells=3, min.features=200))
names(object) <- names(obj)

# ------------------------------------------------------------------------------
# 2. Add Group Metadata to Each Object
# ------------------------------------------------------------------------------
add_group <- function(x){
  obj <- object[[x]]
  obj$group <- x
  return(obj)
}
object <- lapply(names(object), function(x) add_group(x))
names(object) <- names(obj)

# ------------------------------------------------------------------------------
# 3. Quality Control Metrics Calculation
# Calculate mitochondrial percentage for QC filtering
# ------------------------------------------------------------------------------
calqc <- function(x){
  obj <- x
  # Calculate percentage of mitochondrial genes (pattern "^mt-" for mouse)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^mt-")
  return(obj)
}
object <- lapply(object, function(x) calqc(x))

# Visualize QC metrics for each sample
for(i in names(object)){
  VlnPlot(object[[i]], 
          features=c("nCount_RNA", "nFeature_RNA", "percent.mt"),
          ncol=4, pt.size=0)
  ggsave(file=paste0(i, "_qc.pdf"), width=9, height=5)
}

# ------------------------------------------------------------------------------
# 4. Cell Filtering Based on QC Thresholds
# Methods: "Cells with fewer than 200 or more than 6,000 detected genes, 
#          or with mitochondrial content over 5%, were excluded"
# ------------------------------------------------------------------------------
dat <- lapply(object, function(x) 
  subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5))

# ------------------------------------------------------------------------------
# 5. Normalization
# Methods: "Gene expression was normalized using SCTransform"
# First perform standard normalization for cell cycle scoring
# ------------------------------------------------------------------------------
for (i in 1:length(dat)) {
  dat[[i]] <- NormalizeData(dat[[i]], verbose=T)
}

# ------------------------------------------------------------------------------
# 6. Cell Cycle Scoring
# Acquire mouse cell cycle genes from reference database
# ------------------------------------------------------------------------------
cell_cycle_genes <- read.csv(
  "https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv")

library(GenomicFeatures)
library(AnnotationHub)
library(tidyverse)

# Access Ensembl database for mouse annotation
ah <- AnnotationHub()
ahDb <- query(ah, pattern=c("Mus musculus", "EnsDb"), ignore.case=TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n=1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, return.type="data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Map cell cycle genes to current annotation
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, 
                                       by=c("geneID"="gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase=="S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase=="G2/M") %>%
  pull("gene_name")

# Perform cell cycle scoring
for (i in 1:length(dat)) {
  dat[[i]] <- CellCycleScoring(dat[[i]], 
                               s.features=s_genes, 
                               g2m.features=g2m_genes)
}

# ------------------------------------------------------------------------------
# 7. SCTransform Normalization
# Methods: "Gene expression was normalized using SCTransform"
# Regress out mitochondrial percentage and feature count variation
# ------------------------------------------------------------------------------
for (i in 1:length(dat)) {
  dat[[i]] <- SCTransform(dat[[i]], 
                         vars.to.regress=c('percent.mt', 'nFeature_RNA'))
}

# Load sample metadata
group <- read.csv("sample_metadata.csv", row.names=1)

# ------------------------------------------------------------------------------
# 8. Integration with Harmony
# Methods: "count matrices from different groups were integrated with Harmony"
# ------------------------------------------------------------------------------

# Select integration features (5000 most variable genes)
integ_features <- SelectIntegrationFeatures(object.list=dat, nfeatures=5000) 

# Merge normalized samples
merged_seurat <- merge(x=dat[[1]], y=dat[2:length(dat)], merge.data=TRUE)
DefaultAssay(merged_seurat) <- "SCT"

# Manually set variable features of merged Seurat object
VariableFeatures(merged_seurat) <- integ_features

# Calculate PCs using manually set variable features
merged_seurat <- RunPCA(merged_seurat, assay="SCT", npcs=50)

# Run Harmony integration
library(harmony)
sam <- RunHarmony(merged_seurat, 
                  group.by.vars=c("orig.ident"), 
                  reduction="pca", 
                  assay.use="SCT", 
                  reduction.save="harmony")

# ------------------------------------------------------------------------------
# 9. Dimensionality Reduction and Visualization
# Methods: "visualized with uniform manifold approximation and projection (UMAP)"
# ------------------------------------------------------------------------------
sam <- RunUMAP(sam, reduction="harmony", assay="SCT", dims=1:50)

# ------------------------------------------------------------------------------
# 10. Graph-based Clustering
# Methods: "Graph-based clustering in Seurat v4 at a resolution of 0.3"
# ------------------------------------------------------------------------------
sam <- FindNeighbors(sam, reduction="harmony")

# Test multiple resolutions to find optimal clustering
sam <- FindClusters(
  object=sam, 
  dims.use=1:40, 
  resolution=seq(0.1, 2, 0.1), 
  print.output=FALSE, 
  save.SNN=TRUE
)

# Check number of clusters at each resolution
sapply(grep("^SCT_snn_res", colnames(sam@meta.data), value=TRUE),
       function(x) length(unique(sam@meta.data[,x])))

# Select resolution 0.3 based on methods
sam <- FindClusters(sam, resolution=0.3)

# Add group information from metadata
sam$Group <- group[sam@meta.data$orig.ident, "Group"]

# Visualize clusters
library(rcolors)
DimPlot(sam, label=T, split.by="Group", repel=T, cols=distcolor[1:15])
dev.print(pdf, file="umap_Group.pdf")
DimPlot(sam, label=T, cols=distcolor[1:15])
dev.print(pdf, file="umap_cluster.pdf")

# ------------------------------------------------------------------------------
# 11. Marker Gene Identification
# Methods: "marker genes identified using FindAllMarkers"
# ------------------------------------------------------------------------------
DefaultAssay(sam) <- "RNA"
markers <- FindAllMarkers(sam, only.pos=T)

# ------------------------------------------------------------------------------
# 12. Cell Type Annotation
# Methods: "Cell types were annotated using CellMarker 2.0 and PanglaoDB, 
#          with manual adjustments"
# ------------------------------------------------------------------------------
cells <- cellMarker(markers, species="mouse", weight=0)

# Manual annotation based on known markers:
# Schwann cells: S100b, Gpm6b, Prx, Mbp
# Pericytes: Notch3, Acta2, Higd1b, Abcc9, Rgs5, Pdgfrb
# Endothelial: Egfl7, Cldn5, Esam, Sox17, Pecam1, Vwf
# Fibroblasts: Lum, Dcn, Pdgfra, Col6a1
# Macrophages: H2-Ab1, C1qb, Cd74, Mrc1

# Define marker genes for visualization
genes <- c("Mpz","Plp1","Sox10","Cldn19","Cnp","Ntng1","Itgb4","Lmo7",
           "Msln","Slc2a1","Scn7a","Cadm4","Sostdc1","Art3","Mag",
           "Lum","Dcn","Pdgfra","Col6a1","Dpt","Pi16","Osr1","Dpp4",
           "Dhh","Pmp22","H2-Ab1","C1qb","Cd74","Mrc1","Pecam1","Cldn5",
           "Cd300lg","Egfl7","Rgs5","Kcnj8","Notch3","Acta2","Des","Tpm2",
           "Lpl","Nnat","Pparg","Fabp4","Synpo2","Mbp","Gpm6b","Prx")

# Rename clusters to cell types based on marker expression
sam <- RenameIdents(sam, 
                    `0`="mySC",      # Myelinating Schwann cells
                    `1`="Perineurial", 
                    `2`="nmSC",      # Non-myelinating Schwann cells
                    `3`="EndoFib1",  # Endoneurial Fibroblast 1
                    `4`="mySC", 
                    `5`="EpiFib",    # Epineurial Fibroblast
                    `6`="ImmSC",     # Immature Schwann cells
                    `7`="Mac",       # Macrophages
                    `8`="Endo",      # Endothelial cells
                    `9`="Perineurial", 
                    `10`="Pericytes", 
                    `11`="VSMC",     # Vascular Smooth Muscle Cells
                    `12`="VSMC",
                    `13`="SC3",      # Schwann cell subtype 3
                    `14`='EndoFib2') # Endoneurial Fibroblast 2

sam$celltype <- Idents(sam)

# Create major cell type categories
sam$CellType <- ifelse(sam$celltype %in% c("mySC","SC3","nmSC","ImmSC"), 
                       "SC", sam$celltype)
sam$CellType <- ifelse(sam$celltype %in% c("EndoFib1","EpiFib","EndoFib2"), 
                       "Fib", sam$CellType)
sam$CellType <- ifelse(sam$celltype %in% c("Mac"), "Mac", sam$CellType)
sam$CellType <- ifelse(sam$celltype %in% c("Perineurial"), 
                       "Perineurial", sam$CellType)
sam$CellType <- ifelse(sam$celltype %in% c("VSMC"), "SMC", sam$CellType)
sam$CellType <- ifelse(sam$celltype %in% c("Pericytes"), 
                       "Pericytes", sam$CellType)
sam$CellType <- ifelse(sam$celltype %in% c("Endo"), "Endo", sam$CellType)

# Set cell type factor levels for consistent ordering
sam$celltype <- factor(as.character(sam$celltype), 
                       levels=c("ImmSC","nmSC","mySC","SC3","Perineurial",
                               "EndoFib1","EndoFib2","EpiFib","Mac","Endo",
                               "Pericytes","VSMC"))
Idents(sam) <- "celltype"

# Visualize marker genes with dot plot
DotPlot(sam, features=genes) + 
  coord_flip() + 
  scale_color_viridis_c() + 
  theme_light(base_size=12) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + 
  xlab("") + ylab("")
dev.print(pdf, file="celltype_marker_new.pdf")

# Visualize major cell types
sam$CellType <- factor(as.character(sam$CellType), 
                       levels=c("SC","Perineurial","Fib","Mac","Endo",
                               "Pericytes","SMC"))
Idents(sam) <- "CellType"
DotPlot(sam, features=genes) + 
  coord_flip() + 
  scale_color_viridis_c() + 
  theme_light(base_size=12) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + 
  xlab("") + ylab("")
dev.print(pdf, file="celltype_marker_major.pdf")

# Generate feature plots for each marker gene
for(i in genes){
  FeaturePlot(sam, features=i, pt.size=0.1, label=T, repel=T) + 
    scale_color_viridis_c()
  ggsave(file=paste0(i, "_umap.pdf"), height=3.6, width=4.6)
}

# Finalize visualization with custom colors
library(rcolors)
mycols <- distcolor[1:12]
names(mycols) <- levels(sam$celltype)

DimPlot(sam, label=T, cols=mycols, repel=T)
dev.print(pdf, file="umap_celltype.pdf")

DimPlot(sam, label=T, split.by="Group", repel=T, cols=mycols)
dev.print(pdf, file="umap_Group.pdf")

DimPlot(sam, label=T, split.by="Phase", repel=T, cols=mycols)
dev.print(pdf, file="umap_cycle.pdf")

# Create condition labels (Group_CellType)
sam$condition <- paste0(sam$Group, "_", sam$celltype)
Idents(sam) <- "condition"

# Extract metadata
meta <- sam@meta.data
Idents(sam) <- "celltype"

# ------------------------------------------------------------------------------
# 13. Cell Type Proportion Analysis
# Visualize cell type composition for each experimental group
# ------------------------------------------------------------------------------

library(ggrepel)

# Generate pie charts for each group showing cell type proportions

# SD (Standard Diet) group
meta %>%
  dplyr::filter(Group=="SD") %>%
  dplyr::select(Group, celltype) %>%
  group_by(Group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  ggplot(aes(x=2, cell, fill=celltype)) +
  geom_bar(width=1, size=0.1, color="white", stat="identity") +
  scale_fill_manual(values=mycols) +
  geom_text_repel(aes(label=paste0(round(100*cell, 1), "%")),
                  position=position_stack(vjust=0.5), max.overlaps=100) + 
  xlim(0.5, 2.5) +
  coord_polar("y") + 
  guides(fill=guide_legend(ncol=2)) +  
  theme_classic() +
  labs(x=NULL, y=NULL, fill=NULL, title="SD") +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(hjust=0.5, color="#666666"))
dev.print(pdf, file="SD_circle.pdf")

# HFD (High-Fat Diet) group
meta %>%
  dplyr::filter(Group=="HFD") %>%
  dplyr::select(Group, celltype) %>%
  group_by(Group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  ggplot(aes(x=2, cell, fill=celltype)) +
  geom_bar(width=1, size=0.1, color="white", stat="identity") +
  scale_fill_manual(values=mycols) +
  geom_text_repel(aes(label=paste0(round(100*cell, 1), "%")),
                  position=position_stack(vjust=0.5), max.overlaps=100) + 
  xlim(0.5, 2.5) +
  coord_polar("y") + 
  guides(fill=guide_legend(ncol=2)) +  
  theme_classic() +
  labs(x=NULL, y=NULL, fill=NULL, title="HFD") +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(hjust=0.5, color="#666666"))
dev.print(pdf, file="HFD_circle.pdf")

# EX (Exercise) group
meta %>%
  dplyr::filter(Group=="EX") %>%
  dplyr::select(Group, celltype) %>%
  group_by(Group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  ggplot(aes(x=2, cell, fill=celltype)) +
  geom_bar(width=1, size=0.1, color="white", stat="identity") +
  scale_fill_manual(values=mycols) +
  geom_text_repel(aes(label=paste0(round(100*cell, 1), "%")),
                  position=position_stack(vjust=0.5), max.overlaps=100) + 
  xlim(0.5, 2.5) +
  coord_polar("y") + 
  guides(fill=guide_legend(ncol=2)) +  
  theme_classic() +
  labs(x=NULL, y=NULL, fill=NULL, title="EX") +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(hjust=0.5, color="#666666"))
dev.print(pdf, file="EX_circle.pdf")

# DREX (Diet Restriction + Exercise) group
meta %>%
  dplyr::filter(Group=="DREX") %>%
  dplyr::select(Group, celltype) %>%
  group_by(Group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  ggplot(aes(x=2, cell, fill=celltype)) +
  geom_bar(width=1, size=0.1, color="white", stat="identity") +
  scale_fill_manual(values=mycols) +
  geom_text_repel(aes(label=paste0(round(100*cell, 1), "%")),
                  position=position_stack(vjust=0.5), max.overlaps=100) + 
  xlim(0.5, 2.5) +
  coord_polar("y") + 
  guides(fill=guide_legend(ncol=2)) +  
  theme_classic() +
  labs(x=NULL, y=NULL, fill=NULL, title="DREX") +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(hjust=0.5, color="#666666"))
dev.print(pdf, file="DREX_circle.pdf")

# DR (Diet Restriction) group
meta %>%
  dplyr::filter(Group=="DR") %>%
  dplyr::select(Group, celltype) %>%
  group_by(Group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  ggplot(aes(x=2, cell, fill=celltype)) +
  geom_bar(width=1, size=0.1, color="white", stat="identity") +
  scale_fill_manual(values=mycols) +
  geom_text_repel(aes(label=paste0(round(100*cell, 1), "%")),
                  position=position_stack(vjust=0.5), max.overlaps=100) + 
  xlim(0.5, 2.5) +
  coord_polar("y") + 
  guides(fill=guide_legend(ncol=2)) +  
  theme_classic() +
  labs(x=NULL, y=NULL, fill=NULL, title="DR") +
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.title=element_text(hjust=0.5, color="#666666"))
dev.print(pdf, file="DR_circle.pdf")

# ------------------------------------------------------------------------------
# 14. Statistical Testing of Cell Type Proportion Differences
# Wilcoxon rank-sum test for group comparisons
# ------------------------------------------------------------------------------

library(broom)

# HFD vs SD comparison
meta %>%
  group_by(Group, group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  dplyr::filter(Group %in% c("HFD","SD")) %>%
  ungroup %>%
  dplyr::select(Group, celltype, cell) %>%
  group_by(celltype) %>%
  do(tidy(wilcox.test(cell~Group, data=.))) %>%
  write.csv(file="HFDvsSD_wilcox.csv")

# HFD vs DR comparison
meta %>%
  group_by(Group, group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  dplyr::filter(Group %in% c("HFD","DR")) %>%
  ungroup %>%
  dplyr::select(Group, celltype, cell) %>%
  group_by(celltype) %>%
  do(tidy(wilcox.test(cell~Group, data=.))) %>%
  write.csv(file="HFDvsDR_wilcox.csv")

# DR vs SD comparison
meta %>%
  group_by(Group, group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  dplyr::filter(Group %in% c("DR","SD")) %>%
  ungroup %>%
  dplyr::select(Group, celltype, cell) %>%
  group_by(celltype) %>%
  do(tidy(wilcox.test(cell~Group, data=.))) %>%
  write.csv(file="DRvsSD_wilcox.csv")

# DREX vs DR comparison
meta %>%
  group_by(Group, group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  dplyr::filter(Group %in% c("DR","DREX")) %>%
  ungroup %>%
  dplyr::select(Group, celltype, cell) %>%
  group_by(celltype) %>%
  do(tidy(wilcox.test(cell~Group, data=.))) %>%
  write.csv(file="DREXvsDR_wilcox.csv")

# DREX vs EX comparison
meta %>%
  group_by(Group, group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  dplyr::filter(Group %in% c("EX","DREX")) %>%
  ungroup %>%
  dplyr::select(Group, celltype, cell) %>%
  group_by(celltype) %>%
  do(tidy(wilcox.test(cell~Group, data=.))) %>%
  write.csv(file="DREXvsEX_wilcox.csv")

# EX vs SD comparison
meta %>%
  group_by(Group, group, celltype) %>%
  summarise(count=n()) %>%
  mutate(cell=count/sum(count)) %>%
  dplyr::filter(Group %in% c("EX","SD")) %>%
  ungroup %>%
  dplyr::select(Group, celltype, cell) %>%
  group_by(celltype) %>%
  do(tidy(wilcox.test(cell~Group, data=.))) %>%
  write.csv(file="EXvsSD_wilcox.csv")

# ------------------------------------------------------------------------------
# 15. Differential Expression Analysis Between Conditions
# Methods: "Differentially expressed genes (DEGs) between clusters were 
#          identified using Seurat's FindMarkers (adjusted p-value<0.05)"
# For pseudo-bulk analysis, see next section
# ------------------------------------------------------------------------------

Idents(sam) <- "condition"

# HFD vs SD for each cell type
hfdsd <- lapply(unique(sam$celltype), 
                function(x) FindMarkers(sam, 
                                       ident.1=paste0("HFD_", x),
                                       ident.2=paste0("SD_", x),
                                       logfc.threshold=0))
names(hfdsd) <- unique(sam$celltype)
sapply(names(hfdsd), function(x) 
  write.csv(hfdsd[[x]], file=paste0("HFDvsSD_", x, ".csv")))

# DR vs HFD for each cell type
drhfd <- lapply(unique(sam$celltype), 
                function(x) FindMarkers(sam,
                                       ident.1=paste0("DR_", x),
                                       ident.2=paste0("HFD_", x),
                                       logfc.threshold=0))
names(drhfd) <- unique(sam$celltype)
sapply(names(drhfd), function(x) 
  write.csv(drhfd[[x]], file=paste0("DRvsHFD_", x, ".csv")))

# DREX vs DR for each cell type
drexdr <- lapply(unique(sam$celltype), 
                 function(x) FindMarkers(sam,
                                        ident.1=paste0("DREX_", x),
                                        ident.2=paste0("DR_", x),
                                        logfc.threshold=0))
names(drexdr) <- unique(sam$celltype)
sapply(names(drexdr), function(x) 
  write.csv(drexdr[[x]], file=paste0("DREXvsDR_", x, ".csv")))

# EX vs DREX for each cell type
exdrex <- lapply(unique(sam$celltype), 
                 function(x) FindMarkers(sam,
                                        ident.1=paste0("EX_", x),
                                        ident.2=paste0("DREX_", x),
                                        logfc.threshold=0))
names(exdrex) <- unique(sam$celltype)
sapply(names(exdrex), function(x) 
  write.csv(exdrex[[x]], file=paste0("EXvsDREX_", x, ".csv")))

# EX vs HFD for each cell type
exhfd <- lapply(unique(sam$celltype), 
                function(x) FindMarkers(sam,
                                       ident.1=paste0("EX_", x),
                                       ident.2=paste0("HFD_", x),
                                       logfc.threshold=0))
names(exhfd) <- unique(sam$celltype)
sapply(names(exhfd), function(x) 
  write.csv(exhfd[[x]], file=paste0("EXvsHFD_", x, ".csv")))

# DR vs SD for each cell type
drsd <- lapply(unique(sam$celltype), 
               function(x) FindMarkers(sam,
                                      ident.1=paste0("DR_", x),
                                      ident.2=paste0("SD_", x),
                                      logfc.threshold=0))
names(drsd) <- unique(sam$celltype)
sapply(names(drsd), function(x) 
  write.csv(drsd[[x]], file=paste0("DRvsSD_", x, ".csv")))

# Save the final Seurat object
save(sam, file="sc.rdata")

