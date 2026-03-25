# DREX: Multi-omics Reveals How Diet and Exercise Improve Peripheral Neuropathy

This repository contains the analysis code accompanying our manuscript:

> Eid SA, Guo K, Noureldein MH, Jang D-G, Allouch AM, Miller CM, Kiriluk CP, Lentz W, Boutary S, Mule JJ, Pennathur S, Bhatt DK, Feldman EL, Hur J. Multi-omics reveals how diet and exercise improve peripheral neuropathy. *[Journal TBD]*, 2026.

The study investigates how dietary reversal (DR), exercise (EX), and their combination (DREX) reverse obesity-related peripheral neuropathy through Schwann cell metabolic reprogramming, integrating single-cell RNA sequencing, metabolomics, and fluxomics of sciatic nerve tissue.

## Study Overview

**Experimental Groups (5):**

- **SD** — Standard Diet (control)
- **HFD** — High-Fat Diet (obesity/prediabetes model)
- **DR** — Dietary Reversal (switched from HFD to SD at 18 weeks)
- **EX** — Exercise (maintained on HFD with running wheel access)
- **DREX** — Dietary Reversal + Exercise (combined intervention)

**Cell Types Identified (12):**
Myelinating Schwann cells (mySC), non-myelinating Schwann cells (nmSC), immature Schwann cells (ImmSC), SC3, endothelial cells (Endo), endoneurial fibroblasts (EndoFib1, EndoFib2), epineurial fibroblasts (EpiFib), macrophages (Mac), pericytes, perineurial cells, vascular smooth muscle cells (VSMC)

---

## Repository Structure

```
├── README.md
├── .gitignore
├── scRNAseq/
│   └── scripts/                    # 9 R scripts for scRNA-seq analyses
├── metabolomics/
│   ├── scripts/                    # 3 R scripts for metabolomics analyses
│   ├── results/                    # Acylcarnitine heatmaps and ratio plots
│   └── Carnitine_annotation_with_chainClass_07032025.csv
├── fluxomics/
│   └── scripts/                    # 1 R script for isotopologue fluxomics (26 weeks)
└── integration/
    ├── scripts/                    # 4 R scripts for multi-omics integration
    ├── data/                       # Pre-processed intermediate data for integration
    ├── results_clustering/         # Mfuzz soft clustering outputs (k=10)
    └── results_network/            # Pathway-based network visualizations
```

---

## Figure-to-Script Mapping

The table below maps each manuscript figure to the script(s) that generated or contributed to it.

| Figure | Description | Primary Script(s) |
|--------|-------------|-------------------|
| Fig. 3A–B | UMAP of cell types and Schwann cell subtypes | `scRNAseq/scripts/run_step_final.r` |
| Fig. 3C–F | KEGG pathway network plots from SC DEGs | `scRNAseq/scripts/run_enrichment_network.r`, `run_pathway_enrichment.R` |
| Fig. 4A–B | Metabolomics VIP pathway enrichment (polar, untargeted) | `metabolomics/scripts/run_metabolite.R`, `run_pathway_gsva.R` |
| Fig. 5A–D | Transcriptome-metabolome network plots | `integration/scripts/run_network_visualization.R`, `run_pathway_integration.R`, `add_fc_gene_network.R` |
| Fig. 6 | Fluxomics isotopologue analysis (13C-glucose tracing) | `fluxomics/scripts/run_fluxomics_26wk_m0_sum.R` |
| Fig. 7A–C | Multi-omics Mfuzz clustering by SC subtype | `integration/scripts/run_mfuzz_clustering_JH.R` |
| Fig. S3 | Cell markers and composition | `scRNAseq/scripts/run_step_final.r`, `run_cell_fractions.r` |
| Fig. S4 | SC subtype composition and shared DEGs | `scRNAseq/scripts/run_cell_fractions.r`, `run_shared_DEGs_heatmap.r` |
| Fig. S5–S7 | SC subtype DEGs and KEGG enrichment | `scRNAseq/scripts/run_step_final.r`, `run_gsva.r`, `run_pathway_enrichment.R` |
| Fig. S8 | Human–mouse transcriptomic overlap | `scRNAseq/scripts/run_human_comparison.R` |
| Fig. S9 | PLS-DA of polar and untargeted metabolomics | `metabolomics/scripts/run_metabolite.R` |
| Fig. S10 | Metabolomics pathway VIP enrichment | `metabolomics/scripts/run_pathway_gsva.R` |
| Fig. S11–S12 | Fluxomics isotopologue heatmaps | `fluxomics/scripts/run_fluxomics_26wk_m0_sum.R` |
| Fig. S13–S15 | Mfuzz clustering curves per SC subtype | `integration/scripts/run_mfuzz_clustering_JH.R` |
| Fig. S16 | Phenotype-cluster dot plots | `integration/scripts/run_mfuzz_clustering_JH.R` |
| Table S1–S4 | DEG tables for SC subtypes | `scRNAseq/scripts/run_step_final.r` |
| Table S5–S7 | KEGG enrichment for SC subtypes | `scRNAseq/scripts/run_pathway_enrichment.R`, `run_gsva.r` |
| Table S8–S9 | Human–mouse overlap and enrichment | `scRNAseq/scripts/run_human_comparison.R` |
| Table S10 | Mfuzz cluster features | `integration/scripts/run_mfuzz_clustering_JH.R` |

---

## scRNAseq Analysis Scripts

All scripts are located in `scRNAseq/scripts/` (9 scripts).

### Core Pipeline

| Script | Description | Input | Output |
|--------|-------------|-------|--------|
| `run_step_final.r` | Complete scRNA-seq processing: Cell Ranger import, QC filtering (200–6000 genes, <5% mito), SCTransform normalization, Harmony batch integration, UMAP clustering (resolution 0.3), cell type annotation, and differential expression (FindMarkers) for all pairwise comparisons | Cell Ranger filtered_feature_bc_matrix (h5 files), sample metadata | `sc.rdata` (annotated Seurat object), UMAP plots, DEG CSVs per cell type |
| `run_gsva.r` | Gene Set Variation Analysis using KEGG and HALLMARK gene sets. Computes pathway activity scores per cell and tests for significant differences using Wilcoxon tests with BH correction | `sc.rdata`, KEGG mouse annotations, MSigDB HALLMARK gene sets | GSVA scores, significant pathway lists (padj < 0.05), violin plots, heatmaps |
| `run_pseudo.r` | Pseudobulk aggregation and DESeq2 differential expression. Aggregates single-cell counts per sample/cell type, filters low-count genes (>60 total), runs DESeq2 | Pseudobulk count data, sample metadata | `sc_integration_predata.rdata`, DESeq2 results CSVs per cell type |

### Pathway and Enrichment Analyses

| Script | Description |
|--------|-------------|
| `run_pathway_enrichment.R` | KEGG pathway enrichment bubble plots colored by -log10(p-value) and sized by rich factor |
| `run_enrichment_network.r` | Enrichment network visualizations with pathway term trimming and padj filtering |

### Pseudobulk and Cross-Species Comparisons

| Script | Description |
|--------|-------------|
| `run_pseudobulk_DE.R` | Comprehensive pseudobulk DE with integrated GSVA pathway scoring and fGSEA |
| `run_human_comparison.R` | Mouse scRNA-seq DEGs versus human sural nerve transcriptomics via homologene |

### Cell Composition and Shared DEGs

| Script | Description |
|--------|-------------|
| `run_cell_fractions.r` | Cell type composition counts and fractions across groups |
| `run_shared_DEGs_heatmap.r` | Shared DEGs across contrasts with Venn diagrams and heatmaps |

---

## Metabolomics Analysis Scripts

All scripts are located in `metabolomics/scripts/` (3 scripts).

| Script | Description |
|--------|-------------|
| `run_metabolite.R` | Main metabolomics pipeline for 4 datasets (SCN, Acylcarnitines, Polar, Untargeted): log2 transformation, Wilcoxon tests (BH correction), PLS-DA with VIP scoring, KEGG enrichment |
| `run_pathway_gsva.R` | GSVA scoring of metabolic pathways followed by Wilcoxon testing and PLS-DA |
| `omics.r` | Helper function library (~3,200 lines) for PCA, PLS-DA, VIP scoring, and visualization |

**Pairwise Comparisons (5):** HFD_18WK vs WT_18WK, HFD_26WK vs WT_26WK, DR_26WK vs HFD_26WK, EX_26WK vs HFD_26WK, DREX_26WK vs HFD_26WK

---

## Fluxomics Analysis Scripts

All scripts are located in `fluxomics/scripts/` (1 script). These process isotopologue labeling data from [U-13C6]-glucose tracer experiments in sciatic nerve at 26 weeks.

| Script | Description |
|--------|-------------|
| `run_fluxomics_26wk_m0_sum.R` | 26-week isotopologue analysis with M+0 normalization: each isotopologue is normalized to the corresponding M+0 value, followed by statistical testing across all 5 groups with bar plots and heatmaps |

---

## Integration Analysis Scripts

All scripts are located in `integration/scripts/` (4 scripts). Supporting data files are in `integration/data/`.

### Temporal Clustering (Mfuzz)

| Script | Description |
|--------|-------------|
| `run_mfuzz_clustering_JH.R` | Multi-omics Mfuzz soft clustering with `run_mfuzz_batch()`, seed setting, multiple metabolite configurations, and KEGG enrichment across 5 groups (k=10) for 11 cell types |

### Pathway-Based Network Integration

| Script | Description |
|--------|-------------|
| `run_pathway_integration.R` | Gene-metabolite-pathway network construction per cell type via KEGG annotations |
| `add_fc_gene_network.R` | Metabolite fold-change integration from 4 sources into pathway networks |
| `run_network_visualization.R` | Publication-quality network plots (v2.1) with three color scales and Cytoscape export |

### Integration Supporting Data (`integration/data/`)

| File | Description |
|------|-------------|
| `sc_integration_predata.rdata` | Pre-processed pseudobulk single-cell data (averaged by cell type and condition) |
| `metabol_integration_predata.rdata` | Pre-processed metabolomics data for integration |
| `fluxomics_26wk_final_m0_avg.csv` | Averaged M+0-normalized fluxomics data (26 weeks) |
| `sc_deg_bulk_all.csv` | Pseudobulk differential expression results |
| `step_path.csv` | KEGG pathway mapping for metabolites |
| `gene_metabolite_pathway_group.csv` | Gene-metabolite-pathway association table |
| `gene_metabolite_pathway_group_MetaboliteFC_v2.csv` | Gene-metabolite-pathway table with metabolite fold-change values |

---

## R Package Dependencies

**Core:** Seurat (v4+), DESeq2, harmony, SCTransform

**Pathway and Enrichment:** GSVA, scGSVA, fgsea, richR, enrichR

**Metabolomics and Multivariate:** mixOmics, FactoMineR, Mfuzz

**Visualization:** ggplot2, pheatmap, ComplexHeatmap, ggraph, tidygraph, igraph, VennDiagram, VennDetail, UpSetR

**Utilities:** homologene, readxl, dplyr, tidyr, rstatix, data.table

---

## Data Availability

Raw sequencing data (scRNA-seq) and metabolomics/fluxomics datasets are not included in this repository due to size. Raw data will be available through GEO (accession number TBD) upon publication. The `integration/data/` directory contains pre-processed intermediate files required for integration analyses.

**Input data requirements:**

- **scRNA-seq:** Cell Ranger v7.1.0 output (filtered_feature_bc_matrix, mm10 reference genome)
- **Metabolomics:** CSV files for 4 metabolite datasets (SCN, Acylcarnitines, Polar, Untargeted)
- **Fluxomics:** Excel files with isotopologue distributions from [U-13C6]-glucose tracer experiments

---

## Statistical Methods

| Analysis | Method | Threshold | Correction |
|----------|--------|-----------|------------|
| Cell QC filtering | nFeature + mito% | 200–6000 genes, <5% mito | — |
| Graph-based clustering | Seurat Louvain | Resolution 0.3 | — |
| Single-cell DE | Wilcoxon rank-sum | p < 0.05 | BH |
| Pseudobulk DE | DESeq2 | p < 0.01 | BH |
| Pathway analysis (GSVA) | Wilcoxon rank-sum | padj < 0.05 | BH |
| Metabolomics DE | Wilcoxon rank-sum | p < 0.05 | BH |
| PLS-DA VIP | Variable Importance in Projection | VIP > 1 | — |
| Gene filtering (pseudobulk) | Total count | >= 60 counts | — |
| Mfuzz soft clustering | Fuzzy c-means | Membership > 0.5 | — |

---

## Notes on Reproducibility

Scripts contain `setwd()` calls with paths specific to the original analysis environment. Users will need to update these paths to match their local directory structure. All scripts assume the input data files (raw data, Seurat objects, etc.) are available in the working directory or at the specified relative paths.

---

## License

This repository is provided for academic and research use. Please cite the accompanying manuscript if you use these scripts.

---

## Contact

- **Junguk Hur, PhD** — Department of Biomedical Sciences, University of North Dakota, Grand Forks, ND (computational analysis)
- **Eva L. Feldman, MD, PhD** — Department of Neurology, University of Michigan, Ann Arbor, MI (corresponding author)
