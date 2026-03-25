#' ============================================================================
#' Metabolomics Analysis Pipeline
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
#'   Main metabolomics pipeline for 4 datasets: SCN (sciatic nerve, 40 metabolites),
#   Acylcarnitines, Polar, and Untargeted. Performs log2 transformation, Wilcoxon
#   tests with BH correction, fold-change calculations, PLS-DA with VIP scoring
#   (threshold > 1), and KEGG pathway enrichment.
#'
#' Input:
#'   Metabolite CSV files (SCN, Acyl, Polar, Untargeted), Sciatic nerves
#   tissues -untargeted expanded_40metabolites.csv
#'
#' Output:
#'   Wilcoxon result CSVs, heatmaps (p < 0.05), Venn diagrams, PLS-DA plots,
#   VIP score CSVs/PDFs, KEGG enrichment CSVs/PDFs
#'
#' Related figures:
#'   Fig. 4A-B, Fig. S9
#'
#' ============================================================================
# ==============================================================================
# Metabolomics Data Analysis
# Methods: "Polar and untargeted datasets underwent log2 transformation"
# PLS-DA, VIP scores, and pathway enrichment analysis
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Load Multiple Metabolomics Datasets
# ------------------------------------------------------------------------------
filenames <- list.files(pattern="*csv")
d <- lapply(filenames, function(x) read.csv(x, check.names=F, row.names=1))
names(d) <- c("SCN", "Acyl", "Polar", "Untargeted")

# Extract metadata for each dataset
metas <- d$SCN[, 1:2]
metaa <- d$Acyl[, 1:2]
metap <- d$Polar[, 1:2]
metau <- d$Untargeted[, 1:2]

# Extract numerical data (metabolite abundances)
dd <- lapply(d, function(x) x[, 3:ncol(x)])

# ------------------------------------------------------------------------------
# 2. Data Transformation
# Methods: "Polar and untargeted datasets underwent log2 transformation"
# ------------------------------------------------------------------------------

# Log2 transform SCN dataset
dd$SCN <- log2(dd$SCN + 1)

# Load and merge additional SCN metabolites
scn2 <- read.csv("../Sciatic nerves tissues -untargeted expanded_40metabolites.csv",
                 row.names=1, check.names=F)
scn2 <- scn2[rownames(dd$SCN), ]
dd$SCN <- cbind(dd$SCN, log2(scn2[, 3:ncol(scn2)] + 1))

# Log2 transform Untargeted dataset
dd$Untargeted <- log2(dd$Untargeted + 1)

# Note: One metabolite overlaps between SCN and Untargeted 
# (N-Acetylaspartylglutamic acid, correlation=0.41)
# Use Untargeted version and add unique SCN metabolites
dd$Untargeted <- cbind(dd$Untargeted, 
                       dd$SCN[, setdiff(colnames(dd$SCN), 
                                       "N-Acetylaspartylglutamic acid")])

# ------------------------------------------------------------------------------
# 3. Wilcoxon Rank-Sum Test for Group Comparisons
# Statistical testing for all metabolites
# ------------------------------------------------------------------------------

library(rstatix)

# SCN dataset
scnw <- dd$SCN %>%
  rownames_to_column(var="id") %>%
  gather(met, val, -id) %>%
  mutate(Group=metas[id, "Group"]) %>%
  group_by(met) %>%
  wilcox_test(val~Group) %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj")

# Acyl dataset
acyw <- dd$Acyl %>%
  rownames_to_column(var="id") %>%
  gather(met, val, -id) %>%
  mutate(Group=metaa[id, "Group"]) %>%
  group_by(met) %>%
  wilcox_test(val~Group) %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj")

# Untargeted dataset
untarw <- dd$Untargeted %>%
  rownames_to_column(var="id") %>%
  gather(met, val, -id) %>%
  mutate(Group=metau[id, "Group"]) %>%
  group_by(met) %>%
  wilcox_test(val~Group) %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj")

# Polar dataset
polarw <- dd$Polar %>%
  rownames_to_column(var="id") %>%
  gather(met, val, -id) %>%
  mutate(Group=metap[id, "Group"]) %>%
  group_by(met) %>%
  wilcox_test(val~Group) %>%
  adjust_pvalue(method="BH") %>%
  add_significance("p.adj")

# Save statistical results
write.csv(scnw, file="SCN_wilox_all.csv")
write.csv(acyw, file="Acyl_wilox.csv")
write.csv(untarw, file="Untargeted_wilox_all.csv")
write.csv(polarw, file="Polar_wilox.csv")

# ------------------------------------------------------------------------------
# 4. Generate Comparison Matrix
# Define pairwise comparisons for analysis
# ------------------------------------------------------------------------------
make.comparison <- function(group, ref=NULL){
  if(is.null(ref)){
    ref <- group
  }
  tmp <- expand.grid(group, ref, stringsAsFactors=F)
  tmp <- tmp[tmp$Var1 != tmp$Var2, ]
  # Remove duplicate pairs
  tmp <- tmp[!duplicated(t(apply(tmp, 1, sort))), ]
  return(tmp)
}

# Define comparisons against reference groups
comp <- make.comparison(unique(metas$Group), 
                        ref=c("WT_18WK", "HFD_18WK", "WT_26WK", "HFD_26WK"))

# ------------------------------------------------------------------------------
# 5. Calculate Mean Expression and Fold-Changes
# Methods: "Log2Fold-Changes were calculated between groups"
# ------------------------------------------------------------------------------

# Calculate mean for each metabolite per group
scnm <- dd$SCN %>%
  rownames_to_column(var="id") %>%
  gather(met, val, -id) %>%
  mutate(Group=metas[id, "Group"]) %>%
  group_by(met, Group) %>%
  summarise(mu=mean(val)) %>%
  spread(Group, mu)

acym <- dd$Acyl %>%
  rownames_to_column(var="id") %>%
  gather(met, val, -id) %>%
  mutate(Group=metaa[id, "Group"]) %>%
  group_by(met, Group) %>%
  summarise(mu=mean(val)) %>%
  spread(Group, mu)

polarm <- dd$Polar %>%
  rownames_to_column(var="id") %>%
  gather(met, val, -id) %>%
  mutate(Group=metap[id, "Group"]) %>%
  group_by(met, Group) %>%
  summarise(mu=mean(val)) %>%
  spread(Group, mu)

untm <- dd$Untargeted %>%
  rownames_to_column(var="id") %>%
  gather(met, val, -id) %>%
  mutate(Group=metau[id, "Group"]) %>%
  group_by(met, Group) %>%
  summarise(mu=mean(val)) %>%
  spread(Group, mu)

# Calculate fold-changes for each comparison

# SCN: log2-transformed, so use subtraction for fold-change
scnf <- cbind(scnm, do.call(cbind, apply(comp, 1, function(x){
  nam <- paste0(x[1], "vs", x[2])
  res <- scnm[, x[1]] - scnm[, x[2]]  # Difference in log2 space
  names(res) <- nam
  return(res)
})))

# Acyl: not log-transformed, use log2 ratio
acyf <- cbind(acym, do.call(cbind, apply(comp, 1, function(x){
  nam <- paste0(x[1], "vs", x[2])
  res <- log2(acym[, x[1]] / acym[, x[2]])
  names(res) <- nam
  return(res)
})))

# Polar: not log-transformed, use log2 ratio
polarf <- cbind(polarm, do.call(cbind, apply(comp, 1, function(x){
  nam <- paste0(x[1], "vs", x[2])
  res <- log2(polarm[, x[1]] / polarm[, x[2]])
  names(res) <- nam
  return(res)
})))

# Untargeted: log2-transformed, use subtraction
untf <- cbind(untm, do.call(cbind, apply(comp, 1, function(x){
  nam <- paste0(x[1], "vs", x[2])
  res <- untm[, x[1]] - untm[, x[2]]
  names(res) <- nam
  return(res)
})))

# ------------------------------------------------------------------------------
# 6. Merge Statistical Results with Fold-Changes
# ------------------------------------------------------------------------------
scnr <- scnw %>%
  mutate(comparision=paste0(group1, "vs", group2)) %>%
  filter(group2 %in% c("WT_18WK", "HFD_18WK", "WT_26WK", "HFD_26WK")) %>%
  left_join(scnf, by=c("met"="met"))

acyr <- acyw %>%
  mutate(comparision=paste0(group1, "vs", group2)) %>%
  filter(group2 %in% c("WT_18WK", "HFD_18WK", "WT_26WK", "HFD_26WK")) %>%
  left_join(acyf, by=c("met"="met"))

polarr <- polarw %>%
  mutate(comparision=paste0(group1, "vs", group2)) %>%
  filter(group2 %in% c("WT_18WK", "HFD_18WK", "WT_26WK", "HFD_26WK")) %>%
  left_join(polarf, by=c("met"="met"))

untr <- untarw %>%
  mutate(comparision=paste0(group1, "vs", group2)) %>%
  filter(group2 %in% c("WT_18WK", "HFD_18WK", "WT_26WK", "HFD_26WK")) %>%
  left_join(untf, by=c("met"="met"))

# ------------------------------------------------------------------------------
# 7. Heatmap Visualization of Significant Metabolites
# ------------------------------------------------------------------------------
library(pheatmap)

# SCN heatmap (p < 0.05)
pheatmap(t(dd$SCN[, unique(subset(scnw, p<0.05)$met)]),
         scale="row",
         annotation_col=metas[rownames(dd$SCN), "Group", drop=F],
         cluster_cols=F,
         show_colnames=F,
         color=colorRampPalette(c("darkblue", "white", "red"))(128),
         border_color="white",
         cellwidth=8,
         cellheight=10)
dev.print(pdf, file="SCN_wilcox_p_all.pdf")

# Acyl heatmap
pheatmap(t(dd$Acyl[, unique(subset(acyw, p<0.05)$met)]),
         scale="row",
         annotation_col=metaa[rownames(dd$Acyl), "Group", drop=F],
         cluster_cols=F,
         show_colnames=F,
         color=colorRampPalette(c("darkblue", "white", "red"))(128),
         border_color="white",
         cellwidth=8,
         cellheight=10)
dev.print(pdf, file="Acyl_wilcox_p.pdf")

# Polar heatmap
pheatmap(t(dd$Polar[, unique(subset(polarw, p<0.05)$met)]),
         scale="row",
         annotation_col=metap[rownames(dd$Polar), "Group", drop=F],
         cluster_cols=F,
         show_colnames=F,
         color=colorRampPalette(c("darkblue", "white", "red"))(128),
         border_color="white",
         cellwidth=8,
         cellheight=10)
dev.print(pdf, file="Polar_wilcox_p.pdf")

# Untargeted heatmap
pheatmap(t(dd$Untargeted[, unique(subset(untarw, p<0.05)$met)]),
         scale="row",
         annotation_col=metau[rownames(dd$Untargeted), "Group", drop=F],
         cluster_cols=F,
         show_colnames=F,
         color=colorRampPalette(c("darkblue", "white", "red"))(128),
         border_color="white",
         cellwidth=8,
         cellheight=10)
dev.print(pdf, file="Untargeted_wilcox_p.pdf")

# ------------------------------------------------------------------------------
# 8. Venn Diagram Analysis of Significantly Changed Metabolites
# ------------------------------------------------------------------------------

# Split results by comparison
scnrs <- split(scnr, scnr$comparision)
acyrs <- split(acyr, acyr$comparision)
polarrs <- split(polarr, polarr$comparision)
untrs <- split(untr, untr$comparision)

# Convert to data frames
acyrd <- lapply(acyrs, function(x) as.data.frame(x))
polarrd <- lapply(polarrs, function(x) as.data.frame(x))
untrd <- lapply(untrs, function(x) as.data.frame(x))

# Extract significant metabolites (p < 0.05)
scnrg <- lapply(scnrs, function(x) subset(x, p<0.05)$met)
acyrg <- lapply(acyrs, function(x) subset(x, p<0.05)$met)
polarrg <- lapply(polarrs, function(x) subset(x, p<0.05)$met)
untrg <- lapply(untrs, function(x) subset(x, p<0.05)$met)

# Define key comparisons for Venn diagrams
cont <- c("HFD_18WKvsWT_18WK", "HFD_26WKvsWT_26WK", "DR_26wkvsHFD_26WK", 
          "EX_26WKvsHFD_26WK", "DREX_26WKvsHFD_26WK")

# Generate Venn diagrams
plot(venndetail(untrg[cont]))
write.csv(getFeature(venndetail(untrg[cont]), rlist=untrd[cont], 
                     userowname=F, gind=rep(1,5)),
          file="untargeted_venndetail.csv")
dev.print(pdf, file="Untargeted_venn.pdf")

plot(venndetail(scnrg[cont]))
dev.print(pdf, file="SCN_venn.pdf")

plot(venndetail(acyrg[cont]))
dev.print(pdf, file="Acyr_venn.pdf")
write.csv(getFeature(venndetail(acyrg[cont]), rlist=acyrd[cont], 
                     userowname=F, gind=rep(1,5)),
          file="Acyr_venndetail.csv")

plot(venndetail(polarrg[cont]))
dev.print(pdf, file="Polar_venn.pdf")
write.csv(getFeature(venndetail(polarrg[cont]), rlist=polarrd[cont], 
                     userowname=F, gind=rep(1,5)),
          file="Polar_venndetail.csv")

# ------------------------------------------------------------------------------
# 9. PLS-DA Analysis
# Methods: "Multivariate partial least squares-discriminant analysis (PLS-DA) 
#          was performed using mixOmics, with scaled and centered metabolite levels"
# ------------------------------------------------------------------------------

# Untargeted PLS-DA
plun <- make.PLSDA.model(metau$Group, scale(dd$Untargeted))

# Visualize with sample names
plotIndiv(plun$plsda, ellipse=F, style="ggplot2", title="Untargeted",
          group=metau$Group, legend=T, ind.names=T)
dev.print(pdf, file="untarget_plsda_name.pdf")

# Visualize without sample names
plotIndiv(plun$plsda, ellipse=F, style="ggplot2", title="Untargeted",
          group=metau$Group, legend=T, ind.names=F)
dev.print(pdf, file="untarget_plsda_wo_name.pdf")

# Visualize with ellipses
plotIndiv(plun$plsda, ellipse=T, style="ggplot2", title="Untargeted",
          group=metau$Group, legend=T, ind.names=F)
dev.print(pdf, file="untarget_plsda.pdf")

# ------------------------------------------------------------------------------
# 10. VIP Score Analysis
# Methods: "Variable importance in projection (VIP) scores were computed to 
#          rank metabolites contributing to group separation"
# ------------------------------------------------------------------------------

# Extract and plot VIP scores for untargeted metabolomics
vipun <- plun$vip
plot.vip2(vipun)
dev.print(pdf, file="plsda_vip_untargeted.pdf")

# Acyl PLS-DA
plac <- make.PLSDA.model(metaa$Group, scale(dd$Acyl))
plotIndiv(plac$plsda, ellipse=F, style="ggplot2", title="Acyl",
          group=metaa$Group, legend=T, ind.names=T)
dev.print(pdf, file="Acyl_plsda_name.pdf")

plotIndiv(plac$plsda, ellipse=F, style="ggplot2", title="Acyl",
          group=metaa$Group, legend=T, ind.names=F)
dev.print(pdf, file="Acyl_plsda_wo_name.pdf")

plotIndiv(plac$plsda, ellipse=T, style="ggplot2", title="Acyl",
          group=metaa$Group, legend=T, ind.names=F)
dev.print(pdf, file="Acyl_plsda.pdf")

vipac <- plac$vip
plot.vip2(vipac)
dev.print(pdf, file="plsda_vip_Acyl.pdf")

# SCN PLS-DA
plsc <- make.PLSDA.model(metas$Group, scale(dd$SCN))
plotIndiv(plsc$plsda, ellipse=F, style="ggplot2", title="SCN",
          group=metas$Group, legend=T, ind.names=T)
dev.print(pdf, file="SCN_plsda_name.pdf")

plotIndiv(plsc$plsda, ellipse=F, style="ggplot2", title="SCN",
          group=metas$Group, legend=T, ind.names=F)
dev.print(pdf, file="SCN_plsda_wo_name.pdf")

plotIndiv(plsc$plsda, ellipse=T, style="ggplot2", title="SCN",
          group=metas$Group, legend=T, ind.names=F)
dev.print(pdf, file="SCN_plsda.pdf")

vipsc <- plsc$vip
plot.vip2(vipsc)
dev.print(pdf, file="plsda_vip_SCN.pdf")

# Polar PLS-DA
plsp <- make.PLSDA.model(metap$Group, scale(dd$Polar))
plotIndiv(plsp$plsda, ellipse=F, style="ggplot2", title="Polar",
          group=metap$Group, legend=T, ind.names=T)
dev.print(pdf, file="Polar_plsda_name.pdf")

plotIndiv(plsp$plsda, ellipse=F, style="ggplot2", title="Polar",
          group=metap$Group, legend=T, ind.names=F)
dev.print(pdf, file="Polar_plsda_wo_name.pdf")

plotIndiv(plsp$plsda, ellipse=T, style="ggplot2", title="Polar",
          group=metap$Group, legend=T, ind.names=F)
dev.print(pdf, file="Polar_plsda.pdf")

vipp <- plsp$vip
plot.vip2(vipp)
dev.print(pdf, file="plsda_vip_Polar.pdf")

# ------------------------------------------------------------------------------
# 11. Pairwise PLS-DA for Specific Comparisons
# Function to perform PLS-DA for specific group comparisons
# ------------------------------------------------------------------------------
plsdacom <- function(meta, group, data, cv=F, ncomp=4, log=T){
  mm <- subset(meta, Group %in% group)
  id <- rownames(mm)
  dd <- data[id, ]
  plsp <- make.PLSDA.model(mm$Group, scale(dd), cv=cv, ncomp=ncomp, log=log)
  return(list(plsp=plsp, group=mm$Group))
}

# Define key pairwise comparisons
# SD vs HFD, HFD vs DR, HFD vs EX, HFD vs DREX
compp <- comp[c(1, 12, 16:18), ]

# SCN pairwise PLS-DA (log-transformed data)
scnv <- list()
for(i in 1:nrow(compp)){
  group <- compp[i, ]
  tmp <- plsdacom(metas, group=group, dd$SCN, cv=F, ncomp=4, log=T)
  nam <- paste0("SCN_", group[1], "vs", group[2])
  scnv[[nam]] <- tmp
}

# Acyl pairwise PLS-DA (non-log-transformed)
acyv <- list()
for(i in 1:nrow(compp)){
  group <- compp[i, ]
  tmp <- plsdacom(metaa, group=group, dd$Acyl, cv=F, ncomp=4, log=F)
  nam <- paste0("Acyl_", group[1], "vs", group[2])
  acyv[[nam]] <- tmp
}

# Polar pairwise PLS-DA (non-log-transformed)
ploarv <- list()
for(i in 1:nrow(compp)){
  group <- compp[i, ]
  tmp <- plsdacom(metap, group=group, dd$Polar, cv=F, ncomp=4, log=F)
  nam <- paste0("Polar_", group[1], "vs", group[2])
  ploarv[[nam]] <- tmp
}

# Untargeted pairwise PLS-DA (log-transformed)
untv <- list()
for(i in 1:nrow(compp)){
  group <- compp[i, ]
  tmp <- plsdacom(metau, group=group, dd$Untargeted, cv=F, ncomp=4, log=T)
  nam <- paste0("Untarget_", group[1], "vs", group[2])
  untv[[nam]] <- tmp
}

# ------------------------------------------------------------------------------
# 12. Export VIP Scores with Cutoff
# Methods: "Metabolites with VIP scores>1 were considered significant"
# ------------------------------------------------------------------------------
write.vip2 <- function(data, name, cutoff=1){
  data <- as.data.frame(data)
  data$species <- rownames(data)
  dd <- data %>%
    filter(comp1 >= cutoff) %>%
    arrange(desc(comp1))
  write.csv(dd, file=paste(name, "_VIP_cut.csv", sep=""), row.names=F)
}

# Save VIP results for selected SCN comparisons
for(i in names(scnv)[c(1, 4)]){
  tmp <- scnv[[i]]
  
  # Plot with sample names
  plotIndiv(tmp$plsp$plsda, ellipse=F, style="ggplot2", title="SCN",
            group=tmp$group, legend=T, ind.names=T)
  ggsave(file=paste0(i, "_plsda_name.pdf"), width=5.5, height=4.2)
  
  # Plot with ellipses
  plotIndiv(tmp$plsp$plsda, ellipse=T, style="ggplot2", title="SCN",
            group=tmp$group, legend=T, ind.names=F)
  ggsave(file=paste0(i, "_plsda_wo_name.pdf"), width=5.5, height=4.2)
  
  # VIP plot
  vipt <- tmp$plsp$vip
  plot.vip2(vipt)
  ggsave(file=paste0(i, "_plsda_vip.pdf"), width=12, height=5)
  write.vip2(vipt, name=i, cutoff=1)
}

# Save VIP results for all Acyl comparisons
for(i in names(acyv)){
  tmp <- acyv[[i]]
  plotIndiv(tmp$plsp$plsda, ellipse=F, style="ggplot2", title="Acyl",
            group=tmp$group, legend=T, ind.names=T)
  ggsave(file=paste0(i, "_plsda_name.pdf"), width=5.5, height=4.2)
  
  plotIndiv(tmp$plsp$plsda, ellipse=T, style="ggplot2", title="Acyl",
            group=tmp$group, legend=T, ind.names=F)
  ggsave(file=paste0(i, "_plsda_wo_name.pdf"), width=5.5, height=4.2)
  
  vipt <- tmp$plsp$vip
  plot.vip2(vipt, fc=tmp$plsp$fc)
  ggsave(file=paste0(i, "_plsda_vip_fc.pdf"), width=5.77, height=3.5)
  write.vip2(vipt, name=i, cutoff=1)
}

# Save VIP results for all Polar comparisons
for(i in names(ploarv)){
  tmp <- ploarv[[i]]
  plotIndiv(tmp$plsp$plsda, ellipse=F, style="ggplot2", title="Polar",
            group=tmp$group, legend=T, ind.names=T)
  ggsave(file=paste0(i, "_plsda_name.pdf"), width=5.5, height=4.2)
  
  plotIndiv(tmp$plsp$plsda, ellipse=T, style="ggplot2", title="Polar",
            group=tmp$group, legend=T, ind.names=F)
  ggsave(file=paste0(i, "_plsda_wo_name.pdf"), width=5.5, height=4.2)
  
  vipt <- tmp$plsp$vip
  plot.vip2(vipt, fc=tmp$plsp$fc)
  ggsave(file=paste0(i, "_plsda_vip_fc.pdf"), width=5.77, height=4)
  write.vip2(vipt, name=i, cutoff=1)
}

# Save VIP results for all Untargeted comparisons
for(i in names(untv)){
  tmp <- untv[[i]]
  plotIndiv(tmp$plsp$plsda, ellipse=F, style="ggplot2", title="Untargeted",
            group=tmp$group, legend=T, ind.names=T)
  ggsave(file=paste0(i, "_plsda_name.pdf"), width=5.5, height=4.2)
  
  plotIndiv(tmp$plsp$plsda, ellipse=T, style="ggplot2", title="Untargeted",
            group=tmp$group, legend=T, ind.names=F)
  ggsave(file=paste0(i, "_plsda_wo_name.pdf"), width=5.5, height=4.2)
  
  vipt <- tmp$plsp$vip
  plot.vip2(vipt, fc=tmp$plsp$fc)
  ggsave(file=paste0(i, "_plsda_vip_fc.pdf"), width=5.77, height=5.88)
  write.vip2(vipt, name=i, cutoff=1)
}

# ------------------------------------------------------------------------------
# 13. Extract Metabolites with VIP > 1 for Pathway Analysis
# Methods: "Metabolites with VIP scores>1 were considered significant and 
#          visualized using dot plots"
# ------------------------------------------------------------------------------

# Extract metabolites with VIP > 1 for each dataset
untvp <- lapply(untv, function(x) 
  rownames(subset(as.data.frame(x$plsp$vip), comp1 > 1)))
names(untvp) <- sub('Untarget_', '', names(untvp))

ploarvp <- lapply(ploarv, function(x) 
  rownames(subset(as.data.frame(x$plsp$vip), comp1 > 1)))
names(ploarvp) <- sub('Polar_', '', names(ploarvp))

acyvp <- lapply(acyv, function(x) 
  rownames(subset(as.data.frame(x$plsp$vip), comp1 > 1)))
names(acyvp) <- sub('Acyl_', '', names(acyvp))

# ------------------------------------------------------------------------------
# 14. Venn Diagram Analysis of High-VIP Metabolites
# ------------------------------------------------------------------------------
library(VennDetail)

# Create Venn diagrams for comparisons 2-5 (excluding first comparison)
untvpr <- venndetail(untvp[2:5])
ploarvpr <- venndetail(ploarvp[2:5])
acyvpr <- venndetail(acyvp[2:5])

plot(untvpr)
dev.print(pdf, file="Untargeted_venn_vip.pdf")

plot(ploarvpr)
dev.print(pdf, file="Polar_venn_vip.pdf")

plot(acyvpr)
dev.print(pdf, file="Acyl_venn_vip.pdf")

# Save Venn diagram details
write.csv(result(untvpr), file="Untargeted_venndetail_wo18wk.csv")
write.csv(result(ploarvpr), file="Polar_venndetail_wo18wk.csv")
write.csv(result(acyvpr), file="Acyl_venndetail_wo18wk.csv")

# ------------------------------------------------------------------------------
# 15. VIP Plots with Fold-Change Visualization
# Methods: "Log2Fold-Changes were calculated between groups to indicate the 
#          direction and relative magnitude of metabolite differences"
# ------------------------------------------------------------------------------

# Untargeted VIP plots with fold-change
for(i in names(untv)[c(2, 3, 5)]){
  tmp <- untv[[i]]
  vipt <- tmp$plsp$vip
  nm <- sub('Untarget_', '', i)
  fc <- as.data.frame(untf[, c("met", nm)])
  nm <- unlist(strsplit(nm, "vs"))
  fc$dir <- ifelse(fc[, 2] > 0, nm[1], nm[2])
  rownames(fc) <- fc[, 1]
  fc <- fc[rownames(vipt), ]
  colnames(fc)[2] <- "log2FC"
  
  plot.vip2(vipt, fc=fc) + 
    scale_color_manual(values=c("grey", "grey")) + 
    ylim(c(-2.5, 2.5)) + 
    geom_hline(yintercept=0)
  ggsave(file=paste0(i, "_plsda_vip_fc.pdf"), width=14, height=9)
}

# Acyl VIP plots with fold-change
for(i in names(acyv)){
  tmp <- acyv[[i]]
  vipt <- tmp$plsp$vip
  nm <- sub('Acyl_', '', i)
  fc <- as.data.frame(acyf[, c("met", nm)])
  nm <- unlist(strsplit(nm, "vs"))
  fc$dir <- ifelse(fc[, 2] > 0, nm[1], nm[2])
  rownames(fc) <- fc[, 1]
  fc <- fc[rownames(vipt), ]
  colnames(fc)[2] <- "log2FC"
  
  plot.vip2(vipt, fc=fc) + 
    scale_color_manual(values=c("grey", "grey")) + 
    ylim(c(-3, 2.2)) + 
    geom_hline(yintercept=0)
  ggsave(file=paste0(i, "_plsda_vip_fc.pdf"), width=5, height=3)
}

# Polar VIP plots with fold-change
for(i in names(ploarv)){
  tmp <- ploarv[[i]]
  vipt <- tmp$plsp$vip
  nm <- sub('Polar_', '', i)
  fc <- as.data.frame(polarf[, c("met", nm)])
  nm <- unlist(strsplit(nm, "vs"))
  fc$dir <- ifelse(fc[, 2] > 0, nm[1], nm[2])
  rownames(fc) <- fc[, 1]
  fc <- fc[rownames(vipt), ]
  colnames(fc)[2] <- "log2FC"
  
  plot.vip2(vipt, fc=fc) + 
    scale_color_manual(values=c("grey", "grey")) + 
    ylim(c(-2, 2.4)) + 
    geom_hline(yintercept=0)
  ggsave(file=paste0(i, "_plsda_vip_fc.pdf"), width=6, height=6)
}

# ------------------------------------------------------------------------------
# 16. KEGG Pathway Enrichment Analysis
# Methods: "Pathway enrichment analysis was performed on metabolites with 
#          VIP>1 using KEGG annotations via the richR package, mapped through 
#          Human Metabolome Database (HMDB)"
# ------------------------------------------------------------------------------

# Load pathway annotation
colnames(path)[1] <- "metabolite"
path <- path[, c(1, 3)]

library(richR)

# Perform KEGG enrichment for metabolites with VIP > 1
untvps <- lapply(untvp, function(x) enrich(x, path))
acyvps <- lapply(acyvp, function(x) enrich(x, path))
ploarvps <- lapply(ploarvp, function(x) enrich(x, path))

# Visualize enrichment results
comparedot(compareResult(untvps, include.all=T), include.all=T, usePadj=F)
dev.print(pdf, file="Untarget_KEGG_enrich.pdf")

comparedot(compareResult(acyvps, include.all=T), include.all=T, usePadj=F)
dev.print(pdf, file="Acy_KEGG_enrich.pdf")

comparedot(compareResult(ploarvps, include.all=T), include.all=T, usePadj=F)
dev.print(pdf, file="Polar_KEGG_enrich.pdf")

# Save enrichment results
write.csv(compareResult(ploarvps, include.all=T), file="Polar_KEGG_enrich.csv")
write.csv(compareResult(acyvps, include.all=T), file="Acy_KEGG_enrich.csv")
write.csv(compareResult(untvps, include.all=T), file="Untarget_KEGG_enrich.csv")
