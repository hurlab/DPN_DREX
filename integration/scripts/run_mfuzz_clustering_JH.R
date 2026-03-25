#' ============================================================================
#' Multi-Omics Mfuzz Clustering (Refactored)
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
#'   Refactored Mfuzz clustering with reusable run_mfuzz_batch() function, seed setting
#   for reproducibility, and multiple metabolite configurations (noSCN, noSCN_noAcyl
#   variants). Includes KEGG pathway enrichment on cluster gene members.
#'
#' Input:
#'   sc_integration_predata.rdata, metabol_integration_predata.rdata,
#   fluxomics_26wk_final_m0_avg.csv, step_path.csv, DREX Phenotyping Averages.txt
#'
#' Output:
#'   Organized in output_JH/cn{k}_{config}/ subdirectories with PDFs, CSVs,
#   RData per cell type and cluster configuration
#'
#' Related figures:
#'   Fig. 7, Fig. S13-S16, Table S10
#'
#' ============================================================================
################################################################################
## JH 2025/08/21
################################################################################
# Set working directory to the script's location (if using RStudio)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  setwd(dirname(script_path))
  base_dir <- dirname(script_path)
}

################################################################################
## Load required packages
################################################################################
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(Mfuzz, dplyr, tidyr, richR, VennDetail, ggplot2, Biobase)

## Set output directory and create if it doesn't exist
output_dir <- './output_JH'
dir.create(file.path(base_dir, output_dir), showWarnings = FALSE)
################################################################################

################################################################################
## Load Kai's data
################################################################################
load("sc_integration_predata.rdata")
load("metabol_integration_predata.rdata")
################################################################################

## ------------------------------------------------------------------
## Main runner - to replicate what Kai previously did
## ------------------------------------------------------------------
run_mfuzz_batch <- function(
    scd,
    dl,
    flux_file,
    pheno_file = "DREX Phenotyping Averages.txt",
    cn = 12,
    output_dir = "output",
    suffix = NULL,           # if NULL, will derive from flux_file stem
    pdf_width = 12,          # PDF width in inches
    pdf_height = 8,          # PDF height in inches
    seed = 123               # set seed for reproducibility
){
  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages({
    library(Biobase)
    library(Mfuzz)
    library(ggplot2)
    library(tools)  # file_path_sans_ext
  })
  
  # Ensure base output_dir exists (kept exactly as passed in)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Prepare suffix: if not provided, derive from flux_file; strip any leading underscores
  flux_stem <- file_path_sans_ext(basename(flux_file))
  suffix_final <- if (is.null(suffix)) flux_stem else gsub("^_+", "", suffix)
  
  # Subfolder format: cn<cn>_<suffix> (e.g., cn10_sum, cn10_fluxomics_... )
  out_dir <- file.path(output_dir, paste0("cn", cn, "_", suffix_final))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Reproducibility
  set.seed(seed)
  
  # 1) Load and harmonize flux + pheno + dl
  flux <- t(read.csv(flux_file, row.names = 1, check.names = FALSE))
  group_cols <- if ("mySC" %in% names(scd)) colnames(scd$mySC) else colnames(scd[[1]])
  colnames(flux) <- group_cols
  
  renamelist <- function(x){
    colnames(x) <- c("DR","DREX","EX","HFD","SD")
    x
  }
  dl  <- lapply(dl, renamelist)
  ddl <- do.call(rbind, dl)
  
  pheno <- read.delim(pheno_file, sep = "\t", row.names = 1, check.names = FALSE)
  pheno <- pheno[, c(3,5,4,2,1)]
  colnames(pheno) <- c("DR","DREX","EX","HFD","SD")
  
  mname <- c("SD","HFD","DR","EX","DREX")
  
  # 2) Per cell-type mfuzz
  do_mfuzz <- function(x, cn){
    mat <- rbind(ddl, flux, x, pheno)
    mat <- mat[, mname]
    eset <- new("ExpressionSet", exprs = as.matrix(mat))
    eset <- standardise(eset)
    m <- mestimate(eset)
    cl <- mfuzz(eset, c = cn, m = m)
    mat <- as.data.frame(mat)
    mat$cluster <- cl$cluster[rownames(mat)]
    list(m = m, eset = eset, cl = cl, mat = mat)
  }
  
  mfrow_vals <- c(2, ceiling(cn / 2))
  mf <- lapply(scd, function(ct_mat) do_mfuzz(ct_mat, cn = cn))
  names(mf) <- names(scd)
  
  save_mfuzz_plot_pdf <- function(eset, cl, file_base){
    pdf(file.path(out_dir, paste0(file_base, ".pdf")), width = pdf_width, height = pdf_height)
    mfuzz.plot(eset, cl = cl, mfrow = mfrow_vals, time.labels = mname, new.window = FALSE)
    dev.off()
  }
  
  invisible(lapply(names(mf), function(ct){
    write.csv(
      mf[[ct]]$mat,
      file = file.path(out_dir, paste0(ct, "_mfuzz_cluster_", cn, "_", suffix_final, ".csv"))
    )
    save_mfuzz_plot_pdf(
      mf[[ct]]$eset, mf[[ct]]$cl,
      paste0(ct, "_mfuzz_cluster_", cn, "_", suffix_final)
    )
  }))
  
  # 4) Summary phenotype plot (same logic as before)
  clu <- do.call(
    rbind,
    lapply(mf, function(x) x$mat[grep("C2", rownames(x$mat)), "cluster", drop = FALSE])
  )
  clu$group <- sub("\\..*", "", rownames(clu))
  clu$pheno <- sub(".*\\.", "", rownames(clu))
  clu <- clu[grep("\\(", clu$pheno), ]
  clu$cluster <- paste0("cluster", clu$cluster)
  
  p <- ggplot(clu, aes(pheno, cluster)) +
    geom_point() +
    coord_flip() +
    facet_wrap(~group) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  pdf(file.path(out_dir, paste0("pheno_cluster_avg_", cn, "_", suffix_final, ".pdf")),
      width = pdf_width, height = pdf_height)
  print(p)
  dev.off()
  
  # 5) Save objects
  rdata_path <- file.path(out_dir, paste0("mf_clu", cn, "_", suffix_final, ".rdata"))
  save(mf, clu, cn, flux_file, pheno_file, ddl, mname, seed, file = rdata_path, compress = TRUE)
  
  invisible(list(mf = mf, clu = clu, out_dir = out_dir, rdata = rdata_path))
}

## Using avg file; suffix auto-derived from flux filename if not given
run_mfuzz_batch(scd, dl,flux_file = "fluxomics_26wk_final_m0_avg.csv", pdf_height=6,
                cn = 12,output_dir = output_dir, suffix = "m0avg", seed = 123)

run_mfuzz_batch(scd, dl, flux_file = "fluxomics_26wk_final_m0_avg.csv", pdf_height=4.5,
                cn = 10,output_dir = output_dir, suffix = "m0avg", seed = 123)

# Using sum file; explicit suffix "sum" (no leading underscore needed)
run_mfuzz_batch(scd, dl, flux_file = "fluxomics_avg_sum.csv", pdf_height=4,
                cn = 12, output_dir = output_dir, suffix = "sum", seed = 123)

run_mfuzz_batch(scd, dl, flux_file = "fluxomics_avg_sum.csv", pdf_height=4,
                cn = 10, output_dir = output_dir, suffix = "sum", seed = 123)

## The above function works well as expected. However, it should be noted
## that Kai's original code did not include seed setting. Therefore, 
## I cannot recreate the exact same output. 
## In my function, setting seed option has been added to ensure
## reproducibility of results, especially for clustering.
## Kai's results and my results seem to be consistent, except for the 
## order of the clusters.






## ------------------------------------------------------------------
## Remove redundant SCN metabolites
## ------------------------------------------------------------------
## Find the row names that are common to both data frames
common_rows <- intersect(rownames(dl$SCN), rownames(dl$Untar))

## Merge the two data frames using the common row names
## We add suffixes to tell the columns apart if they have the same name
comparison_df <- merge(dl$SCN[common_rows, ], dl$Untar[common_rows, ], 
                       by = "row.names", suffixes = c("_SCN", "_Untar"))

## Note: As seen above, I noticed that SCN metabolites are already
##       included in the Untar group. 
##       Therefore, I will remove the SCN metabolites from the SCN group.
##       This is to avoid redundancy in the analysis.
dl_new <- dl
dl_new$SCN <- NULL

run_mfuzz_batch(scd, dl_new, flux_file = "fluxomics_26wk_final_m0_avg.csv", pdf_height=4.5,
                cn = 10,output_dir = output_dir, suffix = "m0avg_noSCN", seed = 123)

run_mfuzz_batch(scd, dl_new, flux_file = "fluxomics_26wk_final_m0_avg.csv", pdf_height=4.5,
                cn = 12,output_dir = output_dir, suffix = "m0avg_noSCN", seed = 123)

## Note: After reviewing the results, we still seem to have very similar 
##       clustering patterns even after removing SCN metabolites. 

## ------------------------------------------------------------------
## Adjust the metabolite. Remove acyl groups, except for carnitin,
## which needs to be moved to 'untargetted'
## ------------------------------------------------------------------
## Add dl$Acy["carnitine",] to dl_new$Untar, but "Carnitine" (note the case diff)
dl_new_2 <- dl_new
dl_new_2$Untar <- rbind(dl_new_2$Untar, Carnitine = dl$Acy["carnitine", ])
dl_new_2$Acy <- NULL

res1 <- run_mfuzz_batch(scd, dl_new_2, flux_file = "fluxomics_26wk_final_m0_avg.csv", pdf_height=4.5,
                cn = 10,output_dir = output_dir, suffix = "m0avg_noSCNAcyl", seed = 123)

res2 <- run_mfuzz_batch(scd, dl_new_2, flux_file = "fluxomics_26wk_final_m0_avg.csv", pdf_height=4.5,
                       cn = 12,output_dir = output_dir, suffix = "m0avg_noSCNAcyl", seed = 123)


## ------------------------------------------------------------------
## Perform pathway enrichment for metabolite data
## ------------------------------------------------------------------
mpath <- read.csv("./pathway/step_path.csv",
                  row.names = 1, stringsAsFactors = FALSE)

# # mySC, cluster 6
# ids_c6 <- sub('.*\\.', '', rownames(subset(res$mf$mySC$mat, cluster == 6)))
# resk_mySC_cluster6 <- enrich(ids_c6, mpath)
# 
# # mySC, cluster 1
# ids_c1 <- sub('.*\\.', '', rownames(subset(res$mf$mySC$mat, cluster == 1)))
# resk_mySC_cluster1 <- enrich(ids_c1, mpath)
# 
# 
# 
# # resk_mySC_cluster6<-enrich(sub('.*\\.','',rownames(subset(mf$mySC$mat,cluster==6))),mpath)
# 
# 
# get_cluster_ids <- function(res, cell_type, k) {
#   sub('.*\\.', '', rownames(subset(res$mf[[cell_type]]$mat, cluster == k)))
# }
# 
# enrich_by_celltype <- function(res, cell_type, mpath) {
#   ks <- sort(unique(res$mf[[cell_type]]$mat$cluster))
#   setNames(lapply(ks, function(k) enrich(get_cluster_ids(res, cell_type, k), mpath)),
#            paste0("cluster", ks))
# }
# 
# # example
# enrich_mySC <- enrich_by_celltype(res, "mySC", mpath)
# 
# clusters_df <- do.call(rbind, lapply(names(res$mf), function(ct) {
#   df <- res$mf[[ct]]$mat[, "cluster", drop = FALSE]  # column, not row
#   df$feature  <- rownames(df)
#   df$celltype <- ct
#   df[, c("celltype", "feature", "cluster")]
# }))
# rownames(clusters_df) <- NULL
# 
# # clusters_df has columns: cluster, feature, celltype

## Performe KEGG enrichment for metabolites
export_enrichment <- function(res, mpath, db="KEGG", cells=names(res$mf), w=12, h=8, top_n=20){
  suppressPackageStartupMessages(library(ggplot2))
  has_richR <- requireNamespace("richR", quietly = TRUE)
  
  out_db <- file.path(res$out_dir, db)
  dir.create(out_db, showWarnings = FALSE, recursive = TRUE)
  
  pick1 <- function(nms, ...) { hit <- intersect(c(...), nms); if (length(hit)) hit[1] else NA_character_ }
  
  # ID -> pathway name map from mpath
  mp <- as.data.frame(mpath, stringsAsFactors = FALSE)
  stopifnot("pathway" %in% names(mp))
  mp_id <- pick1(names(mp), "Annot","ID","ID.Name","PathwayID","Pathway.ID","pathID","term_id","kegg")
  if (is.na(mp_id)) stop("Could not find an ID column in mpath")
  id2name <- setNames(as.character(mp$pathway), as.character(mp[[mp_id]]))
  
  res$enrich <- setNames(vector("list", length(cells)), cells)
  
  for (ct in cells){
    ks <- sort(unique(res$mf[[ct]]$mat$cluster))
    res$enrich[[ct]] <- setNames(vector("list", length(ks)), paste0("cluster", ks))
    
    for (k in ks){
      ids <- sub('.*\\.', '', rownames(subset(res$mf[[ct]]$mat, cluster == k)))
      if (!length(ids)) next
      
      ek <- try(enrich(ids, mpath), silent = TRUE)
      if (inherits(ek, "try-error") || !NROW(ek)) next
      
      # Convert to df and replace Term with human-readable pathway names
      dd <- as.data.frame(ek, stringsAsFactors = FALSE)
      nms <- names(dd)
      idc <- pick1(nms, "Annot","ID","PathwayID","Pathway.ID","pathID","ID.Name","Term") # Term may be ID in your ek
      if (is.na(idc)) next
      if (!"Term" %in% nms) dd$Term <- dd[[idc]]
      repl <- id2name[ as.character(dd[[idc]]) ]
      dd$Term[!is.na(repl) & nzchar(repl)] <- repl[!is.na(repl) & nzchar(repl)]
      
      # ggdot sometimes expects Count
      if (!"Count" %in% names(dd) && "Significant" %in% names(dd)) dd$Count <- dd$Significant
      
      # Save the modified enrichment table
      csv_path <- file.path(out_db, sprintf("%s_c%d_%s.csv", ct, k, db))
      write.csv(dd, csv_path, row.names = FALSE)
      
      # Plot
      pdf(file.path(out_db, sprintf("%s_c%d_%s.pdf", ct, k, db)), width = w, height = h, useDingbats = FALSE)
      if (has_richR) {
        print(richR::ggdot(dd, top = top_n, usePadj = TRUE))  #, order = TRUE
      } 
      dev.off()
      
      res$enrich[[ct]][[paste0("cluster", k)]] <- list(table = ek, csv = csv_path)
    }
  }
  res
}

# after run_mfuzz_batch(...)
res1 <- export_enrichment(res1, mpath, db="KEGG", w=12, h=8)
res2 <- export_enrichment(res2, mpath, db="KEGG", w=12, h=8)





# #####
# flux<-t(read.csv("fluxomics_26wk_final_m0_avg.csv",row.names = 1))
# ###rename the group to make sure we can combine them with rbind
# colnames(flux)<-colnames(scd$mySC)
# renamelist<-function(x){
#   colnames(x)<-c("DR","DREX","EX","HFD","SD")
#   return(x)
# }
# dl<-lapply(dl, function(x)renamelist(x))
# ###
# pheno<-read.delim("DREX Phenotyping Averages.txt",sep="\t",row.names = 1)
# pheno<-pheno[,c(3,5,4,2,1)]
# colnames(pheno)<-c("DR","DREX","EX","HFD","SD")
# ##
# ddl<-do.call(rbind,dl)
# mname<-c("SD","HFD","DR","EX","DREX")
# do_mfuzz<-function(x,cn=10){
#   mname<-c("SD","HFD","DR","EX","DREX")
#   mat<- rbind(ddl,flux,x,pheno)
#   mat<-mat[,mname]
#   cat(nrow(mat),"\n")
#   eset <- new("ExpressionSet", exprs=as.matrix(mat))
#   eset <- standardise(eset)
#   m <- mestimate(eset)
#   c <- cn  # Number of clusters
#   cl <- mfuzz(eset, c=c, m=m)
#   mat<-as.data.frame(mat)
#   mat$cluster<-cl$cluster[rownames(mat)]
#   return(list(m=m,eset=eset,cl=cl,mat=mat))
# }
# 
# mf<-lapply(scd, function(x)do_mfuzz(x,cn=12))
# mfuzz.plot(mf$Endo$eset, cl=mf$Endo$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Endo_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$EndoFib1$eset, cl=mf$EndoFib1$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EndoFib1_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$EndoFib2$eset, cl=mf$EndoFib2$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EndoFib2_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$EpiFib$eset, cl=mf$EpiFib$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EpiFib_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$ImmSC$eset, cl=mf$ImmSC$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="ImmSC_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$Mac$eset, cl=mf$Mac$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Mac_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$mySC$eset, cl=mf$mySC$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="mySC_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$nmSC$eset, cl=mf$nmSC$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="nmSC_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$Pericytes$eset, cl=mf$Pericytes$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Pericytes_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$Perineurial$eset, cl=mf$Perineurial$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Perineurial_mfuzz_cluster_12.pdf")
# mfuzz.plot(mf$VSMC$eset, cl=mf$VSMC$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="VSMC_mfuzz_cluster_12.pdf")
# ########
# sapply(names(mf), function(x)write.csv(mf[[x]]$mat,file=paste0(x,"_mfuzz_cluster_12.csv")))
# clu<-do.call(rbind,lapply(mf, function(x)x$mat[grep("C2",rownames(x$mat)),"cluster",drop=F]))
# clu$group<-sub('\\..*','',rownames(clu))
# clu$pheno<-sub('.*\\.','',rownames(clu))
# clu<-clu[grep("\\(",clu$pheno),]
# clu$cluster<-paste0("cluster",clu$cluster)
# ggplot(clu,aes(pheno,cluster))+geom_point()+coord_flip()+facet_wrap(~group)+theme_light()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1))
# dev.print(pdf,file="pheno_cluster_avg_12.pdf")
# save(list=ls(),file="mf_clu12.rdata",compress=T)
# ###################
# mf<-lapply(scd, function(x)do_mfuzz(x,cn=10))
# mfuzz.plot(mf$Endo$eset, cl=mf$Endo$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Endo_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$EndoFib1$eset, cl=mf$EndoFib1$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EndoFib1_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$EndoFib2$eset, cl=mf$EndoFib2$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EndoFib2_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$EpiFib$eset, cl=mf$EpiFib$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EpiFib_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$ImmSC$eset, cl=mf$ImmSC$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="ImmSC_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$Mac$eset, cl=mf$Mac$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Mac_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$mySC$eset, cl=mf$mySC$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="mySC_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$nmSC$eset, cl=mf$nmSC$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="nmSC_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$Pericytes$eset, cl=mf$Pericytes$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Pericytes_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$Perineurial$eset, cl=mf$Perineurial$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Perineurial_mfuzz_cluster_10.pdf")
# mfuzz.plot(mf$VSMC$eset, cl=mf$VSMC$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="VSMC_mfuzz_cluster_10.pdf")
# ########
# sapply(names(mf), function(x)write.csv(mf[[x]]$mat,file=paste0(x,"_mfuzz_cluster_10.csv")))
# clu<-do.call(rbind,lapply(mf, function(x)x$mat[grep("C2",rownames(x$mat)),"cluster",drop=F]))
# clu$group<-sub('\\..*','',rownames(clu))
# clu$pheno<-sub('.*\\.','',rownames(clu))
# clu<-clu[grep("\\(",clu$pheno),]
# clu$cluster<-paste0("cluster",clu$cluster)
# ggplot(clu,aes(pheno,cluster))+geom_point()+coord_flip()+facet_wrap(~group)+theme_light()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1))
# dev.print(pdf,file="pheno_cluster_avg_10.pdf")
# save(list=ls(),file="mf_clu10.rdata",compress=T)
# #####
# flux<-t(read.csv("fluxomics_avg_sum.csv",row.names = 1))
# ###rename the group to make sure we can combine them with rbind
# colnames(flux)<-colnames(scd$mySC)
# ##########
# mf<-lapply(scd, function(x)do_mfuzz(x,cn=12))
# mfuzz.plot(mf$Endo$eset, cl=mf$Endo$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Endo_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$EndoFib1$eset, cl=mf$EndoFib1$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EndoFib1_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$EndoFib2$eset, cl=mf$EndoFib2$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EndoFib2_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$EpiFib$eset, cl=mf$EpiFib$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EpiFib_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$ImmSC$eset, cl=mf$ImmSC$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="ImmSC_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$Mac$eset, cl=mf$Mac$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Mac_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$mySC$eset, cl=mf$mySC$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="mySC_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$nmSC$eset, cl=mf$nmSC$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="nmSC_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$Pericytes$eset, cl=mf$Pericytes$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Pericytes_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$Perineurial$eset, cl=mf$Perineurial$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Perineurial_mfuzz_cluster_12_sum.pdf")
# mfuzz.plot(mf$VSMC$eset, cl=mf$VSMC$cl, mfrow=c(2,6), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="VSMC_mfuzz_cluster_12_sum.pdf")
# ########
# sapply(names(mf), function(x)write.csv(mf[[x]]$mat,file=paste0(x,"_mfuzz_cluster_12_sum.csv")))
# clu<-do.call(rbind,lapply(mf, function(x)x$mat[grep("C2",rownames(x$mat)),"cluster",drop=F]))
# clu$group<-sub('\\..*','',rownames(clu))
# clu$pheno<-sub('.*\\.','',rownames(clu))
# clu<-clu[grep("\\(",clu$pheno),]
# clu$cluster<-paste0("cluster",clu$cluster)
# ggplot(clu,aes(pheno,cluster))+geom_point()+coord_flip()+facet_wrap(~group)+theme_light()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1))
# dev.print(pdf,file="pheno_cluster_avg_12_sum.pdf")
# save(list=ls(),file="mf_clu12_sum.rdata",compress=T)
# ###################
# mf<-lapply(scd, function(x)do_mfuzz(x,cn=10))
# mfuzz.plot(mf$Endo$eset, cl=mf$Endo$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Endo_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$EndoFib1$eset, cl=mf$EndoFib1$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EndoFib1_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$EndoFib2$eset, cl=mf$EndoFib2$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EndoFib2_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$EpiFib$eset, cl=mf$EpiFib$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="EpiFib_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$ImmSC$eset, cl=mf$ImmSC$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="ImmSC_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$Mac$eset, cl=mf$Mac$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Mac_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$mySC$eset, cl=mf$mySC$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="mySC_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$nmSC$eset, cl=mf$nmSC$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="nmSC_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$Pericytes$eset, cl=mf$Pericytes$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Pericytes_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$Perineurial$eset, cl=mf$Perineurial$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="Perineurial_mfuzz_cluster_10_sum.pdf")
# mfuzz.plot(mf$VSMC$eset, cl=mf$VSMC$cl, mfrow=c(2,5), time.labels=mname, new.window=FALSE)
# dev.print(pdf,file="VSMC_mfuzz_cluster_10_sum.pdf")
# ########
# sapply(names(mf), function(x)write.csv(mf[[x]]$mat,file=paste0(x,"_mfuzz_cluster_10_sum.csv")))
# ######
# clu<-do.call(rbind,lapply(mf, function(x)x$mat[grep("C2",rownames(x$mat)),"cluster",drop=F]))
# clu$group<-sub('\\..*','',rownames(clu))
# clu$pheno<-sub('.*\\.','',rownames(clu))
# clu<-clu[grep("\\(",clu$pheno),]
# clu$cluster<-paste0("cluster",clu$cluster)
# ggplot(clu,aes(pheno,cluster))+geom_point()+coord_flip()+facet_wrap(~group)+theme_light()+theme(axis.text.x=element_text(angle=90,vjust = 0.5,hjust = 1))
# dev.print(pdf,file="pheno_cluster_avg_10_sum.pdf")
# #######
# save(list=ls(),file="mf_clu10_sum.rdata",compress=T)
# 
# ############
# resk_mySC_cluster6<-enrich(sub('.*\\.','',rownames(subset(mf$mySC$mat,cluster==6))),mpath)
# resk_mySC_cluster1<-enrich(sub('.*\\.','',rownames(subset(mf$mySC$mat,cluster==1))),mpath)
# 


