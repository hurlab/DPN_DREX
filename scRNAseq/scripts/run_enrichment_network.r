#' ============================================================================
#' Enrichment Network Visualization
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
#'   Creates enrichment network plots from KEGG pathway analysis results. Trims long
#   pathway terms (70-character limit), filters by adjusted p-value, and generates
#   network graph visualizations linking enriched pathways.
#'
#' Input:
#'   KEGG enrichment CSVs from major and standard enrichment analyses
#'
#' Output:
#'   Network plots (PDF), enrichment tables, padj-filtered subsets
#'
#' Related figures:
#'   Fig. 3C-F
#'
#' ============================================================================
################################################################################
# Initialize the workspace
################################################################################
# Set working directory to the script's location (if using RStudio)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  setwd(dirname(script_path))
  base_dir <- dirname(script_path)
}

# Install and load required packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, gridExtra, richR, ggrepel)

################################################################################
# Helper: Trim long terms (70 characters)
################################################################################
trim_term <- function(term) {
  if (nchar(term) > 70) {
    words        <- strsplit(term, " ")[[1]]
    out          <- ""
    char_counter <- 0
    for (w in words) {
      sep_len <- ifelse(out == "", 0, 1)
      if (char_counter + nchar(w) + sep_len <= 70) {
        out          <- paste(out, w, sep = ifelse(out == "", "", " "))
        char_counter <- nchar(out)
      } else {
        break
      }
    }
    return(paste0(out, "*"))
  }
  term
}

################################################################################
# Directories for major‐pvalue and pvalue runs & their outputs
################################################################################
major_enrich_dir <- file.path(base_dir, "Enrichment_major_pvalue001")
major_out_dir    <- file.path(base_dir, "250515_EnrichNet-Major_JH")
pval_enrich_dir  <- file.path(base_dir, "Enrichment_Pvalue001")
pval_out_dir     <- file.path(base_dir, "250515_EnrichNet_JH")

dir.create(major_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pval_out_dir,  recursive = TRUE, showWarnings = FALSE)

################################################################################
# Function to process one enrichment directory, with optional “major_” requirement
################################################################################
process_enrichment <- function(enrich_dir, output_dir, require_major = FALSE,
                               padj.cutoff   = NULL) {
  
  # build a filename suffix if cutoff is set
  suffix <- if (!is.null(padj.cutoff)) paste0("_padj", padj.cutoff) else ""
  
  # choose your pattern as before …
  files <- list.files(
    path       = enrich_dir,
    pattern    = if (require_major) 
      "^[^.].*_KEGG_enrichment_major_p(val001)?\\.csv$"
    else
      "^[^.].*_KEGG_enrichment_p(val001)?\\.csv$",
    full.names = TRUE
  )
  
  for (file in files) {
    # ---- print file being processed ----
    message("Processing: ", basename(file), suffix)
    
    # 1) read + trim + compute LogPadj
    df <- read.csv(file, row.names = 1, stringsAsFactors = FALSE) %>%
      mutate(
        TrimmedTerm = sapply(Term, trim_term),
        LogPadj     = -log10(Padj)
      )
    
    # 1.1) optional Padj filtering
    if (!is.null(padj.cutoff)) {
      message("  ↳ Filtering to Padj ≤ ", padj.cutoff)
      df <- df[df$Padj <= padj.cutoff, , drop = FALSE]
    }
    
    # 2) skip if too few terms
    if (nrow(df) <= 2) {
      message("  ↳ Skipped (only ", nrow(df), " term",
              if (nrow(df)==1) "" else "s", " after filtering)")
      next
    }
    
    # 2) coerce GeneID to character if needed
    if (!is.character(df$GeneID)) {
      message("  * Coercing GeneID to character in ", basename(file))
      df$GeneID <- as.character(df$GeneID)
    }
    
    # 3) trim & compute LogPadj
    df <- df %>%
      mutate(
        TrimmedTerm = sapply(Term, trim_term),
        LogPadj     = -log10(Padj)
      )
    
    # 4) safely split GeneID
    gene_names <- tryCatch({
      unique(unlist(strsplit(df$GeneID, ",")))
    }, error = function(e) {
      message("  ! Error splitting GeneID in ", basename(file), ": ", e$message)
      return(character(0))
    })
    if (length(gene_names) == 0) next
    
    # 5) build the network, embedding suffix in the Cytoscape file base name
    network_base <- paste0(tools::file_path_sans_ext(basename(file)),
                           suffix, "_ggnetwork-020")
    
    p <- richR::ggnetwork(
      df, gene_names,
      top             = nrow(df),
      weightcut       = 0.2,
      writeCyt        = TRUE,
      cytoscapeFile   = file.path(output_dir, network_base),
      cytoscapeFormat = "edgelist"
    )
    
    # 6) merge in LogPadj & shrink nodes
    p$data$LogPadj <- df$LogPadj[match(p$data$label, df$TrimmedTerm)]
    p$data$size    <- p$data$size * 0.6
    
    # 7) remap point layer
    pt <- which(sapply(p$layers, function(l) inherits(l$geom, "GeomPoint")))
    p$layers[[pt]]$mapping <- aes(x = x, y = y, colour = LogPadj, size = size)
    
    # 8) add scales & theme, suppress the “adding another scale” warnings
    p <- suppressMessages(
      p +
        scale_colour_gradient(
          name  = expression(-log[10]~"(padj)"),
          low   = "yellow",
          high  = "darkviolet",
          guide = guide_colorbar()
        ) +
        guides(
          colour = guide_colorbar()  # ensure the legend is drawn
        ) +
        scale_size_identity() +
        theme_void() +
        theme(legend.position = "right")
    )
    
    # 9) force labels
    txt <- which(sapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel")))
    p$layers[[txt]]$geom_params$max.overlaps <- Inf
    
    # 10) save PDF, again including suffix
    ggsave(
      file.path(output_dir, paste0(network_base, ".pdf")),
      plot   = p,
      width  = 8,
      height = 6
    )
  }
}

################################################################################
# Run for both directories, toggling major_ as needed
################################################################################
process_enrichment(major_enrich_dir, major_out_dir, require_major = TRUE)
process_enrichment(pval_enrich_dir,  pval_out_dir,  require_major = FALSE)
################################################################################

## Note that on Windows, the filenames maybe too long. 
pval_out_dir     <- file.path(base_dir, "EnrichOut_JH")
dir.create(pval_out_dir,  recursive = TRUE, showWarnings = FALSE)

process_enrichment(pval_enrich_dir,  pval_out_dir,  require_major = FALSE,
                   padj.cutoff = 0.05)



