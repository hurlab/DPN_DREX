#' ============================================================================
#' Pathway-Based Transcriptome-Metabolome Integration
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
#'   Integrates scRNA-seq, metabolomics, and fluxomics data through KEGG pathway
#   annotations. Maps DEGs to enriched KEGG pathways, cross-references with
#   metabolite-pathway mappings, and constructs gene-metabolite-pathway network
#   graphs per cell type.
#'
#' Input:
#'   sc_integration_predata.rdata, sc_deg_bulk_all.csv, step_path.csv,
#   KEGG mouse annotations (via richR)
#'
#' Output:
#'   Network PDFs (ggnet2/tidygraph) per cell type (mySC, ImmSC, nmSC,
#   Pericytes)
#'
#' Related figures:
#'   Fig. 5A-D
#'
#' ============================================================================
setwd("Desktop/Project/Step/integration/")
load("sc_integration_predata.rdata")
library(richR)
library(tidyverse)
library(richR)
mmko<-buildAnnot(species = "mouse",keytype = "SYMBOL",anntype = "KEGG",builtin = F)
gene<-lapply(scd, function(x)rownames(x))
resk<-lapply(gene, function(x)richKEGG(x,mmko,builtin = F))
mpath<-read.csv("step_path.csv",row.names = 1)
reskk<-lapply(resk, function(x)detail(x))
mgene<-lapply(reskk, function(x)left_join(mpath,x,by=c("pathway"="Annot")))
mgene<-lapply(mgene, function(x)na.omit(x))              
mg<-lapply(mgene, function(x)unique(x[,c("X","pathway","GeneID")]))
mmg<-lapply(mg, function(x)x%>%mutate(val=1))
mga<-do.call(rbind,mmg)
mga$cell<-sub('\\..*','',rownames(mga))
generate_net<-function(x){
  g<-graph_from_data_frame(x[,c(1,6,4)])
  meta<-unique(x$X)
  g<-simplify(g)
  V(g)$color<-"darkorange"
  V(g)$color<-ifelse(V(g)$name%in%meta,"darkorange","skyblue")
  ggnet2(g,node.color = V(g)$color,node.alpha = 0.7,
         edge.alpha = 0.8,edge.size = 0.1,node.size = 3)+geom_text_repel(aes(label=V(g)$name),size=3,max.overlaps=1000)
}

generate_net(subset(mga,cell%in%c("mySC","ImmSC","nmSC")))

library(GGally)
library(ggrepel)
generate_net2<-function(x){
  g<-graph_from_data_frame(unique(x[,c(1,3,4)]))
  meta<-unique(x$X)
  g<-simplify(g)
  V(g)$color<-"darkorange"
  V(g)$color<-ifelse(V(g)$name%in%meta,"darkorange","skyblue")
  ggnet2(g,node.color = V(g)$color,node.alpha = 0.7,
         edge.alpha = 0.8,edge.size = 0.1,node.size = 3)+geom_text_repel(aes(label=V(g)$name),size=3)
}
####
generate_net(mmg$mySC)
dev.print(pdf,file="mysc_gene_meta_network.pdf")
####
generate_net(mmg$ImmSC)
dev.print(pdf,file="Immsc_gene_meta_network.pdf")
####
generate_net(mmg$nmSC)
dev.print(pdf,file="nmSC_gene_meta_network.pdf")
####
generate_net(mmg$Pericytes)
dev.print(pdf,file="Pericytes_gene_meta_network.pdf")
##############
sc<-read.csv("sc_deg_bulk_all.csv",row.names = 1)
scs<-subset(sc,cell%in%c("mySC","nmSC","ImmSC"))
scs$group<-paste(scs$cell,scs$comparision,sep="_")
scs$GeneID<-paste0(scs$cell,"_",scs$gene)
####
mga$Group<-paste0(mga$cell,"_",mga$GeneID)
d<-left_join(mga,scs,by=c("Group"="GeneID"))
dd<-na.omit(d)
##only use against HFD
###############################
dds<-subset(dd,comparision%in%c("HFD_vs_SD","DR_vs_HFD","EX_vs_HFD","DREX_vs_HFD"))
###############
library(tidygraph)
data<-dds
####
nodes_X <- data %>% select(X) %>% distinct() %>% mutate(type = "metabolite")
nodes_pathway <- data %>% select(pathway) %>% distinct() %>% mutate(type = "pathway")
nodes_group <- data %>% select(Group, comparision) %>% distinct() %>% 
  rename(name = Group) %>% mutate(type = "gene")

# Combine all nodes
nodes <- bind_rows(
  nodes_X %>% rename(name = X),
  nodes_pathway %>% rename(name = pathway),
  nodes_group
) %>% 
  mutate(id = row_number())

# Step 2: Create edges dataframe
# Create edges from X to pathway
edges_X_pathway <- data %>% 
  select(X, pathway, val) %>% 
  distinct() %>%
  rename(from = X, to = pathway)

# Create edges from pathway to Group
edges_pathway_group <- data %>% 
  select(pathway, Group, val) %>% 
  distinct() %>%
  rename(from = pathway, to = Group)

# Combine all edges
edges <- bind_rows(
  edges_X_pathway,
  edges_pathway_group
) %>%
  # Join with nodes to get the IDs
  left_join(nodes %>% select(name, id), by = c("from" = "name")) %>%
  rename(from_id = id) %>%
  left_join(nodes %>% select(name, id), by = c("to" = "name")) %>%
  rename(to_id = id) %>%
  select(from_id, to_id, val)

# Step 3: Create the graph
graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = val), alpha = 0.5) +
  geom_node_point(aes(color = comparision, shape = type, size = type),alpha=0.75) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_shape_manual(values = c("metabolite" = 15, "pathway" = 17, "gene" = 19)) +
  scale_size_manual(values = c("metabolite" = 6, "pathway" = 4, "gene" = 3)) +
  scale_color_manual(values=distcolor[1:11]) +
  scale_edge_width(range = c(0.2, 0.5), guide = "none") +
  theme_graph() +
  labs(color = "Comparison",
       shape = "Node Type",
       size = "Node Type")
ggsave(file="network_against_HFD_new.pdf")
##############
data<-dd
ggsave(file="network_against_all.pdf")
#####################
# Enhanced Network Visualization Function
# This function creates a multi-level network visualization connecting metabolites, pathways, and genes
# with options for filtering, layout, and appearance customization

# Enhanced Network Visualization Function
# This function creates a multi-level network visualization connecting metabolites, pathways, and genes
# with options for filtering, layout, and appearance customization

create_network_visualization <- function(data, 
                                         selected_comparisons = NULL,
                                         layout_type = "fr", 
                                         node_size_factor = 1,
                                         text_size = 3,
                                         edge_alpha = 0.5,
                                         node_alpha = 0.75,
                                         custom_colors = NULL) {
  
  require(dplyr)
  require(tidygraph)
  require(ggraph)
  
  # Filter by comparisons if specified
  if (!is.null(selected_comparisons)) {
    data <- data %>% filter(comparision %in% selected_comparisons)
  }
  
  # Create node dataframes
  nodes_X <- data %>% 
    select(X) %>% 
    distinct() %>% 
    mutate(type = "metabolite")
  
  nodes_pathway <- data %>% 
    select(pathway) %>% 
    distinct() %>% 
    mutate(type = "pathway")
  
  nodes_group <- data %>% 
    select(Group, comparision) %>% 
    distinct() %>% 
    rename(name = Group) %>% 
    mutate(type = "gene")
  
  # Handle the issue of nodes appearing in multiple comparisons
  # Create a helper function to collapse comparisons for the same node
  collapse_comparisons <- function(df) {
    df %>%
      group_by(name, type) %>%
      summarize(comparision = paste(unique(comparision), collapse = ";"),
                .groups = "drop")
  }
  
  nodes_group <- collapse_comparisons(nodes_group)
  
  # Combine all nodes
  nodes <- bind_rows(
    nodes_X %>% rename(name = X) %>% mutate(comparision = NA_character_),
    nodes_pathway %>% rename(name = pathway) %>% mutate(comparision = NA_character_),
    nodes_group
  ) %>% 
    mutate(id = row_number())
  
  # Create edges dataframe
  # Create edges from X to pathway
  edges_X_pathway <- data %>% 
    select(X, pathway, val, comparision) %>% 
    distinct() %>%
    rename(from = X, to = pathway)
  
  # Create edges from pathway to Group
  edges_pathway_group <- data %>% 
    select(pathway, Group, val, comparision) %>% 
    distinct() %>%
    rename(from = pathway, to = Group)
  
  # Combine all edges
  edges <- bind_rows(
    edges_X_pathway,
    edges_pathway_group
  ) %>%
    # Join with nodes to get the IDs
    left_join(nodes %>% select(name, id), by = c("from" = "name")) %>%
    rename(from_id = id) %>%
    left_join(nodes %>% select(name, id), by = c("to" = "name")) %>%
    rename(to_id = id) %>%
    select(from_id, to_id, val, comparision)
  
  # Step 3: Create the graph
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  # Define colors if not provided
  if (is.null(custom_colors)) {
    # Set of visually distinct colors as default
    default_colors <- c(
      "#377EB8", "#4DAF4A", "#984EA3", 
      "#FF7F00", "#FFFF33", "#A65628", "#F781BF", 
      "#999999", "#66C2A5", "#FC8D62"
    )
    custom_colors <- default_colors
  }
  
  # Create a function to extract the first comparison for each node
  get_node_color <- function(comparison_str) {
    if (is.na(comparison_str)) return(NA)
    # Get the first comparison if multiple are present
    first_comp <- strsplit(comparison_str, ";")[[1]][1]
    return(first_comp)
  }
  
  # Function to extract the first comparison for coloring
  get_node_color <- function(comparison_str) {
    if (is.na(comparison_str)) return(NA)
    # Get the first comparison if multiple are present
    first_comp <- strsplit(comparison_str, ";")[[1]][1]
    return(first_comp)
  }
  
  # Create the visualization
  # Determine if we're showing a single comparison or multiple
  single_comparison <- length(unique(na.omit(edges$comparision))) == 1
  
  # Create the base plot
  p <- ggraph(graph, layout = layout_type)
  
  # Add edges - different styling for single vs multiple comparisons
  if (single_comparison) {
    p <- p + geom_edge_link(aes(width = val), 
                            alpha = edge_alpha,
                            color = "grey50")
  } else {
    p <- p + geom_edge_link(aes(width = val, color = comparision), 
                            alpha = edge_alpha,
                            show.legend = TRUE)
  }
  
  # Add nodes
  p <- p + geom_node_point(aes(color = comparision, 
                               shape = type, 
                               size = type),
                           alpha = node_alpha) +
    geom_node_text(aes(label = name), 
                   repel = TRUE, 
                   size = text_size) +
    scale_shape_manual(values = c("metabolite" = 15, 
                                  "pathway" = 17, 
                                  "gene" = 19)) +
    scale_size_manual(values = c("metabolite" = 6 * node_size_factor, 
                                 "pathway" = 4 * node_size_factor, 
                                 "gene" = 3 * node_size_factor)) +
    scale_edge_width(range = c(0.2, 0.8), guide = "none") +
    theme_graph() +
    labs(color = "Comparison",
         shape = "Node Type",
         size = "Node Type",
         edge_color = "Comparison")
  
  # Apply custom colors if there are enough for the number of comparisons
  unique_comparisons <- na.omit(unique(edges$comparision))
  if (length(custom_colors) >= length(unique_comparisons)) {
    p <- p + 
      scale_color_manual(values = custom_colors, na.value = "grey50")
  }
  
  # Modify legend title based on single/multiple comparisons
  if (single_comparison) {
    p <- p + labs(color = "Node Type",
                  shape = "Node Type",
                  size = "Node Type")
  }
  
  return(p)
}

# Example usage:
# Load required packages
library(dplyr)
library(tidygraph)
library(ggraph)

# Enhanced Network Visualization Function with log2FoldChange support
# This function creates a multi-level network visualization connecting metabolites, pathways, and genes
# with coloring based on up/down regulation (log2FoldChange)

create_network_visualization_color <- function(data, 
                                         selected_comparisons = NULL,
                                         layout_type = "fr", 
                                         node_size_factor = 1,
                                         text_size = 3,
                                         edge_alpha = 0.5,
                                         node_alpha = 0.75,
                                         up_color = "#E41A1C",     # Red for up-regulated
                                         down_color = "#377EB8",   # Blue for down-regulated
                                         custom_colors = NULL) {
  
  require(dplyr)
  require(tidygraph)
  require(ggraph)
  
  # Filter by comparisons if specified
  if (!is.null(selected_comparisons)) {
    data <- data %>% filter(comparision %in% selected_comparisons)
  }
  
  # Check if log2FoldChange column exists
  has_log2fc <- "log2FoldChange" %in% colnames(data)
  
  # Create node dataframes with regulation status if available
  if (has_log2fc) {
    nodes_X <- data %>% 
      select(X, log2FoldChange, comparision) %>% 
      distinct() %>% 
      mutate(type = "metabolite")
    
    nodes_pathway <- data %>% 
      select(pathway, log2FoldChange, comparision) %>% 
      distinct() %>% 
      mutate(type = "pathway")
    
    nodes_group <- data %>% 
      select(Group, log2FoldChange, comparision) %>% 
      distinct() %>% 
      rename(name = Group) %>% 
      mutate(type = "gene")
  } else {
    nodes_X <- data %>% 
      select(X, comparision) %>% 
      distinct() %>% 
      mutate(type = "metabolite")
    
    nodes_pathway <- data %>% 
      select(pathway, comparision) %>% 
      distinct() %>% 
      mutate(type = "pathway")
    
    nodes_group <- data %>% 
      select(Group, comparision) %>% 
      distinct() %>% 
      rename(name = Group) %>% 
      mutate(type = "gene")
  }
  
  # Handle the issue of nodes appearing in multiple comparisons
  # Create a helper function to collapse multiple rows for the same node
  collapse_node_data <- function(df) {
    if (has_log2fc) {
      df %>%
        group_by(name, type) %>%
        summarize(
          comparision = paste(unique(comparision), collapse = ";"),
          log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
          regulation = ifelse(mean(log2FoldChange, na.rm = TRUE) > 0, "Up", "Down"),
          .groups = "drop"
        )
    } else {
      df %>%
        group_by(name, type) %>%
        summarize(
          comparision = paste(unique(comparision), collapse = ";"),
          .groups = "drop"
        )
    }
  }
  
  # Apply the collapse function
  nodes_X <- nodes_X %>% rename(name = X) %>% collapse_node_data()
  nodes_pathway <- nodes_pathway %>% rename(name = pathway) %>% collapse_node_data()
  nodes_group <- collapse_node_data(nodes_group)
  
  # Combine all nodes
  nodes <- bind_rows(
    nodes_X,
    nodes_pathway,
    nodes_group
  ) %>% 
    mutate(id = row_number())
  
  # Create edges dataframe
  # Create edges from X to pathway
  edges_X_pathway <- data %>% 
    select(X, pathway, val, comparision) %>% 
    distinct() %>%
    rename(from = X, to = pathway)
  
  # Create edges from pathway to Group
  edges_pathway_group <- data %>% 
    select(pathway, Group, val, comparision) %>% 
    distinct() %>%
    rename(from = pathway, to = Group)
  
  # Combine all edges
  edges <- bind_rows(
    edges_X_pathway,
    edges_pathway_group
  ) %>%
    # Join with nodes to get the IDs
    left_join(nodes %>% select(name, id), by = c("from" = "name")) %>%
    rename(from_id = id) %>%
    left_join(nodes %>% select(name, id), by = c("to" = "name")) %>%
    rename(to_id = id) %>%
    select(from_id, to_id, val, comparision)
  
  # Create the graph
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  # Determine if we're showing a single comparison or multiple
  single_comparison <- length(unique(na.omit(edges$comparision))) == 1
  
  # Create the base plot
  p <- ggraph(graph, layout = layout_type)
  
  # Add edges - different styling for single vs multiple comparisons
  if (single_comparison) {
    p <- p + geom_edge_link(aes(width = val), 
                            alpha = edge_alpha,
                            color = "grey50")
  } else {
    p <- p + geom_edge_link(aes(width = val, color = comparision), 
                            alpha = edge_alpha,
                            show.legend = TRUE)
  }
  
  # Add nodes - color by regulation if log2FoldChange is available
  if (has_log2fc) {
    p <- p + geom_node_point(aes(color = regulation, 
                                 shape = type, 
                                 size = type),
                             alpha = node_alpha) +
      scale_color_manual(values = c("Up" = up_color, "Down" = down_color))
  } else {
    # Fall back to coloring by comparison
    p <- p + geom_node_point(aes(color = comparision, 
                                 shape = type, 
                                 size = type),
                             alpha = node_alpha)
    
    # Apply custom colors if there are enough for the number of comparisons
    if (!is.null(custom_colors)) {
      unique_comparisons <- na.omit(unique(edges$comparision))
      if (length(custom_colors) >= length(unique_comparisons)) {
        p <- p + scale_color_manual(values = custom_colors, na.value = "grey50")
      }
    }
  }
  
  # Add node labels and styling
  p <- p + 
    geom_node_text(aes(label = name), 
                   repel = TRUE, 
                   size = text_size) +
    scale_shape_manual(values = c("metabolite" = 15, 
                                  "pathway" = 17, 
                                  "gene" = 19)) +
    scale_size_manual(values = c("metabolite" = 6 * node_size_factor, 
                                 "pathway" = 4 * node_size_factor, 
                                 "gene" = 3 * node_size_factor)) +
    scale_edge_width(range = c(0.2, 0.8), guide = "none") +
    theme_graph()
  
  # Set appropriate labels based on what's being shown
  if (has_log2fc) {
    p <- p + labs(color = "Regulation",
                  shape = "Node Type",
                  size = "Node Type")
  } else if (single_comparison) {
    p <- p + labs(color = "Node Type",
                  shape = "Node Type",
                  size = "Node Type")
  } else {
    p <- p + labs(color = "Comparison",
                  shape = "Node Type",
                  size = "Node Type",
                  edge_color = "Comparison")
  }
  
  return(p)
}

# Example usage:
# Load required packages
# library(dplyr)
# library(tidygraph)
# library(ggraph)

# Usage with log2FoldChange data:
# network_plot <- create_network_visualization(dds, 
#                                            up_color = "#D55E00",   # Customizable colors for up-regulation
#                                            down_color = "#0072B2") # Customizable colors for down-regulation
# print(network_plot)

# For a specific comparison only:
# dr_plot <- create_network_visualization(dds, selected_comparisons = "DR_vs_HFD")
# print(dr_plot)

# Save the plot
# ggsave("network_visualization.png", network_plot, width = 12, height = 10, dpi = 300)