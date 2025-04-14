################################################################################
##################### Set Current Working Directory ###########################
################################################################################

if (!requireNamespace("rstudioapi", quietly = TRUE)) install.packages("rstudioapi")
library(rstudioapi)

current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(current_path)
base_dir <- dirname(current_path)

output_dir <- file.path(base_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
data_dir <- file.path(base_dir, "data")

################################################################################
##################### Packages #################################################
################################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c(
  "tidygraph", "tidyverse", "richR", "ggplot2", "dplyr", "ggraph", "igraph"
))

# # Map "Helvetica" to "Arial" on Windows
# windowsFonts(Helvetica = windowsFont("Arial"))


################################################################################
##################### Load and process the data   ##############################
################################################################################

data<-read.csv("input/gene_metabolite_pathway_group.csv",row.names=1)

# Remove any spaces at the beginning of the 'X' column
data$X <- trimws(data$X)

data2<-read.csv("input/gene_metabolite_pathway_group_MetaboliteFC_v2.csv",row.names=1)

# Filter and simplify data2
data2_simple <- data2 %>%
  filter(source != "scnf") %>%
  select(met.x, comparison, val.y) %>%
  distinct() %>%
  mutate(
    join_met = tolower(met.x),
    join_comp = tolower(comparison)
  )

# Add lowercase join keys to data
data <- data %>%
  mutate(
    join_met = tolower(X),
    join_comp = tolower(comparison)
  )

# Perform the join
data_joined <- data %>% left_join(data2_simple, by = c("join_met", "join_comp"))

# Remove temporary join columns if you don’t need them
data <- data_joined %>% select(-join_met, -join_comp, -met.x, -comparison.y)

# Change the' val.y' column name to 'met.name'
data <- data %>% rename(met.Log2FC = val.y, comparison = comparison.x)



################################################################################
# This is a modified verion of the function that creates a network visualization
# The network structure is maintained across different sets. 
################################################################################
# Date: 04/14/2025
# Version: 4

create_network_visualization_JH_4 <- function(data, 
                                                layout_type = "fr", 
                                                node_size_factor = 1,
                                                text_size = 3,
                                                edge_alpha = 0.5,
                                                node_alpha = 0.75,
                                                export_network_file = "cytoscape_network.graphml",
                                                output_dir = ".",
                                                fixed_layout_file = file.path(output_dir, "fixed_layout.rds")) {
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  ######
  # Build node data frames from key identifier columns
  ######
  nodes_X <- data %>% 
    select(X) %>% 
    distinct() %>% 
    mutate(type = "metabolite")
  
  nodes_pathway <- data %>% 
    select(pathway) %>% 
    distinct() %>% 
    mutate(type = "pathway")
  
  nodes_gene <- data %>% 
    select(gene) %>% 
    distinct() %>% 
    mutate(type = "gene")
  
  # Combine nodes and assign unique IDs
  nodes <- bind_rows(
    nodes_X %>% rename(name = X),
    nodes_pathway %>% rename(name = pathway),
    nodes_gene %>% rename(name = gene)
  ) %>% mutate(id = row_number())
  
  ######
  # Build edge data frames
  ######
  edges_X_pathway <- data %>% 
    select(X, pathway, val) %>% 
    distinct() %>%
    rename(from = X, to = pathway)
  
  edges_pathway_gene <- data %>% 
    select(pathway, gene, val) %>% 
    distinct() %>%
    rename(from = pathway, to = gene)
  
  edges <- bind_rows(
    edges_X_pathway,
    edges_pathway_gene
  ) %>%
    left_join(nodes %>% select(name, id), by = c("from" = "name")) %>%
    rename(from_id = id) %>%
    left_join(nodes %>% select(name, id), by = c("to" = "name")) %>%
    rename(to_id = id)
  
  ######
  # Create the base graph and fixed layout
  ######
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  # Check if a saved fixed layout exists; if yes, load it, otherwise compute and save it.
  if (file.exists(fixed_layout_file)) {
    fixed_layout <- readRDS(fixed_layout_file)
  } else {
    fixed_layout <- create_layout(graph, layout = layout_type)
    saveRDS(fixed_layout, fixed_layout_file)
  }
  
  ######
  # Identify unique comparisons (assumes data has a 'comparison' column)
  ######
  unique_comparisons <- na.omit(unique(data$comparison))
  
  # Loop over each unique comparison to generate network images
  for (comp in unique_comparisons) {
    
    # Extract node names associated with current comparison (from X, pathway, gene)
    highlight_nodes <- data %>% 
      filter(comparison == comp) %>% 
      select(X, pathway, gene, log2FoldChange) %>%
      pivot_longer(cols = c(X, pathway, gene), values_to = "node_name") %>%
      distinct(node_name) %>% 
      pull(node_name)
    
    # Extract gene log2FoldChange values for current comparison
    gene_fc <- data %>% 
      filter(comparison == comp) %>% 
      select(gene, log2FoldChange) %>% 
      distinct()
    
    # Extract metabolite-specific met.Log2FC values, using 'X' and 'comparison' as key.
    metabolite_fc <- data %>% 
      filter(comparison == comp) %>% 
      select(X, met.Log2FC) %>% 
      distinct()
    
    # Merge the fixed_layout (with x, y) with the original nodes and the FC info.
    node_data <- fixed_layout %>% 
      left_join(nodes %>% select(name, type), by = "name") %>% 
      left_join(gene_fc, by = c("name" = "gene")) %>% 
      left_join(metabolite_fc, by = c("name" = "X")) %>% 
      # Resolve potential duplicate 'type' columns using coalesce
      mutate(type = coalesce(`type.x`, `type.y`)) %>%
      mutate(
        highlight = name %in% highlight_nodes,
        log2FoldChange_for_plot = ifelse(type == "gene" & highlight, log2FoldChange, NA_real_),
        met_Log2FC_for_plot       = ifelse(type == "metabolite" & highlight, met.Log2FC, NA_real_),
        # For pathway nodes, assign a discrete color (orange when highlighted)
        node_color = case_when(
          type == "pathway" & highlight ~ "orange",
          TRUE ~ "grey80"
        )
      ) %>% 
      select(-`type.x`, -`type.y`)
    
    # Split node_data into gene, metabolite and pathway subsets
    node_data_gene <- node_data %>% filter(type == "gene")
    node_data_met  <- node_data %>% filter(type == "metabolite")
    node_data_path <- node_data %>% filter(type == "pathway")
    
    ######
    # Build the ggplot using fixed layout coordinates
    ######
    # Set seed for reproducibility of repelled text positions
    set.seed(123)
    p <- ggraph(graph, layout = "manual", x = fixed_layout$x, y = fixed_layout$y) +
      # Draw edges
      geom_edge_link(aes(width = val),
                     alpha = edge_alpha,
                     color = "grey50",
                     show.legend = FALSE) +
      
      # --- Pathway nodes (discrete color) ---
      geom_node_point(data = node_data_path,
                      aes(x = x, y = y, color = node_color, shape = type, size = type),
                      alpha = node_alpha) +
      scale_color_manual(values = c("orange" = "orange", "grey80" = "grey80"),
                         guide = "none") +
      
      # Reset color scale for gene nodes
      ggnewscale::new_scale_color() +
      
      # --- Gene nodes (continuous gradient) ---
      geom_node_point(data = node_data_gene,
                      aes(x = x, y = y, color = log2FoldChange_for_plot, shape = type, size = type),
                      alpha = node_alpha) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red",
                            midpoint = 0, na.value = "grey80", 
                            name = "Gene log2FC",
                            guide = guide_colorbar(order = 3)) +
      
      # Reset color scale for metabolite nodes
      ggnewscale::new_scale_color() +
      
      # --- Metabolite nodes (continuous gradient) ---
      geom_node_point(data = node_data_met,
                      aes(x = x, y = y, color = met_Log2FC_for_plot, shape = type, size = type),
                      alpha = node_alpha) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red",
                            midpoint = 0, na.value = "grey80", 
                            name = "Met log2FC",
                            guide = guide_colorbar(order = 4)) +
      
      # Define shape and size scales with specified guide orders:
      scale_shape_manual(values = c("metabolite" = 15, "pathway" = 17, "gene" = 19),
                         guide = guide_legend(order = 2)) +
      scale_size_manual(values = c("metabolite" = 6 * node_size_factor,
                                   "pathway" = 6 * node_size_factor,
                                   "gene" = 6 * node_size_factor),
                        #guide = guide_legend(order = 1)) +
                        guide = "none") +
      
      scale_edge_width(range = c(0.2, 0.8), guide = "none") +
      theme_graph(base_family = "Sans") +
      
      # Add text labels with repelling (reproducibly)
      geom_node_text(data = node_data,
                     aes(x = x, y = y, label = name),
                     size = text_size,
                     repel = TRUE,
                     color = "black") +
      #labs(shape = "Node Type", size = "Node Size") +
      labs(shape = "Node Type") +
      ggtitle(comp)
    
    # Save output images: PNG and PDF (using cairo_pdf for PDF)
    png_file_path <- file.path(output_dir, paste0("network_", comp, "_MetsFC.png"))
    ggsave(filename = png_file_path, plot = p, width = 10, height = 8)
    
    pdf_file_path <- file.path(output_dir, paste0("network_", comp, "_MetsFC.pdf"))
    ggsave(filename = pdf_file_path, plot = p, width = 10, height = 8, device = cairo_pdf)
    
    print(p)
  }
  
  ######
  # Export the complete (uncolored) network as a GraphML file for Cytoscape
  ######
  igraph_obj <- as.igraph(graph)
  igraph::write_graph(igraph_obj, file = file.path(output_dir, export_network_file), format = "graphml")
  
  return(list(graph = graph, fixed_layout = fixed_layout))
}



create_network_visualization_JH_4(data, output_dir=output_dir,
                                  fixed_layout_file=file.path(output_dir, "fixed_layout.rds"))
















################################################################################
# below is Kai's original version
# Note: comparision should be changed to 'comparison'
################################################################################

# # Enhanced Network Visualization Function
# # This function creates a multi-level network visualization connecting metabolites, pathways, and genes
# # with options for filtering, layout, and appearance customization
# 
# create_network_visualization <- function(data, 
#                                          selected_comparisons = NULL,
#                                          layout_type = "fr", 
#                                          node_size_factor = 1,
#                                          text_size = 3,
#                                          edge_alpha = 0.5,
#                                          node_alpha = 0.75,
#                                          custom_colors = NULL) {
# 
#   # Filter by comparisons if specified
#   if (!is.null(selected_comparisons)) {
#     data <- data %>% filter(comparision %in% selected_comparisons)
#   }
#   
#   # Create node dataframes
#   nodes_X <- data %>% 
#     select(X) %>% 
#     distinct() %>% 
#     mutate(type = "metabolite")
#   
#   nodes_pathway <- data %>% 
#     select(pathway) %>% 
#     distinct() %>% 
#     mutate(type = "pathway")
#   
#   nodes_group <- data %>% 
#     select(Group, comparision) %>% 
#     distinct() %>% 
#     rename(name = Group) %>% 
#     mutate(type = "gene")
#   
#   # Handle the issue of nodes appearing in multiple comparisons
#   # Create a helper function to collapse comparisons for the same node
#   collapse_comparisons <- function(df) {
#     df %>%
#       group_by(name, type) %>%
#       summarize(comparision = paste(unique(comparision), collapse = ";"),
#                 .groups = "drop")
#   }
#   
#   nodes_group <- collapse_comparisons(nodes_group)
#   
#   # Combine all nodes
#   nodes <- bind_rows(
#     nodes_X %>% rename(name = X) %>% mutate(comparision = NA_character_),
#     nodes_pathway %>% rename(name = pathway) %>% mutate(comparision = NA_character_),
#     nodes_group
#   ) %>% 
#     mutate(id = row_number())
#   
#   # Create edges dataframe
#   # Create edges from X to pathway
#   edges_X_pathway <- data %>% 
#     select(X, pathway, val, comparision) %>% 
#     distinct() %>%
#     rename(from = X, to = pathway)
#   
#   # Create edges from pathway to Group
#   edges_pathway_group <- data %>% 
#     select(pathway, Group, val, comparision) %>% 
#     distinct() %>%
#     rename(from = pathway, to = Group)
#   
#   # Combine all edges
#   edges <- bind_rows(
#     edges_X_pathway,
#     edges_pathway_group
#   ) %>%
#     # Join with nodes to get the IDs
#     left_join(nodes %>% select(name, id), by = c("from" = "name")) %>%
#     rename(from_id = id) %>%
#     left_join(nodes %>% select(name, id), by = c("to" = "name")) %>%
#     rename(to_id = id) %>%
#     select(from_id, to_id, val, comparision)
#   
#   # Step 3: Create the graph
#   graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
#   
#   # Define colors if not provided
#   if (is.null(custom_colors)) {
#     # Set of visually distinct colors as default
#     default_colors <- c(
#       "#377EB8", "#4DAF4A", "#984EA3", 
#       "#FF7F00", "#FFFF33", "#A65628", "#F781BF", 
#       "#999999", "#66C2A5", "#FC8D62"
#     )
#     custom_colors <- default_colors
#   }
#   
#   # Create a function to extract the first comparison for each node
#   get_node_color <- function(comparison_str) {
#     if (is.na(comparison_str)) return(NA)
#     # Get the first comparison if multiple are present
#     first_comp <- strsplit(comparison_str, ";")[[1]][1]
#     return(first_comp)
#   }
#   
#   
#   # Create the visualization
#   # Determine if we're showing a single comparison or multiple
#   single_comparison <- length(unique(na.omit(edges$comparision))) == 1
#   
#   # Create the base plot
#   p <- ggraph(graph, layout = layout_type)
#   
#   # Add edges - different styling for single vs multiple comparisons
#   if (single_comparison) {
#     p <- p + geom_edge_link(aes(width = val), 
#                             alpha = edge_alpha,
#                             color = "grey50")
#   } else {
#     p <- p + geom_edge_link(aes(width = val), #, color = comparision), 
#                             alpha = edge_alpha,
#                             show.legend = TRUE)
#   }
#   
#   # Add nodes
#   p <- p + geom_node_point(aes(shape = type, 
#                                size = type),
#                            alpha = node_alpha) +
#     geom_node_text(aes(label = name), 
#                    repel = TRUE, 
#                    size = text_size) +
#     scale_shape_manual(values = c("metabolite" = 15, 
#                                   "pathway" = 17, 
#                                   "gene" = 19)) +
#     scale_size_manual(values = c("metabolite" = 6 * node_size_factor, 
#                                  "pathway" = 4 * node_size_factor, 
#                                  "gene" = 3 * node_size_factor)) +
#     theme_graph(base_family = "Arial")  + # Explicitly set a supported font
#     scale_edge_width(range = c(0.2, 0.8), guide = "none") +
#     theme_graph() +
#     labs(shape = "Node Type",
#          size = "Node Type",
#          edge_color = "Comparison")
#   
#   # Apply custom colors if there are enough for the number of comparisons
#   unique_comparisons <- na.omit(unique(edges$comparision))
#   if (length(custom_colors) >= length(unique_comparisons)) {
#     p <- p + 
#       scale_color_manual(values = custom_colors, na.value = "grey50")
#   }
#   
#   # Modify legend title based on single/multiple comparisons
#   if (single_comparison) {
#     p <- p + labs(color = "Node Type",
#                   shape = "Node Type",
#                   size = "Node Type")
#   }
#   
#   return(p)
# }
# #####
# create_network_visualization(data, selected_comparisons = "DR_vs_HFD")
