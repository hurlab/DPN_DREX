library(richR)
library(tidyverse)
library(richR)
data<-read.csv("input/gene_metabolite_pathway_group.csv",row.names=1)
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
#####
create_network_visualization(data, selected_comparisons = "DR_vs_HFD")