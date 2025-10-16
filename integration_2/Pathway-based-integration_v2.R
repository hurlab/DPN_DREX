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
    p <- p + geom_edge_link(aes(width = val), #, color = comparision), 
                            alpha = edge_alpha,
                            show.legend = TRUE)
  }
  
  # Add nodes
  p <- p + geom_node_point(aes(shape = type, 
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
    theme_graph(base_family = "Arial")  + # Explicitly set a supported font
    scale_edge_width(range = c(0.2, 0.8), guide = "none") +
    theme_graph() +
    labs(shape = "Node Type",
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




################################################################################
# This is a modified verion of the function that creates a network visualization
# The network structure is maintained across different sets. 
################################################################################
create_network_visualization_JH_3 <- function(data, 
                                              layout_type = "fr", 
                                              node_size_factor = 1,
                                              text_size = 3,
                                              edge_alpha = 0.5,
                                              node_alpha = 0.75,
                                              export_network_file = "cytoscape_network.graphml",
                                              output_dir = ".") {
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Build node data frames from key identifier columns
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
  
  # Build edge data frames
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
  
  # Create the base graph and compute a fixed layout
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  fixed_layout <- create_layout(graph, layout = layout_type)
  
  # (We're reverting to specifying x and y directly from fixed_layout)
  
  # Identify unique comparisons (assumes data has a 'comparison' column)
  unique_comparisons <- na.omit(unique(data$comparison))
  
  for (comp in unique_comparisons) {
    
    # Extract node names (from columns X, pathway, gene) associated with the current comparison,
    # and get the corresponding log2FoldChange values.
    highlight_nodes <- data %>% 
      filter(comparison == comp) %>% 
      select(X, pathway, gene, log2FoldChange) %>%
      pivot_longer(cols = c(X, pathway, gene), values_to = "node_name") %>%
      distinct(node_name) %>% 
      pull(node_name)
    
    # Extract gene log2FoldChange values for the current comparison
    gene_fc <- data %>% 
      filter(comparison == comp) %>% 
      select(gene, log2FoldChange) %>% 
      distinct()
    
    # Merge the fixed_layout with our original nodes (which contains 'type')
    # and then merge the gene fold-change info.
    node_data <- fixed_layout %>% 
      left_join(nodes %>% select(name, type), by = "name") %>% 
      left_join(gene_fc, by = c("name" = "gene")) %>% 
      mutate(
        type = coalesce(`type.x`, `type.y`),
        highlight = name %in% highlight_nodes,
        log2FoldChange_for_plot = ifelse(type == "gene" & highlight, log2FoldChange, NA_real_),
        node_color = case_when(
          type == "metabolite" & highlight ~ "green",
          type == "pathway" & highlight ~ "orange",
          TRUE ~ "grey80"
        )
      ) %>% select(-type.x, -type.y)
    
    
    # Split node_data into non-gene and gene subsets
    node_data_non_gene <- node_data %>% filter(type != "gene")
    node_data_gene <- node_data %>% filter(type == "gene")
    
    # Build the ggplot using fixed layout coordinates provided directly from fixed_layout
    set.seed(123)
    p <- ggraph(graph, layout = "manual", x = fixed_layout$x, y = fixed_layout$y) +
      # Plot edges with width scaled by 'val'
      geom_edge_link(aes(width = val),
                     alpha = edge_alpha,
                     color = "grey50",
                     show.legend = FALSE) +
      
      # Plot non-gene nodes using discrete colors stored in node_color
      geom_node_point(data = node_data_non_gene,
                      aes(x = x, y = y, color = node_color, shape = type, size = type),
                      alpha = node_alpha) +
      # Define a discrete color scale for non-gene nodes and remove its legend
      scale_color_manual(values = c("green" = "green", "orange" = "orange", "grey80" = "grey80"),
                         guide = "none") +
      
      # Start a new color scale for gene nodes
      ggnewscale::new_scale_color() +
      
      # Plot gene nodes with a continuous gradient mapped from log2FoldChange_for_plot
      geom_node_point(data = node_data_gene,
                      aes(x = x, y = y, color = log2FoldChange_for_plot, shape = type, size = type),
                      alpha = node_alpha) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red",
                            midpoint = 0, na.value = "grey80", name = "log2FC") +
      
      # Define shape and size scales for all nodes
      scale_shape_manual(values = c("metabolite" = 15, "pathway" = 17, "gene" = 19)) +
      scale_size_manual(values = c("metabolite" = 6 * node_size_factor,
                                   "pathway" = 6 * node_size_factor,
                                   "gene" = 6 * node_size_factor)) +
      scale_edge_width(range = c(0.2, 0.8), guide = "none") +
      theme_graph(base_family = "Sans") +
      
      # Add text labels for all nodes
      geom_node_text(data = node_data,
                     aes(x = x, y = y, label = name),
                     size = text_size,
                     repel = TRUE,
                     color = "black") +
      labs(shape = "Node Type", size = "Node Type") +
      
      # Add the comparison as the title of the plot
      ggtitle(comp)

    # Save the plot to a file named with the current comparison value
    file_name <- paste0("network_", comp, ".png")
    file_path <- file.path(output_dir, file_name)
    ggsave(filename = file_path, plot = p, width = 10, height = 8)
    file_name <- paste0("network_", comp, ".pdf")
    file_path <- file.path(output_dir, file_name)
    ggsave(filename = file_path, plot = p, width = 10, height = 8, device = cairo_pdf)
    print(p)
  }
  
  # Export the complete (uncolored) network as GraphML for Cytoscape
  igraph_obj <- as.igraph(graph)
  igraph::write_graph(igraph_obj, file = paste0(output_dir,export_network_file), format = "graphml")
  
  return(list(graph = graph, fixed_layout = fixed_layout))
}

create_network_visualization_JH_3(data, output_dir="output/")


