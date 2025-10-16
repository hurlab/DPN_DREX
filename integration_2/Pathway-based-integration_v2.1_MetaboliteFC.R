################################################################################
##################### Set Current Working Directory ###########################
################################################################################

if (!requireNamespace("rstudioapi", quietly = TRUE)) install.packages("rstudioapi")
library(rstudioapi)

current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(current_path)
base_dir <- dirname(current_path)

output_dir <- file.path(base_dir, "output2.1")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
data_dir <- file.path(base_dir, "data")

################################################################################
##################### Packages #################################################
################################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c(
  "tidygraph", "tidyverse", "richR", "ggplot2", "dplyr", "ggraph", "igraph",
  "patchwork"
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
data <- data_joined %>% select(-join_met, -join_comp, -X, -comparison.y)

# Change the 'val.y' column name to 'met.name'
data <- data %>% rename(met.Log2FC = val.y, comparison = comparison.x,
                        X=met.x)



################################################################################
# This is a modified verion of the function that creates a network visualization
# The network structure is maintained across different sets. 
################################################################################
# Date: 10/16/2025
# Version: 6
# Updates: removed node label, edge color change for active pairs
#          bigger node size
################################################################################

create_network_visualization_JH <- function(
    data,
    layout_type = "fr",
    node_size_factor = 1.2,
    edge_alpha = 0.5,
    node_alpha = 0.85,
    export_network_file = "cytoscape_network.graphml",
    output_dir = ".",
    fixed_layout_file = file.path(output_dir, "fixed_layout.rds"),
    plot_order = NULL
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ###### Build node and edge data ######
  nodes_X <- data %>% select(X) %>% distinct() %>% mutate(type = "metabolite")
  nodes_pathway <- data %>% select(pathway) %>% distinct() %>% mutate(type = "pathway")
  nodes_gene <- data %>% select(gene) %>% distinct() %>% mutate(type = "gene")
  
  nodes <- bind_rows(
    nodes_X %>% rename(name = X),
    nodes_pathway %>% rename(name = pathway),
    nodes_gene %>% rename(name = gene)
  ) %>% mutate(id = row_number())
  
  edges_X_pathway <- data %>%
    select(X, pathway, val) %>% distinct() %>%
    rename(from = X, to = pathway)
  
  edges_pathway_gene <- data %>%
    select(pathway, gene, val) %>% distinct() %>%
    rename(from = pathway, to = gene)
  
  edges <- bind_rows(edges_X_pathway, edges_pathway_gene) %>%
    left_join(nodes %>% select(name, id), by = c("from" = "name")) %>%
    rename(from_id = id) %>%
    left_join(nodes %>% select(name, id), by = c("to" = "name")) %>%
    rename(to_id = id)
  
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  if (file.exists(fixed_layout_file)) {
    fixed_layout <- readRDS(fixed_layout_file)
  } else {
    fixed_layout <- create_layout(graph, layout = layout_type)
    saveRDS(fixed_layout, fixed_layout_file)
  }
  
  if (!"name" %in% names(fixed_layout)) fixed_layout$name <- nodes$name
  
  unique_comparisons <- na.omit(unique(data$comparison))
  plot_list <- list()
  
  for (comp in unique_comparisons) {
    highlight_nodes <- data %>%
      filter(comparison == comp) %>%
      select(X, pathway, gene, log2FoldChange) %>%
      pivot_longer(cols = c(X, pathway, gene), values_to = "node_name") %>%
      distinct(node_name) %>%
      pull(node_name)
    
    gene_fc <- data %>%
      filter(comparison == comp) %>%
      select(gene, log2FoldChange) %>%
      distinct()
    
    metabolite_fc <- data %>%
      filter(comparison == comp) %>%
      select(X, met.Log2FC) %>%
      distinct()
    
    node_data <- fixed_layout %>%
      left_join(nodes %>% select(name, type), by = "name")
    
    if ("type.x" %in% names(node_data) || "type.y" %in% names(node_data)) {
      node_data <- node_data %>%
        mutate(type = dplyr::coalesce(.data[["type.x"]], .data[["type.y"]])) %>%
        select(-tidyselect::any_of(c("type.x", "type.y")))
    }
    
    node_data <- node_data %>%
      left_join(gene_fc, by = c("name" = "gene")) %>%
      left_join(metabolite_fc, by = c("name" = "X")) %>%
      mutate(
        highlight               = name %in% highlight_nodes,
        log2FoldChange_for_plot = ifelse(type == "gene" & highlight, log2FoldChange, NA_real_),
        met_Log2FC_for_plot     = ifelse(type == "metabolite" & highlight, met.Log2FC, NA_real_),
        node_color = dplyr::case_when(
          type == "pathway" & highlight ~ "orange",
          type == "pathway"             ~ "grey80",
          TRUE                          ~ NA_character_
        )
      )
    
    node_data_gene <- node_data %>% filter(type == "gene")
    node_data_met  <- node_data %>% filter(type == "metabolite")
    node_data_path <- node_data %>% filter(type == "pathway")
    
    # Prepare edges with active flag and linewidth
    edge_plot_df <- edges %>%
      transmute(
        from = from_id,
        to   = to_id,
        val  = val,
        x    = fixed_layout$x[from],
        y    = fixed_layout$y[from],
        xend = fixed_layout$x[to],
        yend = fixed_layout$y[to],
        n_from = nodes$name[from],
        n_to   = nodes$name[to]
      ) %>%
      mutate(
        active = n_from %in% highlight_nodes & n_to %in% highlight_nodes,
        linewidth = ifelse(active, 1.1, 0.6)
      )
    
    set.seed(123)
    p <- ggraph(graph, layout = "manual", x = fixed_layout$x, y = fixed_layout$y) +
      
      # --- Edges (darker when active, linewidth fixed, legend suppressed) ---
      geom_segment(
        data = edge_plot_df,
        aes(x = x, y = y, xend = xend, yend = yend, color = active),
        linewidth = edge_plot_df$linewidth,
        alpha = edge_alpha,
        lineend = "round",
        show.legend = FALSE  # removes only edge legend
      ) +
      scale_color_manual(values = c(`TRUE` = "grey10", `FALSE` = "grey70"), guide = "none") +
      
      ggnewscale::new_scale_colour() +
      
      # --- Pathway nodes (discrete color) ---
      geom_node_point(
        data = node_data_path,
        aes(x = x, y = y, color = node_color, shape = type),
        size = 7 * node_size_factor,
        alpha = node_alpha
      ) +
      scale_color_manual(values = c("orange" = "orange", "grey80" = "grey80"), guide = "none") +
      
      ggnewscale::new_scale_colour() +
      
      # --- Gene nodes (continuous colorbar) ---
      geom_node_point(
        data = node_data_gene,
        aes(x = x, y = y, color = log2FoldChange_for_plot, shape = type),
        size = 7 * node_size_factor,
        alpha = node_alpha
      ) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red",
                            midpoint = 0, na.value = "grey80",
                            name = "Gene log2FC",
                            guide = guide_colorbar(order = 3)) +
      
      ggnewscale::new_scale_colour() +
      
      # --- Metabolite nodes (continuous colorbar) ---
      geom_node_point(
        data = node_data_met,
        aes(x = x, y = y, color = met_Log2FC_for_plot, shape = type),
        size = 7 * node_size_factor,
        alpha = node_alpha
      ) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red",
                            midpoint = 0, na.value = "grey80",
                            name = "Met log2FC",
                            guide = guide_colorbar(order = 4)) +
      
      # --- Node type legend kept intact ---
      scale_shape_manual(
        values = c("metabolite" = 15, "pathway" = 17, "gene" = 19),
        guide = guide_legend(order = 2, title = "Node Type")
      ) +
      
      theme_graph(base_family = "Sans") +
      labs(shape = "Node Type") +
      ggtitle(comp)
    
    ggsave(file.path(output_dir, paste0("network_", comp, "_MetsFC.png")), plot = p, width = 12, height = 8)
    ggsave(file.path(output_dir, paste0("network_", comp, "_MetsFC.pdf")), plot = p, width = 12, height = 8, device = cairo_pdf)
    
    print(p)
    plot_list[[comp]] <- p
  }
  
  igraph::write_graph(as.igraph(graph), file = file.path(output_dir, export_network_file), format = "graphml")
  
  if (!is.null(plot_order) && all(plot_order %in% names(plot_list))) {
    plot_list <- plot_list[plot_order]
  }
  
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2, nrow = 2)
  ggsave(file.path(output_dir, "combined_networks.png"), plot = combined_plot, width = 24, height = 20)
  ggsave(file.path(output_dir, "combined_networks.pdf"), plot = combined_plot, width = 24, height = 20, device = cairo_pdf)
  
  return(list(graph = graph, fixed_layout = fixed_layout, plots = plot_list))
}



create_network_visualization_JH(data, output_dir=output_dir,
                                  fixed_layout_file=file.path(output_dir, "fixed_layout.rds"),
                                  plot_order=c("HFD_vs_SD", "DR_vs_HFD", "DREX_vs_HFD"))


# End of script

