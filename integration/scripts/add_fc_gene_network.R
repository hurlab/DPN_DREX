#' ============================================================================
#' Metabolite Fold-Change Network Integration
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
#'   Integrates fold-change values from 4 metabolite sources (acylcarnitines, polar,
#   SCN, untargeted) into the gene-metabolite-pathway network. Standardizes
#   metabolite names and harmonizes comparison labels across datasets.
#'
#' Input:
#'   meta_fc.rdata, gene_metabolite_pathway_group.csv
#'
#' Output:
#'   gene_metabolite_pathway_group_withFC.csv (network table with metabolite
#   fold-changes)
#'
#' Related figures:
#'   Fig. 5A-D (metabolite FC overlay)
#'
#' ============================================================================


library(tidyverse)
load("meta_fc.rdata")
net<-read.csv("gene_metabolite_pathway_group.csv",row.names = 1)

extract_comparison_data <- function(df) {
  # List of comparison columns we want to extract
  comparison_cols <- c(
    "HFD_26WKvsWT_26WK", 
    "DR_26wkvsWT_26WK", 
    "EX_26WKvsWT_26WK", 
    "DREX_26WKvsWT_26WK",
    "DR_26wkvsHFD_26WK", 
    "EX_26WKvsHFD_26WK", 
    "DREX_26WKvsHFD_26WK"
  )
  
  # Find which of these columns exist in the dataframe
  existing_cols <- comparison_cols[comparison_cols %in% colnames(df)]
  
  # Always keep the 'met' column for merging
  cols_to_keep <- c("met", existing_cols)
  
  # Extract only these columns
  df_comp <- df[, cols_to_keep]
  
  # Trim whitespace from metabolite names
  df_comp$met <- trimws(df_comp$met)
  
  # Standardize case: only first letter capitalized
  df_comp$met <- sapply(df_comp$met, standardize_case)
  
  return(df_comp)
}

# Apply the extraction function to each table
acyf_comp <- extract_comparison_data(acyf)
polarf_comp <- extract_comparison_data(polarf)
scnf_comp <- extract_comparison_data(scnf)
untf_comp <- extract_comparison_data(untf)

# Step 2: Add a source column to identify which table each row came from
acyf_comp$source <- "acyf"
polarf_comp$source <- "polarf"
scnf_comp$source <- "scnf"
untf_comp$source <- "untf"

# Step 3: Bind the four tables together
all_comp <- rbind(acyf_comp, polarf_comp, scnf_comp, untf_comp)

# Step 4: Rename the WT to SD in column names
colnames(all_comp) <- gsub("WT_26WK", "SD_26WK", colnames(all_comp))

######
comparison_map <- list(
  "HFD_26WKvsSD_26WK" = "HFD_vs_SD",
  "DR_26wkvsHFD_26WK" = "DR_vs_HFD",
  "EX_26WKvsHFD_26WK" = "EX_vs_HFD",
  "DREX_26WKvsHFD_26WK" = "DREX_vs_HFD"
)

# Columns to remove
columns_to_remove <- c(
  "DR_26wkvsSD_26WK",
  "EX_26WKvsSD_26WK",
  "DREX_26WKvsSD_26WK"
)

# Remove unwanted columns
all_comp <- all_comp[, !(colnames(all_comp) %in% columns_to_remove)]

# Rename the remaining columns according to the mapping
for (old_name in names(comparison_map)) {
  if (old_name %in% colnames(all_comp)) {
    new_name <- comparison_map[[old_name]]
    colnames(all_comp)[colnames(all_comp) == old_name] <- new_name
  }
}

# Display the updated column names
cat("Updated column names after removing unwanted columns:\n")
print(colnames(all_comp))

all_comp_long<-all_comp%>%gather(group,val,-met,-source)
net_renamed <- net
names(net_renamed)[names(net_renamed) == "X"] <- "met"
names(net_renamed)[names(net_renamed) == "comparision"] <- "comparison"
standardize_case <- function(x) {
# First convert everything to lowercase
x <- tolower(x)
# Then capitalize first letter
substr(x, 1, 1) <- toupper(substr(x, 1, 1))
return(x)
}
net_renamed$met<-sub('^ ','',net_renamed$met)
net_renamed$met <- sapply(net_renamed$met, standardize_case)
split_comp_long <- data.frame()

for(i in 1:nrow(all_comp_long)) {
  met_name <- all_comp_long$met[i]
  
  # Check if the metabolite name contains " & " indicating multiple compounds
  if(grepl(" & ", met_name, fixed = TRUE)) {
    # Split the metabolite name by " & "
    met_parts <- unlist(strsplit(met_name, " & ", fixed = TRUE))
    
    # Create a row for each part with the same values
    for(part in met_parts) {
      new_row <- all_comp_long[i,]
      new_row$met <- trimws(part)  # Trim whitespace from the split name
      split_comp_long <- rbind(split_comp_long, new_row)
    }
  } else {
    # If no " & " is found, just add the original row
    split_comp_long <- rbind(split_comp_long, all_comp_long[i,])
  }
}

# Replace the original data with the split version
all_comp_long <- split_comp_long

# Step 2: Standardize case in all_comp_long
all_comp_long$met <- trimws(all_comp_long$met)
all_comp_long$met <- sapply(all_comp_long$met, standardize_case)

net_renamed$x<-paste0(net_renamed$comparison,net_renamed$met)
all_comp_long$x<-paste0(all_comp_long$group,all_comp_long$met)

result<-left_join(net_renamed,all_comp_long,by=c("x"="x"),relationship = "many-to-many")

write.csv(result,file="gene_metabolite_pathway_group_withFC.csv")

