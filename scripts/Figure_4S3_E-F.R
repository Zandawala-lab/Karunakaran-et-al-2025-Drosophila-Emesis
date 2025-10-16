#-----Figure_4S3_E_F.R-----
#-------------------------------------------------------------------------------
# Influence Analysis using Codex for 5HTNs, PPL1s.
#
# From Codex v783 (09/20/25):
# connections.csv, classification.csv, consolidated_cell_types.csv
# From GitHub funkelab/drosophila_neurotransmitters (09/20/25):
# gt_data.csv
#
# Influence score takes ~5 hrs to run on initial run.
# We used a lenient const c of +24, matching BANC, for adjusted_influence that
# corresponds to the minimum accepted influence of 3.78e-11.
#-------------------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(influencer)
library(viridis)

# Configuration variables ------------------------------------------------------

write_data <- TRUE
write_plots <- TRUE
neuron_filename <- "NOI_influence" 

# intialize:--------------------------------------------------------------------

input_path <- "./input/"
output_path <- "./output/"
dataset_path <- paste0(Sys.getenv("R_USER"), "/drosophila_connectome_shared_data/") #Documents/drosophila_connectome_shared_data of current user.
if(!dir.exists(file.path(dataset_path))){ #Uses input_path otherwise.
  dataset_path <- input_path
}


dir.create(paste0(input_path), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_path), showWarnings = FALSE, recursive = TRUE)
#-------------------------------------------------------------------------------
# Load Data
#-------------------------------------------------------------------------------

# gt data from funkelab/drosophila_neurotransmitters -----
gt_data <- read_csv(paste0(dataset_path, "/drosophila_neurotransmitters/gt_data.csv"))


#-------------------------------------------------------------------------------
# Synapse coordinates - Princeton
connections <- read_csv(paste0(dataset_path,"connections.csv"), 
                        col_types = cols(pre_root_id = col_character(),
                                         post_root_id = col_character()))
connections <- connections %>% group_by(pre_root_id, post_root_id) %>% summarise(synapses = sum(syn_count))

#-------------------------------------------------------------------------------
#Load CSV
classification <- read_csv(paste0(dataset_path,"classification.csv"),
                           col_types = cols(root_id = col_character()))

# Use consolidated cell types---------------------------------------------------

consolidated_cell_types <- read_csv(paste0(dataset_path, "consolidated_cell_types.csv"),
                                    col_types= cols(root_id = col_character())) %>%
  rename(cell_type = primary_type) %>%
  select(-`additional_type(s)`)

#classification <- classification %>%
#  rows_update(consolidated_cell_types, by = "root_id")
classification <- classification %>%
  left_join(consolidated_cell_types, by = "root_id")

# Append super class to NA Cell types to help prevent duplicates----------------
classification$cell_type[is.na(classification$cell_type)] <- paste("NA", classification$super_class[is.na(classification$cell_type)],sep = "_")

# Change remaining NAs to chr
classification[is.na(classification)] <- "NA"

# Make all cell types only 1 super class----------------------------------------
# Cell types should only belong to one super_class.

# LC10b optic -> visual_projection
classification <- classification %>%
  mutate(super_class = ifelse(cell_type == "LC10b", "visual_projection", super_class))
# MeMe_e02 optic -> visual_projection
classification <- classification %>%
  mutate(super_class = ifelse(cell_type == "MeMe_e02", "visual_projection", super_class))
# R7 optic -> sensory
classification <- classification %>%
  mutate(super_class = ifelse(cell_type == "R7", "sensory", super_class))

# differentiated more pairs of cell types with >1 super class
classification <- classification %>% 
  group_by(cell_type) %>%
  mutate(cell_type = if(n_distinct(super_class) > 1)
    paste(cell_type, super_class, sep = "_")
    else(cell_type))

rm(consolidated_cell_types)

classification <- classification %>% ungroup()
#-------------------------------------------------------------------------------






# Take the highest confidence serotonin verification.
# If two equal confidence conflicts known and unknown (1/-1 and 0), take the known.
# If it conflicts positive and negative (1 and -1), send a warning.
ser_types <- gt_data %>%
  filter(neurotransmitter_verified_confidence >= 5) %>% # Optional filter
  group_by(cell_type) %>%
  filter(neurotransmitter_verified_confidence == max(neurotransmitter_verified_confidence)) %>%
  # in case of ties, rank serotonin values: 1 > -1 > 0
  arrange(desc(abs(serotonin)), desc(serotonin)) %>%
  summarise(
    serotonin_val = first(serotonin),
    serotonin_all = list(serotonin),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    conflict = (1 %in% serotonin_all) & (-1 %in% serotonin_all)
  ) %>%
  ungroup() %>%
  filter(serotonin_val == 1)

# warn on conflicts only (-1 and 1 together)
if (any(ser_types$conflict)) {
  warning("Conflicting serotonin assignments for cell_types: ",
          paste(ser_types$cell_type[ser_types$conflict], collapse = ", "))
}

# final set of serotonin-positive types
ser_types <- ser_types %>%
  filter(!conflict) %>%
  pull(cell_type)

# Standardize text
print(ser_types)
# Split by comma
ser_types <- strsplit(ser_types, ",") %>% unlist()

#-------------------------------------------------------------------------------
# Get IDs
#-------------------------------------------------------------------------------

classification <- classification %>% mutate(class = ifelse(cell_type %in% ser_types, "SER", class))

# Serotonin from funkelab
ser <- classification %>%
  filter(cell_type %in% ser_types)
ser_no_bSEL <- ser %>%
  filter(cell_type != "DNg28")
# PPL1's
dan <- classification %>%
  filter(root_id %in% c("720575940621040737", "720575940617691170"))
# KC alphabeta and gamma
kc <- classification %>%
  filter(cell_type %in% c("KCab", "KCab-p", "KCg-m", "KCg-d", "KCg-s1", "KCg-s2", "KCg-s3", "KCg-s"))
# KC alphabeta
kcab <- classification %>%
  filter(cell_type %in% c("KCab", "KCab-p"))
# KC gamma
kcg <- classification %>%
  filter(cell_type %in% c("KCg-m", "KCg-d", "KCg-s1", "KCg-s2", "KCg-s3", "KCg-s"))
#-------------------------------------------------------------------------------
# Data Preparation
#-------------------------------------------------------------------------------

# connections should have columns: pre_root_id, post_root_id, synapses
# classification should have columns: root_id, cell_type, super_class, class

# Create edge list with synapse counts and postsynaptic-normalisation
# This normalisation step helps calibrate for neuron size and connection number
# i.e. a neuron with more inputs likely needs more of them 'activated' to have an effect
edges.table <- connections %>%
  select(pre = pre_root_id, post = post_root_id, count = synapses) %>%
  group_by(post) %>%
  mutate(post_count = sum(count, na.rm = TRUE)) %>%
  mutate(norm = round(count/post_count, 4)) %>%
  ungroup() %>%
  mutate(pre = as.character(pre),
         post = as.character(post)) %>%
  filter(count > 5)  # Apply synapse count threshold

cat("Loaded", nrow(edges.table), "edges from local connectome data\n")
head(edges.table)

# Mutate column names if necessary
classification_mapped <- classification %>%
  #mutate(region = class) %>%
  select(root_id, cell_type, super_class, class)


# Update class column for plotting. Distinguish PPL1s and KC types.
classification_mapped <- classification_mapped %>%
  mutate(class = ifelse(class == "NA", "other", class),
         class = ifelse(cell_type %in% c("KCab", "KCab-p"), "KCab", class),
         class = ifelse(cell_type %in% c("KCg-m", "KCg-d", "KCg-s1", "KCg-s2", "KCg-s3", "KCg-s"), "KCg", class),
         class = ifelse(class == "Kenyon_Cell", "other KC", class),
         class = ifelse(cell_type %in% ser_types, "5HTN", class),
         class = ifelse(cell_type == "PPL101", "PPL1", class),
         class = ifelse(class == "DAN", "other DAN", class)
  )


#-------------------------------------------------------------------------------

# The influence calculator.
# Only run this once. Subsequent runs will clear the cache and result in long runtimes.

if (!exists("ic")){
  ic <- influence_calculator_r(edgelist_simple = edges.table, meta = classification_mapped)
}

#-------------------------------------------------------------------------------








analyze_neuron_influence <- function(noi.ids, edges.table, ic, classification_mapped) {
  
  # Get direct targets of the input neurons
  noi.post.ids <- edges.table %>%
    filter(pre %in% noi.ids) %>%
    pull(post) %>%
    unique() %>%
    as.character()
  
  cat("Direct targets of NOI:", length(noi.post.ids), "neurons\n")
  
  # Calculate influence (may be slow the first time)
  system.time({
    noi.influence <- ic$calculate_influence(noi.ids)
  })
  cat("Calculated influence for", nrow(noi.influence), "neurons\n")
  
  # Add metadata to the influence table
  noi.influence$directly_connected <- noi.influence$id %in% noi.post.ids
  noi.influence_meta <- left_join(noi.influence, classification_mapped, by = c("id" = "root_id"))
  
  # Summarize influence by class
  infl_class <- noi.influence_meta %>%
    filter(!is_seed) %>% # Remove self-connections
    group_by(class) %>%
    summarise(mean = mean(adjusted_influence, na.rm = TRUE), .groups = 'drop') %>%
    arrange(desc(mean)) %>%
    mutate(
      class = coalesce(class, "unknown"),
      class = gsub("_", " ", class)
    )
  
  # Summarize influence by superclass
  infl_sc <- noi.influence_meta %>%
    filter(!is_seed) %>% # Remove self-connections
    group_by(super_class) %>%
    summarise(mean = mean(adjusted_influence, na.rm = TRUE), .groups = 'drop') %>%
    arrange(desc(mean)) %>%
    mutate(
      super_class = coalesce(super_class, "unknown"),
      super_class = gsub("_", " ", super_class)
    )
  
  # Return all results in a named list
  return(list(
    meta = noi.influence_meta,
    class = infl_class,
    super_class = infl_sc
  ))
}

#-------------------------------------------------------------------------------
# 2. Set Up and Run the Loop
# This loop processes each neuron group using the function defined above.
#-------------------------------------------------------------------------------

# Create a list that defines the neuron groups to analyze
# Each item has the root IDs and a filename stem.
neuron_groups_to_process <- list(
  SER = list(ids = ser$root_id, filename = "SER_CONF5"),
  DAN = list(ids = dan$root_id, filename = "DAN_CONF5"),
  SER_no_bSEL  = list(ids = ser_no_bSEL$root_id,  filename = "SER_no_bSEL_CONF5")
)

# Loop through each neuron group
for (group_name in names(neuron_groups_to_process)) {
  
  cat("\n=======================================================\n")
  cat("Processing Neuron Group:", group_name, "\n")
  cat("=======================================================\n")
  
  current_group <- neuron_groups_to_process[[group_name]]
  
  # Run the main analysis function
  results <- analyze_neuron_influence(
    noi.ids = current_group$ids,
    edges.table = edges.table,
    ic = ic,
    classification_mapped = classification_mapped
  )
  
  # Write the detailed metadata file if required
  if (exists("write_data") && write_data) {
    output_filename <- paste0(output_path, current_group$filename, "_influence_meta.csv")
    write_csv(results$meta, file = output_filename)
    cat("Metadata saved to:", output_filename, "\n")
  }
  
  # Create the final variables (e.g., infl_class_SER, infl_sc_SER) in your global environment
  assign(paste0("infl_class_", group_name), results$class, envir = .GlobalEnv)
  assign(paste0("infl_sc_", group_name), results$super_class, envir = .GlobalEnv)
}




#-------------------------------------------------------------------------------
# Plotting Figure E
#-------------------------------------------------------------------------------

infl_class_SER$neuron_type <- "5HTN"
infl_class_DAN$neuron_type <- "PPL1"

infl_class_combined <- rbind(infl_class_SER, infl_class_DAN)

# Order
infl_class_combined$neuron_type <- factor(
  infl_class_combined$neuron_type,
  levels = rev(c("5HTN", "PPL1"))
)

# Fill in blanks
infl_class_combined <- infl_class_combined %>%
  complete(class, neuron_type)


# Order by highest influence, independent of from which source.
class_order <- infl_class_combined %>%
  group_by(class) %>%
  summarise(max_mean = max(mean, na.rm = TRUE)) %>%
  arrange(max_mean) %>%
  pull(class)


# Heatmap by class for all three neuron types
class_plot <- ggplot(infl_class_combined, aes(x = factor(class, levels = class_order), y = neuron_type, fill = mean)) +
  geom_raster() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis(
    option = "cividis",
    direction = 1,
    breaks = c(min(infl_class_combined$mean, na.rm = TRUE), 
               (min(infl_class_combined$mean, na.rm = TRUE) + max(infl_class_combined$mean, na.rm = TRUE)) / 2, 
               max(infl_class_combined$mean, na.rm = TRUE)),
    labels = round(c(min(infl_class_combined$mean, na.rm = TRUE), 
                     (min(infl_class_combined$mean, na.rm = TRUE) + max(infl_class_combined$mean, na.rm = TRUE)) / 2, 
                     max(infl_class_combined$mean, na.rm = TRUE)), 1),
    na.value = "white"
  ) + 
  labs(title = NULL,
       x = NULL,
       y = NULL,
       fill = "Mean influence") +
  theme_minimal() +
  coord_fixed() +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5 
    )
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        legend.title.position = "bottom",
        legend.key.width = unit(1, "null"),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal")

print(class_plot)

if(write_plots){
  ggsave(paste0(output_path, "Figure_4S3_F.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
  ggsave(paste0(output_path, "Figure_4S3_F.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
}








#-------------------------------------------------------------------------------
# Plotting Figure F
#-------------------------------------------------------------------------------

# Remove autocrine
infl_class_SER_no_bSEL <- infl_class_SER_no_bSEL %>% filter(class != "5HTN")

# Combine
infl_class_SER$neuron_type <- "5HTN"
infl_class_SER_no_bSEL$neuron_type <- "No Bitter SEL"

infl_class_combined <- rbind(infl_class_SER, infl_class_SER_no_bSEL)


# Order
infl_class_combined$neuron_type <- factor(
  infl_class_combined$neuron_type,
  levels = rev(c("5HTN", "No Bitter SEL"))
)


# Calculate overall mean by class across all neuron types for ordering
class_order <- infl_class_combined %>%
  group_by(class) %>%
  summarise(overall_mean = mean(mean, na.rm = TRUE)) %>%
  arrange(overall_mean) %>%
  pull(class)

# Heatmap by class for all three neuron types
class_plot <- ggplot(infl_class_combined, aes(x = factor(class, levels = class_order), y = neuron_type, fill = mean)) +
  geom_raster() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis(
    option = "cividis",
    direction = 1,
    breaks = c(min(infl_class_combined$mean, na.rm = TRUE), 
               (min(infl_class_combined$mean, na.rm = TRUE) + max(infl_class_combined$mean, na.rm = TRUE)) / 2, 
               max(infl_class_combined$mean, na.rm = TRUE)),
    labels = round(c(min(infl_class_combined$mean, na.rm = TRUE), 
                     (min(infl_class_combined$mean, na.rm = TRUE) + max(infl_class_combined$mean, na.rm = TRUE)) / 2, 
                     max(infl_class_combined$mean, na.rm = TRUE)), 1),
    na.value = "white"
  ) + 
  labs(title = NULL,
       x = NULL,
       y = NULL,
       fill = "Mean influence") +
  theme_minimal() +
  coord_fixed() +
  guides(
    fill = guide_colorbar(
      title.position = "top",  
      title.hjust = 0.5 
    )
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        legend.title.position = "bottom",
        legend.key.width = unit(1, "null"),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal")

print(class_plot)

if(write_plots){
  ggsave(paste0(output_path, "Figure_4S3_G.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
  ggsave(paste0(output_path, "Figure_4S3_G.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
}

