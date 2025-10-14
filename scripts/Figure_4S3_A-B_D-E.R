#----Figure_4S3_A-B_D-E.R----
#-------------------------------------------------------------------------------
# FAFB v783 brainmesh of 5HTNs, PPL1 DANs
# Princeton Synapses (July 2025). Data from Codex.
# Serotonin neurons taken from funkelab/drosophila_neurotransmitters on 09-22-25
# and filtered for confidence >= 5.
# load packages:----------------------------------------------------------------
library(tidyverse)
library(data.table)
library(fafbseg)
library(elmr)
library(zoo)
# set var:----------------------------------------------------------------------

read_data = FALSE
write_data = TRUE
write_plots = TRUE

options(scipen = 30)
options(dplyr.summarise.inform = FALSE)

# set colors:-------------------------------------------------------------------

# intialize:--------------------------------------------------------------------

input_path <- "./input/"
output_path <- "./output/"
dataset_path <- paste0(Sys.getenv("R_USER"), "/drosophila_connectome_shared_data/") #Documents/drosophila_connectome_shared_data of current user.
if(!dir.exists(file.path(dataset_path))){ #Uses input_path otherwise.
  dataset_path <- input_path
}


dir.create(paste0(input_path, "meshes/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_path, "brainmesh/"), showWarnings = FALSE, recursive = TRUE)
#-------------------------------------------------------------------------------
# Load Data
#-------------------------------------------------------------------------------
# Synapse coordinates - Princeton
if(!exists("coords")){
  coords <- fread(paste0(dataset_path, "fafb_v783_princeton_synapse_table.csv"),
                  colClasses = c("pre_root_id_720575940" = "character",
                                 "post_root_id_720575940" = "character"))
}

# Cell Types from FAFB v783 Codex -----
consolidated_cell_types <- read_csv(paste0(dataset_path, "consolidated_cell_types.csv"),
                                    col_types= cols(root_id = col_character())) %>%
  rename(cell_type = primary_type) %>%
  select(-`additional_type(s)`)

# gt data from funkelab/drosophila_neurotransmitters -----
gt_data <- read_csv(paste0(dataset_path, "/drosophila_neurotransmitters/gt_data.csv"))
# Take the highest confidence serotonin verification.
# If two equal confidence conflicts known and unknown (1/-1 and 0), take the known.
# If it conflicts positive and negative (1 and -1), send a warning.
gt_ser <- gt_data %>%
  filter(neurotransmitter_verified_confidence >= 5) %>% # Optional filter
  group_by(cell_type) %>%
  filter(neurotransmitter_verified_confidence == max(neurotransmitter_verified_confidence)) %>%
  # in case of ties, rank serotonin values: 1 > -1 > 0
  arrange(desc(abs(serotonin)), desc(serotonin)) %>%
  summarise(
    serotonin_val = first(serotonin),
    serotonin_all = list(serotonin),
    source = list(neurotransmitter_verified_source),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    conflict = (1 %in% serotonin_all) & (-1 %in% serotonin_all)
  ) %>%
  ungroup() %>%
  filter(serotonin_val == 1)

# warn on conflicts only (-1 and 1 together)
if (any(gt_ser$conflict)) {
  warning("Conflicting serotonin assignments for cell_types: ",
          paste(gt_ser$cell_type[gt_ser$conflict], collapse = ", "))
}

# final set of serotonin-positive types
ser_types <- gt_ser %>%
  filter(!conflict) %>%
  pull(cell_type)

# Standardize text
print(ser_types)
# Split by comma
ser_types <- strsplit(ser_types, ",") %>% unlist()

#-------------------------------------------------------------------------------
# Get IDs
#-------------------------------------------------------------------------------
# Serotonin from funkelab
ser <- consolidated_cell_types %>%
  filter(cell_type %in% ser_types)# %>% pull(root_id) %>% paste0(collapse = ",")
# PPL1's
dan <- consolidated_cell_types %>%
  filter(root_id %in% c("720575940621040737", "720575940617691170"))
#-------------------------------------------------------------------------------
# Plot SER
#-------------------------------------------------------------------------------
NOI_list <- ser$root_id
neuron_filename <- "5HT_CONF5"
input_coords <- coords[post_root_id_720575940 %in% substring(NOI_list, 10),
                       .(x = ctr_x, y = ctr_y, z = ctr_z)]
output_coords <- coords[pre_root_id_720575940 %in% substring(NOI_list, 10),
                        .(x = ctr_x, y = ctr_y, z = ctr_z)]
# Download meshes if necessary
existing_meshes <- file.exists(paste0(dataset_path, "meshes/", NOI_list, ".obj"))
if(!all(existing_meshes)){
  fafbseg::save_cloudvolume_meshes(NOI_list[!existing_meshes], savedir = paste0(input_path, "meshes"))
}

# Read meshes
NEURON_MESHES <- nat::read.neurons(paste0(dataset_path, "meshes/", NOI_list, ".obj"))

# Plot neurons ---
open3d()
par3d(skipRedraw = TRUE) # Do not render
view3d(userMatrix = rotationMatrix(90*pi/90, 1, 0, 0),zoom=0.5)
par3d(windowRect = c(0, 0, 2560, 1440))
#plot brain mesh
plot3d(FAFB14.surf,add=T, alpha=0.1,col="grey")
#plot neurons
plot3d(NEURON_MESHES, col="black", add=T, WithNodes=F, soma=T)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_neurons.png'))

# Plot input ---
plot3d(input_coords$x,input_coords$y,input_coords$z, add=T, col="magenta3", size=3, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_postsynapses.png'))

# Plot output ---
pop3d()
plot3d(output_coords$x,output_coords$y,output_coords$z, add=T, col="limegreen", size=4, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_presynapses.png'))

# Plot both ---
pop3d()
plot3d(input_coords$x,input_coords$y,input_coords$z, add=T, col="magenta3", size=3, alpha = 1)
plot3d(output_coords$x,output_coords$y,output_coords$z, add=T, col="limegreen", size=4, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_synapses.png'))
close3d()

# Plot input ---
open3d()
par3d(skipRedraw = TRUE) # Do not render
view3d(userMatrix = rotationMatrix(90*pi/90, 1, 0, 0),zoom=0.5)
par3d(windowRect = c(0, 0, 2560, 1440))
#plot brain mesh
plot3d(FAFB14.surf,add=T, alpha=0.1,col="grey")
#plot synapses
plot3d(input_coords$x,input_coords$y,input_coords$z, add=T, col="magenta3", size=3, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_postsynapses_only.png'))

# Plot output ---
pop3d()
plot3d(output_coords$x,output_coords$y,output_coords$z, add=T, col="limegreen", size=4, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_presynapses_only.png'))

# Plot both ---
pop3d()
plot3d(input_coords$x,input_coords$y,input_coords$z, add=T, col="magenta3", size=3, alpha = 0.5)
plot3d(output_coords$x,output_coords$y,output_coords$z, add=T, col="limegreen", size=4, alpha = 0.5)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_synapses_only.png'))
close3d()

#-------------------------------------------------------------------------------
# Plot DAN
#-------------------------------------------------------------------------------
NOI_list <- dan$root_id
neuron_filename <- "DAN"
input_coords <- coords[post_root_id_720575940 %in% substring(NOI_list, 10),
                       .(x = ctr_x, y = ctr_y, z = ctr_z)]
output_coords <- coords[pre_root_id_720575940 %in% substring(NOI_list, 10),
                        .(x = ctr_x, y = ctr_y, z = ctr_z)]
# Download meshes if necessary
existing_meshes <- file.exists(paste0(dataset_path, "meshes/", NOI_list, ".obj"))
if(!all(existing_meshes)){
  fafbseg::save_cloudvolume_meshes(NOI_list[!existing_meshes], savedir = paste0(input_path, "meshes"))
}

# Read meshes
NEURON_MESHES <- nat::read.neurons(paste0(dataset_path, "meshes/", NOI_list, ".obj"))

# Plot neurons ---
open3d()
par3d(skipRedraw = TRUE) # Do not render
view3d(userMatrix = rotationMatrix(90*pi/90, 1, 0, 0),zoom=0.5)
par3d(windowRect = c(0, 0, 2560, 1440))
#plot brain mesh
plot3d(FAFB14.surf,add=T, alpha=0.1,col="grey")
#plot neurons
plot3d(NEURON_MESHES, col="black", add=T, WithNodes=F, soma=T)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_neurons.png'))

# Plot input ---
plot3d(input_coords$x,input_coords$y,input_coords$z, add=T, col="magenta3", size=3, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_postsynapses.png'))

# Plot output ---
pop3d()
plot3d(output_coords$x,output_coords$y,output_coords$z, add=T, col="limegreen", size=4, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_presynapses.png'))

# Plot both ---
pop3d()
plot3d(input_coords$x,input_coords$y,input_coords$z, add=T, col="magenta3", size=3, alpha = 1)
plot3d(output_coords$x,output_coords$y,output_coords$z, add=T, col="limegreen", size=4, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_synapses.png'))
close3d()

# Plot input ---
open3d()
par3d(skipRedraw = TRUE) # Do not render
view3d(userMatrix = rotationMatrix(90*pi/90, 1, 0, 0),zoom=0.5)
par3d(windowRect = c(0, 0, 2560, 1440))
#plot brain mesh
plot3d(FAFB14.surf,add=T, alpha=0.1,col="grey")
#plot synapses
plot3d(input_coords$x,input_coords$y,input_coords$z, add=T, col="magenta3", size=3, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_postsynapses_only.png'))

# Plot output ---
pop3d()
plot3d(output_coords$x,output_coords$y,output_coords$z, add=T, col="limegreen", size=4, alpha = 1)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_presynapses_only.png'))

# Plot both ---
pop3d()
plot3d(input_coords$x,input_coords$y,input_coords$z, add=T, col="magenta3", size=3, alpha = 0.5)
plot3d(output_coords$x,output_coords$y,output_coords$z, add=T, col="limegreen", size=4, alpha = 0.5)
rgl.snapshot(paste0(output_path, "brainmesh/", neuron_filename, '_synapses_only.png'))
close3d()





