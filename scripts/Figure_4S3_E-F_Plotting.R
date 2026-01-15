#-----Figure_4S3_E_F_Plotting.R-----
#-------------------------------------------------------------------------------
# Continued from Figure_4S3_E_F.R
#
# Loads influence .csv's generated and plots as heatmap.
#-------------------------------------------------------------------------------

library(tidyverse)
library(influencer)
library(viridis)

scaled_heatmap_palette <- colorRampPalette(c('#2166AC', '#4393C3', '#92C5DE', '#F4A582', '#D6604D', '#B2182B'))(100)

dan <- read_csv(paste0("output/", "DAN_CONF5_influence_meta.csv"),
                col_types= cols(id = col_character()))
dan_full <- dan %>%
  rename(target = class) %>%
  adjust_influence()

ser <- read_csv(paste0("output/", "SER_CONF5_influence_meta.csv"),
                col_types= cols(id = col_character()))
ser_full <- ser %>%
  rename(target = class) %>%
  adjust_influence()

ser_alt <- read_csv(paste0("output/", "SER_no_bSEL_CONF5_influence_meta.csv"),
                    col_types= cols(id = col_character()))
ser_alt_full <- ser_alt %>%
  rename(target = class) %>%
  adjust_influence()


#-------------------------------------------------------------------------------
# Plotting Figure E
#-------------------------------------------------------------------------------
process_list <- list(PPL1 = dan_full,
                     `5HTN` = ser_full,
                     `No Bitter SEL` = ser_alt_full)

process_influence <- function(df, neuron_type_label) {
  df %>%
    mutate(neuron_type = neuron_type_label) %>%
    
    # set names to original
    rename(
      class = target,
      mean  = adjusted_influence_norm_by_sources_and_targets
    ) %>%
    
    mutate(mean = ifelse(is_seed, NA, mean)) %>%
    
    # Normalize scale 0–1 within neuron_type
    group_by(neuron_type) %>%
    mutate(mean = (mean - min(mean, na.rm = TRUE)) /
             (max(mean, na.rm = TRUE) - min(mean, na.rm = TRUE))) %>%
    ungroup() %>%
    
    # Fill in blanks
    complete(class, neuron_type) %>%
    
    # remove seed (always max influence to self)
    filter(!is_seed)
}

infl_class <- process_list %>% imap(~ process_influence(.x, .y))


# Heatmap by class for PPL1 - labels on bottom
class_plot <- ggplot(infl_class$`PPL1`, aes(x = fct_reorder(class, mean), y = neuron_type, fill = mean)) +
  geom_raster() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(
    colours = scaled_heatmap_palette,
    breaks = c(
      min(infl_class$`PPL1`$mean, na.rm = TRUE),
      max(infl_class$`PPL1`$mean, na.rm = TRUE)
    ),
    labels = round(c(
      min(infl_class$`PPL1`$mean, na.rm = TRUE),
      max(infl_class$`PPL1`$mean, na.rm = TRUE)
    ), 1),
    na.value = "white"
  ) +
  labs(title = NULL,
       x = NULL,
       y = NULL,
       fill = "mean influence (normalized)") +
  theme_minimal() +
  coord_fixed() +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5 
    )
  )+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.33),
        axis.text.y = element_blank(), # Make text in illustrator so ggsave is aligned
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        legend.title = element_text(angle = 90),
        legend.title.position = "top",
        legend.key.height = unit(1.02, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.position = "right",
        legend.direction = "vertical")

print(class_plot)

if(write_plots){
  ggsave(paste0(output_path, "Figure_4S3_PPL1.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
  ggsave(paste0(output_path, "Figure_4S3_PPL1.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
}


# Heatmap by class for 5HTN - labels on top
class_plot <- ggplot(infl_class$`5HTN`, aes(x = fct_reorder(class, mean), y = neuron_type, fill = mean)) +
  geom_raster() +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(
    colours = scaled_heatmap_palette,
    breaks = c(
      min(infl_class$`5HTN`$mean, na.rm = TRUE),
      max(infl_class$`5HTN`$mean, na.rm = TRUE)
    ),
    labels = round(c(
      min(infl_class$`5HTN`$mean, na.rm = TRUE),
      max(infl_class$`5HTN`$mean, na.rm = TRUE)
    ), 1),
    na.value = "white"
  ) +
  labs(title = NULL,
       x = NULL,
       y = NULL,
       fill = "mean influence (normalized)") +
  theme_minimal() +
  coord_fixed() +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5 
    )
  )+
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.33),
        axis.text.y = element_blank(), # Make text in illustrator so ggsave is aligned
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        legend.title = element_text(angle = 90),
        legend.title.position = "top",
        legend.key.height = unit(1.02, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.position = "right",
        legend.direction = "vertical")

print(class_plot)

if(write_plots){
  ggsave(paste0(output_path, "Figure_4S3_5HTN.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
  ggsave(paste0(output_path, "Figure_4S3_5HTN.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
}

#-------------------------------------------------------------------------------
# Plotting Figure F
#-------------------------------------------------------------------------------



# Heatmap by class for 5HTN without Bitter SEL - labels on top
class_plot <- ggplot(infl_class$`No Bitter SEL`, aes(x = fct_reorder(class, mean), y = neuron_type, fill = mean)) +
  geom_raster() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradientn(
    colours = scaled_heatmap_palette,
    breaks = c(
      min(infl_class$`No Bitter SEL`$mean, na.rm = TRUE),
      max(infl_class$`No Bitter SEL`$mean, na.rm = TRUE)
    ),
    labels = round(c(
      min(infl_class$`No Bitter SEL`$mean, na.rm = TRUE),
      max(infl_class$`No Bitter SEL`$mean, na.rm = TRUE)
    ), 1),
    na.value = "white"
  ) +
  labs(title = NULL,
       x = NULL,
       y = NULL,
       fill = "mean influence (normalized)") +
  theme_minimal() +
  coord_fixed() +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5 
    )
  )+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.33),
        axis.text.y = element_blank(), # Make text in illustrator so ggsave is aligned
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        legend.title = element_text(angle = 90),
        legend.title.position = "top",
        legend.key.height = unit(1.02, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.position = "right",
        legend.direction = "vertical")

print(class_plot)

if(write_plots){
  ggsave(paste0(output_path, "Figure_4S3_5HTN_no_bSEL.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
  ggsave(paste0(output_path, "Figure_4S3_5HTN_no_bSEL.pdf"), 
         plot = class_plot, width = 20, height = 16, units = "cm")
}






