#-----Figure_2-6_scRNA.R-----
#-------------------------------------------------------------------------------
# Karunakaran et al. - single cell / single nucleus analysis 
# A gut-brain axis for aversive interoception drives innate and anticipatory 
# emesis in Drosophila
#-------------------------------------------------------------------------------
# This script reproduces the figures used in the manuscript. It:
#  1) loads Seurat objects (gut, aminergic neurons, Kenyon cells),
#  2) prepares subsets used for plotting (enterocytes, enteroendocrine cells,
#     Kenyon cell subtypes),
#  3) creates DimPlot and DotPlot figures and saves them to disk.
#-------------------------------------------------------------------------------

# Load packages ----------------------------------------------------------------
library(Seurat)
library(ggplot2)

# General variables -----------------------------------------------------------
input_path <- "./input/"
output_path <- "./output/"
set.seed(1234)

# Load datasets ----------------------------------------------------------------
#load Gut single nucleus data
load(file.path(input_path, "Gut_seurat.rda"))

#load serotonin and dopamine neuron clusters
amine_seurat <- readRDS(paste0(input_path,"amine_seurat.rds"))

#load Kenyon cell data
brain_waddell_seurat <- readRDS(file = paste0(input_path,"Brain_Waddell_seurat.rds"))

# Subset data for plotting ---------------------------------------------------
# Define annotations and capture cell IDs for subsets used in figures.

eec_annotation <- "enteroendocrine cell"
ec_annotation <- c(
  "enterocyte of anterior adult midgut epithelium",
  "adult differentiating enterocyte",
  "enterocyte of posterior adult midgut epithelium",
  "adult midgut enterocyte",
  "enterocyte-like"
)

# Subset cell IDs 
eec_ids <- Cells(subset(x = Gut_seurat, subset = annotation == eec_annotation))
ec_ids <- Cells(subset(x = Gut_seurat, subset = annotation %in% ec_annotation))

# Subset seurat
ec <- subset(x = Gut_seurat, subset = annotation %in% ec_annotation)
message(sprintf("Number of EC cells: %d", length(Cells(ec))))

# Set cell type
ec@meta.data$celltype <- "EC"

eec <- subset(x = Gut_seurat, subset = annotation == eec_annotation)
message(sprintf("Number of EEC cells: %d", length(Cells(eec))))

# Figure 2D -------------------------------------------------------------------
# t-SNE plot highlighting EC and EEC clusters. 
plot_a <- DimPlot(
  Gut_seurat,
  reduction = "tsne",
  group.by = "annotation",
  pt.size = 0.2,
  label = FALSE,
  cells.highlight = list(eec_ids, ec_ids),
  cols = "darkgray",
  cols.highlight = c("blue", "red"),
  repel = TRUE
) +
  NoLegend() +
  labs(title = NULL)

ggsave(
  filename = file.path(output_path, "Figure_2D_tSNE.pdf"),
  plot = plot_a,
  width = 10,
  height = 10,
  units = "cm",
  dpi = 300
)

# DotPlot showing genes expressed in ECs -------------------------------------
genes_to_plot_ec <- c(
  "Nox", "Duox", "AstA", "Mip", "AstC", "CCHa1",
  "CCHa2", "Dh31", "NPF", "Orcokinin", "Tk"
)

plot_b <- DotPlot(
  object = ec,
  features = genes_to_plot_ec,
  group.by = "celltype",
  dot.min = 0.05,
  cols = c("#FFFFFF", "blue")
) +
  xlab(NULL) +
  ylab(NULL) +
  coord_flip() +
  scale_size(range = c(1, 8)) +
  theme(panel.grid.major = element_line(color = "black", size = 0.1, linetype = 2)) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "italic")
  )

ggsave(
  filename = file.path(output_path, "Figure_2D_dotplot.pdf"),
  plot = plot_b,
  width = 8.5,
  height = 8,
  units = "cm",
  dpi = 300
)

# Figure 2G ------------------------------------------------------------------
# DotPlot showing TRP genes expressed in EECs
trp_genes_eec <- c(
  "Trpm","TrpA1","nompC","pain","trp","trpl","Trpγ",
  "pyx","wtrw","iav","nan","Pkd2","Trpml","CG42638"
)

plot_c <- DotPlot(
  object = eec,
  features = trp_genes_eec,
  group.by = "annotation",
  dot.min = 0.05,
  cols = c("#FFFFFF", "#FF0000")
) +
  xlab(NULL) +
  ylab(NULL) +
  coord_flip() +
  scale_size(range = c(1, 8)) +
  theme(panel.grid.major = element_line(color = "black", size = 0.1, linetype = 2)) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "italic")
  )

ggsave(
  filename = file.path(output_path, "Figure_2G.pdf"),
  plot = plot_c,
  width = 8.5,
  height = 9,
  units = "cm",
  dpi = 300
)

# Figure 3A ------------------------------------------------------------------
# DotPlot showing neuropeptides expressed in EECs
neuropeptides_eec <- c(
  "pros", "amon", "svr", "Phm", "AstA", "AstB", "AstC", "CCHa1",
  "CCHa2", "Dh31", "NPF", "Orcokinin", "Tk"
)

plot_d <- DotPlot(
  object = eec,
  features = neuropeptides_eec,
  group.by = "annotation",
  dot.min = 0.05,
  cols = c("#FFFFFF", "#FF0000")
) +
  xlab(NULL) +
  ylab(NULL) +
  coord_flip() +
  scale_size(range = c(1, 8)) +
  theme(panel.grid.major = element_line(color = "black", size = 0.1, linetype = 2)) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "italic")
  )

ggsave(
  filename = file.path(output_path, "Figure_3A.pdf"),
  plot = plot_d,
  width = 8.5,
  height = 9,
  units = "cm",
  dpi = 300
)

# Figure 4D ------------------------------------------------------------------
# DotPlot showing dopamine and serotonin pathway genes in aminergic neurons
receptors_amine <- c(
  "SerT", "Trh", "Vmat", "Ddc", "DAT", "ple", "Dop1R1", "Dop1R2",
  "Dop2R", "DopEcR", "5-HT1A", "5-HT1B", "5-HT2A", "5-HT2B", "5-HT7"
)

plot_e <- DotPlot(
  object = amine_seurat,
  features = receptors_amine,
  dot.min = 0.05,
  scale = FALSE
) +
  xlab(NULL) +
  ylab(NULL) +
  coord_flip() +
  scale_size(range = c(1, 8)) +
  scale_colour_gradient2(low = "grey", mid = "darkgrey", high = "black") +
  theme(panel.grid.major = element_line(colour = "black", size = 0.1, linetype = 2)) +
  theme(axis.text.y = element_text(face = "italic")) +
  RotatedAxis()

ggsave(
  filename = file.path(output_path, "Figure_4D.pdf"),
  plot = plot_e,
  width = 10,
  height = 14,
  units = "cm",
  dpi = 300
)

# Figure 5A ------------------------------------------------------------------
# DotPlot showing neuropeptide receptors in aminergic neurons
neuropeptide_receptors_amine <- c(
  "SerT", "Trh", "Vmat", "Ddc", "DAT", "ple", "AstA-R1", "AstA-R2",
  "SPR", "AstC-R1", "AstC-R2", "CCHa1-R", "CCHa2-R", "Dh31-R", "NPFR",
  "Tkr86C", "TkR99D"
)

plot_f <- DotPlot(
  object = amine_seurat,
  features = neuropeptide_receptors_amine,
  dot.min = 0.05,
  scale = FALSE
) +
  xlab(NULL) +
  ylab(NULL) +
  coord_flip() +
  scale_size(range = c(1, 8)) +
  scale_colour_gradient2(low = "grey", mid = "darkgrey", high = "black") +
  theme(panel.grid.major = element_line(colour = "black", size = 0.1, linetype = 2)) +
  theme(axis.text.y = element_text(face = "italic")) +
  RotatedAxis()

ggsave(
  filename = file.path(output_path, "Figure_5A.pdf"),
  plot = plot_f,
  width = 10,
  height = 14,
  units = "cm",
  dpi = 300
)

# Supplementary Figure S6-2E -------------------------------------------------
# Prepare Kenyon cell subset and relabel clusters for readability.

# Preserve old identities in metadata for traceability
brain_waddell_seurat$old_ident <- Idents(object = brain_waddell_seurat)

# Rename identity classes 
brain_waddell_seurat <- RenameIdents(
  object = brain_waddell_seurat,
  `9` = "KCalphabeta",
  `10` = "KCgamma",
  `13` = "KCa'b'"
)

# Save new cluster labels in metadata
brain_waddell_seurat@meta.data$new_clusters <- Idents(brain_waddell_seurat)

# Add an 'assay' metadata column to identify the source of these data
brain_waddell_seurat@meta.data$assay <- "Waddell"

# Subset Kenyon cells (only include the named subclusters)
kc_annotation <- c("KCalphabeta", "KCa'b'", "KCgamma")
kc_subset <- subset(x = brain_waddell_seurat, subset = new_clusters %in% kc_annotation)

# DotPlot showing dopamine and serotonin receptors expressed in KCs
kc_genes <- c(
  "ab", "sNPF", "ey", "Dop1R1", "Dop1R2", "Dop2R", "DopEcR",
  "5-HT1A", "5-HT1B", "5-HT2A", "5-HT2B", "5-HT7"
)

plot_g <- DotPlot(
  object = kc_subset,
  features = kc_genes,
  dot.min = 0.05,
  scale = FALSE
) +
  xlab(NULL) +
  ylab(NULL) +
  coord_flip() +
  scale_size(range = c(1, 8)) +
  scale_colour_gradient2(low = "grey", mid = "darkgrey", high = "black") +
  theme(panel.grid.major = element_line(colour = "black", size = 0.1, linetype = 2)) +
  theme(axis.text.y = element_text(face = "italic")) +
  RotatedAxis()

ggsave(
  filename = file.path(output_path, "Figure_S6_2E.pdf"),
  plot = plot_g,
  width = 10,
  height = 12,
  units = "cm",
  dpi = 300
)

# End of script ----------------------------------------------------------------
message("All plots generated and saved to: ", normalizePath(output_path, mustWork = FALSE))
