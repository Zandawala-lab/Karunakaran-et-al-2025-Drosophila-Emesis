# Karunakaran-et-al-2025-Drosophila-Emesis
Code to replicate scRNA analysis and Figure S4-3 connectomic analysis. https://doi.org/10.1101/2025.10.01.679515

## Requirements
Files corresponding to the figures numbers may be ran in their entirety to generate plots. The outputs of these files have been included in the /output/ folder.

### Figure_4S3_A-D.R
From [Codex v783](https://codex.flywire.ai/api/download?dataset=fafb) (09/20/25), download and place in /input/ (folder must be created):    
- connections.csv, classification.csv, consolidated_cell_types.csv  
- fafb_v783_princeton_synapse_table.csv
  
From [GitHub funkelab/drosophila_neurotransmitters](https://github.com/funkelab/drosophila_neurotransmitters) (09/20/25), download and place in /input/:  
- gt_data.csv

### Figure_4S3_E-F.R, Figure_4S3_E-F_Plotting.R
From [Codex v783](https://codex.flywire.ai/api/download?dataset=fafb) (09/20/25), download and place in /input/:  
- connections.csv, classification.csv, consolidated_cell_types.csv
  
From [GitHub funkelab/drosophila_neurotransmitters](https://github.com/funkelab/drosophila_neurotransmitters) (09/20/25), download and place in /input/:  
- gt_data.csv

### Figure_2-6_scRNA.R
From [Zenodo](https://zenodo.org/uploads/18261450), download and place in /input/:  
- Gut_seurat.rda, amine_seurat.rds, Brain_Waddell_seurat.rds

## Questions:
Address questions or requests for files to mzandawala[at]unr.edu or zeyuc[at]unr.edu.
