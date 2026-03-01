# README document for manuscript "Heterogeneous gene network estimation for single-cell transcriptomic data via a joint regularized deep neural network"






##  Data

###  Directory Structure Description
The data used in this study is stored in the `data` folder of the repository, which is divided into two core subdirectories:
  - `RealData`: Contains the five publicly available single-cell transcriptomic datasets described below, including preprocessed expression matrices, cell type annotations, and downstream analysis results.
- `Application-Based Simulation data`: Includes three biologically realistic simulation datasets (GSD, scMultiSim-T3, SERGIO-DS1) derived from real biological networks, used for method benchmarking under controlled conditions.

###  Abstract
The single-cell transcriptomic datasets used in this study include five publicly available datasets:
  - Human lung adenocarcinoma (LUAD) cell lines, derived from three cell lines (HCC827, H1975, H2228) with sequencing data from multiple platforms.
- Human peripheral blood mononuclear cells (PBMC), consisting of purified immune cell subtypes (CD56+ natural killer cells, CD19+ B cells, CD4+/CD25+ regulatory T cells).
- Mouse embryonic stem cells (mESCs), sequenced under three culture conditions (2i, a2i, lif).
- Mouse liver cells from the Mouse Cell Atlas (MCA), including diverse cell types involved in immune regulation and tissue function.
- Mouse uterus cells from MCA, containing cell subtypes such as endothelial, glandular, and immune cells.

These datasets vary in sample size (from 704 to 12,418 cells), gene count (from 6,237 to 15,618 genes), and zero-inflation rates (from 30.77% to 91.12%), enabling comprehensive evaluation of the proposed JRDNN-KM method.

###  Availability
All datasets are publicly available for download via the following links:
  
  - **LUAD cell lines**: Available at https://github.com/LuyiTian/sc_mixology (Tian et al., 2019).
- **PBMC**: Available at https://github.com/10XGenomics/single-cell-3prime-paper (Zheng et al., 2017).
- **mESCs**: Available in the ArrayExpress database under accession number E-MTAB-2600 (Kolodziejczyk et al., 2015).
- **Mouse liver and uterus cells**: Available at https://figshare.com/articles/dataset/MCA_DGE_Data/5435866/8 (Han et al., 2018).

###  Real Dataset Description

#### Human Lung Adenocarcinoma (LUAD) Cell Lines
- **Contributors**: Tian, L., Dong, X., Freytag, S., et al.
- **Citation**: Tian, L., et al. (2019). Benchmarking single cell RNA-sequencing analysis pipelines using mixture control experiments. *Nature Methods*, 16(6), 479–487.
- **File format**: Gene expression matrices (count data), metadata with cell line labels.
- **Preprocessing**:
  - Cells expressing fewer than 500 genes and genes expressed in fewer than 100 cells were removed.
- Batch effect correction was performed using MNN (Haghverdi et al., 2018).
- Library size normalization was applied using total count normalization (TC) (Cole et al., 2019).
- 100 highly variable genes (HVGs) were selected using the variance stabilizing transformation (vst) method in the Seurat R package (Hafemeister & Satija, 2019).

#### Human Peripheral Blood Mononuclear Cells (PBMC)
- **Contributors**: Zheng, G. X., Terry, J. M., Belgrader, P., et al.
- **Citation**: Zheng, G. X., et al. (2017). Massively parallel digital transcriptional profiling of single cells. *Nature Communications*, 8(1), 14049.
- **File format**: Filtered barcode matrix (HDF5), spatial metadata, and cell type annotations.
- **Preprocessing**:
  - Focused on three cell types: CD56+ natural killer cells, CD19+ B cells, and CD4+/CD25+ regulatory T cells.
- Quality control: Cells with <500 expressed genes and genes with <100 expressing cells were filtered out.
- Normalization and batch correction as described for LUAD, followed by selection of 100 HVGs.

#### Mouse Embryonic Stem Cells (mESCs)
- **Contributors**: Kolodziejczyk, A. A., Kim, J. K., Tsang, J. C., et al.
- **Citation**: Kolodziejczyk, A. A., et al. (2015). Single cell RNA-sequencing of pluripotent states unlocks modular transcriptional variation. *Cell Stem Cell*, 17(4), 471–485.
- **File format**: Gene expression matrices, metadata with culture condition labels (2i, a2i, lif).
- **Preprocessing**:
  - Quality control filters: Cells with <500 expressed genes and genes with <100 expressing cells were excluded.
- Normalized using Seurat’s `LogNormalize` and batch-corrected with MNN.
- 100 HVGs were selected using vst.

#### Mouse Liver Cells (MCA)
- **Contributors**: Han, X., Wang, R., Zhou, Y., et al.
- **Citation**: Han, X., et al. (2018). Mapping the mouse cell atlas by microwell-seq. *Cell*, 172(5), 1091–1107.
- **File format**: DGE (digital gene expression) matrices, cell type annotations.
- **Preprocessing**:
  - Cell types representing <5% of the total population were excluded, resulting in 5 cell subgroups.
- Quality control and normalization steps identical to above, with 100 HVGs selected.

#### Mouse Uterus Cells (MCA)
- **Contributors**: Han, X., Wang, R., Zhou, Y., et al.
- **Citation**: Han, X., et al. (2018). Mapping the mouse cell atlas by microwell-seq. *Cell*, 172(5), 1091–1107.
- **File format**: DGE matrices, cell type annotations (endothelial, glandular, macrophage, etc.).
- **Preprocessing**:
  - Cell types with <5% representation were excluded, retaining 6 cell subgroups.
- Quality control, normalization, and HVG selection as described for other datasets.

###  Application-Based Simulation Data Description
To further validate the performance of JRDNN-KM under biologically realistic conditions, we employed three application-based simulation datasets derived from real biological networks, which are stored in the `Application-Based Simulation data` subdirectory. These datasets incorporate practical biological processes and provide a challenging benchmark for network inference and cell subgroup identification algorithms.

| Dataset | Core Characteristics | Citation |
  |---------|---------------------|----------|
  | **Gonadal Sex Determination (GSD)** | Based on a Boolean network model of gonadal sex determination, capturing the bipotential differentiation of gonads into male (Sertoli cells) or female (Granulosa cells) cell fates. | Pratapa et al. (2020), *Nature Methods* |
  | **scMultiSim-T3** | Generated by scMultiSim, a multi-modality single-cell data simulator integrating key biological factors including gene regulatory networks (GRNs) and cell-cell interactions (CCIs). | Li et al. (2025), *Nature Methods* |
  | **SERGIO-DS1** | Constructed by the SERGIO simulator based on an *E. coli* derived GRN, simulating single-cell gene expression profiles with realistic biological noise. | Dibaeinia & Sinha (2020), *Cell Systems* |
  
  
  
  
  ##  Reproduce the results
  
  ###  Abstract  
  The `code_and_data` directory contains all core resources to reproduce results:  
  - `Simulation Generate/`: Standalone code to generate custom synthetic single-cell expression data and network topologies (SBM, scale-free, star-chain hybrid networks) with tunable biological parameters.  
- `RealDemo/`: A streamlined, self-contained pipeline for rapid validation of JRDNN-KM and competing methods on the LUAD dataset (includes preprocessed input data and method implementation scripts).  
- `Simulations/`: Scripts and results for simulation studies.  
- `RealData/`: Scripts and results for real single-cell data analyses (including full and subsampling-based robustness validation).  



###  Instructions for Use  
All results in the manuscript (simulations and real data analyses) are fully reproducible. The `code_and_data` directory is structured hierarchically to facilitate step-by-step reproduction, with strict path management via the `here` package (**critical: set the working directory root to `code_and_data` before running any scripts**).  

#### Simulation Generate  

This directory contains standalone code to **generate custom simulation frameworks** (network topologies + synthetic single-cell expression data) from scratch, with full control over key simulation parameters:  
  
  | Exact File List | Content Description |
  |-----------------|---------------------|
  | `generate_SBM.R` | Generates synthetic networks based on the Stochastic Block Model (SBM) (a modular network topology with predefined community structures): <br> - Outputs: Adjacency matrices defining gene-gene interaction networks (SBM-based) <br> - Derives cell × gene expression matrices aligned with SBM network topology <br> - Tunable parameters: <br>   ✔ Zero-inflation rate  <br>   ✔ Non-linear complexity of gene-gene interactions <br>   ✔ Number of network communities  <br>   ✔ Edge density within/between blocks |
  | `generate_scalefree.R` | Generates scale-free networks (power-law degree distribution, mimicking biological gene networks): <br> - Outputs: Adjacency matrices for scale-free gene networks <br> - Derives cell × gene expression matrices consistent with scale-free topology <br> - Tunable parameters: <br>   ✔ Zero-inflation rate <br>   ✔ Non-linear complexity of gene regulation <br>   ✔ Power-law exponent (degree distribution) <br>   ✔ Total number of genes/nodes |
  | `generate_starchain.R` | Generates star-chain hybrid networks (combining star-shaped hub networks and linear chain sub-networks): <br> - Outputs: Adjacency matrices for hybrid gene networks <br> - Derives cell × gene expression matrices matching hybrid topology <br> - Tunable parameters: <br>   ✔ Zero-inflation rate <br>   ✔ Non-linear complexity of hub-gene regulation <br>   ✔ Proportion of star vs. chain sub-networks <br>   ✔ Number of hub genes |
  
  To generate custom simulation data:  
  1. Set working directory to `code_and_data/Simulation Generate/`;  
2. Modify parameter values (zero-inflation, non-linearity) at the top of the `.R` script;  

#### RealDemo  
This directory offers a streamlined, self-contained pipeline to reproduce analysis results (for the LUAD dataset) of the proposed JRDNN-KM method and competing methods, designed for quick validation:  
  
  It includes **input data files** and **method implementation scripts**:  
  - Input data:  
  - `luad.csv` (preprocessed cell×gene expression matrix) + `label.txt` (ground-truth cell line labels for benchmarking) – a paired set of expression and annotation files;  
- `luad.Rdata` – a consolidated R data file containing both the preprocessed cell×gene expression matrix and corresponding cell line labels.  
- Method implementation scripts:  
  - The `JRDNN-KM` folder houses the core code for the proposed method.  
- Competing methods are implemented via dedicated scripts: `CSCORE+SPQN.R` (for the CSCORE+SPQN pipeline), `locCSN.py` (for locCSN), `Normalisr.py` (for Normalisr), and `other_competing_methods.R` (a wrapper for additional comparative methods).  



####  Simulations  
Updated with exact file mappings for simulation studies:  
  
  | Subdirectory | Exact File List | Content Description |
  |--------------|-----------------|---------------------|
  | `code/` | `ComplexGenerative_draw_Figure5.R` | Generates plots for simulation results under complex generative mechanisms (main Figure 5) |
  |          | `Simulation_draw_FigureS23-S27.R` | Generates plots for standard simulation scenarios (Supplementary Figures S23–S27) |
  | `result_data/` | `Complex generative mechanisms/` (directory) | Simulation results under complex generative settings (e.g., non-linear gene interactions, dynamic subgroup structures) |
  |               | `simulation result/` (directory) | Results from standard simulation scenarios (balanced/imbalanced subgroups, shared network information, linear interactions, high dropout, signal-noise variation) |
  
  
  
  ####  RealData
  #####  Total sample  
  This subdirectory contains full-scale analysis results and plotting scripts for real single-cell datasets, with the exact file breakdown as follows:  
  
  | Subdirectory | Exact File List | Content Description |
  |--------------|-----------------|---------------------|
  | `code/` | `ARI_draw_Figure2B,S4.R` | Generates ARI (Adjusted Rand Index) plots for main Figure 2B and Supplementary Figure S4 |
  |          | `ARI_NMI_draw_Figure2A.R` | Generates ARI vs. NMI (Normalized Mutual Information) performance plots for main Figure 2A |
  |          | `Convergence_draw_FigureS1.r` | Generates model convergence curves for Supplementary Figure S1 |
  |          | `Convergence_draw_FigureS2.R` | Generates model convergence curves for Supplementary Figure S2 |
  |          | `Graphreal_draw_Figure3,S5,S6.R` | Generates network structure visualizations for real datasets (main Figure 3, Supplementary Figures S5–S6) |
  |          | `ModelCheck_FigureS3.R` | Generates model validation/ diagnostic plots for Supplementary Figure S3 |
  |          | `Network_draw_Figure4.R` | Generates detailed gene network visualizations for main Figure 4 |
  |          | `Network_draw_FigureS7.R`- `Network_draw_FigureS10.R`| Generates detailed gene network visualizations for Supplementary Figure S7-S10 |
  |          | `Sensitivity_draw_FigureS17.R` | Generates sensitivity analysis plots for Supplementary Figure S17 |
  |          | `Upset_draw_FigureS11.R` | Generates Upset plots for cross-method comparison (Supplementary Figure S11) |
  | `result_data/` | `all_experiments.RData`, `centers_evolution.RData`, `Estimated_networks/` (directory), `loss_data.RData`, `model checking/` (directory), `plot_data.RData`, `Result_ground_truth/` (directory), `Sensitive/` (directory), `subgroup_results_BLGGM.csv`, `subgroup_results_JRDNN_KM.csv`, `subgroup_results_SC3.csv`, `subgroup_results_Seurat.csv` | Precomputed results to support plotting scripts in the `code/` directory, including aggregated experiment data, model training metrics, inferred gene networks, subgroup assignment results, and diagnostic data for all 5 real datasets (LUAD, PBMC, mESCs, mouse liver, mouse uterus). |
  
  **Key Note**: All R scripts in `code/` use the `here` package to reference files in `result_data/` (e.g., `here("RealData", "Total sample", "result_data", "Estimated_networks")`). Ensure the working directory is set to the `code_and_data` root before execution.  


#####  Sub Sample  
Used for subsampling-based robustness validation (50 subsamples per dataset to evaluate method stability), with the following structure:  
  
  | Subdirectory | Exact File List | Content Description |
  |--------------|-----------------|---------------------|
  | `code/` | `sub-sample_draw_FigureS18-S21.R` | Generates robustness plots (ARI/NMI distributions across subsamples) for Supplementary Figures S18–S21 |
  | `result_data/` | `Cluster/` (directory) | Subsampling results for clustering performance (ARI/NMI scores, cluster consistency metrics) across 5 real datasets |
  |               | `Network estimation/` (directory) | Subsampling results for network inference (edge consistency, modularity stability, edge weight variability) across 5 real datasets |
  
  
  
  
  
  
  
  
  
  