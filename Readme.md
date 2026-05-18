---
output:
  pdf_document:
    latex_engine: xelatex
---

# README document for manuscript ” Integrated Bayesian non-parametric spatial modeling for cross-sample identification of spatially variable genes”

## Data

All datasets (real spatial transcriptomic data and Simulation data) used in
this study are publicly available and stored in the GitHub repository zhoumeng123456/IBaySVG-analysis.

### Directory Structure Description

The data used in this study is stored in the `data/` folder of the
repository, which is divided into two core subdirectories:

- `Realdataset`: Contains the two publicly available spatial transcriptomic datasets described below, 
  including both raw and preprocessed gene expression matrices, benchmark scRNA/snRNA-seq reference 
  matrices, spatial coordinates, and cell-type deconvolution proportions.
- `Simulation data`: Includes simulation datasets spanning three
  canonical spatial patterns and Three additional heterogeneous spatial
  structures. For each canonical pattern, multiple data-integration
  scenarios are constructed, with and without covariates, to
  systematically evaluate model performance under varying structural
  complexities.

### Abstract

The spatial transcriptomic datasets used in this study include two
publicly available datasets:

- Human dorsolateral prefrontal cortex (DLPFC), comprising laminar
  cortical structures and evaluated under within-donor and cross-donor
  section integration scenarios.
- Human squamous cell carcinoma (SCC), a representative tumor tissue
  dataset characterized by spatial heterogeneity.

These datasets differ in spot resolution (~1,000–4,000 spots), tissue
architecture, and sparsity levels, providing complementary benchmarks
for evaluating model performance under diverse spatial transcriptomic
settings.

### Availability

All datasets are publicly available for download via the following
links:

- **DLPFC**: The raw datasets used in this study include both spatial transcriptomics (ST) and single-nucleus RNA-seq (snRNA-seq) references (Maynard et al., 2021). URLs: 
  - ST Dataset: Available at the GitHub repository <https://github.com/LieberInstitute/HumanPilot> 
  - snRNA-seq Reference: <https://www.dropbox.com/s/5919zt00vm1ht8e/sce_DLPFC_annotated.zip?dl=1>.
- **SCC**: Available at the Gene Expression Omnibus (GEO) under accession
  number GSE144240 and GSE144236(Ji et al., 2020). URLs: 
  - ST Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144240
  - scRNA-seq Reference: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144236
    
### Real Dataset Description

#### Human dorsolateral prefrontal cortex (DLPFC)

- **Contributors**: Maynard, K.R. and Collado-Torres, L. , et al.
- **Citation**: Maynard, K., Collado-Torres, L. et al. (2021),
  ‘Transcriptome-scale spatial gene expression in the human dorsolateral
  prefrontal cortex’, Nature Neuroscience 24, 425–436.
- **File format**: Gene expression matrix (HDF5), spatial coordinates
  (TXT), snRNA-seq count matrix (HDF5), and snRNA-seq cell-type and gene
  annotations (RDS).
- **Preprocessing**:
  - Quality control: Cells with \<500 expressed genes and genes with
    \<100 expressing cells were filtered out.
  - Selecting the top 8,000 highly variable genes (HVGs) per slice using
    the Seurat package and identifying the consensus HVGs across all
    four slices.
  - Focused on seven cell types (Astro, EndoMural, MicroOligo, Oligo,
    OPC, Excit, and Inhib), with spot-wise proportions estimated using
    Redeconve (Zhou et al. 2023).

#### Human squamous cell carcinoma (SCC)

- **Contributors**: Ji, A. L., Rubin, A. J., Thrane, K., Jiang, S.,
  Reynolds, D. L., Meyers, R. M., Guo, M., George, J., Mollbrink, A. et
  al.
- **Citation**: Ji, A. L., Rubin, A. J., Thrane, K., Jiang, S.,
  Reynolds, D. L., Meyers, R. M., Guo, M., George, J., Mollbrink, A. et
  al. (2020), ‘Multimodal analysis of composition and spatial
  architecture in human squamous cell carcinoma’, Cell 182(2),
  497–514.e22.
- **File format**: Gene expression matrix (TSV), spatial coordinates
  (TSV), scRNA-seq reference counts (TXT), and scRNA-seq reference
  cell-type annotations (TXT).
- **Preprocessing**:
  - Quality control: Cells with \<100 expressed genes and genes with
    \<10 expressing cells were filtered out.
  - Selecting the top 8,000 highly variable genes (HVGs) per slice using
    the Seurat package and identifying the consensus HVGs across all
    four slices.
  - Focused on three major cell types (Keratinocyte, TSK, and
    Melanocyte), with spot-wise proportions estimated using Redeconve.

### Simulation Data Description

To further validate the performance of IBaySVG under controlled
conditions, we employed simulation datasets spanning three canonical
spatial patterns and nine additional heterogeneous spatial structures,
which are stored in the Simulation data subdirectory. For each canonical
pattern, multiple data-integration scenarios were constructed, with and
without covariates, enabling systematic evaluation of the method under
varying spatial architectures and integration complexities. These
simulations provide a structured and challenging benchmark for assessing
the accuracy and robustness of IBaySVG in identifying spatially variable
genes.

| Dataset | Core Characteristics |
|----|----|
| **Three canonical spatial patterns** | Each pattern includes 12 integration scenarios with varying signal strengths (high to low) and zero-inflation rates (high, medium, low), plus additional scenarios without covariates. |
| **Three additional heterogeneous spatial structures** | Comprised of three classes of complex spatial architectures: (i) zero-inflated nearest-neighbor Gaussian process (ZINNGP) capturing covariance-based spatial dependencies, (ii) ZINB-NonSpa hybrid models with linear-focal, linear-periodic, and focal-periodic fusion patterns, and (iii) enhanced ZINB-NonSpa patterns including sigmoidal activation and four polynomial functional forms. |

All the simulated data are generated by ZINB distributions or Gaussian
process framework with covariates from Dirichlet sampling to evaluate
model performance under diverse and non-canonical spatial
configurations.

## Reproduce the results

### Abstract

The `code_and_data` directory contains all core resources to reproduce
results:

- `Demo/`: the main IBaySVG implementation, usage instructions for 
  both IBaySVG and comparison methods, and a small example dataset 
  for testing the full analysis workflow.
- `Simulation Generate/`: Standalone code to generate synthetic spatial
  expression data across multiple canonical (linear, focal, periodic)
  and heterogeneous spatial patterns, with tunable signal strengths,
  zero-inflation rates, and covariate configurations.
- `Simulations/`: Scripts and results for simulation studies(including
  parameter sensitivity analyses).  
- `RealData/`: Scripts and results for real spatial transcriptomic
  datasets (including stability assessments, model checking and
  parameter sensitivity analysis).

### Instructions for Use

All results in the manuscript (simulations and real data analyses) are
fully reproducible. The `code_and_data` directory is structured
hierarchically to facilitate step-by-step reproduction, with strict path
management via the `here` package (**critical: set the working directory
root to `code_and_data` before running any scripts**).

#### Simulation Generate

This directory contains standalone code to **generate custom simulation
frameworks** (network topologies + synthetic single-cell expression
data) from scratch, with full control over key simulation parameters:

| Exact File List | Content Description |
|----|----|
| generate_simu.R | Generates synthetic spatial transcriptomic datasets under a zero-inflated negative binomial (ZINB) model with configurable spatial effect structures: <br> - Outputs: <br> ✔ Spot × gene count matrix (N × G) <br> ✔ Spot × 2 spatial coordinate matrix (N × 2) <br> ✔ Optional spot × covariate matrix (N × J) <br> - Supports canonical and heterogeneous spatial patterns (e.g., linear, focal, periodic, and extended spatial structures) <br> - Tunable parameters: <br> ✔ Proportion of spatially variable (SV) genes <br> ✔ Spatial signal strength for SV/non-SV genes <br> ✔ Zero-inflation rate (dropout level) <br> ✔ Dispersion parameter of the ZINB distribution <br> ✔ Baseline mean expression level <br> ✔ Covariate inclusion and distribution scenarios |

To generate simulation data:  
1. Set working directory to `code_and_data/Simulation Generate/`;  
2. choose the spatial pattern and modify parameter values
(zero-inflation, spatial signal strength, Covariate inclusion) according
to the function documentation;

#### Simulations

Updated with exact file mappings for simulation studies:

| Subdirectory | Exact File List | Content Description |
|---|----------|------|
| `code/` | `simulation_spotregion_draw_Figure_S20.R` | Generates plots for partitioned spot regions considered in the basic simulations.(Supplementary Figures S20) |
|  | `simulation_draw_Figure_S21-S29.R` | Generates plots for basic simulation results under three canonical spatial patterns(Supplementary Figures S21–S29) |
|  | `additional_expression_pattern_Figure_S30.R` | Generates plots for examples of spatial expressions under additional spatial structures(Supplementary Figures S30) |
|  | `Simulation_additional_draw_Figure_S31-S33.R` | Generates plots for simulation results under nine additional spatial structures (Supplementary Figures S31–S33) |
|  | `simulation_nocov_draw_Figure_S34-S36.R` | Generates plots for simulation results under three canonical spatial patterns without covariates (Supplementary Figures S34–S36) |
|  | `sensitive_simulation_draw_Figure_S40-S42.R` | Generates plots for sensitivity analysis of hyperparameters in the simulated dataset under three canonical spatial patterns(Supplementary Figures S40–S42) |
|  | `comparison_between_ZINB_NB_Table_S35-S36.R` | Illustrates the comparison between ZINB and NB Models under simulation scenarios for Table S35-S36 |
|  | `comparison_between_two_statistics_Table_S1-S3.R` | Illustrates the comparison of the performance between two statistics under simulation scenarios for Table S1-S3 |
|  | `Computer_time_S4.R` | Illustrates the average execution time per gene under different spots numbers for Table S4 |
| `result_data/` | `sensitivity analysis of hyperparameters/` (directory) | Contains results from sensitivity analyses of key hyperparameters evaluated on simulated datasets under three canonical spatial patterns (linear, focal, and periodic) |
|  | `simulation result/` (directory) | Contains results from baseline simulation scenarios (linear, focal, and periodic patterns), evaluated both with and without covariate effects, as well as additional simulation settings under nine alternative spatial structures. |
|  | `additional spatial pattern/` (directory) | Contains plots of spatial expression of additional spatial structures selected |
|  | `noinflation/`(directory) | Contains the code and results of IBaySVG without zero inflation structure |
|  | `compare statistics/` (directory) | Contains results of using different statistics in IBaySVG |
||`compute time.csv`|Contains the results of computing time of IBaySVG|
#### Demo 

The `Demo/` folder provides a minimal working example for reproducing the main analysis pipeline and benchmarking procedures.

It contains the following components:

- **`Implement_of_IBaySVG.pdf`**：Provides detailed documentation of the IBaySVG algorithm and instructions for its implementation.
- **`IBaySVG_main.R`**：The main function implementing the proposed IBaySVG method.
- **`Comparison_methods_R.R`**: Implementation of competing methods in the R environment.
- **`Comparison_methods_python.py`**: Implementation of competing methods in the Python environment.
- **Example data files**: The remaining files provide example datasets in both R and CSV formats:
  - `data_example.RData` — R-format example dataset.
  - `matrix*_count_example.csv` — UMI count matrices (CSV format).
  - `matrix*_position_example.csv` — Spatial coordinate matrices (CSV format).

The R-formatted files are intended for use in the R environment, while the CSV files are provided to ensure compatibility with Python-based analysis.



#### RealData

This directory contains full-scale analysis results and plotting scripts
for real spatial transcriptomic datasets, with the exact file breakdown
as follows:

| Subdirectory | Exact File List | Content Description |
|---|----------|------|
| `code/` | `spatial_expression_draw_Figure_S1.R` | Generates plots of Examples of linear, focal, periodic, and one more complex spatial patterns for Supplementary Figure S1 |
|  | `venn_plot_draw_Figure_S2-S4.R` | Generates venn plots of SV genes identified in real datasets for Supplementary Figure S2-S4 |
|  | `celltype_plot_draw_Figure_S5-S7.R` | Generates plots of distribution of cellular composition in real datasets for Supplementary Figure S5-S7 |
|  | `Model_diagnostics_draw_Figure_S8_Table_S5-S7.R` | Generates plots of likelihood ratio test comparison of ZINB and NB models and results of voung test comparing linear/focal/period model in real datasets for Supplementary Figure S8 |
|  | `upset_draw_Figure_S9-S11.R` | Generates upset plots of SV gene identification through multiple integration strategie in real datasets for Supplementary Figure S9-S11 |
|  | `gene_cluster_draw_Figure_S14-S19_Table_S17-S19.R` | Generates expression and heatmap plots of dominant gene clusters of SV genes identified by IBaySVG in real datasets for Supplementary Figure S14-S19 and the enrichment result for supplementary Table S17-S19 |
|  | `sensitive_realdata_draw_Figure_S37-S39.R` | Generates plots for sensitivity analysis of hyperparameters in real datasets (Supplementary Figures S37–S39) |
|  | `spatial_domain_draw_Figure3.R` | Generates spatial domains derived from SV genes identified by IBaySVG for main Figure3 |
|  | `Stability_analysis_Table_S8-S13.R` | Illustrates the results of subsampling validation of SV gene detection in real dataset for Table S8-S13 |
|  | `marker_gene_identify_Table_S15-S16.R` | Illustrates the performance comparison of different methods in identifying layer-specific marker for Table S15 and S16 |
|  | `comparison_spotcluster_index_Table_S20-S34.R` | Illustrates the results of comparative evaluation of spot clustering performance in real dataset for Table S20-S34 |
|  | `comparative_analysis_Figure2_S12_S13.R`|Illustrates the process of comparative analysis of the identified SV genes and generates the plots of spatial expression for Figure2,S12 and S13|
|  | `data_preprocess.R`|Illustrates the comprehensive data cleaning, quality control, and preprocessing pipeline starting from the raw datasets. |
| `result_data/` | `plot of figure_S1/`(directory)| Contains prepared plots of figure S1|
|                |`realdata svgene/`(directory)| Contains identified SV genes of each method |
||`genecluster/`(directory)|Contains the results of inferred gene cluster|
||`spotcluster/`(directory)|Contains the results of inferred spot cluster|
||`stability/`(directory)|Contains the results of stability analysis of each method |
||`realdata dataset/`(directory)|Contains the preprocessed dataset used in realdata analysis|
||`model check/`(directory)|Contains the results of model checking|
||`sensitive analysis/`(directory)|Contains the results of sensitive analysis of hyperparameters|
||`upset plot/`(directory)|Contains the preprocessed upset plot|
||`benchmark/`(directory)|Contains the benchmark results obtained from Zeng et al.(2012) and Maynard et al.(2021) |

#### Downloading Large Reproducibility Files

Some intermediate outputs and large data files (e.g.,  simulation and realdata analysis outputs, and processed count matrices) are excluded from this repository due to their size. To fully reproduce the results in the paper:

1. Download the required files from: https://github.com/zhoumeng123456/IBaySVG-analysis

2. Specifically, download and replace the following directories:
  
   - `data/`
   - `RealData/result_data/benchmark/`
   - `RealData/result_data/genecluster/`
   - `RealData/result_data/spotcluster/`
   - `RealData/result_data/plot of figure_S1/`
   - `Simulations/result_data/additional spatial pattern/`
   - `Simulations/result_data/noinflation/example_linear_inf5.RData`

3. Place the downloaded folders in the same relative directory structure as provided in this repository.

4. Do not modify folder names, as the scripts rely on relative paths.

All analysis scripts will run without modification once the directory structure is preserved.


### Reproducibility Workflow

The recommended order for reproducing the results in this manuscript is:

#### Step 1: Download Required Files from GitHub before running any scripts

#### Step 2: Set the Working Directory

All R scripts rely on the `here` package for path management. Set the working directory to the `code_and_data` root directory before running any scripts.

#### Step 3: Run the Demo Example (Quick Start)

Users are encouraged to begin with the example provided in the `Demo/` directory. This minimal working example demonstrates:

- Implementation of the IBaySVG method by reading the `Implement_of_IBaySVG.pdf` and running `IBaySVG_main.R`
- Implementation of the  Comparative  method by running: 
  - `Comparison_methods_R.R`
  - `Comparison_methods_python.py`
- Expected input and output formats.

#### Step 4: Reproduce Simulation Studies

Simulation experiments can be reproduced using the `Simulation Generate/` and `Simulations/` directory. Users may either regenerate the simulation datasets or directly reproduce results using the provided data files.

- Custom simulation data can be generated using `Simulation Generate/generate_simu.R`;
- or load the specific simulation data from `data/Simulation data/`;
- All simulation figures and tables can be regenerated using scripts in:
  - `Simulations/code/`
  - `Simulations/result_data/`


#### Step 5: Reproduce Real Data Analyses

All real-data analyses are contained in the `RealData/` directory. To fully reproduce the results, we recommend the following execution order:

**(1) Model Diagnostics**

Begin with model checking and diagnostic analyses to validate the distributional assumptions (e.g., ZINB vs. NB models) and spatial functional forms.

Scripts are located in:

- `RealData/code/Model_diagnostics_draw_Figure_S8_Table_S5-S7.R`
- `RealData/result_data/model check/`

**(2) Identification of Spatially Variable (SV) Genes**

Running the IBaySVG and competing methods produced in `Demo/` to identify SV genes based on the preprocessed datasets provided in:

- Option A (Quick Start): Users can directly utilize the preprocessed benchmarks provided in `RealData/result_data/realdata dataset/` for an immediate run.
- Option B (Full Reproduction): Alternatively, reproduction users can completely regenerate these preprocessed datasets from the raw matrices by running the provided data cleaning pipeline `RealData/code/data_preprocess.R`.

The resulting SV gene lists are stored in:

- `RealData/result_data/realdata svgene/`

**(3) Downstream Visualization and Biological Interpretation**

After identifying SV genes, reproduce downstream analyses and visualization, including:

- Venn and UpSet plots: 
  - `RealData/code/venn_plot_draw_Figure_S2-S4.R`
  - `RealData/code/upset_draw_Figure_S9-S11.R`
- Cell-type composition analysis: 
  - `RealData/code/celltype_plot_draw_Figure_S5-S7.R`
- comparative analysis: 
  - `RealData/code/comparative_analysis_Figure2_S12_S13.R`

**(4) Sensitive Analysis, Stability Analysis and Comparative Benchmarking**

Next, reproduce robustness validation and benchmarking analyses, including:

- Sensitive Analysis of hyperparameters: 
  - `RealData/code/sensitive_realdata_draw_Figure_S37-S39.R`
- Subsampling-based stability assessment: 
  - `RealData/code/Stability_analysis_Table_S8-S13.R`
- Comparative performence in identifying marker genes: 
  - `RealData/code/marker_gene_identify_Table_S15-S16.R`

**(5) Gene Clustering and Spot Clustering**

Finally, reproduce clustering analyses based on identified SV genes:

- Gene cluster inference and visualization: 
  - `RealData/code/gene_cluster_draw_Figure_S14-S19_Table_S17-S19.R`
- Spot clustering and spatial domain identification: 
  - `RealData/code/comparison_spotcluster_index_Table_S20-S34.R`
  - `RealData/code/spatial_domain_draw_Figure3.R`

**Key Note**: All R scripts in `code/` use the `here` package to
reference files in `result_data/`. Ensure the working directory is set
to the `code_and_data` root before execution.





