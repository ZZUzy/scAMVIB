# Adaptive Multi-view Information Bottleneck for Multi-omics Data Clustering
## Introduction
This repository introduces scAMVIB, an adaptive multi-view clustering method based on the variational information bottleneck framework. Our approach uses similarity network fusion to construct an omics-fused sample matrix that captures intrinsic inter-omics relationships, while integrating it with distinct omics profiles to build a unified multi-view feature space—effectively redefining the analytical framework for single-cell multi-omics data.
## Framework
[Framework.pdf](https://github.com/user-attachments/files/23602916/Framework.pdf)

## System Requirements
The required python dependencies are given below. 

h5py
numpy==1.17.2
numba==0.51.1
pandas==0.25.1
scikit-learn==0.23
tqdm==4.33.0
colorlog==4.0.2
colored==1.3.93
seaborn==0.9.1
python-highcharts==0.4.2
umap-learn==0.3.10
snfpy==0.2.2
vit-pytorch

## Datasets
The seven real multi-omics datasets used in this study can be categorized into two types based on their omics composition: the first type combines scRNA-seq and scATAC-seq data (containing mRNA and ATAC information), including the human cell line mixture CellLine dataset (GSE126074) and the human peripheral blood mononuclear cell PBMC3k dataset; the second type integrates scRNA-seq with scADT-seq data (containing mRNA and ADT information), comprising the umbilical cord blood mononuclear cell CBMN dataset (GSE100866) and two batches of mouse spleen lymph node data (SLN208 and SLN111 from GSE150599). All datasets were obtained from the [GEO database](https://www.ncbi.nlm.nih.gov/geo/) except for the human PBMC3k dataset, which was downloaded from the 10x [Genomics website](https://www.10xgenomics.com/).

## Run scAMVIB
To execute the provided MATLAB script (`run_script.m`), follow these steps:

1. **Preparation**: Ensure your file structure includes the script along with three required subdirectories (`MainFunction`, `BasicFunction`, and `clustering performance`) in the same parent directory. Verify that all input data files exist at their specified paths:
   - RNA data: `D:/Dataset/Tesedemo/from scMDC/CITEseq_sln111/Batch01/rna_cite_sln111_Batch1.tsv`
   - ADT data: `D:/Dataset/Tesedemo/from scMDC/CITEseq_sln111/Batch01/adt_cite_sln111_Batch1.tsv`
   - Label data: `D:/Dataset/Tesedemo/from scMDC/CITEseq_sln111/Batch01/label_cite_sln111_Batch1.tsv`
   - Information-enhanced data: `D:/Experiment/Code/Python Code/scRNAseq code package/MoINER-master/results_preprocessing/SLN111D1/4.Multi-Omics_Data(Information Enhanced).csv`

2. **Execution**: Open MATLAB and navigate to the script's directory in the command window. Run the script using either:
   - Direct call: `run_script`
   - Or using the run function: `run('run_script.m')`

3. **Process**: The script will automatically:
   - Load and preprocess multi-view data (RNA, ADT, and information-enhanced matrices)
   - Perform random subset selection of cells
   - Execute auto-weighted multi-view information bottleneck clustering
   - Calculate and display clustering performance metrics (ARI and NMI)
   - Save results to `Result/multi_omics/SLN111D1/` directory

4. **Output**: Monitor the command window for real-time progress updates and final clustering performance statistics. Results will be automatically saved while temporary variables are cleared upon completion.
