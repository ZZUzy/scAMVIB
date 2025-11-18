# Adaptive Multi-view Information Bottleneck for Multi-omics Data Clustering
## Introduction
This repository introduces scAMVIB, an adaptive multi-view clustering method based on the variational information bottleneck framework. Our approach uses similarity network fusion to construct an omics-fused sample matrix that captures intrinsic inter-omics relationships, while integrating it with distinct omics profiles to build a unified multi-view feature space—effectively redefining the analytical framework for single-cell multi-omics data.
## Framework
<img width="2000" height="2074" alt="image" src="https://github.com/user-attachments/assets/4b56e999-1f50-4d69-92a4-d71fd91559db" />

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
