## Project name : Autism Spectrum Disorder single cell transcriptome profiles

#### Project time : 2022.03.08 - 2022.06.14

#### Description :
* The project was conducted as part of "single-cell genomics" class in the Department of Bioinformatics at Soongsil University, Korea.
* The purpose of this project is to analyze Autism Spectrum Disorder (ASD) in terms of single-cell RNA sequencing.
* According to this project, I want to check gene expression of ASD (Autism Spectrum Disorder) for each cell type to check E/I imbalance hypothesis.

#### Software & Package version
* R >= 4.1.1
* cellranger 6.1.2
* Seurat 4.0

#### Our Dataset :
* GSE67835 (Healthy human brain samples, raw count matrix)
* SRR9262932 (ASD prefrontal cortex tissue, Fastq files)

#### We uploaded these files:
* control_c1.rds : The healthy brain cortex Seurat object data. (Fluidigm C1 platform)
* ASD_10x.rds : The ASD prefrontal cortex tissue Seurat object data. (10x Genomics platform)
* C1_control.R : scRNA-seq analysis code for healthy human cortex.
* 10x_ASD.R : scRNA-seq analysis code for ASD prefrontal cortex tissue.
* Integration.R : scRNA-seq analysis code for integrated data. (healthy + ASD)
* filtered_feature_bc_matrix : ASD scRNA-seq analysis raw count matrix. (using cellranger pipeline)
* control_counts.csv : Healthy brain cortex raw count matrix.
