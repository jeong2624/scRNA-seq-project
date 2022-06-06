# Single_cell_Genomics_Class

This is my scRNA-seq analysis project for Single cell Genomics class in 1st semester, 2022.

My project is that analyzing Autism Spectrum Disorder (ASD) in terms of single-cell RNA sequencing.

### We used the analysis tools for my project.
* R >= 4.1.1
* cellranger 7.0.0
* python >= 3.7.13

### We used the analysis R package for my project.
* Seurat 4.0

### We uploaded these files:
* control_c1.rds : The healthy brain cortex Seurat object data. (Fluidigm C1 platform)
* ASD_10x.rds : The ASD prefrontal cortex tissue Seurat object data. (10x Genomics platform)
* C1_control.R : scRNA-seq analysis code for healthy human cortex.
* 10x_ASD.R : scRNA-seq analysis code for ASD prefrontal cortex tissue.
* Integration.R : scRNA-seq analysis code for integrated data. (healthy + ASD)
