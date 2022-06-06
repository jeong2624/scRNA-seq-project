"
Created on Tue May 31 2022

@author: Jeong-Woon, Park

"

### Setup the Seurat Object
Load_packages =  c("dplyr","Seurat","patchwork","readr")
lapply(Load_packages, library, character.only = TRUE)

# Load the human cortex dataset
control = readRDS("control_c1.rds")
control

# The [[ operator can add columns to object metadata. 
# This is a great place to stash QC stats
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object,
# i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# we visualize QC metrics, and use these to filter cells.
# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have > 5% mitochondrial counts
control <- subset(control, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)

### Normalizing the data
control <- NormalizeData(control, normalization.method = "LogNormalize", scale.factor = 10000)

### Identification of highly variable features (feature selection)
control <- FindVariableFeatures(control, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(control), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(control)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(list(plot1, plot2))

### Scaling the data
all.genes <- rownames(control)
control <- ScaleData(control, features = all.genes)

### Perform linear dimensional reduction
control <- RunPCA(control, features = VariableFeatures(object = control))

# visualizing both cells and features that define the PCA
# Examine and visualize PCA results a few different ways
print(control[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(control, dims = 1:2, reduction = "pca")
DimPlot(control, reduction = "pca")

# when trying to decide which PCs to include for further downstream analyses
DimHeatmap(control, dims = 1:15, balanced = TRUE)

### Determine the ‘dimensionality’ of the dataset
control <- JackStraw(control, num.replicate = 100)
control <- ScoreJackStraw(control, dims = 1:20)

JackStrawPlot(control, dims = 1:15)
ElbowPlot(control)

### Cluster the cells
control <- FindNeighbors(control, dims = 1:15)
control <- FindClusters(control, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(control), 5)

### Run non-linear dimensional reduction (UMAP)
control <- RunUMAP(control, dims = 1:15)
DimPlot(control, reduction = "umap")

### Finding differentially expressed features (cluster biomarkers)
if (FALSE){
VlnPlot(control, features = c("AQP4", "GJA1", "SLC39A12", "SLC4A4"), 
        slot = "counts", log = TRUE) ## cluster 2 --> Astrocyte

VlnPlot(control, features = c("ANO3", "NEUROD6"), 
        slot = "counts", log = TRUE) ## cluster 1 --> Excitatory neurons

VlnPlot(control, features = c("GAD1", "VIP", "GAD2"), 
        slot = "counts", log = TRUE) ## cluster 0--> Inhibitory neurons

VlnPlot(control, features = c("OLIG2", "OLIG1", "SOX10", "VCAN"),
        slot = "counts", log = TRUE) ## cluster 4 --> OPC

VlnPlot(control, features = c("MOG", "MBP", "BCAS1", "CNP"),
        slot = "counts", log = TRUE) ## cluster 3 --> Oligodendrocyte

VlnPlot(control, features = c("CCL3", "GPR183", "CX3CR1", "CSF1R"),
        slot = "counts", log = TRUE) ## cluster 5 --> Microglia

VlnPlot(control, features = c("PALMD", "APOLD1"),
        slot = "counts", log = TRUE) ## cluster 6 --> Endothelial
}

# find markers for every cluster compared to all remaining cells, 
# report only the positive ones
control.markers <- FindAllMarkers(control, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)
control.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

### Assigning cell type identity to clusters
new.cluster.ids <- c("Inhibitory neurons", "Excitatory neurons", "Astrocyte", 
                     "Oligodendrocytes", "OPC", "Microglia", "Endothelial cell")
names(new.cluster.ids) <- levels(control)
control <- RenameIdents(control, new.cluster.ids)
DimPlot(control, reduction = "umap", label = TRUE, pt.size = 0.5)

# generates an expression heatmap for given cells and features.
control.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(control, features = top10$gene)

# visualizing marker expression. 
VlnPlot(control, features = c("ARID1B", "ASH1L", "CHD2", 
                          "MECP2", "CACNA1C", "SHANK2"), 
        slot = "counts", log = TRUE)
