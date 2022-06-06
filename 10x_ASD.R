"
Created on  May 31 2022

@author: Jeong-Woon, Park

"

### Setup the Seurat Object
Load_packages =  c("dplyr","Seurat","patchwork","readr")
lapply(Load_packages, library, character.only = TRUE)

# Load the human cortex dataset
ASD = readRDS("ASD_10x.rds")
ASD

# The [[ operator can add columns to object metadata. 
# This is a great place to stash QC stats
ASD[["percent.mt"]] <- PercentageFeatureSet(ASD, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(ASD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object,
# i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(ASD, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ASD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(list(plot1, plot2))

# we visualize QC metrics, and use these to filter cells.
# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have > 5% mitochondrial counts
ASD <- subset(ASD, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

### Normalizing the data
ASD <- NormalizeData(ASD, normalization.method = "LogNormalize", scale.factor = 10000)

### Identification of highly variable features (feature selection)
ASD <- FindVariableFeatures(ASD, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ASD), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ASD)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(list(plot1, plot2))

### Scaling the data
all.genes <- rownames(ASD)
ASD <- ScaleData(ASD, features = all.genes)

### Perform linear dimensional reduction
ASD <- RunPCA(ASD, features = VariableFeatures(object = ASD))

# visualizing both cells and features that define the PCA
# Examine and visualize PCA results a few different ways
print(ASD[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ASD, dims = 1:2, reduction = "pca")
DimPlot(ASD, reduction = "pca")

# when trying to decide which PCs to include for further downstream analyses
DimHeatmap(ASD, dims = 1:15, balanced = TRUE)

### Determine the ‘dimensionality’ of the dataset
ASD <- JackStraw(ASD, num.replicate = 100)
ASD <- ScoreJackStraw(ASD, dims = 1:20)

JackStrawPlot(ASD, dims = 1:15)
ElbowPlot(ASD)

### Cluster the cells
ASD <- FindNeighbors(ASD, dims = 1:15)
ASD <- FindClusters(ASD, resolution = 1.0)

# Look at cluster IDs of the first 5 cells
head(Idents(ASD), 5)

### Run non-linear dimensional reduction (UMAP)
ASD <- RunUMAP(ASD, dims = 1:15)
DimPlot(ASD, reduction = "umap")

### Finding differentially expressed features (cluster biomarkers)
if (FALSE){
VlnPlot(ASD, features = c("ANO3", "GAD1", "SLC17A7", "STMN2")) # cluster 0, 2, 3, 7, 8, 9 --> Excitatory neurons
VlnPlot(ASD, features = c("VIP", "TAC1")) # cluster 4, 5, 6 --> inhibitory neurons
VlnPlot(ASD, features = c("GFAP")) # cluster 1 --> Astrocytes
VlnPlot(ASD, features = c("CD34", "VWF")) # cluster 13 --> Endothelial
VlnPlot(ASD, features = c("PECAM1")) # cluster 13 --> microglia
VlnPlot(ASD, features = c("BCAS1", "CNP", "PLP1")) # cluster 11 --> Oligodendrocyte
VlnPlot(ASD, features = c("OLIG1", "OLIG2", "SOX10")) # cluster 10 --> OPC
}

# find markers for every cluster compared to all remaining cells, 
# report only the positive ones
ASD.markers <- FindAllMarkers(ASD, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)
ASD.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

### Assigning cell type identity to clusters
if (FALSE){
new.cluster.ids <- c("Excitatory neurons 1", "Astrocytes" ,"Excitatory neurons 2", 
                     "Excitatory neurons 3", "Inhibitory neurons 1", "Inhibitory neurons 2",
                     "Inhibitory neurons 3", "Excitatory neurons 4", "Excitatory neurons 5",
                     "Excitatory neurons 6", "OPC", "Oligodendrocyte", "Endothelial cell",
                     "Microglia")
}

new.cluster.ids <- c("Excitatory neurons", "Astrocytes" ,"Excitatory neurons", 
                     "Excitatory neurons", "Inhibitory neurons", "Inhibitory neurons",
                     "Inhibitory neurons", "Excitatory neurons", "Excitatory neurons",
                     "Excitatory neurons", "OPC", "Oligodendrocyte", "Endothelial cell",
                     "Microglia")

names(new.cluster.ids) <- levels(ASD)
ASD <- RenameIdents(ASD, new.cluster.ids)
DimPlot(ASD, reduction = "umap", label = TRUE, pt.size = 0.5)

# generates an expression heatmap for given cells and features.
ASD.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(ASD, features = top5$gene)

# visualizing marker expression. 
VlnPlot(ASD, features = c("ARID1B", "ASH1L", "CHD2", 
                          "MECP2", "CACNA1C", "SHANK2"), 
        slot = "counts", log = TRUE)


