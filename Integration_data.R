"
Created on Thr Jun 2 2022

@author: Jeong-Woon, Park

"

### Setup the Seurat Object
Load_packages =  c("dplyr","Seurat","patchwork","readr", "stringr" ,"sjmisc")
lapply(Load_packages, library, character.only = TRUE)

### Loda datasets
control = readRDS("control_c1.rds")
ASD = readRDS("ASD_10x.rds")

# The [[ operator can add columns to object metadata. 
# This is a great place to stash QC stats
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^MT-")
ASD[["percent.mt"]] <- PercentageFeatureSet(ASD, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(ASD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object,
# i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot3 <- FeatureScatter(ASD, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(ASD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4

# we visualize QC metrics, and use these to filter cells.
# We filter cells that have unique feature counts over 2,500 or less than 200
# We filter cells that have > 5% mitochondrial counts
control <- subset(control, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)
ASD <- subset(ASD, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

### Merge datasets
data = merge(control, ASD, merge.data = TRUE)
data

# split the dataset into a list of two seurat objects (CTRL and ASD)
data.list <- SplitObject(data, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list)

### Perform integration
data.anchors <- FindIntegrationAnchors(object.list = data.list, 
                                       anchor.features = features)

# this command creates an 'integrated' data assay
data.combined <- IntegrateData(anchorset = data.anchors)

### Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)

### Perform linear dimensional reduction
data.combined <- RunPCA(data.combined, features = VariableFeatures(object = data.combined))

# Examine and visualize PCA results a few different ways
print(data.combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data.combined, dims = 1:2, reduction = "pca")
DimPlot(data.combined, reduction = "pca")
DimHeatmap(data.combined, dims = 1:15, balanced = TRUE)

### Determine the ‘dimensionality’ of the dataset
data.combined <- JackStraw(data.combined, num.replicate = 100)
data.combined <- ScoreJackStraw(data.combined, dims = 1:20)

JackStrawPlot(data.combined, dims = 1:15)
ElbowPlot(data.combined)

### Cluster the cells
data.combined <- FindNeighbors(data.combined, dims = 1:15)
data.combined <- FindClusters(data.combined, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(data.combined), 5)

### Run non-linear dimensional reduction (UMAP/tSNE)
data.combined <- RunUMAP(data.combined, dims = 1:15)
DimPlot(data.combined, reduction = "umap")

# To visualize the two conditions side-by-side, 
# we can use the split.by argument to show each condition colored by cluster.
DimPlot(data.combined, reduction = "umap", split.by = "orig.ident")

### Finding differentially expressed features (cluster biomarkers)
DefaultAssay(data.combined) <- "RNA"

if (FALSE){
VlnPlot(data.combined, features = c("VIP", "GAD2"), 
        slot = "counts", log = TRUE) ## cluster 2, 6, 7, 9, 11 --> Inhibitory Neurons

VlnPlot(data.combined, features = c("AQP4", "GJA1", "SLC39A12", "SLC4A4"), 
        slot = "counts", log = TRUE) ## cluster 3 --> Astrocyte

VlnPlot(data.combined, features = c("ANO3", "NEUROD6", "SLC17A7", "STMN2"), 
        slot = "counts", log = TRUE) ## cluster 0, 1, 4, 5 --> Excitatory neurons

VlnPlot(data.combined, features = c("OLIG2", "OLIG1", "SOX10", "VCAN"),
        slot = "counts", log = TRUE) ## cluster 8 --> OPC

VlnPlot(data.combined, features = c("MOG", "MBP", "BCAS1", "CNP"),
        slot = "counts", log = TRUE) ## cluster 11 --> Oligodendrocyte

VlnPlot(data.combined, features = c("CCL3", "GPR183", "CX3CR1", "CSF1R"),
        slot = "counts", log = TRUE) ## cluster 12 --> Microglia

VlnPlot(data.combined, features = c("PALMD", "APOLD1"),
        slot = "counts", log = TRUE) ## cluster 10 --> Endothelial
}

# find markers for every cluster compared to all remaining cells, 
# report only the positive ones
markers <- FindAllMarkers(data.combined, only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

### Assigning cell type identity to clusters
if (FALSE){
new.cluster.ids <- c("Excitatory neurons", "Excitatory neurons", "Inhibitory neurons", 
                     "Astrocyte", "Excitatory neurons", "Excitatory neurons", 
                     "Inhibitory neurons", "Inhibitory neurons", "OPC",
                     "Inhibitory neurons", "Endothelial cell", "Oligodendrocyte",
                     "Microglia")

names(new.cluster.ids) <- levels(data.combined)
data.combined <- RenameIdents(data.combined, new.cluster.ids)
DimPlot(data.combined, reduction = "umap", label = TRUE, pt.size = 0.5)
}

new.cluster.ids <- c("Excitatory neurons 1", "Excitatory neurons 2", "Inhibitory neurons 1", 
                     "Astrocyte", "Excitatory neurons 3", "Excitatory neurons 4", 
                     "Inhibitory neurons 2", "Inhibitory neurons 3", "OPC",
                     "Inhibitory neurons 4", "Endothelial cell", "Oligodendrocyte",
                     "Microglia")

names(new.cluster.ids) <- levels(data.combined)
data.combined <- RenameIdents(data.combined, new.cluster.ids)
DimPlot(data.combined, reduction = "umap", label = TRUE, pt.size = 0.5)

### generates an expression heatmap for given cells and features.
DefaultAssay(data.combined) <- "RNA"
data.combined <- ScaleData(data.combined, verbose = FALSE)

markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(data.combined, features = top10$gene)

### visualizing marker expression. (Dot plot)
DefaultAssay(data.combined) <- "RNA"

markers.to.plot <- c("ARID1B", "ASH1L", "CHD2", "AFF2", 
                     "MECP2", "CACNA1C", "SHANK2", "ATRX")
DotPlot(data.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, 
        split.by = "orig.ident") + RotatedAxis()

### visualizing marker expression. (feature plot)
data.combined@meta.data$group = Idents(data.combined)

FeaturePlot(data.combined, features = markers.to.plot, split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

plots <- VlnPlot(data.combined, features = c("ARID1B", "ASH1L", "CHD2"), split.by = "orig.ident", group.by = "group",
                 pt.size = 0, combine = FALSE, split.plot = TRUE)
wrap_plots(plots = plots, ncol = 1)
