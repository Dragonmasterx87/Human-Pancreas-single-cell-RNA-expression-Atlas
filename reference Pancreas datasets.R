# The following code allows for the analysis of 5 single cell RNAseq datasets of the human pancreas
# Information on these datasets can be found in the following locations:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81076
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86469
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5061/
# This code was written by Fahd Qadir PhD. on 06/03/2020 email: mqadir@tulane.edu

# 1. installation and loading of packages
# Devtools
install.packages('devtools')
library(devtools)

# Seuratdata
devtools::install_github('satijalab/seurat-data')

# Seurat wrappers
devtools::install_github('satijalab/seurat-wrappers')

# Load packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(SeuratWrappers)

options(future.globals.maxSize = 4000 * 1024^2)

# Loading of datasets
GSE81076 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/GSE81076.csv", header = TRUE, sep = ",", row.names = 1)
GSE85241 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/GSE85241.csv", header = TRUE, sep = ",", row.names = 1)
GSE86469 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/GSE86469.csv", header = TRUE, sep = ",", row.names = 1)
GSE84133 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/GSE84133.csv", header = TRUE, sep = ",", row.names = 1)
EMTAB5061 <- read.csv("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/pancreas/EMTAB5061.csv", header = TRUE, sep = ",", row.names = 1)

# Create Seurat objects
GSE81076 <- CreateSeuratObject(counts = GSE81076, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
GSE85241 <- CreateSeuratObject(counts = GSE85241, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
GSE86469 <- CreateSeuratObject(counts = GSE86469, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
GSE84133 <- CreateSeuratObject(counts = GSE84133, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)
EMTAB5061 <- CreateSeuratObject(counts = EMTAB5061, project = "SeuratProject", assay = "RNA", min.cells = 3, min.features = 200)

# Sample specific Metadata addition
GSE81076$sample <- "GSE81076"
GSE85241$sample <- "GSE85241"
GSE86469$sample <- "GSE86469"
GSE84133$sample <- "GSE84133"
EMTAB5061$sample <- "EMTAB5061"

# Create a list of datasets containing seurat objects
pancreas.list <- list("GSE81076" = GSE81076, "GSE85241" =GSE85241, "GSE86469" = GSE86469, "GSE84133" = GSE84133, "EMTAB5061" = EMTAB5061)

# Run SCtransform (computationally intensive step, takes time to perform)
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = TRUE)
}

# Select features downstream
pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = TRUE)

# Identify anchors
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = TRUE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = TRUE)

# Normalize based on RNA
pancreas.integrated <- NormalizeData(pancreas.integrated, normalization.method = "LogNormalize", assay = "RNA", scale.factor = 1e6, 
                                     verbose = TRUE)

#Clustering
pancreas.integrated <- RunPCA(pancreas.integrated, verbose = TRUE)
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 1.2)

# Visualization Clustering
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)
plots <- DimPlot(pancreas.integrated, group.by = c("sample", "CellType"))
plots & theme(legend.position = "right") & guides(color = guide_legend(nrow = 14, byrow = TRUE,
                                                                     override.aes = list(size = 5)))
DimPlot(pancreas.integrated)

# Organize clusters
plot <- DimPlot(pancreas.integrated, reduction = "umap")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "beta")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "alpha")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "delta")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "epsilon")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "gamma")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "ductal")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "acinar")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "macrophage")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "Tcell")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "endothelial")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "quiescent stellate")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "activated stellate")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "schwann")
pancreas.integrated <- CellSelector(plot = plot, object = pancreas.integrated, ident = "mast")
levels(pancreas.integrated)

# Saving this information in the metadata slot
head(Idents(pancreas.integrated))
pancreas.integrated$CellType <- Idents(pancreas.integrated)
head(pancreas.integrated@meta.data)

# Run find variable features again running this is questionable, as only the var features from integrated data is useful
# But Seurat recommends re-running this
DefaultAssay(object = pancreas.integrated) <- "RNA"
pancreas.integrated <- FindVariableFeatures(pancreas.integrated, selection.method = "vst", nfeatures = 3000)

# Define an order of cluster identities remember after this step-
# cluster re-assignment occurs, which re-assigns clustering in my_levels
my_levels <- c("beta", "alpha", "delta", "gamma", "epsilon", "ductal", "acinar", "quiescent stellate", "activated stellate", "schwann", "endothelial", "macrophage", "Tcell", "mast")
head(pancreas.integrated@meta.data$CellType)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
pancreas.integrated@meta.data$CellType <- factor(x = pancreas.integrated@meta.data$CellType, levels = my_levels)

save(pancreas.integrated, file="C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/Workspace/pancreas.integrated.Robj")
load("C:/Users/mqadir/Box/Lab 2301/sCell Analysis Project/Refrence Human Pancreas/scRNAseq datasets/Workspace/pancreas.integrated.Robj")

# Check metadata
head(pancreas.integrated@meta.data)
table(pancreas.integrated$sample)

# Check activeidents
head(Idents(pancreas.integrated))

# Change active idents to CellType
Idents(pancreas.integrated) <- "CellType"

# For UMAP visualization
DefaultAssay(object = pancreas.integrated) <- "RNA"
FeaturePlot(object = pancreas.integrated, 
            features = c("ACE2"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 100,
            order = TRUE)

# Visualize information
table(pancreas.integrated$sample)
DefaultAssay(object = pancreas.integrated) <- "RNA"
VlnPlot(pancreas.integrated, c("ACE2", "TMPRSS2"), group.by = "CellType", assay = "RNA", slot = "data")
VlnPlot(pancreas.integrated, c("GABRA1"), group.by = "CellType", assay = "RNA", slot = "data")

# Set cell identity to sample identity so that you can extraxt cell type information for plotting
Idents(object = pancreas.integrated) <- pancreas.integrated@meta.data$celltype

# How can I extract expression matrix for all beta cells
betacells <- subset(pancreas.integrated, idents = c("beta"))

# Violin plot
DefaultAssay(object = betacells) <- "RNA"
VlnPlot(object = betacells, features = c("ACE2", "TMPRSS2"), group.by = "sample", slot = "data")

# How can I extract expression matrix for all beta cells
alphacells <- subset(pancreas.integrated, idents = c("alpha"))

# Violin plot
DefaultAssay(object = alphacells) <- "RNA"
VlnPlot(object = alphacells, features = c("ACE2", "TMPRSS2"), group.by = "sample", slot = "data")

# Set cell identity to sample identity
Idents(object = pancreas.integrated) <- pancreas.integrated@meta.data$celltype

# Find if SRD genes are differentially expressed
beta.integrated.markers <- FindAllMarkers(object = pancreas.integrated, slot = 'data', test.use = 'wilcox')

# How can I calculate the average expression of all cells within a cluster?
cluster.averages <- AverageExpression(pancreas.integrated, assay= "RNA", slot = "data")
head(cluster.averages[["RNA"]][c("ACE2", "TMPRSS2"), 1:14])
