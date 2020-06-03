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

# Run SCtransform
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
pancreas.integratedx <- pancreas.integrated
pancreas.integrated <- pancreas.integratedx 

# Normalize based on RNA
pancreas.integrated <- NormalizeData(pancreas.integrated, normalization.method = "LogNormalize", assay = "RNA", verbose = TRUE)

# Normalize based on SCT:
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = TRUE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = TRUE)

# Visualization Clustering
pancreas.integrated <- RunPCA(pancreas.integrated, verbose = TRUE)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)
plots <- DimPlot(pancreas.integrated, group.by = c("tech", "celltype"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))

# Normalize RNA assay for plotting
DefaultAssay(object = pancreas.integrated) <- "RNA"
pancreas.integrated <- NormalizeData(pancreas.integrated, normalization.method = "LogNormalize", scale.factor = 1e6)

# Visualize information
table(pancreas.integrated$dataset)
DefaultAssay(object = pancreas.integrated) <- "RNA"
VlnPlot(pancreas.integrated, c("ACE2", "TMPRSS2"), group.by = "celltype", assay = "RNA", slot = "data")

# For UMAP visualization
DefaultAssay(object = pancreas.integrated) <- "SCT"
FeaturePlot(object = pancreas.integrated, 
            features = c("ACE2"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            #max.cutoff = 100,
            order = TRUE)

# Set cell identity to sample identity so that you can extraxt cell type information for plotting
Idents(object = pancreas.integrated) <- pancreas.integrated@meta.data$celltype

# How can I extract expression matrix for all beta cells
betacells <- subset(pancreas.integrated, idents = c("beta"))

# Violin plot
DefaultAssay(object = betacells) <- "RNA"
VlnPlot(object = betacells, features = c("ACE2", "TMPRSS2"), group.by = "tech")

# How can I extract expression matrix for all beta cells
alphacells <- subset(pancreas.integrated, idents = c("alpha"))

# Violin plot
DefaultAssay(object = alphacells) <- "RNA"
VlnPlot(object = alphacells, features = c("ACE2", "TMPRSS2"), group.by = "tech", slot = "data")

# Set cell identity to sample identity
Idents(object = pancreas.integrated) <- pancreas.integrated@meta.data$celltype

# Find if SRD genes are differentially expressed
beta.integrated.markers <- FindAllMarkers(object = pancreas.integrated, slot = 'data', test.use = 'wilcox')

# How can I calculate the average expression of all cells within a cluster?
cluster.averages <- AverageExpression(pancreas.integrated, assay= "RNA", slot = "data")
head(cluster.averages[["RNA"]][c("ACE2", "TMPRSS2"), 1:13])
