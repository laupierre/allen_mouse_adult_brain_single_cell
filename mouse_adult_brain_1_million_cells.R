# see https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x
# See jeremy's way to subset the big data, https://community.brain-map.org/t/gene-expression-matrix-cvs-is-too-large-to-load-it/566/17


library(rhdf5)       # For reading hdf5; https://www.bioconductor.org/packages/devel/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
#library(HDF5Array)  # Alternatively option for reading the data matrix with less memory: https://rdrr.io/github/Bioconductor/HDF5Array/
library(data.table)  # For fast reading of csv files
library(Seurat)
library (dplyr)

genes   <- h5read("expression_matrix.hdf5","/data/gene")
samples <- h5read("expression_matrix.hdf5","/data/samples")


## Read in the metadata using fread (fast!)
metadata <- fread("metadata.csv")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_name
# Note that the order of metadata and counts and the number of cells are DIFFERENT !

metadata <- metadata [ ,c("sample_name", "cluster_label", "subclass_label", "neighborhood_label")]

# Remove Meis2 cells (less than 25 cells)
metadata <- metadata [metadata$subclass_label != "Meis2", ]
sort (table (metadata$subclass_label))


cell_types <- unique (metadata$subclass_label)

resl <- list ()

for (i in 1:length (cell_types)) {
use_samples  <- metadata$sample_name[metadata$subclass_label == cell_types[i] ]
print (cell_types[i])
set.seed(42)
num_samples <- length (use_samples)
num_s <- ifelse (num_samples < 200, num_samples, 200)
idx <- match (intersect(sample (use_samples, num_s), samples), samples)
read_samples <- sort (idx)


## Read in only a relevant subset of data using h5read. This is the method that probably works best in most situations. Use HD5Array if there are memory issues  
counts1 <- h5read("expression_matrix.hdf5", "/data/counts", index = list(read_samples, NULL))
# Note that COLUMNS are genes and need a transposition
counts1 <- t(counts1)
subcounts1 <- as(counts1, "dgCMatrix")
rownames(subcounts1) <- as.character(genes)
colnames(subcounts1) <- as.character(samples) [read_samples]
resl[[i]] <- subcounts1
}


res <- do.call ("cbind", resl)

seurat <- CreateSeuratObject(counts=res)         
met <- as.data.frame(metadata[colnames(res),]) # Add the metadata
seurat <- AddMetaData(seurat, met)  

saveRDS(object = seurat, file = "seurat_allen_small_mouse_brain.RDS")


# run SCT normalization and dimensionality reduction
allen_reference <- SCTransform(seurat, ncells = 3000, verbose = FALSE, method = "poisson") %>%
                   RunPCA(verbose = FALSE) %>%
                   RunUMAP(dims = 1:30)

# the annotation is stored in the 'subclass_label' column of object metadata
Idents (allen_reference) <- "subclass_label"
DimPlot(allen_reference, label = TRUE)

saveRDS(object = allen_reference, file = "seurat_allen_small_mouse_brain_with_UMAP.RDS")



