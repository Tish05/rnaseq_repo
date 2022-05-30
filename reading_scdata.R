if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

## script to demonstrate reading single cell matrices in various format 
# and converting to seurat object 

#loading libraries

library(Seurat)
library(Seuratdisk)

# .RDS format

rds_obj <- readRDS('ependymal_cells.rds')

# 10X CellRanger .HDF5 format 

hdf5_obj <- Read10X_h5(filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)

#converting seurat object:

seurat_hdf5 <- CreateSeuratObject(counts = hdf5_obj)


# to look at your created seurat object:
# type in the console : str (name of the object-in this case - seurat_hdf5) click enter


# .mtx file
mtx_obj <- ReadMtx(mtx = "raw_feature_bc_matrix/matrix.mtx.gz",
                   features = "raw_feature_bc_matrix/features.tsv.gz",
                   cells = "raw_feature_bc_matrix/barcodes.tsv.gz")

seurat_mtx <- CreateSeuratObject(counts = mtx_obj)

# to look at your created seurat object: 1st 10 rows and 1st 10 columns
# type in the console : mtx_obj[1:10,1:10]

# .loom files

loom_oj <- Connect(filename = "adult-hem-organs-10X-bone-marrow.loom", mode = 'r')

seurat_loom <- as.Seurat(loom_oj)

# .h5ad format 

# step 1: convert AnnData object to an h5Seurat file

Convert("adata_SS2_for_download.h5ad", dest = "h5seurat", overwrite = TRUE)

# step 2: Load h5Seurat file into a Seurat object 

seurat_anndata <- LoadH5Seurat("adata_SS2_for_download.h5seurat")