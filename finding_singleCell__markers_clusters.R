library(Seurat)
library(tidyverse)
install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install('multtest')
install.packages('metap')
BiocManager::install('limma')


# load data
ifnb_harmony <-readRDS("F:/Biomagician_tutorials_workshop/work_space_rnaseq/rna_seq_projects/ifnb_harmony.rds")

# to see the data
str(ifnb_harmony)
View(ifnb_harmony@meta.data)

# visualize data
clusters <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'stim')
condition|clusters

# find the markers -----------------

# "find conserve markers" are more applicable to a scenario like this data type where
# we have data from two conditions and cells present in both the conditions 

# " find All markers" is more appropriate to use for data where we have data from just one conditional one group

# findAll markers -----------------

FindAllMarkers(ifnb_harmony,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')

# the first parameter would be the name of the object 
# the next parameter is log fold change threshold so this value is the minimal log two fold change average expression
# of  the gene in the cluster 
## you're trying to compare relative to the average expression of the gene and all
## the other clusters combined so the default value is 0.25 

# the next uh parameter is minimum percentage- this value indicates it will only test those genes that are
## detected at a minimum fraction of this number of cells in either of the two populations
### so let's say if we set it as 0.5 so it will only detect those genes and test those genes 
# rather that are detected at 50 frequency in either of the two clusters or the two populations that we are trying to compare
### so 50 is a little stringent the default is 0.1

# so we'll set it we'll go with the default or we can change it according to our prior knowledge
# by setting it up more stringent it will speed up the process of finding the markers but 
## at the same time you lose a lot of markers that are not detected in that many  percentage of the or the
## fraction of cells so let's go by default 

# the next parameter is only positive and true so this will return only positive markers 
# that is markers that are up-regulated  

## if we wish to have both the up regulator and the double regulated markers then we can set it to false

# the next parameter is test use and this function provides you with a lot
# of options for the test to use so you can choose from wilcox to poisson to negative binomial to t test 

# you can also use external packages like 'DESeq2' 
# will provide the slot as counts because it uses the raw counts
# it automatically switches the slot to the counts otherwise the default slot is the data slot which stores the normalized counts data and
# if you do not even specify the test use it goes with the default test

## findConserved markers -------------

# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
# DefaultAssay(ifnb_harmony) <- 'RNA'

DefaultAssay(ifnb_harmony)

markers_cluster3 <- FindConservedMarkers(ifnb_harmony,
                                         ident.1 = 3,
                                         grouping.var = 'stim')

head(markers_cluster3)

# running the "find conserved markers" function i.e- if to identify what cluster 3 - this is going to be my ident 1 
# because this is the first cluster that if you want to identify and to compare the cluster 3 with all the other clusters
# because to essentially get the markers that are up regulated or expressed in cluster 3 
# which will further help me to identify what are the type of the cells that form this cluster

# the first parameter would be the name of the data 
# the second parameter is ident one
# the third parameter is the grouping variable ; to group the cells by stim which has the information--
# on the control or stimulated condition 

## let's save this to a variable that is called markers for cluster three 
## To take a look at our markers from for cluster three

# find conserved markers function internally separates out cells by condition 
# for each condition it compared the cluster 3 to all the other clusters 

# let's visualize top features
FeaturePlot(ifnb_harmony, features = c('FCGR3A'), min.cutoff = 'q10')

# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))

# rename cluster 3 ident
Idents(ifnb_harmony)
ifnb_harmony <- RenameIdents(ifnb_harmony, `3` = 'CD16 Mono')

Idents(ifnb_harmony)

#visualize the data now
DimPlot(ifnb_harmony, reduction = 'umap', label = T)

# cells already have annotations provided in the metadata
View(ifnb_harmony@meta.data)

# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual anotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# setting Idents as Seurat annotations provided (also a sanity check!)
Idents(ifnb_harmony) <- ifnb_harmony@meta.data$seurat_annotations
Idents(ifnb_harmony)

#visualize the data now
DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# findMarkers between conditions ---------------------

ifnb_harmony$celltype.cnd <- paste0(ifnb_harmony$seurat_annotations,'_', ifnb_harmony$stim)

View(ifnb_harmony@meta.data)
Idents(ifnb_harmony) <- ifnb_harmony$celltype.cnd

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(ifnb_harmony, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_cluster3)


FeaturePlot(ifnb_harmony, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')