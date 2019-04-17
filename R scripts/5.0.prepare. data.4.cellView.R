
##  module load r/3.5.0-py2-qqwf6c6
##  source  ~/bin/system.py3.6.5_env/bin/activate


### CellView

## The primary structure of the .Rds file comprises three dataframes, with the following object names and column names

## log2cpm - Your (N x M) genes vs. cells expression matrix. Gene names need to be in Ensembl gene id format (e.g. ENSG, ENSM)

## tsne.data - Your dimensionality reduction and sample clustering information. This dataframe contains M rows and 4 columns: V1, V2, V3, dbCluster
   ##  V1, V2, V3 store the 3 dimensional representation of your data, e.g. from t-SNE, PCA, etc.
   ##  dbCluster contains numerical cluster assignments. The row names of this data frame must correspond to the column names of your expression matrix.


## featuredata - A dataframe representing gene annotations with row names in ensembl gene id format. The following 2 columns are required:
    ##???Chromosome.Name - integers representing chromosome numbers.
??????## Associated.Gene.Name - Gene symbol
???????????????## Other columns are allowed but not utilized. The number of rows can be larger than the number of genes (N) in your expression matrix, but beware of duplicate gene names with unique ENSGIDs. 


# cowplot enables side-by-side ggplots
library("cowplot")
library("Seurat")
library("scran")
library("rsvd")
library(ggplot2)
library(gridExtra)


setwd("/home/haibol/haibol_work/scRNA-seq")

out_objects_dir <- "./results/R.out/data/Robjects"
out_plot_dir <- "./results/R.out/ALRA"
if(!dir.exists(out_plot_dir)) dir.create(out_plot_dir)


library_id <- c("A", "B", "C", "CT2-1NOV", "CT2-30OCT") 
set.seed(1234)

load(file=file.path(out_objects_dir, "ALRA.imputed.RData"))
pbmc.integrated <- RunTSNE(object = pbmc.integrated, reduction = "pca", k.seed = 2,
                           dims = 1:22,  dim.embed = 3)



###############################################################################

options(stringsAsFactors = FALSE, row.names = 1, as.is = T)
log2cpm <- as.matrix(pbmc.integrated@assays$integrated@data)


## rename the genes
pig_genes <- read.delim("All.pig.gene.plus.HGNC.names.txt", header = TRUE, as.is = TRUE)

log2cpm <- merge(log2cpm, pig_genes, by.x ="row.names", by.y = "Gene.name", all.x = TRUE)
rownames(log2cpm) <- log2cpm[, ncol(log2cpm)]

featuredata <- log2cpm[, 1, drop =FALSE]
log2cpm <- log2cpm[, -1]
log2cpm <- log2cpm[, -ncol(log2cpm)]


## chromosome information
pig_chr_genes <- read.delim("pig_genes_v95.txt",header = TRUE, sep = "\t")

featuredata <- merge(featuredata, pig_chr_genes[, c(1,3)], 
                     by.x= "row.names", by.y = "Gene.stable.ID", 
                     sort = TRUE, all.x= TRUE)

rownames(featuredata) <- featuredata[,1]
featuredata <- featuredata[, c(3, 2)]
colnames(featuredata) <- c("Chromosome.Name", "Associated.Gene.Name")

## reorder
featuredata<- featuredata[rownames(log2cpm), ]


tsne.data <- pbmc.integrated@reductions$tsne@cell.embeddings
tsne.data <- as.data.frame(cbind(tsne.data, 
                                 as.numeric(as.character(pbmc.integrated@active.ident))))

colnames(tsne.data) <- c(paste0("V", 1:3), "dbCluster")

## sanity check
all(rownames(tsne.data) == colnames(log2cpm))
all(rownames(featuredata) == rownames(log2cpm))

## output data for cellview
save(log2cpm, featuredata, tsne.data, file = file.path(out_objects_dir, 'ALRA.imputed.PMBC.integrated.scRNA-seq.Rds'))


##
packages <- c("shiny", "plotly", "shinythemes", "DT",
              "threejs", "sm", "RColorBrewer", "mclust", "reshape")
lapply(packages, function(.x){
    BiocManager::install(.x)
})

library(shiny)
library(plotly)
library(shinythemes)
library(ggplot2)
library(DT)
library(pheatmap)
library(threejs)
library(sm)
library(RColorBrewer)
library(mclust)
library(reshape)


shiny::runGitHub("CellView", 'mohanbolisetty')



