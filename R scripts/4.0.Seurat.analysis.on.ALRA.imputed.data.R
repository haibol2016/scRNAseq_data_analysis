##  module load r/3.5.0-py2-qqwf6c6
##  source  ~/bin/system.py3.6.5_env/bin/activate

##  install Seurat Realease 3.0
#install.packages('devtools')
#devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')


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

dataList <- readRDS(file=file.path(out_objects_dir, "ExpressionList_QC.rds"))

m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
fD$keep[is.na(fD$keep)] <- FALSE
rm(dataList)

# Gene and cell filtering

m <- m[fD$keep, pD$PassAll]   ## 11707 genes X 15847 cells
pD <- pD[pD$PassAll, ]
rownames(pD)  <- pD[, 1]
fD <- fD[fD$keep, ]

## pig gene names

pig_genes <- read.delim("All.pig.gene.plus.HGNC.names.txt", header = TRUE, as.is = TRUE)
pig_genes_HGNC <- pig_genes$Gene.name
names(pig_genes_HGNC) <- pig_genes$Gene.stable.ID
rownames(m) <- pig_genes_HGNC[rownames(m)]

# subset data

pbmc <- CreateSeuratObject(counts = m, meta.data = pD)
pbmc.list <- SplitObject(object = pbmc, split.by = "SampleID")


# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
pbmc.list <- lapply(pbmc.list, function(.x){
    temp <- NormalizeData(object = .x)
    temp <- FindVariableFeatures(object = temp, do.plot = FALSE)
    temp
})

# ALRA is undeterministic
# apply ALRA to impute non-technical zeroes
#  A: k=18,  7.49% nonzero to 47.66% nonzero
#  B: k=18,  6.97% nonzero to 46.35% nonzero
#  C: k=20,  7.10% nonzero to 43.39% nonzero
#  CT2-1NOV: k=23, 9.31% nonzero to 44.40% nonzero
#  CT2-30OCT: k=24, 9.48% nonzero to 46.32% nonzero

# Choose k.

choose_k <- function (A_norm,K=100, pval_thresh=1E-10, noise_start=80,q=2) {
    #  Heuristic for choosing rank k for the low rank approximation based on
    #  statistics of the spacings between consecutive singular values. Finds
    #  the smallest singular value \sigma_i such that $\sigma_i - \sigma_{i-1}
    #  is significantly different than spacings in the tail of the singular values.
    #
    #
    # Args:
    #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
    #   K: Number of singular values to compute. Must be less than the smallest dimension of the matrix.
    #   pval_thresh : The threshold for ``significance''
    #   noise_start : Index for which all smaller singular values are considered noise
    #   q : Number of additional power iterations
    #
    # Returns:
    #   A list with three items
    #       1) Chosen k
    #       2) P values of each possible k
    #       3) Singular values of the matrix A_norm

    if (K > min(dim(A_norm))) {
        stop("For an m by n matrix, K must be smaller than the min(m,n).\n")
    }
    if (noise_start > K-5) {
        stop("There need to be at least 5 singular values considered noise.\n")
    }
    noise_svals <- noise_start:K
    rsvd_out <- rsvd(A_norm,K,q=q)
    diffs <- diff(rsvd_out$d)
    pvals <- pnorm(diffs,mean(diffs[noise_svals-1]),sd(diffs[noise_svals-1]))
    k <- max(which( pvals  <pval_thresh))
    return (list( k=k, pvals =pvals, d=rsvd_out$d))
}

k_choice <- lapply(pbmc.list, function(.x) {
    choose_k(as.matrix(.x@assays$RNA@data))
})
# For the results in the paper, automatically chosen k worked quite well, but in
# some cases you might want to take a closer look, as we do here. The k is
# chosen based on the spacings between the singular values, as it can be quite
# hard to identify the ``beginning of noise'' from just looking at the spectrum
# itself. Uncomment the code below to plot them


for (i in 1:5)
{
    pdf(file.path(out_plot_dir, paste0("k-value.chose", i, ".pdf")), width = 25 , height = 8)
    df <- data.frame(x=1:100,y=k_choice[[i]]$d)
    g1<-ggplot(df,aes(x=x,y=y)) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice[[i]]$k)   + theme( axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_i') + ggtitle('Singular values')
    df <- data.frame(x=2:100,y=diff(k_choice[[i]]$d))[3:99,]
    g2<-ggplot(df,aes(x=x,y=y)) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice[[i]]$k+1)   + theme(axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_{i} - s_{i-1}') + ggtitle('Singular value spacings')
    grid.arrange(g1,g2,nrow=1)

    dev.off()
}

ks <- c(18, 18, 20, 23, 24)
pbmc.list <- mapply(function(.x, .y){
    temp <- RunALRA(object=.x, k = .y, q = 10, assay = NULL,
            slot = "data", setDefaultAssay = TRUE, genes.use = NULL,
            K = NULL, p.val.th = 1e-10, noise.start = NULL, q.k = 2,
            k.only = FALSE)
    temp
}, pbmc.list, ks)


save(file=file.path(out_objects_dir, "ALRA.imputed.RData"))





###  Apply Seurat to imputed data

load(file=file.path(out_objects_dir, "ALRA.imputed.RData"))
 
# ###### reconstruct Seurat object using imputed expression
# 
counts <- do.call(cbind, lapply(pbmc.list, function(.x){
    as.matrix(.x@assays$alra@data)
}))
pbmc <- CreateSeuratObject(counts = counts, meta.data = pD)
pbmc.list <- SplitObject(object = pbmc, split.by = "SampleID")

pbmc.list <- lapply(pbmc.list, function(.x){
    temp <- FindVariableFeatures(object = .x, do.plot = FALSE)
    temp
})

pbmc_int <- FindIntegrationAnchors(object.list = pbmc.list, scale = TRUE, dims = 1:30)



## get integrated expression values for all genes
pbmc.integrated <- IntegrateData(anchorset = pbmc_int,
                                 features.to.integrate = rownames(m), dims = 1:30)


## integrated analysis
DefaultAssay(object = pbmc.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
regression <- TRUE

if (regression){
    pbmc.integrated <- ScaleData(object = pbmc.integrated,
                                 vars.to.regress = c("UmiSums", "prcntMito"), verbose = FALSE)
} else {
    ## no regression out
    pbmc.integrated <- ScaleData(object = pbmc.integrated, verbose = FALSE)
}

# # pdf(file.path(out_plot_dir, "Highly variable gene plot.pdf"), width =15, height =5)
# # pbmc.integrated <- FindVariableFeatures(object = pbmc.integrated, do.plot = TRUE)
# # dev.off()
# 
pbmc.integrated <- RunPCA(object = pbmc.integrated, features = pbmc.integrated$integrated@var.features, npcs = 50, verbose = FALSE)

# ## plot variance
# sd <- pbmc.integrated@reductions$pca@stdev
# var <- sd^2/(sum(sd^2))*100
# 
# pdf(file.path(out_plot_dir, "vairance for PCA-ALRA.pdf"), height = 10, width = 10)
# plot(x=1:50, y=var, pch = 16, type= "b", ylab= "Variance (%)", xlab = "Principle component")
# dev.off()
# 
# ## UMAP: This depends on python package umap-learn
# 
pbmc.integrated <- RunUMAP(object = pbmc.integrated,  seed.use = 423,
                           reduction = "pca", dims = 1:22)
## TSNE
pbmc.integrated <- RunTSNE(object = pbmc.integrated, reduction = "pca", k.seed = 2,
                           dims = 1:22)
pbmc.integrated <- FindNeighbors(object = pbmc.integrated, reduction = "pca",
                                 dims = 1:22 )
pbmc.integrated <- FindClusters(object = pbmc.integrated, reduction = "pca",
                                dims = 1:22, save.SNN = TRUE)




pdf(file.path(out_plot_dir, "18 PCA-Tsne and Umap plot of cell clusters-ALRA.pdf"), width =15, height = 12)
p1 <- DimPlot(object = pbmc.integrated, reduction = "tsne", group.by = "SampleID", pt.size =0.5)
p2 <- DimPlot(object = pbmc.integrated, reduction = "tsne", do.return = TRUE, label = TRUE,  pt.size = 0.5)
p3 <- DimPlot(object = pbmc.integrated, reduction = "umap", group.by = "SampleID", pt.size =0.5)
p4 <- DimPlot(object = pbmc.integrated, reduction = "umap", do.return = TRUE, label = TRUE,  pt.size = 0.5)

plot_grid(p1, p2, p3, p4, nrow =2)
dev.off()


save.image(file=file.path(out_objects_dir, "Seurat.integrated.all.features-ALRA.RData"))

overlay <- function(file_name, width, height, features, reduction)
{
    pdf(file.path(out_plot_dir,file_name), height = height, width = width)
    print(FeaturePlot(object = pbmc.integrated, features = features,
                      dims = c(1, 2), cells = NULL,
                      cols = c("lightgrey", "red"), pt.size = 1, min.cutoff = "q9",
                      max.cutoff = NA, reduction = reduction, split.by = NULL,
                      shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
                      order = NULL, label = FALSE, label.size = 4, ncol = NULL,
                      combine = TRUE, coord.fixed = FALSE))
    dev.off()
}

signature_genes <- c("GNLY",  "CD5", "CD79B",  "MS4A1",  "CD69", "LYZ", "CD3E",
                     "CD8A",  "CD79A",  "TRDC",  "MS4A2",  "JCHAIN", "NCR1",
                     "FCER1G", "CD5",  "CST3", "ENSSSCG00000007978")

overlay(file_name = "overlay of markers on tSNE-ALRA.pdf",
        height =25, width = 25, reduction = "tsne",
        features = signature_genes)
overlay(file_name = "overlay of markers on UMAP.pdf",
        height =25, width = 25, reduction = "umap",
        features = signature_genes)



### plot all naostring detected genes
nanostring_genes <- read.delim("nanostring target ensembl gene ID.txt", header = TRUE, as.is = TRUE)
nanostring_genes <- nanostring_genes$Gene.name[nanostring_genes$Gene.stable.ID != ""]
nanostring_genes <- nanostring_genes[nanostring_genes %in% rownames(m)]

for (i in seq(1, length(nanostring_genes), 24))
{
    if( i <= 160)
    {
        overlay(file_name= paste0("Genes.expr.overlay.on.clusters-ALRA", i, ".pdf"),
                width = 25, height =25, reduction = "tsne",
               features = nanostring_genes[i:(i + 23)])
        next
    }else{
        overlay(file_name= paste0("Genes.expr.overlay.on.clusters-ALRA", i, ".pdf"),
                width = 25, height = 20, reduction = "tsne",
                features = nanostring_genes[i:length(nanostring_genes)])
    }
}

save.image(file=file.path(out_objects_dir, "ALRA.imputed.RData"))
##
all_markers <- FindAllMarkers(object = pbmc.integrated, test.use = "wilcox")
write.table(all_markers, file.path(out_plot_dir, "cluster-specific.markers.genes.across.all.clusters-ALRA.txt"), sep ="\t", quote = FALSE, row.names =TRUE)


#all_markers <- read.delim(file.path(out_plot_dir, "cluster-specific.markers.genes.across.all.clusters-ALRA.txt"), as.is = TRUE)
markers.use  <- subset(all_markers, avg_logFC >= 1 & p_val_adj <= 0.05)$gene

markers.use <- markers.use[markers.use %in% rownames(pbmc.integrated@assays$integrated@scale.data) ]

pdf(file.path(out_plot_dir,"clusterwise.markers.heatmap-ALRA.pdf"), height= 25, width = 15)
DoHeatmap(object = pbmc.integrated, features = markers.use, cells = NULL, 
          group.by = "ident", size = 2.5,
          group.bar = TRUE, disp.min = -3, disp.max = 3,
          slot = "scale.data", assay = "integrated", label = TRUE,
          hjust = 0, angle = 90, combine = TRUE)

dev.off()
