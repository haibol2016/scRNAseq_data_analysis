##  module load r/3.5.0-py2-qqwf6c6
##  source  ~/bin/system.py3.6.5_env/bin/activate

##  install Seurat Realease 3.0
#install.packages('devtools')
#devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')


# cowplot enables side-by-side ggplots
library(cowplot)
library(Seurat)


setwd("/home/haibol/haibol_work/scRNA-seq")

out_objects_dir <- "./results/R.out/data/Robjects"
out_plot_dir <- "./results/R.out"
library_id <- c("A", "B", "C", "CT2-1NOV", "CT2-30OCT")    
    
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

### Integration of 3 pancreatic islet cell datasets
pbmc_int <- FindIntegrationAnchors(object.list = pbmc.list, dims = 1:30)
pbmc.integrated <- IntegrateData(anchorset = pbmc_int, dims = 1:30)


## integrated analysis
DefaultAssay(object = pbmc.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pbmc.integrated <- ScaleData(object = pbmc.integrated, vars.to.regress = c("UmiSums", "prcntMito"), verbose = FALSE)

# pdf(file.path(out_plot_dir, "Highly variable gene plot.pdf"), width =15, height =5)
# pbmc.integrated <- FindVariableFeatures(object = pbmc.integrated, do.plot = TRUE)
# dev.off()

pbmc.integrated <- RunPCA(object = pbmc.integrated, features = pbmc.integrated$integrated@var.features, npcs = 50, verbose = FALSE)

## plot variance
sd <- pbmc.integrated@reductions$pca@stdev
var <- sd^2/(sum(sd^2))*100

pdf(file.path(out_plot_dir, "vairance for PCA.pdf"), height = 10, width = 10)
plot(x=1:50, y=var, pch = 16, type= "b", ylab= "Variance (%)", xlab = "Principle component")
dev.off()

## UMAP: This depends on python package umap-learn

pbmc.integrated <- RunUMAP(object = pbmc.integrated, 
                           reduction = "pca", dims = 1:18)
## TSNE
pbmc.integrated <- RunTSNE(object = pbmc.integrated, reduction = "pca", 
                           dims = 1:18)
pbmc.integrated <- FindNeighbors(object = pbmc.integrated, reduction = "pca", 
              dims = 1:18 )
pbmc.integrated <- FindClusters(object = pbmc.integrated, reduction = "pca", 
                                dims = 1:18, save.SNN = TRUE)

pdf(file.path(out_plot_dir, "18 PCA-Tsne and Umap plot of cell clusters.pdf"), width =15, height =12)
p1 <- DimPlot(object = pbmc.integrated, reduction = "tsne", group.by = "SampleID", pt.size =0.5)
p2 <- DimPlot(object = pbmc.integrated, reduction = "tsne", do.return = TRUE, label = TRUE,  pt.size = 0.5)
p3 <- DimPlot(object = pbmc.integrated, reduction = "umap", group.by = "SampleID", pt.size =0.5)
p4 <- DimPlot(object = pbmc.integrated, reduction = "umap", do.return = TRUE, label = TRUE,  pt.size = 0.5)

plot_grid(p1, p2, p3, p4, nrow =2)
dev.off()

all_markers <- FindAllMarkers(object = pbmc.integrated, test.use = "wilcox")

markers.use=subset(all_markers, avg_logFC >= 1)$gene

pdf("Markers.plot.pdf", height= 20, width =15)
DoHeatmap(object = pbmc.integrated, features = markers.use, cells = NULL, 
          group.by = "ident", size =1.5,
          group.bar = TRUE, disp.min = -2.5, disp.max = NULL,
          slot = "scale.data", assay = NULL, label = TRUE,
          hjust = 0, angle = 90, combine = TRUE)

dev.off()



pdf(file.path(out_plot_dir,"overlay of markers.pdf"), height =15, width = 25)
FeaturePlot(object = pbmc.integrated, features = c("ENSSSCG00000008228", ## GNLY
                         "ENSSSCG00000013115", ## CD5
                         "ENSSSCG00000017283", ## CD79B
                         "ENSSSCG00000021812", ## MS4A1
                         "ENSSSCG00000000653", ## CD69
                         "ENSSSCG00000000492", ## LYZ
                         "ENSSSCG00000040140", ## CD3E
                         "ENSSSCG00000008217", ## CD8A
                         "ENSSSCG00000033684", ## CD79A
                         "ENSSSCG00000022512", ## TRDC 
                         "ENSSSCG00000034506", ## MS4A2
                         "ENSSSCG00000035379", ## JCHAIN
                         "ENSSSCG00000022675", ## NCR1
                         "ENSSSCG00000006357", ## FCER1G
                         "ENSSSCG00000013115", ## CD5
                         "ENSSSCG00000037360", ## CST3
                         "ENSSSCG00000007978"), dims = c(1, 2), cells = NULL,
            cols = c("lightgrey", "red"), pt.size = 1, min.cutoff = "q9",
            max.cutoff = NA, reduction = "tsne", split.by = NULL,
            shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
            order = NULL, label = TRUE, label.size = 4, ncol = NULL,
            combine = TRUE, coord.fixed = FALSE)
dev.off()


pdf(file.path(out_plot_dir,"overlay of markers on umap.pdf"), height =15, width = 25)
FeaturePlot(object = pbmc.integrated, features = c("ENSSSCG00000008228", ## GNLY
                         "ENSSSCG00000013115", ## CD5
                         "ENSSSCG00000017283", ## CD79B
                         "ENSSSCG00000021812", ## MS4A1
                         "ENSSSCG00000000653", ## CD69
                         "ENSSSCG00000000492", ## LYZ
                         "ENSSSCG00000040140", ## CD3E
                         "ENSSSCG00000008217", ## CD8A
                         "ENSSSCG00000033684", ## CD79A
                         "ENSSSCG00000022512", ## TRDC 
                         "ENSSSCG00000034506", ## MS4A2
                         "ENSSSCG00000035379", ## JCHAIN
                         "ENSSSCG00000022675", ## NCR1
                         "ENSSSCG00000006357", ## FCER1G
                         "ENSSSCG00000013115", ## CD5
                         "ENSSSCG00000037360", ## CST3
                         "ENSSSCG00000007978"), dims = c(1, 2), cells = NULL,
            cols = c("lightgrey", "red"), pt.size = 1, min.cutoff = "q9",
            max.cutoff = NA, reduction = "umap", split.by = NULL,
            shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
            order = NULL, label = TRUE, label.size = 4, ncol = NULL,
            combine = TRUE, coord.fixed = FALSE)
dev.off()

save.image(file = "Seurat.RData")

#####
markers <- read.delim("C:\\Users\\Haibo\\Desktop\\Tuggle_lab-help\\Nanostring\\10-30-2018 PBMC scRNAseq/all.markers.between.clusters.txt", header = TRUE)

pig_gene <- read.delim("C:\\Users\\Haibo\\Desktop\\Tuggle_lab-help\\Huber\\10-30-2018 macrophage LPS RNA-seq/pig_genes.Ensembl.92.txt")

markers <- merge(markers, pig_gene, by.x= "gene", by.y = "EnsemblID", all.x=T)
write.table(markers, file = "C:\\Users\\Haibo\\Desktop\\Tuggle_lab-help\\Nanostring\\10-30-2018 PBMC scRNAseq/all.markers.between.clusters.txt", sep ="\t", quote = FALSE, row.names = FALSE)

