##  module load r/3.5.0-py2-qqwf6c6
##  source  ~/bin/system.py3.6.5_env/bin/activate

##  install Seurat Realease 3.0
#install.packages('devtools')
#devtools::install_github(repo = 'satijalab/seurat', ref = 'release/3.0')


# cowplot enables side-by-side ggplots
library(cowplot)
library(Seurat)


out_objects_dir <- "./results3.0/R.out/data/Robjects"
out_plot_dir <- "./results3.0/R.out/plots"

out_document_dir <- "./results3.0/R.out/results"
if(!dir.exists(out_document_dir)) 
{
    dir.create(out_document_dir, recursive = TRUE)
}


library_id <- c("A", "B", "C", "CT2-1NOV", "CT2-30OCT")    

dataList <- readRDS(file=file.path(out_objects_dir, "ExpressionList_QC.rds"))

m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
fD$keep[is.na(fD$keep)] <- FALSE
rm(dataList)

# Gene and cell filtering

m <- m[fD$keep, pD$PassAll]   ## 11707 genes X 15847 cells  # CR3.0 11442 X 16668
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
    temp <- FindVariableFeatures(object = temp)
    temp
})

### Integration of 5 PBMC cell datasets
pbmc_int <- FindIntegrationAnchors(object.list = pbmc.list, dims = 1:30)
pbmc.integrated <- IntegrateData(anchorset = pbmc_int, dims = 1:30)


## integrated analysis
DefaultAssay(object = pbmc.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pbmc.integrated <- ScaleData(object = pbmc.integrated, 
                             vars.to.regress = c("UmiSums", "prcntMito"), 
                             verbose = FALSE)

pbmc.integrated <- RunPCA(object = pbmc.integrated, 
                          features = pbmc.integrated$integrated@var.features, 
                          npcs = 50, verbose = FALSE)

## plot variance
sd <- pbmc.integrated@reductions$pca@stdev
var <- sd^2/(sum(sd^2))*100

pdf(file.path(out_plot_dir, "1.7.Scree plot of vairance for PCA.pdf"), 
    height = 10, width = 10)
plot(x=1:50, y=var, pch = 16, type= "b", 
     ylab= "Variance (%)", xlab = "Principle component")
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

pdf(file.path(out_plot_dir, "1.8.18 PCA-Tsne and Umap plot of cell clusters.pdf"), 
    width = 15, height = 12)
p1 <- DimPlot(object = pbmc.integrated, reduction = "tsne", 
              group.by = "SampleID", pt.size =0.5)
p2 <- DimPlot(object = pbmc.integrated, reduction = "tsne", 
              do.return = TRUE, label = TRUE,  pt.size = 0.5)
p3 <- DimPlot(object = pbmc.integrated, reduction = "umap", 
              group.by = "SampleID", pt.size =0.5)
p4 <- DimPlot(object = pbmc.integrated, reduction = "umap", 
              do.return = TRUE, label = TRUE,  pt.size = 0.5)

plot_grid(p1, p2, p3, p4, nrow =2)
dev.off()

all_markers <- FindAllMarkers(object = pbmc.integrated, test.use = "wilcox")
write.table(all_markers, file = file.path(out_document_dir, 
                                          "1.0.All.marker.genes.no.imputation.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

markers.use <- subset(all_markers, avg_logFC >= 1)$gene

pdf(file.path(out_plot_dir,"1.8.Markers.plot.pdf"), height = 40, width = 15)
DoHeatmap(object = pbmc.integrated, 
          features = markers.use, 
          cells = NULL, 
          group.by = "ident", size =1.5,
          group.bar = TRUE, disp.min = -2.5, disp.max = NULL,
          slot = "scale.data", assay = NULL, label = TRUE,
          hjust = 0, angle = 90, combine = TRUE)
dev.off()

markers <- c("CD3E", "CD4","CD5", "CD8A", "CD8B", "TRDC",  
             "GZMB", "IFNG", "CD79A", "CD79B", "CD19", 
             "CD69", "MS4A1", "FCER1G", "MS4A2", "JCHAIN",  
             "ITGAM","FCGR1A", "CD14", "SERPING1", "MX1",  
             "IL1RAP", "IFNGR1", "CST3", "TLR4",
             "NCR1", "KLRB1",  "GNLY", "LYZ", "MCM2",
             "MCM3", "TOP2A", "CCNB1", "PCNA")
features <- do.call("c", lapply(markers, function(.x) {
    rownames(fD)[grepl(paste0("-", .x, "$"), rownames(fD), perl = TRUE)]
}))

## add hemoglobin alpha gene, ""ENSSSCG00000007978"
features <- c(features, "ENSSSCG00000007978")

pdf(file.path(out_plot_dir,"1.9.Overlay of markers.pdf"), 
    height = 42, width = 31)
FeaturePlot(object = pbmc.integrated, features = features, dims = c(1, 2), 
            cells = NULL, cols = c("lightgrey", "red"), 
            pt.size = 1, min.cutoff = "q9",
            max.cutoff = NA, reduction = "tsne", split.by = NULL,
            shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
            order = NULL, label = TRUE, label.size = 4, ncol = 5,
            combine = TRUE, coord.fixed = TRUE, sort.cell = TRUE)
dev.off()

    
pdf(file.path(out_plot_dir,"2.0.Overlay of markers on umap.pdf"),
    height = 42, width = 31)
FeaturePlot(object = pbmc.integrated, features = features, 
            dims = c(1, 2), cells = NULL,
            cols = c("lightgrey", "red"), pt.size = 1, min.cutoff = "q9",
            max.cutoff = NA, reduction = "umap", split.by = NULL,
            shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
            order = NULL, label = TRUE, label.size = 4, ncol = 5,
            combine = TRUE, coord.fixed = TRUE, sort.cell = TRUE)
dev.off()

save.image(file = file.path(out_objects_dir, 
                            "Seurat.integrated.without.imputation.RData"))
