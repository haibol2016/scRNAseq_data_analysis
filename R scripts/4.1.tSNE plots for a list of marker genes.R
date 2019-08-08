
## Do this after ALRA-based imputation, Seurat alignment, integration, 
## normalization and tSNE dimension reduction

library("Seurat")
library("ggplot2")
library("gridExtra")


## load the R data containing imputed, integrated single cell expression data
load("ALRA.imputed.RData")

signature_genes <- c("CD3E", "CD4", "CD8A", "CD8B", "GZMB", "IFNG",
                     "CD79A", "CD79B", "CD19", "FCER1G", "ITGAM", "FCGR1A",
                     "CD14", "SERPING1", "MX1", "NOD1", "NLRP3", "IFNGR1",
                     "TRDC", "CST3", "NCR1", "KLRB1", "CD5", "GNLY")

overlay.p2 <- FeaturePlot(object = pbmc.integrated, features = signature_genes,
                          dims = c(1, 2), cells = NULL,
                          cols = c("lightgrey", "red"), 
                          pt.size = 0.3, min.cutoff = "q9",
                          max.cutoff = NA, reduction = "tsne", split.by = NULL,
                          shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
                          order = NULL, label = FALSE, label.size = 4, ncol = 6,
                          combine = FALSE, coord.fixed = FALSE)

## apply ggplot comfiguration setting to customize each ggplot object
p2 <- lapply(overlay.p2, function(.x){
    .x + theme(text = element_text(size = 10),
               legend.position="bottom",
               legend.box = "horizontal",
               legend.key.width = unit(1.0, "cm"),
               legend.key.height = unit(0.3,"cm"),
               plot.title = element_text(hjust = 0.5, size = 10, face = "plain"))
})


## ouput each figures
pdf("08-06-2019.tSNE.plots.showing.expression.for.18.marker.genes.pdf", width = 18, height = 13)

grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], p2[[5]], p2[[6]],
             p2[[7]],  p2[[8]], p2[[9]], p2[[10]], p2[[11]], p2[[12]],
             p2[[13]],  p2[[14]], p2[[15]], p2[[16]], p2[[17]], p2[[18]],
             p2[[19]],  p2[[20]], p2[[21]], p2[[22]], p2[[23]], p2[[24]],
             layout_matrix = matrix(1:24, byrow = TRUE, nrow =4))
dev.off()