## load R 3.5 on nova 
## module load r/3.5.0-py2-qqwf6c6

## install core packages from Bioconductor
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()

## install all dependent packages
packages2install <- c("dplyr", "cellrangerRkit", "cowplot", 
                      "knitr", "scran", "limma", "Seurat")
uninstalled_packages <- packages2install[!packages2install 
                                         %in% rownames(installed.packages())]
lapply(uninstalled_packages, BiocManager::install)

## load packages
packages2load <- c("cellrangerRkit", "dplyr", "Seurat")
lapply(packages2load, library, character.only = TRUE)

# Script to prepare cellranger data for downstream analysis
# Read in output from cell ranger using cellrangerRkit
samples <- c("./A", "./B", "./C", "./CT2-1NOV", "./CT2-30OCT")
library_id <- gsub(".+\\/", "", samples, perl = TRUE)

## For CellRanger 3.x or 2.x, use a named vector specifying the path to 10X output 
## for further analysis Within each folder, there are three .gz files for CellRanger V2.x, 
## and four uncompressed files for CellRanger V3.x. The name will be used as prefix 
## of cell barcodes in Seurat object

data_dir <- file.path(getwd(),library_id, "outs/filtered_feature_bc_matrix")
names(data_dir) <- library_id

## Whether the output is from CellRanger v2.x or CellRanger v3.x
pre_version3.0 <- FALSE  # CellRanger v2.x

##### For output from CellRanger v2.x
if (pre_version3.0) 
{
    gene_bc_matrix <- lapply(samples, load_cellranger_matrix, genome="ssc3")
    
    ## Cell Barcode
    pDat <- do.call(rbind, lapply( gene_bc_matrix, function(.x){
        df <- data.frame(pData(.x))
        print(dim(df))
        df
    }))
    
    ## number of detected cells
    numberOfCellperLib <- do.call(c, lapply(gene_bc_matrix, function(.x){
        nrow(pData(.x))
    }))
    
    ## count table
    cDat <- do.call(cbind, lapply(gene_bc_matrix, function(.x){
        df <-as.matrix(exprs(.x))
        df <- df[order(rownames(df)), ]
        df
    }))
    
    # reduce size of matrix
    keep <- rowSums(cDat) > 1
    cDat <- cDat[keep,]   ##  13064 genes,  8767 cells
    
    ## gene symbol mapping
    fDat <- data.frame(fData(gene_bc_matrix[[1]]))[rownames(cDat),]
    
    # Add more info to phenotype Data
    pDat <- mutate(pDat, barcode=as.character(barcode)) %>%
        mutate(SampleID= do.call(c, mapply(function(.x, .y){
            rep(.x, .y)
        }, library_id, numberOfCellperLib, SIMPLIFY = FALSE)))
    
} else {
    lapply(data_dir, dir) ## Should show barcodes.tsv, genes.tsv, and matrix.mtx
    scRNA_data <- Read10X(data.dir = data_dir)
    seurat_object = CreateSeuratObject(counts = scRNA_data)
    
    ## Cell Barcode
    pDat <-data.frame(barcode = colnames(seurat_object))
    pDat$SampleID <- gsub("([^_]+).+", "\\1", pDat$barcode, perl = TRUE)
    
    ## count table
    cDat <-  as.matrix(GetAssayData(object = seurat_object, slot = 'counts'))
    # reduce size of matrix
    keep <- rowSums(cDat) > 1
    cDat <- cDat[keep,]   ##  13064 genes,  8767 cells
    ## gene symbol mapping
    fDat <- data.frame(ID = rownames(cDat))
    rownames(fDat) <- fDat$ID 
    
    ## augmented fDat by adding gene symbols
    con <- gzfile(file.path(data_dir[1], "features.tsv.gz"))
    ssc_genes <- read.delim(con, sep = "\t", header = FALSE, as.is = TRUE)[, 1:2]
    colnames(ssc_genes) <- c("EnsemblID", "Symbol")
    ssc_genes$Symbol<- gsub("_", "-", ssc_genes$Symbol, perl = TRUE)
    
    fDat <- merge(fDat, ssc_genes, by.x ="row.names", 
                  by.y = "Symbol", all.x =TRUE, sort = FALSE)
    rownames(fDat) <- fDat[, 1]
    fDat <- fDat[, -1]

}

# Add more info to the feature Data
## extract mitochondrial genes
gtf <- "/home/haibol/haibol_work/haibo_genome/sus_scrofa/Sus_scrofa.Sscrofa11.1.97.gtf"

mitoGenes <-  system2("grep", 
                    args = c('^MT', gtf, "| grep -o 'ENSSSCG[0-9]*' | uniq"), 
                    stdout = TRUE)

fDat$Mitochondrial <- fDat$EnsemblID %in% mitoGenes


# Save data
stopifnot(identical(rownames(fDat),rownames(cDat)) & 
          identical(colnames(cDat),pDat$barcode))
DataList <- list("phenoData"= pDat, "featureData"=fDat, "counts"=cDat)

out_objects_dir <- "./results3.0/R.out/data/Robjects"
if(!dir.exists(out_objects_dir)) {dir.create(out_objects_dir, recursive = TRUE)}

saveRDS(DataList,file=file.path(out_objects_dir, "ExpressionList.rds"))
