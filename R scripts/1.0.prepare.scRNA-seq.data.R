
## load R 3.5 on nova 
## module load r/3.5.0-py2-qqwf6c6


## install core packages from Bioconductor
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()

## install all dependent packages
packages2install <- c("dplyr", "cellrangerRkit", "cowplot", "knitr", "scran", "limma", "Seurat")
uninstalled_packages <- packages2install[!packages2install %in% rownames(installed.packages())]
lapply(uninstalled_packages, BiocManager::install)

## load packages
packages2load <- c("cellrangerRkit", "dplyr", "Seurat")
lapply(packages2load, library, character.only = TRUE)

# Script to prepare cellranger data for downstream analysis

# ---- ReadData ----
# Read in output from cell ranger using cellrangerRkit
samples <- c("./A", "./B", "./C", "./CT2-1NOV", "./CT2-30OCT")
library_id <- gsub(".+\\/", "", samples, perl = TRUE)


## for CellRanger 3.x or 2.x
## use a named vector specifying the path to 10X output for further analysis
## Within each folder, there are three .gz files for CellRanger V2.x, 
## and four uncompressed files for CellRanger V3.x

data_dir <- c(PBMC1="C:\\Users\\Haibo\\Downloads\\PBMC1\\outs\\filtered_feature_bc_matrix",
              PBMC2="C:\\Users\\Haibo\\Downloads\\PBMC2\\outs\\filtered_feature_bc_matrix")

numberOfCellperLib <- vector("numeric", length(samples))

## Whether the output is from CellRanger v2.x or CellRanger v3.x
pre_version3.0 <- TRUE  # CellRanger v2.x


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
} else {
    lapply(data_dir, dir) ## Should show barcodes.tsv, genes.tsv, and matrix.mtx
    scRNA_data <- Read10X(data.dir = data_dir)
    seurat_object = CreateSeuratObject(counts = scRNA_data)
    
    ## Cell Barcode
    pDat <-data.frame(barcode = colnames(seurat_object))
    ## number of cells detected from each dataset
    numberOfCellperLib <- do.call(c, lapply(names(data_dir), function(.x) {
        sum(grepl(pattern = .x, colnames(seurat_object), fixed = TRUE))
    }))
    ## count table
    cDat <-  as.matrix(GetAssayData(object = seurat_object, slot = 'counts'))
    # reduce size of matrix
    keep <- rowSums(cDat) > 1
    cDat <- cDat[keep,]   ##  13064 genes,  8767 cells
    ## gene symbol mapping
    fDat <- data.frame(feature = rownames(cDat))
}    


# ---- Formatting ----

# Add more info to phenotype Data
pDat <- mutate(pDat, barcode=as.character(barcode)) %>%
    mutate(SampleID= do.call(c, mapply(function(.x, .y){
        rep(.x, .y)
    }, library_id, numberOfCellperLib, SIMPLIFY = FALSE)))
 

# Add more info to the feature Data
## extract mitochondrial genes
## mkdir -p 
## grep '^MT' Sus_scrofa.Sscrofa11.1.92.gtf | grep -o 'ENSSSCG[0-9]*' | uniq > ./results/R.out/data/miscData/MitoGenes.txt

mitoGenes <- read.table("./results/R.out/data/miscData/MitoGenes.txt")

# tfCheck <- read.table("../data/miscData/TFcheckpoint_WithENSID.tsv",
#                       header=TRUE, sep="\t")

fDat$Mitochondrial <- fDat$id %in% mitoGenes$V1


# Save data
stopifnot(identical(rownames(fDat),rownames(cDat)) & identical(colnames(cDat),pDat$barcode))
DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)

out_objects_dir <- "./results/R.out/data/Robjects"
if(!dir.exists(out_objects_dir)) {dir.create(out_objects_dir)}

saveRDS(DataList,file=file.path(out_objects_dir, "ExpressionList.rds"))
