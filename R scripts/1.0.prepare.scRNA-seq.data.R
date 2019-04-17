
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
packages2load <- c("cellrangerRkit", "dplyr")
lapply(packages2load, library, character.only = TRUE)

# Script to prepare cellranger data for downstream analysis

# ---- ReadData ----
# Read in output from cell ranger using cellrangerRkit
samples <- c("./A", "./B", "./C", "./CT2-1NOV", "./CT2-30OCT")
library_id <- gsub(".+\\/", "", samples, perl = TRUE)
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
    df}))

# reduce size of matrix
keep <- rowSums(cDat) > 1
cDat <- cDat[keep,]   ##  13064 genes,  8767 cells

## gene symbol mapping
fDat <- data.frame(fData(gene_bc_matrix[[1]]))[rownames(cDat),]

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
