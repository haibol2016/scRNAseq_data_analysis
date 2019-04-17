# QC-Analysis

## load packages
packages2load <- c("scran", "dplyr", "knitr", "ggplot2", "Rtsne", "cowplot")
lapply(packages2load, library, character.only = TRUE)

source("./functions.R")

# Load Data
out_objects_dir <- "./results/R.out/data/Robjects"
out_plot_dir <- "./results/R.out"

dataList <- readRDS(file.path(out_objects_dir, "ExpressionList.rds"))
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

# ---- QCOverview ----

# Sequencing Depth and Genes detected
pD$UmiSums<- colSums(m)
pD$GenesDetected <- colSums(m!=0)
pdf(file.path(out_plot_dir, "Sequence depth and number of gene detected real scale.pdf"), width = 7, height = 4)
genesDetected <- ggplot(pD, aes(x=SampleID, y=GenesDetected, fill= SampleID)) +
    geom_violin(draw_quantiles=0.5)+
    #scale_y_log10() +
    ylab("Total number of genes detected")
print(genesDetected)    
LibrarySize <- ggplot(pD, aes(x=SampleID,y=UmiSums, fill=SampleID)) +
    geom_violin(draw_quantiles=0.5)+
    #scale_y_log10() + 
    ylab("Total number of molecules(depth)")
print(LibrarySize)    
dev.off()





########################### Genewise QC ###################################
ntop <- 50 
mRel <- t(t(m)/colSums(m))
rownames(mRel)  <- fD$symbol
topExpressed <- rowMedians(mRel) ## library size normalized expression
names(topExpressed) <- rownames(mRel)
topExpressed <- topExpressed %>% sort(.,decreasing=TRUE) %>% names
plotData <- t(mRel)[,topExpressed[1:ntop]] %>%
    reshape2::melt() %>%
    dplyr::rename(Cell=Var1,
                  Gene=Var2,
                  RelativeExpression=value)
## Add mitochondrial gene information
plotData <- merge(plotData, fD, by.x = "Gene", by.y = "symbol", all.x = TRUE)

pdf(file.path(out_plot_dir,"Top 50 gene expression.pdf"), width = 7, height =7)
topGenes <- ggplot(plotData, aes(x=Gene, y=RelativeExpression, color= Mitochondrial)) +
    geom_boxplot() +
    coord_flip() +
    theme_bw()
print(topGenes)


freqOfExp <- m!=0
rownames(freqOfExp) <- fD$symbol
freqOfExp <- sort(rowSums(freqOfExp)/ncol(freqOfExp),decreasing=TRUE)
plotData <- data.frame("Gene"=names(freqOfExp),"Frequency"=freqOfExp)
plotData <- merge(plotData, fD, by.x = "Gene", by.y = "symbol", all.x = TRUE, sorted =FALSE)
plotData <- plotData[order(plotData$Frequency, decreasing= TRUE), ]

## most frequently detected genes across all cells
topFreq <- ggplot(plotData[1:ntop,], aes(x=factor(Gene,levels=Gene), y=Frequency, color= Mitochondrial)) +
    geom_bar(stat="identity", fill ="white") +
    coord_flip() +
    xlab("Gene") +
    theme_bw()
print(topFreq)
dev.off()




#########################  Cell Viability ##################################
theme_set(theme_grey())
pdf(file.path(out_plot_dir, "Cell viability plot.pdf"), width = 12, height = 4)
mMito <- m[fD$Mitochondrial,]
idtop <- fD[fD$symbol %in% names(freqOfExp)[1:ntop],"id"]
mTop <- m[idtop,]!=0
pD$prcntTop <- colSums(mTop)/ntop
pD$prcntMito <- colSums(mMito)/colSums(m)
cellViability <- ggplot(pD, aes(x=prcntMito, y=GenesDetected, color = SampleID))+
    geom_point() + facet_wrap(~SampleID, nrow =1)+theme_get() + ylab("#Genes detected per cell")
cellViability
prcntTop <- ggplot(pD, aes(x=prcntTop, y=GenesDetected, color = SampleID))+
    geom_point() + facet_wrap(~SampleID, nrow =1) +theme_get()
prcntTop
dev.off()

rm(mMito)
rm(mTop)



# ---- Index-Swapping ----
library_id <- c("A", "B", "C", "CT2-1NOV", "CT2-30OCT")
for (i in seq_along(library_id))
{
    pD$barcode <- ifelse(pD$SampleID== library_id[i], gsub("-\\d+", paste0("-",i), pD$barcode), pD$barcode)
}

colnames(m) <- pD$barcode

barcodes <- as.character(pD$barcode)
spt <- strsplit(barcodes, split = "-", fixed = T)
pD$sample <- sapply(spt, function(x) x[2])
pD$bcs <- sapply(spt, function(x) x[1])
pD.add <- data.frame(bcs=names(table(pD$bcs)),
                     bc.obs=(table(pD$bcs)))
pD.add <- pD.add[,-2]
pD <- dplyr::left_join(pD,pD.add)

pdf(file.path(out_plot_dir, "Index swapping.pdf"), width = 7, height =7)
# P1
index.p1 <- ggplot(pD, aes(x=SampleID, fill=as.factor(bc.obs.Freq))) +
    geom_bar() +
    #     ggtitle("Samples contain many shared barcodes") +
    scale_fill_discrete(name="Times barcode observed") +
    theme(legend.position="bottom",
          legend.direction="horizontal")

index.p1

pD$SampleID <- factor(pD$SampleID)

compare <- function(barcodes, samples) {
    out <- NULL
    ids <- levels(samples)
    combs <- combn(ids,m=2, simplify=FALSE)
    for (i in seq_along(combs)) {
        comb <- combs[[i]]
        s1 <- comb[1]
        s2 <- comb[2]
        bc1 <- as.character(barcodes[samples==s1])
        bc2 <- as.character(barcodes[samples==s2])
        x <- length(intersect(bc1,bc2))
        m <- length(bc1)
        n <- 750000
        k <- length(bc2)
        p.val <- phyper(q=x-1,m=m,n=n,k=k,lower.tail=FALSE)
        tmp <- data.frame(s1=s1,
                          s2=s2,
                          n1=m,
                          n2=k,
                          shared=x,
                          p.val=p.val)
        out <- rbind(out,tmp)
    }
    return(out)
}

# P2
compDf <- compare(pD$bcs, pD$SampleID)
index.p2 <- ggplot(compDf, aes(x=shared, y=-log10(p.val))) +
    geom_point() +
    xlab("# shared Barcodes") +
    ylab("-log10(P)") +
    geom_hline(yintercept=2,lty="dashed",color="red") 

index.p2
dev.off()




# ---- QCThresholding ----

# Fixed threshold 0.05 seems too stringent
MitoCutOff <- 0.1

#Thresholding per Sample
smryByGroup <- group_by(pD, SampleID) %>%
    summarize(
        mPrcntMito=median(prcntMito),
        madPrcntMito=mad(prcntMito),
        threshold_PrcntMito=min(mPrcntMito + 3* madPrcntMito, MitoCutOff),
        
        mGenesDetected=median(log10(GenesDetected)),
        madGenesDetected=mad(log10(GenesDetected)),
        
        ## median - 3*MAD
        threshold_GenesDetected=max(mGenesDetected-3*madGenesDetected,log10(500)),
        
        mUmiSums=median(log10(UmiSums)),
        madUmiSums=mad(log10(UmiSums)),
        
        ## median depth - 3*MAD
        threshold_UmiSums=max(mUmiSums-3*madUmiSums,log10(1000))) %>%
    
    select(SampleID,starts_with("threshold")) 
kable(smryByGroup)


#  |SampleID  | threshold_PrcntMito| threshold_GenesDetected| threshold_UmiSums|
#  |:---------|-------------------:|-----------------------:|-----------------:|
#  |A         |           0.0842837|                 2.69897|                 3|
#  |B         |           0.0808965|                 2.69897|                 3|
#  |C         |           0.0752058|                 2.69897|                 3|
#  |CT2-1NOV  |           0.0833021|                 2.69897|                 3|
#  |CT2-30OCT |           0.0721746|                 2.69897|                 3|
    

pD <- mutate(pD,
             ThresholdViability = 0.1,  ## based on plotting, set to 0.1
             ThresholdGenesDet = 0,
             ThresholdLibSize = 0)

grps <- as.character(unique(pD$SampleID))
for (grp in grps) {
    thrs <- filter(smryByGroup, SampleID==grp) %>% select(-SampleID) %>% t() %>% as.vector()
    names(thrs) <- filter(smryByGroup, SampleID==grp) %>% select(-SampleID) %>% t() %>% rownames()
    pD <- mutate(pD,
                 #ThresholdViability= ifelse(SampleID==grp, MitoCutOff, ThresholdViability),
                 ThresholdGenesDet= ifelse(SampleID==grp, 10^thrs["threshold_GenesDetected"],ThresholdGenesDet),
                 ThresholdLibSize= ifelse(SampleID==grp, 10^thrs["threshold_UmiSums"],ThresholdLibSize))
}

pD <- mutate(pD,
             PassViability=prcntMito < ThresholdViability,
             PassGenesDet=GenesDetected > ThresholdGenesDet,
             PassLibSize=UmiSums > ThresholdLibSize,
             PassAll= PassViability & PassGenesDet & PassLibSize & bc.obs.Freq==1)

# Illustrate thresholds
gdHist <- ggplot(pD, aes(x=GenesDetected,y=..density..)) +
    geom_histogram(fill="white",color="black",bins=100) +
    geom_vline(data=smryByGroup,aes(xintercept=10^threshold_GenesDetected),color="red",lty="longdash") +
    scale_x_log10() +
    xlab("Total number of genes detected") +
    facet_wrap(~ SampleID) 

libSizeHist <- ggplot(pD, aes(x=UmiSums,y=..density..)) +
    geom_histogram(fill="white",color="black",bins=100) +
    geom_vline(data=smryByGroup,aes(xintercept=10^threshold_UmiSums),color="red",lty="longdash") +
    scale_x_log10() +
    facet_wrap(~SampleID) +
    xlab("Total number of unique molecules") 

cellViability <- cellViability %+% pD
cellViability <- cellViability + 
    annotate("rect",ymin=-Inf, ymax=Inf, xmax=Inf, xmin=MitoCutOff,
             fill="grey", alpha=0.3) 


# Overview over cells removed
table(pD$SampleID,pD$PassGenesDet)
#           FALSE TRUE
# A           498 3109
# B           407 3249
# C           183 3291
# CT2-1NOV      6 2946
# CT2-30OCT    25 4052


table(pD$SampleID,pD$PassLibSize)
#           FALSE TRUE
# A             0 3607
# B           300 3356
# C             0 3474
# CT2-1NOV      0 2952
# CT2-30OCT     0 4077


table(pD$SampleID,pD$PassViability)
#           FALSE TRUE
# A            18 3589
# B            23 3633
# C            15 3459
# CT2-1NOV    189 2763
# CT2-30OCT   100 3977



table(pD$SampleID,pD$PassAll)
#           FALSE TRUE
# A           584 3023
# B           611 3045
# C           274 3200
# CT2-1NOV    250 2702
# CT2-30OCT   200 3877




pdf(file.path(out_plot_dir, "Thresholding library size and cell viability.pdf"), height =7 ,width =7)
  gdHist
  libSizeHist
  cellViability
dev.off()


nCells <- as.data.frame(table(pD$SampleID))
colnames(nCells) <- c("replicate", "count")

sample_id <- as.character(unique(pD$SampleID))
ave.counts <- vector("list", length = length(sample_id))
for (i in seq_along(sample_id))
{
    ave.counts[[i]] <- rowMeans(m[, which(pD$SampleID == sample_id[i])])
}
names(ave.counts) <- sample_id


pdf(file.path(out_plot_dir, "average counts in each library.pdf"), width =12, height =8)
op <- par(mfrow = c(2,3))

plot_ave <- function(x){hist(log10(x), breaks=100, main="", col="grey80",
                             xlab=expression(log[10]~"average count"))}
lapply(ave.counts, plot_ave)

dev.off()
par(op)

## This is not very good criterion
#fD$keep <- rowMeans(m) > 10^(-2.5)

## based on number of cell (passed filtering), which express a gene, to determine
## filtering of genes (0.1%)
m_filter <- m[, pD$PassAll]
fD$keep <- rowSums(matrix(m_filter !=0, nrow = nrow(m_filter))) > ncol(m_filter)* 1/1000


sum(fD$keep)  ## 11707

# Save Data
stopifnot(identical(as.character(pD$barcode),colnames(m)))
stopifnot(identical(as.character(fD$id),rownames(m)))
out <- list()
out[["counts"]] <- m
out[["phenoData"]] <- pD
out[["featureData"]] <- fD
saveRDS(out,file=file.path(out_objects_dir, "ExpressionList_QC.rds"))