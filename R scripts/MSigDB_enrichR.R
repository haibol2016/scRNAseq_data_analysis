library(msigdbr)
library("clusterProfiler")
msigdbr_species()
all_gene_sets = msigdbr(species = "Mus musculus")
all_gene_sets <- as.data.frame(unique(all_gene_sets[, 1:2]))
all_gene_sets$set_name <- gsub(":", "_", paste0(all_gene_sets[,1], all_gene_sets[,2]))

out <- "MsigDB.out"
if (!dir.exists(out)){
    dir.create(out)
}
deg_files <- dir("./", ".control$")
names(deg_files) <- gsub("(.+?).DE.genes.Mutant.versus.control", 
                         "\\1", deg_files, perl = TRUE)

id_map <- read.delim("../GRCm39.ensemblID2EntrezID.txt",
                     header = TRUE, as.is = TRUE)

#expressed_genes <- read.delim("../Mesemchymal.expressed.genes.symbol.id.txt",
#                              header = TRUE, as.is = TRUE)

expressed_genes <- read.delim("../Epithelium.expressed.genes.symbol.id.txt",
                              header = TRUE, as.is = TRUE)

expressed_genes <- merge(expressed_genes, id_map,
                         by = "EnsemblID", all.x = TRUE)
expressed_genes <- expressed_genes[!is.na(expressed_genes$EntrezgeneID), ]

degs <- lapply(deg_files, function(f){
    dat <- read.delim(f, header = TRUE, as.is = TRUE)
    dat <- dat[!grepl("(bGal|GFP)", dat$Symbol, perl = TRUE), ]
    dat <- merge(dat, expressed_genes,
                 by = "EnsemblID", 
                 all.x = TRUE)
    dat <- dat[dat$p_val_adj <= 0.05, ]
    up <- dat[dat$avg_log2FC >= 0.25, ]
    down <- dat[dat$avg_log2FC <= -0.25, ]
    list(up = up, down = down, all = dat)
})

for (c in seq_along(all_gene_sets$gs_cat)){
    if (all_gene_sets$gs_subcat[c] != "") {
        m_t2g <- msigdbr(species = "Mus musculus", 
                         category = all_gene_sets$gs_cat[c], 
                         subcategory = all_gene_sets$gs_subcat[c]) %>% 
            dplyr::select(gs_name, entrez_gene)
        
    } else {
        m_t2g <- msigdbr(species = "Mus musculus", category =  all_gene_sets$gs_cat[c]) %>% 
            dplyr::select(gs_name, entrez_gene) 
    }

    for (d in seq_along(degs)) {
        for (direction in c("up", "down", "all")) {
            em <- enricher(gene = as.character(degs[[d]][[direction]]$EntrezgeneID), 
                           pvalueCutoff = 0.05, 
                           pAdjustMethod = "BH", 
                           universe = as.character(expressed_genes$EntrezgeneID), 
                           minGSSize = 10, 
                           maxGSSize = 500, 
                           qvalueCutoff = 0.1,
                           TERM2GENE=m_t2g)
            if (!is.null(em) && nrow(em) >= 1){
                write.table(em, file = file.path(out, paste0(names(degs)[d], "-", 
                                                             all_gene_sets$set_name[c], "-",
                                             paste0(direction,
                                                    "-regulated.txt"))),
                            sep = "\t", quote = FALSE, row.names = FALSE)
                tryCatch({pdf(file.path(out, paste0(names(degs)[d], "-", 
                                                    all_gene_sets$set_name[c], "-",
                                     paste0(direction,
                                            "-regulated.pdf"))),
                              width = 8, height = {if (nrow(em)< 10){5} else{min(nrow(em)*0.5, 10)}})
                    print(dotplot(em, showCategory = 20))},
                    error=function(cond) {return(NA)},
                    finally= {dev.off()})
            }
        }
    }
}




## enrichR
install.packages("enrichR")

library(enrichR)
library(cowplot)
library(patchwork)
library(ggplot2)
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()


#[1] "TRRUST_Transcription_Factors_2019"        
#[2] "ChEA_2016"                                
#[3] "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"
#[4] "ENCODE_TF_ChIP-seq_2015"                  
#[5] "Genome_Browser_PWMs"                      
#[6] "TRANSFAC_and_JASPAR_PWMs"
#[7] "TF-LOF_Expression_from_GEO"

dbs <- dbs$libraryName[c(147, 104, 96, 50, 1, 2, 12)]


deg_files <- dir("./", ".control$")
names(deg_files) <- gsub("(.+?).DE.genes.Mutant.versus.control", 
                         "\\1", deg_files, perl = TRUE)


up_dow_genes <- mapply(function(.x, .y){
    dat <- read.delim(.x, header = TRUE, as.is = TRUE)
    dat <- dat[!grepl("(bGal|GFP)", dat$Symbol, perl = TRUE), ]
    dat <- dat[dat$p_val_adj <= 0.05, ]
    up <- dat[dat$avg_log2FC >= 0.25, "Symbol"]
    down <- dat[dat$avg_log2FC <= -0.25, "Symbol"]
    gene_list <- list(up = up, down = down, all = dat[, "Symbol"])
    
    null <- mapply(function(x, y) {
        
        enriched <- enrichr(x, dbs)
        
        out <- "Enrichr.out"
        sub_out <- file.path(out, paste0(.y, y, "-genes.TF.target.enrichment.out"))
        if (!dir.exists(sub_out)) {
            dir.create(sub_out, recursive = TRUE)
        }
        null <-  mapply(function(i, j){
            write.table(i, file = file.path(sub_out, paste0(j, "_table.txt")),
                        sep = "\t", quote = FALSE, row.names = FALSE)
        }, enriched, names(enriched), SIMPLIFY = FALSE)
        
    }, gene_list, names(gene_list), SIMPLIFY = FALSE)
}, deg_files, names(deg_files), SIMPLIFY = FALSE)


## visualize the EnrichR result
folders <- dir("Enrichr.out", "target.enrichment.out$", full.names = TRUE)
out <- file.path("Enrichr.out", "visualization_top10")
if (!dir.exists(out)){
    dir.create(out)
}

for (f in folders){
    
    files <- dir(f, full.names = TRUE)
    out_names <- gsub(".+/(.+?)-genes.TF.target.enrichment.out", "\\1", f, perl = TRUE)
    set_names <- gsub(".+/(.+?)_table.txt", "\\1", files, perl = TRUE)
    tf <- mapply(function(.x, .y){
        tf_enrichment <- read.delim(.x, header = TRUE, as.is = TRUE)[, c(1, 4, 8)]
        
        ## only significant ones
        tf_enrichment <- tf_enrichment[tf_enrichment$Adjusted.P.value < 0.2, ]
        tf_enrichment <- tf_enrichment[order(tf_enrichment$Combined.Score, 
                                             decreasing = TRUE), ]
        tf_enrichment <- tf_enrichment[1:min(10, nrow(tf_enrichment)), ]
        if (nrow(tf_enrichment) > 0) {
            tf_enrichment$set <- .y
        }
        
        tf_enrichment <- tf_enrichment[order(tf_enrichment$Combined.Score), ]
        tf_enrichment$Term <- factor(tf_enrichment$Term, 
                                     levels = tf_enrichment$Term)
        tf_enrichment
    }, files, set_names, SIMPLIFY = FALSE, USE.NAMES = TRUE)
    
    pdf(file.path(out, paste0(out_names, "enriched.transcription.factors.pdf")), 
        height = 8, width = 11)
    p <- lapply(tf, function(x) {
        ggplot(x, aes(x = Term, y = Combined.Score, fill =set)) + 
            geom_bar(stat = "identity",  width= 0.5, 
                     position = position_dodge(width=0.1)) + 
            coord_flip() + 
            ylab("Combined Score")  +
            xlab("Term") + theme_cowplot(font_size = 8) +
            theme(legend.position="top", 
                  axis.text = element_text(size = 6)) 
    })
    
    print(wrap_plots(p, nrow = 4, ncol = 2))
    dev.off()
}
