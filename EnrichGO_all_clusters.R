# Code for GO tree plot
library(clusterProfiler) #clusterProfiler_4.10.0
library(org.Mm.eg.db) #org.Mm.eg.db_3.18.0
library(enrichplot) #enrichplot_1.22.0
library(tidyverse) #tidyverse_2.0.0
library(msigdbr) #msigdbr_7.5.1
library(paletteer) #paletteer_1.6.0

setwd("C:/Users/bbasu/Documents/Mouse_electrical_stim")
#Import DEGs of all clusters
DEG <- readRDS("DEG_subcluster_function_list.rds")

# Make an empty list
plots <- list()
GSEA <- list()


for (df in names(DEG)) {
  tryCatch({
    data = DEG[[df]]
    message("Doing analysis for ", df)
    # Gene Set Enrichment Analysis (GSEA)
    gene.list <- data %>%
      dplyr::filter(abs(avg_log2FC) > 0.20 & p_val_adj < 0.05)%>%
      tibble::rownames_to_column(var = "gene")
    #===============================================================
    # Perform GO enrichment
    #===============================================================
    ego <- enrichGO(gene          = gene.list$gene,
                    OrgDb         = org.Mm.eg.db, # or Org.Hs.eg.db
                    ont           = "MF", 
                    #one of “BP”, “MF”, “CC” or “ALL”
                    pAdjustMethod = "fdr", 
                    #one of “bonferroni”, “BH”, “BY”, “fdr”, “none”
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    keyType = "SYMBOL", 
                    #“ENSEMBL”, “ENTREZID”, “SYMBOL”
                    readable      = TRUE)

    GSEA[[df]] <- ego
    # Make a Tree Plot
    edox2 <- pairwise_termsim(ego)
    
    p1 <- treeplot(edox2, showCategory = 10,
                   hclust_method = "average",  label_format = 30,
                   color = "p.adjust", nWords = 0,
                   hilight = TRUE, hextend = 1, nCluster = 5,
                   offset = 30)+
      scale_color_paletteer_c("grDevices::Reds 3","p.adjust",direction = 1)+
      ggtitle(df)
    
    p1$layers[[7]]$aes_params$size <- 2
    # Store the plots in list
    plots[[df]] <- p1
  }, 
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}

saveRDS(plots, file = "Treeplots_mouse.rds")
saveRDS(GSEA, file = "GO_MF_mouse_all_cluster.rds")