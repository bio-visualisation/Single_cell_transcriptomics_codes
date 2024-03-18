# Code for GSEA tree plot
library(clusterProfiler) #clusterProfiler_4.10.0
library(org.Hs.eg.db) #org.Mm.eg.db_3.18.0
library(enrichplot) #enrichplot_1.22.0
library(tidyverse) #tidyverse_2.0.0
library(msigdbr) #msigdbr_7.5.1
library(paletteer) #paletteer_1.6.0

# Download GO molecular function gene set from MSigDB 
collection <- msigdbr_collections()
msigdbr_show_species()
m_ont <- msigdbr(species = "Homo sapiens", category = "C5", 
                 subcategory = "GO:MF") %>% 
  select(gs_name, gene_symbol)
head(m_ont)

# Remove the GOMF string from the gs_name column
m_ont <- m_ont %>%
  mutate_at("gs_name", str_replace, "GOMF_", "")
head(m_ont)

# Make an empty list
plots <- list()
GSEA <- list()

for (df in names(DEG)) {
  tryCatch({
  data = DEG[[df]]
  message("Doing analysis for ", df)
  # Gene Set Enrichment Analysis (GSEA)
  gene.list <- data %>%
    dplyr::mutate(Score = p_val_adj * sign(avg_log2FC))%>%
    tibble::rownames_to_column(var = "gene")%>%
    dplyr::select(gene, Score)
  # Make the rank file
  ranks <- deframe(gene.list)
  head(ranks)
  # Set decreasing order
  geneList = sort(ranks, decreasing = TRUE)
  #===============================================================
  # Perform GSEA
  #===============================================================
  em2 <- GSEA(geneList, TERM2GENE = m_ont)
  head(em2)
  gsea_result <- em2@result
  
  GSEA[[df]] <- gsea_result
  # Make a Tree Plot
  edox2 <- pairwise_termsim(em2)
  p1 <- treeplot(edox2, showCategory = 10,
                 hclust_method = "average",  label_format = 30,
                 color = "NES", nWords = 0,
                 hilight = TRUE, hextend = 1, nCluster = 5,
                 offset = 30)+
    scale_color_paletteer_c("grDevices::Reds 3", "NES",direction = -1)+
    ggtitle(df)
  
  p1$layers[[7]]$aes_params$size <- 2
  # Store the plots in list
  plots[[df]] <- p1
  }, 
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) 
}

saveRDS(plots, file = "Treeplots.rds")
saveRDS(GSEA, file = "GSEA_MF_human.rds")
