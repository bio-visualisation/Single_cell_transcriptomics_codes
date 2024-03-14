# Code for cnet plot
library(clusterProfiler) #clusterProfiler_4.10.0
library(org.Hs.eg.db) #org.Mm.eg.db_3.18.0
library(enrichplot) #enrichplot_1.22.0
library(tidyverse) #tidyverse_2.0.0
library(msigdbr) #msigdbr_7.5.1
library(paletteer) #paletteer_1.6.0

# Import the data
res <- read.csv(file = "Human_Bulk_Stimulated.csv", 
                header = T)
data <- data.frame(gene = res$gene,
                   pval = res$FDR, 
                   lfc = res$logFC)

genelist <- data %>% dplyr::filter(pval < 0.05)

ego <- enrichGO(gene          = genelist$gene,
                OrgDb         = org.Hs.eg.db, # or Org.Hs.eg.db
                ont           = "MF", 
                #one of “BP”, “MF”, “CC” or “ALL”
                pAdjustMethod = "fdr", 
                #one of “bonferroni”, “BH”, “BY”, “fdr”, “none”
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                keyType = "SYMBOL", 
                #“ENSEMBL”, “ENTREZID”, “SYMBOL”
                readable      = TRUE)

result <- ego@result

# Make CNET plot
fc <- data%>%
  select(gene, lfc)%>%
  deframe()
head(fc)



p1 <- cnetplot(ego, showCategory = 10, circular = TRUE, colorEdge = TRUE, foldChange = fc)+
  paletteer::scale_color_paletteer_c("ggthemes::Red-Blue Diverging", direction = 1)+
  labs(color = "log2FC")

p1

pdf(file = "Human_CNET_plot_logfc.pdf", height =12,
    width = 17)
p1
dev.off()

