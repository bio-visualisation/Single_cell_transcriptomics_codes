library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(msigdbr)
library(fgsea)
library(org.Mm.eg.db)

# Load diff expressed gene Anterior cerebellum

Ant <- read.csv("Ant_WT_vs_KO.csv",
                 header = T, sep = "\t")%>%
  dplyr::select(Gene_name, logFC, adj.P.Val)

gene.list <- Ant %>%
  arrange(adj.P.Val)%>%
  dplyr::mutate(Score = -log10(adj.P.Val + 0.0001)* sign(logFC))%>%
  dplyr::select(Gene_name, Score)

# Make the rank file
ranks <- deframe(gene.list)
head(ranks)

# ##we need to have the annotated gene set first. 
# 
# m_df <- msigdbr(species = "Mus musculus", category = "C5")
# 
# ##filtering and formating gene sets
# fgsea_sets <- m_df %>% split(x = .$gene_symbol, 
#                              f = .$gs_name)
# ###Perform GSEA
# fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
# 
# ###data wrangling
# fgseaResTidy <- fgseaRes %>%
#   as_tibble() %>%
#   arrange(desc(NES))
# 
# head(fgseaResTidy)
# 
# ###Write result as table
# fwrite(fgseaResTidy, file="Plots/Neuroblast_GSEA_identity.csv",
#        sep="\t", sep2=c("", " ", ""))
# 
# #Save the file
# png(filename = "Plots/Neuroblast_signature.png", 
#     height = 3, width = 6,
#     units = "in", res = 300)
# plotEnrichment(fgsea_sets[["FAN_EMBRYONIC_CTX_NSC_2"]],
#                ranks) +
#   labs(title="FAN_EMBRYONIC_CTX_NSC_2 signatures")
# dev.off()

# Set decreasing order
geneList = sort(ranks, decreasing = TRUE)
head(geneList)
#===============================================================
# Perform GSEA
#===============================================================
m_ont <- msigdbr(species = "Mus musculus", 
                 category = "C8") %>% 
  dplyr::select(gs_name, gene_symbol)

head(m_ont)

em2 <- GSEA(geneList, TERM2GENE = m_ont)
head(em2)
gsea_result <- em2@result

# Save the GSEA result
write.csv(gsea_result, 
          file = "Ant_GSEA_cell_type.csv", 
          row.names = F)

# Save the GSEA Plot
png(filename = "GOBP_CYTOPLASMIC_TRANSLATION.png", 
    width = 8, height = 6, units = "in",
    res = 600)
gseaplot2(em2, geneSetID = "GOBP_CYTOPLASMIC_TRANSLATION", 
          title = em2$Description["GOBP_CYTOPLASMIC_TRANSLATION"])+
  ggtitle("GOBP_CYTOPLASMIC_TRANSLATION")
dev.off()

#
png(filename = "GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY.png", 
    width = 8, height = 6, units = "in",
    res = 600)
gseaplot2(em2, geneSetID = "GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY", 
          title = em2$Description["GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY"])+
  ggtitle("GOMF_VOLTAGE_GATED_CHANNEL_ACTIVITY")
dev.off()

# Divergent Plot

Ont <- read.csv("Ant_GSEA_Ontology.csv",
                header = T)
cell_type <- read.csv("Ant_GSEA_cell_type.csv",
                      header = T)

Ont_top <- Ont %>%
  clusterProfiler::slice(1,5,7,10,13,15,17,25,160)%>% 
  mutate(gene_set = "Gene Ontology")

cell_top <- cell_type %>%
  clusterProfiler::slice(1,3,6,7, 70, 108, 124)%>%
  mutate(gene_set = "Cell Type Signature")

df <- rbind(Ont_top, cell_top)

p <- ggplot(df, aes(x = reorder(ID, NES), y=NES))+
  geom_bar(stat='identity',  aes(fill=p.adjust), width=.8, color ="black",
           size = .1)+
  coord_flip()+
  theme_bw()+
  facet_grid(gene_set ~.,space="free", scales="free")+
  labs(x = "", y = "Normalized Enrichment Score", fill = "Adjusted P value")+
  paletteer::scale_fill_paletteer_c("grDevices::Cyan-Magenta",
                                    direction = 1)+
  theme(axis.text.x =element_text(size=8, colour = "black"),
        axis.text.y = element_text(size=7),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        legend.title=element_text(size=6),
        legend.text=element_text(size=5),
        legend.key.width=unit(0.25,"cm"),
        legend.key.height=unit(0.25,"cm"),
        legend.position = c(0.7,0.1))+
  theme(strip.text.y = element_text(size = 8))
p

ggsave(filename = "Tlx3_KO_Ant_divergent.png",
       p,
       height = 5,
       width = 6,
       units = "in",
       dpi = 600)
