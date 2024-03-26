library(tidyverse)
library(paletteer)
library(ComplexHeatmap)

# Load the data
data <- read.csv("Cerebellum_genes_final_normalized_scores.csv",
                 header = T)

head(data)

data <- data %>%
  separate(ID, into = c("ENS_ID", "Gene"), 
           sep = "_", remove = TRUE)

# Load diff expressed gene Anterior cerebellum

Ant <- read.csv("Ant_WT_vs_KO.csv",
                 header = T, sep = "\t")%>%
  select(Gene_name, logFC, adj.P.Val)

de_genes <- Ant %>%
  filter(adj.P.Val < 0.05)

# Filter DE genes for heatmap

df <- data %>%
  filter(data$Gene %in% de_genes$Gene_name)%>%
  select(Gene, Ant_Cere_WT_N1, Ant_Cere_WT_N2,
         Ant_Cere_KO_N1, Ant_Cere_KO_N2,
         Post_Cere_WT_N1, Post_Cere_WT_N2,
         Post_Cere_KO_N1, Post_Cere_KO_N2)%>%
  column_to_rownames(var = "Gene")



mat <- as.matrix(df[,c(1:4)])
# Scale the expression value
scaled_mat <- t(scale(t(mat)))

my_palette <- paletteer::paletteer_c("grDevices::Cyan-Magenta", 30)

ha = HeatmapAnnotation(df = data.frame(Samples = c("Control","Control",
                                                   "cKO","cKO")),
                       col = list(Samples = c("Control" = "khaki3",
                                              "cKO" = "darkmagenta")),
                       annotation_legend_param = list(title = "Ant. Cerebellum",
                                                      title_gp = gpar(fontsize = 10,
                                                                      fontface = "bold")),
                       annotation_name_side = "left",
                       simple_anno_size = unit(2, "mm"))

ht1 <- Heatmap(scaled_mat,cluster_columns = T,
               width = unit(3, "cm"), 
               height = unit(8, "cm"),
               show_column_names = F,
               show_row_names = F,
               col = my_palette,
               bottom_annotation = ha,
               name = "Row Z score",
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10,
                                 fontface = "bold")
               ))

ht1

ann <- df %>%
  filter(rownames(df) %in% c("Rpl38", "Rpl36a",
                             "Kcna1", "Lrrc26",
                             "Cacna2d1", "Ano3",
                             "Hmgb1", "Hmgb2", "Stmn1",
                             "Cenpe", "Birc5", "Pcsk9"))
vrn <- rownames(df) %in% rownames(ann)

ht2 <- ht1 + 
  rowAnnotation(link = anno_mark(at = which(vrn), 
                                 labels = row.names(scaled_mat)[vrn],
                                 labels_gp = gpar(fontsize = 10), 
                                 padding = unit(1, "mm")))
ht2

png("Anterior/Ant.Cerebellum_heatmap.png",
    height = 5, width = 6, units = "in",
    res = 600)
ht2
dev.off()