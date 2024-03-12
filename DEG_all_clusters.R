# Code for calculating DEGs for all clusters 
# Make volcano plots for all clusters
#===============================================================================

library(Seurat)
library(Signac)
library(tidyverse)
library(scCustomize)
library(ggrepel)
library(ggthemes)
library(RColorBrewer)

merged <- readRDS("merged_stimulation_subcluster.rds")

DEG <- list()
plots <- list()

color <- brewer.pal(12, "Paired")
color

Idents(merged) <- merged$celltype

# Vector of colors
library(RColorBrewer)
color <- brewer.pal(12, "Paired")
color


colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", 
            "#33A02C", "#FB9A99", "#E31A1C", 
            "#FDBF6F", "#FF7F00",
            "#CAB2D6")


for(i in sort(levels(merged))){
  cluster = i
  message("Calculating DEG for ",cluster)
  cluster_index = which(sort(levels(merged)) == cluster)
  cluster.marker <- FindMarkers(merged, ident.1 = "Stimulated", 
                                ident.2 = "Baseline", 
                                subset.ident = cluster, 
                                group.by = "condition",
                                verbose = FALSE,
                                test.use = "wilcox")
  DEG[[cluster]] <- cluster.marker
  #Data wrangling
  data <- data.frame(gene = row.names(cluster.marker),
                     pval = sign(cluster.marker$avg_log2FC)*(-log10(cluster.marker$p_val_adj)), 
                     lfc = cluster.marker$avg_log2FC)
  
  
  data <- mutate(data, color = case_when(data$lfc > 0 & abs(data$pval) > 1.3 ~ "Increased",
                                         data$lfc < 0 & abs(data$pval) > 1.3 ~ "Decreased",
                                         data$pval < 1.3 ~ "nonsignificant"))
  
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(data, aes(x = pval, y = lfc, color = color))
  
  
  # Add ggplot2 layers
  p1 <- vol + 
    geom_point(size = 1, na.rm = T) +
    scale_color_manual(name = "Directionality",
                       values = c(Increased = colors[cluster_index], 
                                  Decreased = colors[cluster_index], 
                                  nonsignificant = "gray90")) +
    theme_base() + # change overall theme
    theme(legend.position = "none") + # change the legend
    ylab("") + # Change Y-Axis label
    xlab("")+  # Change X-Axis label
    ylim(-10, 11)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
  p2 <- p1 + geom_text_repel(data= data %>%
                               filter(abs(data$pval) > 1.3)%>%
                               head(9), aes(label=gene),
                             box.padding = unit(.5, "lines"),hjust= 0.30,
                             segment.color = 'black', max.overlaps = Inf,
                             colour = 'black', size = 3)+
    ggtitle(cluster)
  
  plots[[cluster]] <- p2
  
}

saveRDS(DEG, file = "DEG_function_list.rds")
saveRDS(plots, file = "Plots_function_list.rds")

DEG <- readRDS("DEG_function_list.rds")
plots <- readRDS("Plots_function_list.rds")