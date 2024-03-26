library(tidyverse)
library(ggrepel)
df1 <- read.table(file = "HSD17B8_TARGET_GENES.tsv", 
                  sep="\t", header = TRUE)

df1$SYMBOL <- stringr::str_to_title(df1$SYMBOL)

genes <- c("Mki67", "Top2a", "Aurkb", "Foxm1")


plot <- ggplot(data=df1, aes(x=RANK.IN.GENE.LIST, 
                     y= RUNNING.ES)) +
  geom_line(color = "green")+
  geom_point(data = df1 %>%
               filter(SYMBOL %in% genes),
             color = "red") +
  geom_rug(col="black",alpha=0.5, size=0.3,
           sides = 'b')+
  theme_minimal()+
  ggrepel::geom_text_repel(data = df1 %>%
                             filter(SYMBOL %in% genes),
                           aes(label = SYMBOL),
                           box.padding = unit(1, "lines"),
                           hjust= 0.30,
                           segment.color = 'black',
                           colour = 'black')+
  labs(title = "Hsd17b8 Target Genes",
       x = "Ranked gene list",
       y = "Enrichment Score",
       color = "")+
  theme(axis.text.x = element_blank())

plot

png("hsd17b8_target_genes.png",
    height = 3, width = 6, units = "in",
    res = 600)
plot
dev.off()



