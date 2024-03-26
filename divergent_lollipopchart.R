library(tidyverse)

ko_TF <- read.csv("TF_pos.tsv", 
                  sep = "\t", header = T) |> 
  select(NAME, NES, FDR.q.val) |> 
  slice(c(1:2), )###Slicing different rows

cntrl_TF <- read.csv("TF_neg.tsv", 
                     sep = "\t", header = T) |> 
  select(NAME, NES, FDR.q.val) |> 
  slice(c(1,4,6,7,8,9,29), )


TF_merge <- rbind(cntrl_TF, ko_TF)

my_data <- TF_merge |>
  mutate(NAME = fct_reorder(NAME, NES))

my_data$NES <- round(my_data$NES, 3) #Round up NES values upto 3 digits

head(my_data)

#Divergent lollipop chart
p <- ggplot(my_data,
            aes(x = NAME,
                y = NES))+
  geom_segment(aes(y = 0,
                   x = NAME,
                   xend = NAME,
                   yend = NES),
               color = "skyblue")+
  geom_point(stat='identity', 
             aes(size = abs(NES)), 
             col = ifelse(my_data$NES > 0, "deeppink", "blue"))+
  coord_flip()+
  theme_light()+
  labs(title = "Comparison of NIHes1 KO vs Control",
       subtitle = "Gene Set Enrichment Analysis",
       size = "NES",
       x = "")+
  theme(axis.text.y = element_text(color = "black",
                                   size = 12))
p


#Save the file
png("NIHes1_KO_GSEA_chart.png", 
     height = 4, width = 6,
     units = "in", res = 600)
p
dev.off()