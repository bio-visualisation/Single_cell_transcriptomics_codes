library(tidyverse)
# Importing dataset
df <- read.csv(file = "data.csv", header = T)

df$Tlx3_chip <- str_to_title(df$Tlx3_chip)
# data for overlapping
Tlx3_chip <- df %>%
  select(Tlx3_chip)%>%
  deframe()
Gn <- df %>%
  select(scRNA_GN)%>%
  deframe()
Gnp <- df %>%
  select(scRNA_GNP)%>%
  deframe()

# Plot
library(VennDiagram)
x <- list(A = Tlx3_chip, B = Gn, C = Gnp)

overlap <- calculate.overlap(x)
library(data.table)
overlap_list <- data.frame(overlap)

# Function display_venn
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

venn.diagram(x, lwd = 0.8, resolution = 600, 
             fill = c("#d7bde2", "#1abc9c", "deeppink"),
             lty = 'blank',
             filename = "Tlx3_chip_venn.png",
             category.names = c("Tlx3 ChIP",
                                "DEG in GN",
                                "DEG in GNP"),
             cat.cex = 1.5,
             cat.default.pos = 'outer',
             cat.pos = c(-20, 20, 180),
             scaled = T,
             cat.col = c("#6c3483", 
                         "#0b5345", 
                         "deeppink"))

#cat.pos = 0 at 12'O clock