
library(tidyverse)
library(ggridges)
library(ggh4x)

data <- read.csv("VZ_SIZE.csv", header = T)

# We need density plot

df <- data %>%
  gather(key = "Sample",
         value = "Thickness")

df$Tissue <- ifelse(str_starts(df$Sample, "C"),
                    "Control", "KO")
df$Type <- "Cortex"

plot <- ggplot(df, aes(x = Thickness, y = Type,
                       fill = Tissue,
                       point_color = Tissue))+
  geom_density_ridges(jittered_points = TRUE,
                      position = position_raincloud(width = 0.3, 
                                                    height = 0.15),
                      alpha = 0.6, scale = 0.25, point_size = 1,
                      linetype = "blank")+
  scale_fill_manual(values = c("lawngreen", "purple"))+
  scale_discrete_manual("point_color",
                        values = c("lawngreen", "purple"))+
  theme_classic()+ xlim(0,3)+
  guides(x = "axis_truncated", y = "axis_truncated")+
  xlab("Thickness (um)")+ ylab("")+ labs(fill = "",
                                         point_color = "",
                                         color = "")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

plot
png("VZ_thickness.png",
    height = 4, width = 7, units = "in", res = 600)
plot
dev.off()

wilcox.test(df$Thickness ~ df$Tissue)
