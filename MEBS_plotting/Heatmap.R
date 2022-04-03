# Heatmap script adapted from:
# Heatmap of mine pit microbe diversity
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2017-06-05

#rm(list = ls())

# Load dependancies
library("tidyr")
library("ggplot2")
library("data.table")
library("dplyr")
library(stringr)

# Read data and format for heatmap
mine.data <- read.delim(file = "Input/sulfur_input2.txt",
                        stringsAsFactors = FALSE)
#mine.long <- data_trans<-transpose(mine.data)
mine.long <- pivot_longer(data = mine.data,
                          cols = -c(1:4),
                          names_to = "Genes", 
                          values_to = "Sum")

# Remove some genes
mine.long<-mine.long %>% 
  filter(!str_detect(Genes, c('dmsB_2|dmsC')))

# Plot MEBS completeness
mine.heatmap <- ggplot(data = mine.long, mapping = aes(x = Month,
                                                       y = Genes,
                                                       fill = Sum)) +
  geom_tile() +
  xlab(label = "Month") +
  # Facet on depth and drop empty columns
  facet_grid(~ Site, switch = "x", scales = "free_x", space = "free_x") + 
  # Set colors different from default
  scale_fill_gradient(name = "Num of Genes", 
                      low = "#FFFFFF",
                      high = "#012345") +
  theme_bw() +
  
  theme(strip.placement = "outside", # Move depth boxes to bottom of plot
        plot.title = element_text(hjust = 0.5), # Center-justify plot title
        axis.title.y = element_blank(), # Remove y-axis title
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  #ggtitle(label = "MEBS completeness") +
  scale_y_discrete(limits = rev(levels(as.factor(mine.long$Genes))))

mine.heatmap

################################

mine.data <- read.delim(file = "../../MEBS/sfb_arcbac_mebs_2_metadata.tsv_completenes.txt",
                        stringsAsFactors = FALSE)
#mine.long <- data_trans<-transpose(mine.data)
mine.long <- pivot_longer(data = mine.data,
                          cols = -c(1:4),
                          names_to = "Genes", 
                          values_to = "Completeness")

# Transform abundance data for better visualization
#mine.long$Sqrt.abundance <- sqrt(mine.long$Abundance)

#drop some pathways
mine.long<-mine.long %>% 
  filter(!str_detect(Genes, c('aminobutanoate|allantoin|Caffeine|mcrABC|Urea|GABA|coB.coM|L.glutamine|Methanogenesis')))

# Plot MEBS completeness
mine.heatmap <- ggplot(data = mine.long, mapping = aes(x = X,
                                                       y = Genes,
                                                       fill = Completeness)) +
  geom_tile() +
  xlab(label = "Site") +
  # Facet on depth and drop empty columns
  facet_grid(~ Site, switch = "x", scales = "free_x", space = "free_x") + 
  # Set colors different from default
  scale_fill_gradient(name = "Completeness", 
                      low = "#FFFFFF",
                      high = "#012345") +
  theme_bw() +
  
  theme(strip.placement = "outside", # Move depth boxes to bottom of plot
        plot.title = element_text(hjust = 0.5), # Center-justify plot title
        axis.title.y = element_blank(), # Remove y-axis title
        axis.text.x = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
  #ggtitle(label = "MEBS completeness") +
  scale_y_discrete(limits = rev(levels(as.factor(mine.long$Genes))))

mine.heatmap

