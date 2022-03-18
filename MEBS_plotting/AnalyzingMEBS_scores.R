library(reshape2)
library(ggplot2)
library(viridis)
library(tidyverse)

data <- read.delim(file = "../sfbArcBac_mebs_Norm_mod.txt", header = TRUE, sep = "\t")

#Subset the data for just S scores
data_S <- data[1:5]

#Filter the data so just S scores at site 4.1
data_S_4_1 <- data_S%>%filter(str_detect(Site,"4_1"))
#Filter that so just month of July
S_4_1_july <- data_S_4_1[grepl("July", data_S_4_1$Genome),]
