############################ Reads Stacked Bar Plot ##########################################
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)

data <- read.delim("input/Reads_inMAGs_totals.txt")

#rename sites
data <- data %>% separate(Site, c("Month", "Site"), 
                                          sep= "_USGS_")
data$Month<-gsub("SF_*", "", data$Month)
data$Month<-gsub("12|11", "", data$Month)
data$Month<-gsub("Jul", "July", data$Month)

data$Site <- paste0(data$Site, "_", data$Month)

#set order of x axis 
data$Site <- factor(data$Site, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
                                                  "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                                  "13_July", "13_Oct", "13_Jan", "13_May",
                                                  "21_July", "21_Oct", "21_Jan", "21_May",
                                                  "24_July", "24_Oct", "24_Jan", "24_May")) 
  
#melt the data
data_melt<-data[, c("Site", "NumReads", "NumMAGReads")]
data_melt <- melt(data_melt, id = "Site") 

dev.off()
#plot
p<-ggplot(data_melt, aes(fill=variable, y=value, x=Site, label = value)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#bababa", "#d6604d")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5),
        #panel.background = element_blank(),
        axis.ticks = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_text(size = 1.5, position = position_stack(vjust = 0.5))
p

p + annotate(geom="text", x="4_1_July", y=30000000, label="19.1%", size = 2,
               color="black") +
  annotate(geom="text", x="4_1_Oct", y=37000000, label="22.3%", size = 2,
             color="black") +
  annotate(geom="text", x="4_1_Jan", y=40000000, label="25.1%", size = 2,
             color="black") +
  annotate(geom="text", x="4_1_May", y=38000000, label="18.9%", size = 2,
             color="black") +
  annotate(geom="text", x="8_1_July", y=30000000, label="18%", size = 2,
             color="black") +
  annotate(geom="text", x="8_1_Oct", y=38000000, label="20.7%", size = 2,
             color="black") +
  annotate(geom="text", x="8_1_Jan", y=38000000, label="22.5%", size = 2,
             color="black") +
  annotate(geom="text", x="8_1_May", y=10000000, label="2.3%", size = 2,
             color="black") +
  annotate(geom="text", x="13_July", y=24000000, label="14.8%", size = 2,
             color="black") +
  annotate(geom="text", x="13_Oct", y=14000000, label="5.7%", size = 2,
             color="black") +
  annotate(geom="text", x="13_Jan", y=19000000, label="7.8%", size = 2,
             color="black") +
  annotate(geom="text", x="13_May", y=43000000, label="24.6%", size = 2,
             color="black") +
  annotate(geom="text", x="21_July", y=13000000, label="5.3%", size = 2,
             color="black") +
  annotate(geom="text", x="21_Oct", y=15000000, label="7%", size = 2,
             color="black") +
  annotate(geom="text", x="21_Jan", y=48000000, label="19.8%", size = 2,
             color="black") +
  annotate(geom="text", x="21_May", y=21000000, label="9.2%", size = 2,
             color="black") +
  annotate(geom="text", x="24_July", y=11000000, label="4.1%", size = 2,
             color="black") +
  annotate(geom="text", x="24_Oct", y=21000000, label="8.6%", size = 2,
             color="black") +
  annotate(geom="text", x="24_Jan", y=15000000, label="6%", size = 2,
             color="black") +
  annotate(geom="text", x="24_May", y=12000000, label="4.6%", size = 2,
             color="black")

#example input data
# specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
# condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
# value <- abs(rnorm(12 , 0 , 15))
# test <- data.frame(specie,condition,value)





