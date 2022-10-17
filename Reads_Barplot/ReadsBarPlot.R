############################ Reads Stacked Bar Plot ##########################################
setwd("~/Google Drive/My Drive/SF_Bay/SanFranciscoBay/Reads_Barplot/")

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

############################ Coverm Input Data ##########################################

data2 <- read.delim("input/SFBMAGS-vs-Reads-CoverM.tsv")

#rename sites
colnames(data2) = gsub("_filter.METAGENOME.fastq.gz.Relative.Abundance....", "", colnames(data2))
#colnames(data2) = gsub("_sed_USGS", "", colnames(data2))
#colnames(data2) = gsub("SF_", "", colnames(data2))

#transform from wide to long
data2_long <- melt(data2, id = c("Genome")) 

#rename variable column to Site
data2_long<-rename(data2_long, "Site" = "variable")

#clean up Site names
data2_long <- data2_long %>% separate(Site, c("Month", "Site"), 
                          sep= "_USGS_")
data2_long$Month<-gsub("SF_*", "", data2_long$Month)
data2_long$Month<-gsub("12|11", "", data2_long$Month)
data2_long$Month<-gsub("Jul", "July", data2_long$Month)
data2_long$Month<-gsub("_sed", "", data2_long$Month)

data2_long$Site <- paste0(data2_long$Site, "_", data2_long$Month)

#set order of x axis 
data2_long$Site <- factor(data2_long$Site, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
                                        "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                        "13_July", "13_Oct", "13_Jan", "13_May",
                                        "21_July", "21_Oct", "21_Jan", "21_May",
                                        "24_July", "24_Oct", "24_Jan", "24_May")) 

#subset only unmapped
data2_long_plotting<-data2_long[data2_long$Genome=="unmapped",]

#get rid of unmapped so can sum by site for percent mapped
data2_long_mapped<-data2_long[!grepl("unmapped", data2_long$Genome),]

#sum values for percent mapped
data2_long_mapped <- data2_long_mapped %>%
  group_by(Site) %>% 
  summarise(value=sum(value)) #summing the group

#combine mapped and unmapped into 1 data frame
data2_long_plotting<-rename(data2_long_plotting, "Unmapped Reads" = "value")
data2_long_plotting<-select(data2_long_plotting, -Genome)

data2_long_mapped<-rename(data2_long_mapped, "Mapped Reads" = "value")

##vlookup
data2_long_plotting <- data2_long_mapped %>%
  dplyr::select("Site", "Mapped Reads") %>%
  right_join(data2_long_plotting, by = c("Site" = "Site"))

#melt the data
data2_melt<-data2_long_plotting[, c("Site", "Mapped Reads", "Unmapped Reads")]
data2_melt <- melt(data2_melt, id = "Site") 

#add site without month column as option for faceting
data2_melt$Site2 = data2_melt$Site
data2_melt <- data2_melt %>% separate(Site2, c("SiteOnly", "Month"), sep= "(?=_[A-Z])")

#set order of x axis 
data2_melt$SiteOnly <- factor(data2_melt$SiteOnly, levels=c("4_1", "8_1", "13", "21", "24")) 

#plot
dev.off()
#plot
p<-ggplot(data2_melt, aes(fill=forcats::fct_rev(variable), y=value, x=Site, label = round(value, digits = 2))) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#9ecae1", "#3182bd")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA),
        axis.ticks = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Percent Reads",
      title = "CoverM Percent of Reads Mapped to SFB MAGs") +
  guides(fill=guide_legend(title="")) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  facet_grid(.~SiteOnly, scales = "free") 
  
p

ggsave("output/CoverM_ReadsMapped_Barplot.png", width = 10)
ggsave("output/CoverM_ReadsMapped_Barplot_Facet.png", width = 10)



