############################ Reads Stacked Bar Plot ##########################################

library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)

############################ Coverm Input Data ##########################################

data <- read.delim("input/SFBMAGS-vs-Reads-CoverM.tsv")

#rename sites
colnames(data) = gsub("_filter.METAGENOME.fastq.gz.Relative.Abundance....", "", colnames(data))
#colnames(data2) = gsub("_sed_USGS", "", colnames(data2))
#colnames(data2) = gsub("SF_", "", colnames(data2))

#transform from wide to long
data_long <- melt(data, id = c("Genome")) 

#rename variable column to Site
data_long<-rename(data_long, "Site" = "variable")

#clean up Site names
data_long <- data_long %>% separate(Site, c("Month", "Site"), 
                          sep= "_USGS_")
data_long$Month<-gsub("SF_*", "", data_long$Month)
data_long$Month<-gsub("12|11", "", data_long$Month)
data_long$Month<-gsub("Jul", "July", data_long$Month)
data_long$Month<-gsub("_sed", "", data_long$Month)

data_long$Site <- paste0(data_long$Site, "_", data_long$Month)

#set order of x axis 
data_long$Site <- factor(data_long$Site, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
                                        "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                        "13_July", "13_Oct", "13_Jan", "13_May",
                                        "21_July", "21_Oct", "21_Jan", "21_May",
                                        "24_July", "24_Oct", "24_Jan", "24_May")) 

#subset only unmapped
data_long_plotting<-data_long[data_long$Genome=="unmapped",]

#get rid of unmapped so can sum by site for percent mapped
data_long_mapped<-data_long[!grepl("unmapped", data_long$Genome),]

#sum values for percent mapped
data_long_mapped <- data_long_mapped %>%
  group_by(Site) %>% 
  summarise(value=sum(value)) #summing the group

#combine mapped and unmapped into 1 data frame
data_long_plotting<-rename(data_long_plotting, "Unmapped Reads" = "value")
data_long_plotting<-select(data_long_plotting, -Genome)

data_long_mapped<-rename(data_long_mapped, "Mapped Reads" = "value")

##vlookup
data_long_plotting <- data_long_mapped %>%
  dplyr::select("Site", "Mapped Reads") %>%
  right_join(data_long_plotting, by = c("Site" = "Site"))

#melt the data
data_melt<-data_long_plotting[, c("Site", "Mapped Reads", "Unmapped Reads")]
data_melt <- melt(data_melt, id = "Site") 

#add site without month column as option for faceting
data_melt$Site2 = data_melt$Site
data_melt <- data_melt %>% separate(Site2, c("SiteOnly", "Month"), sep= "(?=_[A-Z])")

#set order of x axis 
data_melt$SiteOnly <- factor(data_melt$SiteOnly, levels=c("4_1", "8_1", "13", "21", "24")) 

#plot
dev.off()
#plot
p<-ggplot(data_melt, aes(fill=forcats::fct_rev(variable), y=value, x=Site, label = round(value, digits = 2))) + 
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



