setwd("~/Google Drive/My Drive/SF_Bay/SanFranciscoBay/Abundance_Barplot/")

############################ Abundance Stacked Bar Plot ##########################################
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)

############################ Read in Inputs #############################

#read in coverm input
coverm<-read.delim2(file = "Input/SFBMAGS-vs-Reads-CoverM.tsv")
gtdbtk<-read.delim2(file = "Input/gtdbtk.r207.bacarc.summary.tsv")

############################ Clean up the Coverm Input #############################
colnames(coverm) = gsub("_filter.METAGENOME.fastq.gz.Relative.Abundance....", "", colnames(coverm))
#transform from wide to long
coverm_long <- melt(coverm, id = c("Genome")) 
#rename variable column to Site
coverm_long<-rename(coverm_long, "Site" = "variable")
#clean up Site names
coverm_long <- coverm_long %>% separate(Site, c("Month", "Site"), 
                                        sep= "_USGS_")
coverm_long$Month<-gsub("SF_*", "", coverm_long$Month)
coverm_long$Month<-gsub("12|11", "", coverm_long$Month)
coverm_long$Month<-gsub("Jul", "July", coverm_long$Month)
coverm_long$Month<-gsub("_sed", "", coverm_long$Month)
#make a site column with full site name
coverm_long$SiteFull <- paste0(coverm_long$Site, "_", coverm_long$Month)
#drop unmapped
coverm_long<-coverm_long[!grepl("unmapped", coverm_long$Genome),]

############################ Add GTDBtk taxonomy #############################
coverm_long <- gtdbtk %>%
  dplyr::select("user_genome", "classification") %>%
  right_join(coverm_long, by = c("user_genome" = "Genome"))
#only keep phylum and class of the tax string
coverm_long <- coverm_long %>% separate(classification, c("d", "p", "c", "o", "f", "g", "s"), 
                            sep= ";")
coverm_long <- coverm_long %>% select(c("user_genome","p","c","Month","Site","value","SiteFull"))
#for Proteobacteria and Thermoproteota keep class, everything else, keep phylum
coverm_long_proteo <- coverm_long %>% filter(grepl("p__Proteobacteria", p))
coverm_long_proteo <- coverm_long_proteo %>% select(-c("p"))
coverm_long_proteo <- coverm_long_proteo %>% rename("Taxa" = "c")

coverm_long_thermo <- coverm_long %>% filter(grepl("p__Thermoproteota", p))
coverm_long_thermo <- coverm_long_thermo %>% select(-c("p"))
coverm_long_thermo <- coverm_long_thermo %>% rename("Taxa" = "c")

coverm_long <- coverm_long %>% filter(!grepl("p__Proteobacteria|p__Thermoproteota", p))
coverm_long <- coverm_long %>% select(-c("c"))
coverm_long <- coverm_long %>% rename("Taxa" = "p")
#put the data frames back together
coverm_long_final <- rbind(coverm_long, coverm_long_proteo, coverm_long_thermo)

############################ Set order of sites for plotting #############################
#set order of x axis 
coverm_long_final$Site <- factor(coverm_long_final$Site, levels=c("4_1", "8_1", "13", "21", "24"))
#set order of x axis 
coverm_long_final$SiteFull <- factor(coverm_long_final$SiteFull, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
                                                              "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                                              "13_July", "13_Oct", "13_Jan", "13_May",
                                                              "21_July", "21_Oct", "21_Jan", "21_May",
                                                              "24_July", "24_Oct", "24_Jan", "24_May")) 

############################ Sum abundance by site and taxa, filter for top 10 #############################

coverm_long_10p <- coverm_long_final %>%
  group_by(Taxa, SiteFull, Site) %>% #grouping taxa together for the same site
  #filter(value >= quantile(value, .95)) #take values top 95th percentile of each group (in group_by)
  #.95 to .90 for top 10%
  summarise(value=sum(as.numeric(value))) %>% #summing the group
  ungroup() %>% #don't group by Taxa anymore
  group_by(SiteFull) %>% #now group by site
  arrange(desc(value)) %>%
  mutate(rank=row_number()) %>% #make a new variable called rank where rank values
  filter(rank <= 10) #CHANGE THIS NUMBER TO GET TOP X OF COMMUNITY

############################ Plot the top 10 most abundant per site #############################

#see how many colors you'll need
length(unique(coverm_long_10p$Taxa))

#get most colors automatically, overshooting so I have options
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

#do some customization
col_vector <- c("#7FC97F","#BEAED4","#d9d9d9","#A6CEE3","#8dd3c7","#fb9a99","#BF5B17",
                "#FDC086","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D",
                "#666666","#FFFF99","#1F78B4","#B2DF8A")

####The following produces the bar plot in Figure 1, which was cropped and rotated in Biorender
dev.off()
plot <- coverm_long_10p %>%
  ggplot(aes(x = SiteFull, y = as.numeric(value), fill = Taxa)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  labs(x = "Site", y = "CoverM Percent Relative Abundance") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank()) + #turn this off to get the outline back)
  scale_y_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  ggtitle("Abundant SFB Taxa per Site") + #Change for top X grabbed
  facet_grid(.~Site, scales = "free") 
#scale_y_continuous(breaks = seq(0, 45, by=5))
#coord_flip()
plot

ggsave("Output/CoverM_Abundance_Top10perSite_GTDBv2.1_biorender_trans.png", plot, width = 10, height = 5, dpi = 500,
       bg = "transparent")
ggsave("Output/CoverM_Abundance_Top10perSite_GTDBv2.1_biorender_trans.pdf", plot, width = 10, height = 5, dpi = 500,
       bg = "transparent")



############################ Sum abundance by site and taxa, filter for bottom 10 #############################
############################ This is to understand the rare biosphere #############################


#drop 0s in coverm_long
coverm_long_final <- filter(coverm_long_final, value > 0)

coverm_long_10p_bottom <- coverm_long_final %>%
  group_by(Taxa, SiteFull, Site) %>% #grouping taxa together for the same site
  #filter(value >= quantile(value, .95)) #take values top 95th percentile of each group (in group_by)
  #.95 to .90 for top 10%
  summarise(value=sum(as.numeric(value))) %>% #summing the group
  ungroup() %>% #don't group by Taxa anymore
  group_by(SiteFull) %>% #now group by site
  arrange(value) %>% #Here is the change to get bottom, do not put desc, which gives default, ascend
  mutate(rank=row_number()) %>% #make a new variable called rank where rank values
  filter(rank <= 10) #CHANGE THIS NUMBER TO GET TOP X OF COMMUNITY

############################ Plot the top 10 least abundant per site #############################

#see how many colors you'll need
length(unique(coverm_long_10p_bottom$Taxa))

#get most colors automatically, overshooting so I have options
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))

#do some customization
col_vector <- c("#7FC97F","#BEAED4","#d9d9d9","#A6CEE3","#8dd3c7","#fb9a99","#BF5B17",
                "#FDC086","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D",
                "#666666","#FFFF99","#1F78B4","#B2DF8A", "#fb8072", "#ffed6f")

####The following produces the bar plot in Figure 1, which was cropped and rotated in Biorender
dev.off()
plot <- coverm_long_10p_bottom %>%
  ggplot(aes(x = SiteFull, y = as.numeric(value), fill = Taxa)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  labs(x = "Site", y = "CoverM Percent Relative Abundance") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        #panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
        #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_y_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  ggtitle("Low Abundance SFB Taxa per Site") + #Change for top X grabbed
  facet_grid(.~Site, scales = "free") 
#scale_y_continuous(breaks = seq(0, 45, by=5))
#coord_flip()
plot

ggsave("Output/CoverM_Abundance_Bottom10perSite_GTDBv2.1_trans.png", plot, width = 10, height = 5, dpi = 500,
       bg = "transparent")
ggsave("Output/CoverM_Abundance_Bottom10perSite_GTDBv2.1_trans.pdf", plot, width = 10, height = 5, dpi = 500,
       bg = "transparent")

