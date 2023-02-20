setwd("~/Google Drive/My Drive/SF_Bay/SanFranciscoBay/Abundance_Barplot/")

############################ Abundance Stacked Bar Plot ##########################################
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)

data <- read.delim("Input/input_log.tsv") #read in data log values
data2 <- read.delim("Input/input.tsv") #non log value input
data_renamed <- read.delim("Input/input_log_renamed.tsv")
data2_renamed <- read.delim("Input/input_renamed.tsv")
data_gtdbv2.1 <- read.delim("Input/input_GTDBTk_v2.1.0.txt")


#colnames(data2) #get names of columns
#ncol(data2) #get number of columns
data_long <- melt(data_gtdbv2.1, id = c("Taxa", "Bin")) #transform from wide to long #CHANGE DATA TO DATA2 FOR NON-LOG
colnames(data_long)[3] <- "site" #rename column
#rename sites
data_long$site<-gsub("SF_Jan12_4_1", "4_1_Jan", data_long$site)
data_long$site<-gsub("SF_May12_4_1", "4_1_May", data_long$site)
data_long$site<-gsub("SF_Jul11_4_1", "4_1_July", data_long$site)
data_long$site<-gsub("SF_Oct11_4_1", "4_1_Oct", data_long$site)
data_long$site<-gsub("SF_Jan12_8_1", "8_1_Jan", data_long$site)
data_long$site<-gsub("SF_May12_8_1", "8_1_May", data_long$site)
data_long$site<-gsub("SF_Jul11_8_1", "8_1_July", data_long$site)
data_long$site<-gsub("SF_Oct11_8_1", "8_1_Oct", data_long$site)
data_long$site<-gsub("SF_Jan12_13", "13_Jan", data_long$site)
data_long$site<-gsub("SF_May12_13", "13_May", data_long$site)
data_long$site<-gsub("SF_Jul11_13", "13_July", data_long$site)
data_long$site<-gsub("SF_Oct11_13", "13_Oct", data_long$site)
data_long$site<-gsub("SF_Jan12_21", "21_Jan", data_long$site)
data_long$site<-gsub("SF_May12_21", "21_May", data_long$site)
data_long$site<-gsub("SF_Jul11_21", "21_July", data_long$site)
data_long$site<-gsub("SF_Oct11_21", "21_Oct", data_long$site)
data_long$site<-gsub("SF_Jan12_24", "24_Jan", data_long$site)
data_long$site<-gsub("SF_May12_24", "24_May", data_long$site)
data_long$site<-gsub("SF_Jul11_24", "24_July", data_long$site)
data_long$site<-gsub("SF_Oct11_24", "24_Oct", data_long$site)


data_long$site <- factor(data_long$site, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
                                                  "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                                  "13_July", "13_Oct", "13_Jan", "13_May",
                                                  "21_July", "21_Oct", "21_Jan", "21_May",
                                                  "24_July", "24_Oct", "24_Jan", "24_May")) #set order of x axis 

data_long_5p <- data_long %>%
  group_by(Taxa, site) %>% #grouping taxa together for the same site
  #filter(value >= quantile(value, .95)) #take values top 95th percentile of each group (in group_by)
                                        #.95 to .90 for top 10%
  summarise(value=sum(value)) %>% #summing the group
  ungroup() %>% #don't group by Taxa anymore
  group_by(site) %>% #now group by site
  arrange(desc(value)) %>%
  mutate(rank=row_number()) %>% #make a new variable called rank where rank values
  filter(rank <= 10) #CHANGE THIS NUMBER TO GET TOP X OF COMMUNITY
  

# plot <- ggplot(data_long, aes(x = site, y = as.numeric(value), fill = Taxa)) + 
#   geom_bar(stat = "identity") +
#   labs(x = "Site", y = "Abundance") + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
#         panel.background = element_blank()) +
#   scale_y_continuous(expand = c(0, 0))
# plot
# 
# ggsave("Abundance_StackedBarPlot.pdf", plot, dpi = 500)

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

plot2 <- data_long_5p %>%
  ggplot(aes(x = site, y = as.numeric(value), fill = Taxa)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_vector) +
  labs(x = "Site", y = "Abundance") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Ten most abundant SFB taxa per site") #Change for top X grabbed
plot2


#ggsave("Abundance_StackedBarPlotTop5perSite_GTDB.pdf", plot2, width = 8, dpi = 500)
#ggsave("Abundance_StackedBarPlotTop5perSite_GTDB.png", plot2, width = 8, height = 5, dpi = 500)

ggsave("Output/Abundance_BarPlotTop5perSite_GTDBv2.1.pdf", plot2, width = 8, dpi = 500)
ggsave("Output/Abundance_BarPlotTop5perSite_GTDBv2.1.png", plot2, width = 8, height = 5, dpi = 500)

ggsave("Output/Abundance_BarPlotTop10perSite_GTDBv2.1.pdf", plot2, width = 8, dpi = 500)
ggsave("Output/Abundance_BarPlotTop10perSite_GTDBv2.1.png", plot2, width = 10, height = 5, dpi = 500)


##### Subset abundance table for specific community members/ecophysiologies

#Sulfate reducers aka dsrA encoding that are reductive in phylogeny
#SRB <- read.delim(file = "Input/reductive_dsrAs_sort.txt", header = FALSE)

library(stringr)
library(dplyr)

### Histogram of all abundance
all_Abund <- data_long$value
hist(all_Abund,
     main = "SFB Abundance Distribution",
     xlab = "Abundance")
#line for lowest high abundance of reductive dsrA-encoding Delta
abline(v = .4635, col = "black", lwd=3, lty=2)


############  subset data_long for Deltaproteobacteria to see the distribution
data_long_deltas<-data_long%>%filter(str_detect(Taxa,"Deltaproteobacteria"))
data_long_log_deltas<-data_long%>%filter(str_detect(Taxa,"Deltaproteobacteria"))

#histogram non-log transformed  
Delta_Abund <- data_long_deltas$value
hist(Delta_Abund,
     xlim = c(.01,15),
     breaks = 500,
     #breaks = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,2,3,4,5,6,7,8,9,10,15,20,30),
     main = "Deltaproteobacteria Abundance Distribution",
     xlab = "Abundance")

#histogram log transformed
hist(Delta_Abund,
     main = "Deltaproteobacteria Abundance Distribution",
     xlab = "Abundance")
abline(v = .4635, col = "black", lwd=3, lty=2)

############  Gammas
data_long_gammas<-data_long%>%filter(str_detect(Taxa,"c__Gammaproteobacteria"))
Gamma_Abund <- data_long_gammas$value

#histogram non-log transformed 
hist(Gamma_Abund,
     xlim = c(.01,6),
     breaks = 500,
     main = "Gammaproteobacteria Abundance Distribution",
     xlab = "Abundance")

#histogram log transformed 
hist(Gamma_Abund,
     main = "Gammaproteobacteria Abundance Distribution",
     xlab = "Abundance")
abline(v = .7095, col = "black", lwd=3, lty=2)


############  AmoA-encoding
data_long_thaums<-data_long%>%filter(str_detect(Taxa,"c__Nitrososphaeria"))
Thaums_Abund <- data_long_thaums$value

hist(Thaums_Abund,
     main = "Thaumarchaeota Abundance Distribution",
     xlab = "Abundance")
abline(v = .778, col = "black", lwd=3, lty=2)

###

data_long_comammox<-data_long%>%filter(str_detect(Taxa,"p__Nitrospirota"))
Comammox_Abund <- data_long_comammox$value

hist(Comammox_Abund,
     main = "Nitrospirota Abundance Distribution",
     xlab = "Abundance")
abline(v = .764, col = "black", lwd=3, lty=2)


write_tsv(data_long, file = "data_long_allAbundance.tsv")





############################ CoverM Input Abundance Stacked Bar Plot ##########################################

#read in coverm input
coverm<-read.delim(file = "Input/SFBMAGS-vs-Reads-CoverM.tsv")

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

#set order of x axis 
coverm_long$SiteFull <- factor(coverm_long$SiteFull, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
                                                    "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                                    "13_July", "13_Oct", "13_Jan", "13_May",
                                                    "21_July", "21_Oct", "21_Jan", "21_May",
                                                    "24_July", "24_Oct", "24_Jan", "24_May")) 

#vlookup to add taxonomy
coverm_long <- data_gtdbv2.1 %>%
  dplyr::select("Bin", "Taxa") %>%
  right_join(coverm_long, by = c("Bin" = "Genome"))

#set order of x axis 
coverm_long$Site <- factor(coverm_long$Site, levels=c("4_1", "8_1", "13", "21", "24"))

############################# sum abundance by site and taxonomy
coverm_long_5p <- coverm_long %>%
  group_by(Taxa, SiteFull, Site) %>% #grouping taxa together for the same site
  #filter(value >= quantile(value, .95)) #take values top 95th percentile of each group (in group_by)
  #.95 to .90 for top 10%
  summarise(value=sum(value)) %>% #summing the group
  ungroup() %>% #don't group by Taxa anymore
  group_by(SiteFull) %>% #now group by site
  arrange(desc(value)) %>%
  mutate(rank=row_number()) %>% #make a new variable called rank where rank values
  filter(rank <= 5) #CHANGE THIS NUMBER TO GET TOP X OF COMMUNITY

coverm_long_10p <- coverm_long %>%
  group_by(Taxa, SiteFull, Site) %>% #grouping taxa together for the same site
  #filter(value >= quantile(value, .95)) #take values top 95th percentile of each group (in group_by)
  #.95 to .90 for top 10%
  summarise(value=sum(value)) %>% #summing the group
  ungroup() %>% #don't group by Taxa anymore
  group_by(SiteFull) %>% #now group by site
  arrange(desc(value)) %>%
  mutate(rank=row_number()) %>% #make a new variable called rank where rank values
  filter(rank <= 10) #CHANGE THIS NUMBER TO GET TOP X OF COMMUNITY


#drop 0s in coverm_long
coverm_long <- filter(coverm_long, value > 0)

coverm_long_10p_bottom <- coverm_long %>%
  group_by(Taxa, SiteFull, Site) %>% #grouping taxa together for the same site
  #filter(value >= quantile(value, .95)) #take values top 95th percentile of each group (in group_by)
  #.95 to .90 for top 10%
  summarise(value=sum(value)) %>% #summing the group
  ungroup() %>% #don't group by Taxa anymore
  group_by(SiteFull) %>% #now group by site
  arrange(value) %>%
  mutate(rank=row_number()) %>% #make a new variable called rank where rank values
  filter(rank <= 10) #CHANGE THIS NUMBER TO GET TOP X OF COMMUNITY

############################# plot

#get most colors automatically
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

#do some customization
col_vector <- c("#7FC97F","#BEAED4","#d9d9d9","#A6CEE3","#8dd3c7","#fb9a99","#BF5B17",
                "#FDC086","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D",
                "#666666","#FFFF99","#1F78B4","#B2DF8A")

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
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank()) + #turn this off to get the outline back)
  scale_y_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  ggtitle("CoverM ten least abundant SFB taxa per site") + #Change for top X grabbed
  facet_grid(.~Site, scales = "free") 
  #scale_y_continuous(breaks = seq(0, 45, by=5))
#coord_flip()
plot


ggsave("Output/CoverM_Abundance_Top5perSite_GTDBv2.1_facet.pdf", plot, width = 8, dpi = 500)
ggsave("Output/CoverM_Abundance_Top5perSite_GTDBv2.1_facet.png", plot, width = 10, height = 5, dpi = 500)

ggsave("Output/CoverM_Abundance_Top10perSite_GTDBv2.1_facet.pdf", plot, width = 8, dpi = 500)
ggsave("Output/CoverM_Abundance_Top10perSite_GTDBv2.1_facet.png", plot, width = 10, height = 5, dpi = 500)
ggsave("Output/CoverM_Abundance_Top10perSite_GTDBv2.1_facet_wide.png", plot, width = 15, height = 5, dpi = 500)


ggsave("Output/CoverM_Abundance_Top10perSite_GTDBv2.1_biorender.png", plot, width = 10, height = 5, dpi = 500,
       bg = "transparent")

ggsave("Output/CoverM_Abundance_Top10perSite_GTDBv2.1_biorender_trans.png", plot, width = 10, height = 5, dpi = 500,
       bg = "transparent")


############################ CoverM Input Abundance Stacked Area Plot ##########################################

plot <- coverm_long_10p %>%
  ggplot(aes(x = SiteFull, y = as.numeric(value), fill = Taxa)) + 
  geom_area() #+
  # scale_fill_manual(values = col_vector) +
  # labs(x = "Site", y = "CoverM Percent Relative Abundance") +
  # guides(fill=guide_legend(override.aes = list(size=3))) +
  # theme_bw() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
  #       panel.background = element_blank()) +
  # scale_y_continuous(expand = c(0, 0)) +
  # ggtitle("CoverM ten most abundant SFB taxa per site") + #Change for top X grabbed
  # facet_grid(.~Site, scales = "free")
#coord_flip()
plot

############################ Write data outputs ##########################################

write.table(coverm_long_10p, file = "Output/CoverM_long_top10_Abund.tsv", quote = FALSE, col.names = TRUE,
          row.names = FALSE, sep = "\t")


############################ Plotting rare biosphere ##########################################

#############create input
coverm_long_all <- coverm_long %>%
  group_by(Taxa, Site) %>% #grouping taxa together for the same site
  #filter(value >= quantile(value, .95)) #take values top 95th percentile of each group (in group_by)
  #.95 to .90 for top 10%
  summarise(value=sum(value)) %>% #summing the group
  ungroup() %>% #don't group by Taxa anymore
  group_by(Site) #%>% #now group by site
  #arrange(desc(value))

#drop orgs that have 0s for abundance
coverm_long_all <- filter(coverm_long_all, value > 0)

#separate into sites
coverm_long_all_4_1<-coverm_long_all %>% filter(grepl('4_1', Site))
#sort
coverm_long_all_4_1<-arrange(coverm_long_all_4_1, desc(value))

coverm_long_all_4_1$Taxa <- factor(coverm_long_all_4_1$Taxa, levels=coverm_long_all_4_1$Taxa)

#order the input data and factor
input_gs<-arrange(input_gs, desc(value))
input_gs$Class <- factor(input_gs$Class, levels=input_gs$Class)

###########plot
dev.off()
plot <- coverm_long_all_4_1 %>%
  ggplot(aes(x = Taxa, y = value, fill = Taxa)) + 
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  #scale_fill_manual(values = col_vector) +
  labs(x = "Taxa", y = "CoverM % Relative Abundance") +
  #guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  ggtitle("Site 4.1 CoverM Abundance") + #Change for top X grabbed
  geom_hline(yintercept=1,linetype=2,color="grey") +
  guides(fill="none")
  #facet_grid(.~SiteFull, scales = "free")
  #scale_y_continuous(breaks = seq(0, 45, by=5))
#coord_flip()
plot





