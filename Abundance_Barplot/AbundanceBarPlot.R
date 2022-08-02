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


