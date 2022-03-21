############################ Abundance Stacked Bar Plot ##########################################
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)

data <- read.delim("Input/input_log.tsv") #read in data
data2 <- read.delim("Input/input.tsv")


colnames(data2) #get names of columns
ncol(data2) #get number of columns
data_long <- melt(data2, id = c("Taxa", "Bin")) #transform from wide to long
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
  filter(rank <= 5) #grab the ones ranked top number inserted
  

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

plot2 <- data_long_5p %>%
  ggplot(aes(x = site, y = as.numeric(value), fill = Taxa)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = as.vector(rev(kelly(n = 13)))) +
  labs(x = "Site", y = "Abundance") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Five most abundant SFB taxa per site")
plot2


ggsave("Abundance_StackedBarPlotTop5perSite.pdf", plot2, width = 8, dpi = 500)
ggsave("Abundance_StackedBarPlotTop5perSite.png", plot2, width = 8, height = 5, dpi = 500)

