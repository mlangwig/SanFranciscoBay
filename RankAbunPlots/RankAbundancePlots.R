############################ Abundance Stacked Bar Plot ##########################################
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)

############################ Read in Inputs #############################

#coverm input
coverm <- read.delim2(file = "../Abundance_Barplot/Output/coverm_long_final.tsv")

############################ Get rid of month info #############################
coverm <- coverm %>%
  group_by(Taxa, Site) %>% #grouping taxa together for the same site
  #filter(value >= quantile(value, .95)) #take values top 95th percentile of each group (in group_by)
  #.95 to .90 for top 10%
  summarise(value=sum(as.numeric(value))) %>% #summing the group
  ungroup() %>% #don't group by Taxa anymore
  group_by(Site) #%>% #now group by site

############################ Generate input per site #############################
#generate separate data frames
coverm_4_1 <- coverm %>% 
  filter(grepl("4_1", Site)) %>%
  arrange(desc(value))
coverm_4_1$Taxa<-factor(coverm_4_1$Taxa, levels=coverm_4_1$Taxa)

coverm_8_1 <- coverm %>% 
  filter(grepl("8_1", Site)) %>%
  arrange(desc(value))
coverm_8_1$Taxa<-factor(coverm_8_1$Taxa, levels=coverm_8_1$Taxa)

coverm_13 <- coverm %>% 
  filter(grepl("13", Site)) %>%
  arrange(desc(value))
coverm_13$Taxa<-factor(coverm_13$Taxa, levels=coverm_13$Taxa)

coverm_21 <- coverm %>% 
  filter(grepl("21", Site)) %>%
  arrange(desc(value))
coverm_21$Taxa<-factor(coverm_21$Taxa, levels=coverm_21$Taxa)

coverm_24 <- coverm %>% 
  filter(grepl("24", Site)) %>%
  arrange(desc(value))
coverm_24$Taxa<-factor(coverm_24$Taxa, levels=coverm_24$Taxa)

############################ Plotting #############################

############################ 4_1 #############################
dev.off()
p4_1 <- ggplot(coverm_4_1, aes(x = Taxa, y = value, group=1)) + 
  geom_point() +
  geom_line() +
  labs(x = "Taxa", y = "CoverM % Relative Abundance") +
  ggtitle("Site 4.1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
        #panel.border = element_blank()) +
  #scale_fill_viridis_d() +
  geom_hline(yintercept=1,linetype=2,color="black") +
  guides(fill="none") +
  #facet_grid(.~SiteFull, scales = "free")
  scale_x_discrete(limits=rev) +
  scale_y_continuous(expand = c(0,0), limits = c(0,60)) +
  coord_flip()
p4_1

############################ 8_1 #############################

dev.off()
p8_1 <- ggplot(coverm_8_1, aes(x = Taxa, y = value, fill = Taxa, group=1)) + 
  geom_point() +
  geom_line() +
  labs(x = "Taxa", y = "CoverM % Relative Abundance") +
  ggtitle("Site 8.1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) +
  scale_fill_viridis_d() +
  geom_hline(yintercept=1,linetype=2,color="black") +
  guides(fill="none") +
  #facet_grid(.~SiteFull, scales = "free")
  scale_x_discrete(limits=rev) +
  scale_y_continuous(expand = c(0,0), limits = c(0,50)) +
  coord_flip()
p8_1

############################ 13 #############################

dev.off()
p13 <- ggplot(coverm_13, aes(x = Taxa, y = value, fill = Taxa, group=1)) + 
  geom_point() +
  geom_line() +
  labs(x = "Taxa", y = "CoverM % Relative Abundance") +
  ggtitle("Site 13") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) +
  scale_fill_viridis_d() +
  geom_hline(yintercept=1,linetype=2,color="black") +
  guides(fill="none") +
  #facet_grid(.~SiteFull, scales = "free")
  scale_x_discrete(limits=rev) +
  scale_y_continuous(expand = c(0,0), limits = c(0,60)) +
  coord_flip()
p13

############################ 21 #############################

dev.off()
p21 <- ggplot(coverm_21, aes(x = Taxa, y = value, fill = Taxa, group = 1)) + 
  geom_point() +
  geom_line() +
  labs(x = "Taxa", y = "CoverM % Relative Abundance") +
  ggtitle("Site 21") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) +
  scale_fill_viridis_d() +
  geom_hline(yintercept=1,linetype=2,color="black") +
  guides(fill="none") +
  #facet_grid(.~SiteFull, scales = "free")
  scale_x_discrete(limits=rev) +
  scale_y_continuous(expand = c(0,0), limits = c(0,50)) +
  coord_flip()
p21

############################ 24 #############################

dev.off()
p24 <- ggplot(coverm_24, aes(x = Taxa, y = value, fill = Taxa, group = 1)) + 
  geom_point() +
  geom_line() +
  labs(x = "Taxa", y = "CoverM % Relative Abundance") +
  ggtitle("Site 24") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) +
  scale_fill_viridis_d() +
  geom_hline(yintercept=1,linetype=2,color="black") +
  guides(fill="none") +
  #facet_grid(.~SiteFull, scales = "free")
  scale_x_discrete(limits=rev) +
  scale_y_continuous(expand = c(0,0), limits = c(0,60)) +
  coord_flip()
p24

############################ cow plot to put them together #############################

library(patchwork)
plot_patch <- p4_1 + p8_1 + p13 + p21 + p24 +
  plot_layout(ncol = 2) +
  plot_annotation(title = "SFB Rank Abundances")
plot_patch

ggsave("Output/RankAbundances.png", plot_patch, width = 10, height = 16, dpi = 600)
ggsave("Output/RankAbundances.pdf", plot_patch, width = 10, height = 16, dpi = 600)

