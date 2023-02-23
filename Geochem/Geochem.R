############################ Geochem Plotting ##########################################

library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)

############################ read in and process input ##########################################

geochem <- read.delim("Input/GeochemData.tsv")

#drop some data
geochem <- select(geochem, c(SiteFull, Site, Month, Temp, Sal, 
                             NO3, NH4, C_N))

#melt to long format
geochem_long <- melt(geochem, id = c("SiteFull", "Site", "Month"))

#change Jul to July
geochem_long$Month <- gsub("Jul", "July", geochem_long$Month)
geochem_long$Month <- factor(geochem_long$Month, levels = c("July", "Oct", "Jan", "May"))

############################ Set order of sites for plotting #############################
#set order of x axis 
geochem_long$Site <- factor(geochem_long$Site, levels=c("4_1", "8_1", "13", "21", "24"))
#set order of x axis 
geochem_long$SiteFull <- factor(geochem_long$SiteFull, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
                                                                          "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                                                          "13_July", "13_Oct", "13_Jan", "13_May",
                                                                          "21_July", "21_Oct", "21_Jan", "21_May",
                                                                          "24_July", "24_Oct", "24_Jan", "24_May"))

####################################### Plotting Supplemental ##########################################

library(grid)
text_2011 <- textGrob("2011", gp=gpar(fontsize=15, fontface="bold"))
text_2012 <- textGrob("2012", gp=gpar(fontsize=15, fontface="bold"))

dev.off()
plot <- geochem_long %>%
  ggplot(aes(x = Month, y = as.numeric(value), group = variable)) + #group to make this behave
  geom_line(aes(color = variable, linetype = variable), size = 2) +
  #scale_fill_manual(values = col_vector) +
  labs(x = "Month", y = "Value") +
  #guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        #legend.background = element_rect(color = "white"),
        #legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
        #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_y_continuous(expand = c(0, 0), breaks = c(0,5,10,15,20,25,30)) + #turn this on to make it look aligned with ticks
  scale_color_discrete(name = "",labels = c("Temp (ºC)", "Salinity (PSU)", "NO3 (mM)",
                                 "NH4 (mM)", "C/N")) +
  scale_linetype_discrete(name = "",labels = c("Temp (ºC)", "Salinity (PSU)", "NO3 (mM)",
                                             "NH4 (mM)", "C/N")) + #need both color and linetype to only get 1 legend/merge them
  ggtitle("SFB Geochemistry") + #Change for top X grabbed
  facet_wrap(.~Site, scales = "free") +
  annotation_custom(text_2011, xmin=2.1,xmax=1,ymin=-37) +
  annotation_custom(text_2012, xmin=6,xmax=1,ymin=-37) +
  coord_cartesian(ylim = c(0,31), xlim = c(1.5,3.5), clip="off") #messed with xlim to get less space between plot and margin
#scale_y_continuous(breaks = seq(0, 45, by=5))
  #coord_flip()
plot

# , width = 15, height = 7,
ggsave("Output/Geochem_Supplemental.png", plot, width = 15, height = 7, dpi = 500)




############################## Input for main figure geochem, salinity only ##########################################

geochem_sal <- geochem_long %>% filter(grepl(c("Sal"), variable))

############################ Set order of sites for plotting #############################
#set order of x axis 
geochem_sal$Site <- factor(geochem_sal$Site, levels=c("4_1", "8_1", "13", "21", "24"))
#set order of x axis 
geochem_sal$SiteFull <- factor(geochem_sal$SiteFull, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
                                                                "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                                                "13_July", "13_Oct", "13_Jan", "13_May",
                                                                "21_July", "21_Oct", "21_Jan", "21_May",
                                                                "24_July", "24_Oct", "24_Jan", "24_May"))


############### plotting #############

text_2011 <- textGrob("2011", gp=gpar(fontsize=22, fontface="bold"))
text_2012 <- textGrob("2012", gp=gpar(fontsize=22, fontface="bold"))

dev.off()
p2 <- geochem_sal %>%
  ggplot(aes(x = Month, y = as.numeric(value), group = Site)) + #group to make this behave
  geom_line(aes(color = Site), size = 2) +
  #scale_fill_manual(values = col_vector) +
  labs(x = "Month", y = "Salinity (PSU)") +
  #guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
    text = element_text(size = 25),
    legend.key.size = unit(1, 'cm'),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) + #turn this off to get the outline back)
  #scale_y_continuous(expand = c(0, 0), breaks = c(7,9,11,13,15,17,19,21)) + #turn this on to make it look aligned with ticks
  scale_color_manual(values = c("#4F508C", "#B56478", "#CF9A28", "#28827A", "#3F78C1"),
                     name = "",
                     labels = c("4.1", "8.1", "13", "21", "24")) +
  ggtitle("") +
  annotation_custom(text_2011, xmin=2.1,xmax=1,ymin=-45) +
  annotation_custom(text_2012, xmin=6,xmax=1,ymin=-45) +
  coord_cartesian(ylim = c(0,31), xlim = c(1.5,3.5), clip="off") #messed with xlim to get less space between plot and margin
#scale_y_continuous(breaks = seq(0, 45, by=5))
#coord_flip()
p2


#save
ggsave("Output/Salinity.png", p2, dpi = 500, bg = "transparent")



