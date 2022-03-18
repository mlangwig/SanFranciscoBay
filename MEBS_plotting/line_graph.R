library(reshape2)
library(ggplot2)
library(viridis)
library(tidyverse)
library(wesanderson)
#devtools::install_github("katiejolly/nationalparkcolors")
library(nationalparkcolors)

####################################### Real Data Scores by Site ###############################################
#read MEBS scores input
data_long <- read.delim(file = "Input/sfbArcBac_mebs_Norm_mod.txt", header = TRUE, sep = "\t")

#read in geochemical data
geochem_data<-read.delim("Input/Geochemical_Data.txt", header = TRUE)

#collapse the MEBS scores into 1 column instead of 5
data_longer<-pivot_longer(data_long, cols = S:N, values_to = "Scores")
#change the name of the 5th column to Cycle
colnames(data_longer)[5] <- "Cycle"

#sum values for each cycle at month at a specific site for each month
plot_input <- data_longer %>%
  group_by(Cycle, Month, Site) %>%
  summarise(mean(Scores), sd(Scores))
#change the annoying column name
colnames(plot_input)
colnames(plot_input)[4] <- "mean"
colnames(plot_input)[5] <- "STD"
#factor sites so can plot in correct order of sites
plot_input$Site_f = factor(plot_input$Site, levels = c('4_1', '8_1', '13', '21', '24'))

#drop C, Fe, and O
plot_input<-plot_input[!grepl("C", plot_input$Cycle),]
plot_input<-plot_input[!grepl("Fe", plot_input$Cycle),]
plot_input<-plot_input[!grepl("O", plot_input$Cycle),]

#add a unique column in plot_input so can map geochem data
new_col5<-paste(plot_input$Site,plot_input$Month,sep="_")
plot_input$Site_Month<-new_col5
#map salinity and temperature from geochem_data to the plot_input
plot_input <- geochem_data %>%
  dplyr::select(Site_Month, Sal, Temp) %>%
  right_join(plot_input, by = c("Site_Month" = "Site_Month"))

################## generate the plot ################
##create the text objects to plot below x axis
library(grid)
text_2011 <- textGrob("2011", gp=gpar(fontsize=6, fontface="bold"))
text_2012 <- textGrob("2012", gp=gpar(fontsize=6, fontface="bold"))

##plot
dev.off()
all_plot <- plot_input %>%
  ggplot(aes(x= factor(Month, levels = c("7", "10", "1", "5")), y=mean, group=Cycle, color=Cycle)) +
  #geom_ribbon(aes(ymin = mean - STD, ymax = mean + STD), alpha = 0.2, color = NA) +
  #geom_errorbar()
  geom_line(size = 1) +
  #scale_color_viridis(discrete = TRUE) +
  ggtitle("MEBS scores") +
  #theme_ipsum() +
  theme(legend.key.size = unit(.5, 'cm')) +
  theme(legend.text = element_text(size=8), axis.title.x = element_text(vjust = -2.8)) +
  scale_color_manual(values = rev(park_palette("Everglades"))) +
  ylab("MEBS Scores\n(average score per month)") +
  xlab("Month") +
  ylim(.2,.8) +
  scale_x_discrete(labels=c("1" = "Jan", "5" = "May", "7" = "July", "10" = "Oct")) +
  geom_line(aes(y=Sal/39), linetype = "dashed", color = "black") + #change to color = "black" and move linetype out of aes and change to linetype = "dashed"
  #geom_line(aes(y=Sal/39, linetype = "Salinity"), color = "red") + #for red line
  scale_y_continuous(sec.axis = sec_axis(~.*39, name = "Salinity (ppt)", breaks = seq(0, 35, 5))) +
  scale_linetype('') +
  theme(text = element_text(size = 9), axis.line.y.right = element_line(color = "black", linetype = "dashed")) + 
  #theme(text = element_text(size = 9), axis.line.y.right = element_line(color = "red")) + #for red line
  annotation_custom(text_2011, xmin=2.1,xmax=1,ymin=-0.25,ymax=-0.07) +
  annotation_custom(text_2012, xmin=6.1,xmax=1,ymin=-0.25,ymax=-0.07) +
  coord_cartesian(ylim=c(0,0.8), clip="off")
line_graph_facet<-all_plot + facet_grid(~Site_f)
line_graph_facet

ggsave(line_graph_facet, filename = "Output/MEBS_SFB_lineGraph_salinity.png", height = 2.5, width = 8, device = "png", dpi = 2500)
#ggsave(line_graph_facet, filename = "MEBS_SFB_lineGraph.pdf", height = 2, width = 8, dpi = 800)

####################################### Making python input for Val ###############################################

#modify data_longer for Val:
#add Site_Month column for geochem mapping
new_col<-paste(data_longer$Site,data_longer$Month,sep="_")
data_longer$Site_Month<-new_col

#map geochem data
data_longer_geo <- geochem_data %>%
  dplyr::select(Site_Month, Sal, Temp) %>%
  right_join(data_longer, by = c("Site_Month" = "Site_Month"))

#map abundance data
##read input
abundance<-read.delim("../Abundance_Barplot/Input/input_log.tsv")
##pivot long
abundance_long<-pivot_longer(abundance, cols = SF_Jan12_13:SF_Oct11_8_1, values_to = "Abundance_log")
colnames(abundance_long)[3] <- "Site"
##add Site_Month column for mapping
abundance_long$Site2 = abundance_long$Site
abundance_long$Site2<-gsub("4_1", "4.1", abundance_long$Site2)
abundance_long$Site2<-gsub("8_1", "8.1", abundance_long$Site2)
abundance_long_sep<-abundance_long %>% separate(Site2, into = c("Month", "Site"), sep = "_(?=[^_]+$)")
abundance_long_sep$Month<-gsub("SF_Oct11", "10", abundance_long_sep$Month)
abundance_long_sep$Month<-gsub("SF_May12", "5", abundance_long_sep$Month)
abundance_long_sep$Month<-gsub("SF_Jan12", "1", abundance_long_sep$Month)
abundance_long_sep$Month<-gsub("SF_Jul11", "7", abundance_long_sep$Month)

abundance_long_sep$Site<-gsub("4.1", "4_1", abundance_long_sep$Site)
abundance_long_sep$Site<-gsub("8.1", "8_1", abundance_long_sep$Site)

new_col2<-paste(abundance_long_sep$Site,abundance_long_sep$Month,sep="_")
abundance_long_sep$Site_Month<-new_col2

##create unique bin_site column to map abundance
new_col3<-paste(abundance_long_sep$Site_Month,abundance_long_sep$Bin,sep="_")
abundance_long_sep$Site_Month_Bin<-new_col3

new_col4<-paste(data_longer_geo$Site_Month,data_longer_geo$Genome,sep="_")
data_longer_geo$Site_Month_Bin<-new_col4

##FINALLY HERE --> MAP ABUNDANCE INFORMATION ONTO THE WHOLE DATA FRAME
data_longer_geo_abun <- abundance_long_sep %>%
  dplyr::select(Site_Month_Bin, Abundance_log) %>%
  right_join(data_longer_geo, by = c("Site_Month_Bin" = "Site_Month_Bin"))

#write data_longer_geo as is
write_delim(data_longer_geo_abun, file = "Output/Val_output_allMEBS.txt")

#make a version just sulfur and nitrogen
data_longer_geo_abun<-data_longer_geo_abun[!grepl("C", data_longer_geo_abun$Cycle),]
data_longer_geo_abun<-data_longer_geo_abun[!grepl("Fe", data_longer_geo_abun$Cycle),]
data_longer_geo_abun<-data_longer_geo_abun[!grepl("O", data_longer_geo_abun$Cycle),]
#write it
write_delim(data_longer_geo_abun, file = "Output/Val_output_MEBS_SandN.txt")



####################################### Real Data Sulfur ###############################################

#read input
data_long <- read.delim(file = "../sfbArcBac_mebs_Norm_mod.txt", header = TRUE, sep = "\t")

# #melt to change from wide to long
# data_long <- melt(data, id = c("Genome", "Taxonomy", "Site", "Month"))
# write_delim(data_long, "SFB_MEBS_data_long.tsv", delim = "\t")

#Get only sulfur values
#sulfur_data <- subset(data_long, variable=='sulfur')
sulfur_data <- data_long[1:5]

#sum sulfur values for each month at a specific site
sulfur_data_grouped <- sulfur_data %>%
  group_by(Site, Month) %>%
  summarise(mean(S), sd(S))
#change the annoying column name
colnames(sulfur_data_grouped)
colnames(sulfur_data_grouped)[3] <- "mean"
colnames(sulfur_data_grouped)[4] <- "STD"

#plot
sulf_plot <- sulfur_data_grouped %>%
  ggplot( aes(x= factor(Month, levels = c("1", "5", "7", "10")), 
              y=mean, group=Site, color=Site)) +
  geom_ribbon(aes(ymin = mean - STD, ymax = mean + STD), alpha = 0.2, color = NA) +
  #geom_errorbar()
  geom_line(size = 1) +
  #scale_color_viridis(discrete = TRUE) +
  ggtitle("MEBS sulfur score across sites") +
  #theme_ipsum() +
  theme(legend.key.size = unit(1, 'cm')) +
  theme(legend.text = element_text(size=10)) +
  scale_color_manual(values = wes_palette("Rushmore1")) +
  ylab("MEBS Sulfur Score") +
  xlab("Month")
sulf_plot + facet_grid(~Site)



####################################### Real Data Nitrogen ###############################################

#Get only sulfur values
nitrogen_data <- subset(data_long, variable=='nitrogen')

#sum sulfur values for each month at a specific site
nitrogen_data_grouped <- nitrogen_data %>%
  group_by(Site, Month) %>%
  summarise(sum(value))
#change the annoying column name
colnames(nitrogen_data_grouped)
colnames(nitrogen_data_grouped)[3] <- "sum"

#plot
nitr_plot <- nitrogen_data_grouped %>%
  ggplot( aes(x= factor(Month, levels = c("Jan", "May", "July", "Oct")), 
              y=sum, group=Site, color=Site)) +
  geom_line(size = 1) +
  #scale_color_viridis(discrete = TRUE) +
  ggtitle("MEBS nitrogen score across sites") +
  theme_ipsum(base_family = "sans") +
  scale_color_manual(values = wes_palette("Rushmore1")) +
  ylab("MEBS Nitrogen Score") +
  xlab("Month")
nitr_plot

####################################### Real Data Iron ###############################################

#Get only sulfur values
iron_data <- subset(data_long, variable=='iron')

#sum sulfur values for each month at a specific site
iron_data_grouped <- iron_data %>%
  group_by(Site, Month) %>%
  summarise(sum(value))
#change the annoying column name
colnames(iron_data_grouped)
colnames(iron_data_grouped)[3] <- "sum"

#plot
iron_plot <- iron_data_grouped %>%
  ggplot( aes(x= factor(Month, levels = c("Jan", "May", "July", "Oct")), 
              y=sum, group=Site, color=Site)) +
  geom_line(size = 1) +
  #scale_color_viridis(discrete = TRUE) +
  ggtitle("MEBS iron score across sites") +
  theme_ipsum(base_family = "sans") +
  scale_color_manual(values = wes_palette("Rushmore1")) +
  ylab("MEBS Iron Score") +
  xlab("Month")
iron_plot

####################################### Real Data Oxygen ###############################################

#Get only sulfur values
oxy_data <- subset(data_long, variable=='oxygen')

#sum sulfur values for each month at a specific site
oxy_data_grouped <- oxy_data %>%
  group_by(Site, Month) %>%
  summarise(sum(value))
#change the annoying column name
colnames(oxy_data_grouped)
colnames(oxy_data_grouped)[3] <- "sum"

#plot
oxy_plot <- oxy_data_grouped %>%
  ggplot( aes(x= factor(Month, levels = c("Jan", "May", "July", "Oct")), 
              y=sum, group=Site, color=Site)) +
  geom_line(size = 1) +
  #scale_color_viridis(discrete = TRUE) +
  ggtitle("MEBS oxygen score across sites") +
  theme_ipsum(base_family = "sans") +
  scale_color_manual(values = wes_palette("Rushmore1")) +
  ylab("MEBS Oxygen Score") +
  xlab("Month")
oxy_plot

####################################### Facet Plots ###############################################

library(gridExtra)
library(cowplot)

combo_plot <- plot_grid(sulf_plot + theme(legend.position="none"), nitr_plot + theme(legend.position="none"),
                        iron_plot + theme(legend.position="none"), oxy_plot + theme(legend.position="none"),
                        labels = c('A', 'B', 'C', 'D'),
                        ncol = 2)
combo_plot

legend <- get_legend(
  oxy_plot +
    theme(legend.key.size = unit(1, 'cm'),
          legend.text = element_text(size=15),
          legend.title = element_text(size = 17)
          ) #legend.position = "bottom"
  )
  # create some space to the left of the legend
  #sulf_plot + theme(legend.box.margin = margin(0, 0, 0, -5))
#)

facet <- plot_grid(combo_plot, legend, rel_widths = c(3, .4))
facet

#test <- combo_plot + draw_grob(legend, 2/3.3, 0, .3/3.3, 1, vjust = .47, hjust = 1.5)
#test

ggsave("MEBS_scores_SFBcombo_plot_Norm.pdf", facet, width = 15, height = 10, dpi = 500)




####################################### Example ###############################################

# # Libraries
# library(ggplot2)
# library(babynames) # provide the dataset: a dataframe called babynames
# library(dplyr)
# library(hrbrthemes)
# library(viridis)
# 
# # Keep only 3 names
# don <- babynames %>% 
#   filter(name %in% c("Ashley", "Patricia", "Helen")) %>%
#   filter(sex=="F")
# 
# # Plot
# don %>%
#   ggplot( aes(x=year, y=n, group=name, color=name)) +
#   geom_line() +
#   scale_color_viridis(discrete = TRUE) +
#   ggtitle("Popularity of American names in the previous 30 years") +
#   theme_ipsum() +
#   ylab("Number of babies born")


