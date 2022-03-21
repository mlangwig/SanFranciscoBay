
rm(list=ls())


###################################### nitrogen ###############################################################

library(tidyr)
library(devtools)
library(ggplot2)
library(RColorBrewer)

#read in data
data <- read.delim(file = "Input/SFB_mebs_normScore_test.txt", header = TRUE)
#unique(data$site) 
data$site2 <- factor(data$site2, levels = c("4_1","8_1", "13", "21", "24"))
data$month <- factor(data$month, levels = c("July","Oct", "Jan", "May"))



#data<-data%>%separate(col=site,into = c("site2","month"), sep = "_")


dev.off()
bp_n <- ggplot(data, aes(x = site2, y = nitrogen, 
                       color=month)) + geom_boxplot(alpha = 1/2) + 
  labs(title = "MEBS Normalized Nitrogen Scores", y = "MEBS Normalized Nitrogen Score", 
       x = "Site") +
  theme_bw()
  #theme(legend.position = "none")
  #scale_y_continuous(expand=c(0,0), limits=c(0,1))

bp_n

bp_n_pretty <- bp_n + scale_color_manual(breaks = c("July", "Oct", "Jan",
                                                "May"),
                                     values=c("#88CCEE", "#44AA99", "#3F578F",
                                              "#BAB8B9"), #"#882255"), #, "#CAB2D6" #, "#6A3D9A"
                                     name = "Month",
                                    labels=c("July", "Oct", "Jan", "May")
                                    )
bp_n_pretty

#ggsave("Nitrogen_MEBS_sfb_boxplot.pdf", bp_n_pretty, width = 10, dpi = 500)

###################################### sulfur ###############################################################

dev.off()
bp_s <- ggplot(data, aes(x = site2, y = sulfur, 
                       color=month)) + geom_boxplot(alpha = 1/2) + 
  labs(title = "MEBS Normalized Sulfur Scores", y = "MEBS Normalized Sulfur Score", 
       x = "Site") +
  theme_bw()
#theme(legend.position = "none")
#scale_y_continuous(expand=c(0,0), limits=c(0,1))

bp_s

bp_s_pretty <- bp_s + scale_color_manual(breaks = c("July", "Oct", "Jan",
                                                    "May"),
                                     values=c("#88CCEE", "#44AA99", "#3F578F",
                                              "#BAB8B9"), #"#882255"), #, "#CAB2D6" #, "#6A3D9A"
                                     name = "Month",
                                     labels=c("July", "Oct", "Jan", "May")
)
bp_s_pretty

#ggsave("Sulfur_MEBS_sfb_boxplot.pdf", bp_s_pretty, width = 10, dpi = 500)

###################################### carbon ###############################################################

dev.off()
bp_c <- ggplot(data, aes(x = site2, y = carbon, 
                       color=month)) + geom_boxplot(alpha = 1/2) + 
  labs(title = "MEBS Normalized Carbon Scores", y = "MEBS Normalized Carbon Score", 
       x = "Site") +
  theme_bw()
#theme(legend.position = "none")
#scale_y_continuous(expand=c(0,0), limits=c(0,1))

bp_c

bp_c_pretty <- bp_c + scale_color_manual(breaks = c("Jan", "May", "July",
                                                "Oct"),
                                     values=c("#88CCEE", "#44AA99", "#3F578F",
                                              "#BAB8B9"), #"#882255"), #, "#CAB2D6" #, "#6A3D9A"
                                     name = "Month",
                                     labels=c("Jan", "May", "July", "Oct")
)
bp_c_pretty

ggsave("Carbon_MEBS_sfb_boxplot.pdf", bpc__pretty, width = 10, dpi = 500)

###################################### oxygen ###############################################################

dev.off()
bp_o <- ggplot(data, aes(x = site2, y = oxygen, 
                       color=month)) + geom_boxplot(alpha = 1/2) + 
  labs(title = "MEBS Normalized Oxygen Scores", y = "MEBS Normalized Oxygen Score", 
       x = "Site") +
  theme_bw()
#theme(legend.position = "none")
#scale_y_continuous(expand=c(0,0), limits=c(0,1))

bp_o

bp_o_pretty <- bp_o + scale_color_manual(breaks = c("Jan", "May", "July",
                                                "Oct"),
                                     values=c("#88CCEE", "#44AA99", "#3F578F",
                                              "#BAB8B9"), #"#882255"), #, "#CAB2D6" #, "#6A3D9A"
                                     name = "Month",
                                     labels=c("Jan", "May", "July", "Oct")
)
bp_o_pretty

ggsave("Oxygen_MEBS_sfb_boxplot.pdf", bp_o_pretty, width = 10, dpi = 500)


###################################### iron ###############################################################

dev.off()
bp_i <- ggplot(data, aes(x = site2, y = iron, 
                       color=month)) + geom_boxplot(alpha = 1/2) + 
  labs(title = "MEBS Normalized Iron Scores", y = "MEBS Normalized Iron Score", 
       x = "Site") +
  theme_bw()
#theme(legend.position = "none")
#scale_y_continuous(expand=c(0,0), limits=c(0,1))

bp_i

bp_i_pretty <- bp_i + scale_color_manual(breaks = c("Jan", "May", "July",
                                                "Oct"),
                                     values=c("#88CCEE", "#44AA99", "#3F578F",
                                              "#BAB8B9"), #"#882255"), #, "#CAB2D6" #, "#6A3D9A"
                                     name = "Month",
                                     labels=c("Jan", "May", "July", "Oct")
)
bp_i_pretty

ggsave("Iron_MEBS_sfb_boxplot.pdf", bp_i_pretty, width = 10, dpi = 500)

################################### combo ###################################

library(gridExtra)
library(cowplot)

combo_plot <- plot_grid(bp_n_pretty, bp_s_pretty, labels = c('A', 'B'),
                        ncol = 2)
combo_plot


ggsave("Output/combo_plot_NS.png", combo_plot, width = 10, height = 5, dpi = 500)


###################################

combo_plot2 <- plot_grid(bp_i_pretty, bp_n_pretty, 
                        bp_o_pretty, bp_s_pretty, labels = c('A', 'B', 'C',
                                                             'D'),
                        ncol = 2)
combo_plot2


ggsave("combo_plot2.png", combo_plot2, width = 15, height = 10, dpi = 500)






