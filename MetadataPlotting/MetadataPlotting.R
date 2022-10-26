#setwd
setwd("~/Google Drive/My Drive/SF_Bay/SanFranciscoBay/MetadataPlotting/")

#libraries
library(dplyr)
library(tidyverse)
library(wesanderson)

#read metadata input
input<-read_tsv(file = "input/SFB_MAGs_metadata.txt")

#separate taxonomy for math and plotting
input <- input %>% separate(GTDB, c("d", "p", "c", "o", "f", "g", "s"), 
                            sep= ";")

#count number of genomes in each class
c_count <- input %>%
  group_by(c) %>% 
  summarise(count_c=length(c)) %>%
  arrange(desc(count_c))

#obtain average genome size per class
input_gs <- input %>%
  group_by(c) %>% 
  summarise(value=mean(GenomeSize_MB)) %>%
  arrange(desc(value))

#map the class count to the input_gs table
input_gs <- c_count %>%
  dplyr::select("c", "count_c") %>%
  right_join(input_gs, by = c("c" = "c"))

#add parentheses around numbers in column
input_gs <- input_gs %>%
  mutate(count_c = paste0("(", count_c, ")"))

#merge the columns to show value when plotting
input_gs$Class <- paste0(input_gs$c, " ", input_gs$count_c)

#drop the genome with no class
input_gs<-filter(input_gs, c != "c__")

#order the input data and factor
input_gs<-arrange(input_gs, desc(value))
input_gs$Class <- factor(input_gs$Class, levels=input_gs$Class)

dev.off()
#plot
p1 <- ggplot(input_gs, aes(x = value, y = Class, fill = Class)) + 
  geom_bar(stat = "identity") + 
  ylab("Class (# of MAGs)") +
  xlab("Genome Size (MB)") +
  ggtitle("Average Genome Size per Class of SFB MAGs") +
  scale_x_continuous(limits = c(0,8), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(0,101), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) +
  guides(fill="none") +
  scale_fill_hue(l=50) +
  scale_y_discrete(limits=rev)
p1

ggsave("output/GenomeSize_SFB_MAGs.pdf", p1, dpi = 500)
ggsave("output/GenomeSize_SFB_MAGs.png", p1, dpi = 500)





