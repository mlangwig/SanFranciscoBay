############### Plotting sfb misc metadata #####################

#libraries
library(dplyr)
library(tidyverse)
library(wesanderson)
library(ggplot2)

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
        panel.border = element_rect(colour = "grey", fill=NA)) +
  guides(fill="none") +
  scale_fill_viridis_d() +
  scale_y_discrete(limits=rev)
p1

ggsave("output/AvgGenomeSize_SFB_MAGs.pdf", p1, dpi = 500, width = 8, height = 8)
ggsave("output/AvgGenomeSize_SFB_MAGs.png", p1, dpi = 500, width = 8, height = 8)


############### GC content #####################
#read input
gc <- read.delim2(file = "input/gc_out.txt", sep = "\t", header = TRUE)
nxra_mags <- read.delim2(file = "input/nxrA_mags.txt", sep = "\t", header = FALSE)
nxra_scafs <- read.delim2(file = "input/nxrA_scaffs_fna.txt", sep = "\t", header = FALSE)
nxra_new <- read.delim2(file = "input/nxrA_mags_unknown.txt", sep = "\t", header = FALSE)

#copy column for splitting
gc$Bin = gc$ID
gc <- gc %>%
  mutate(Bin = str_replace_all(Bin, ">", ""))
#split so have a MAG name column
gc <- gc %>% separate(Bin, c("Bin", NA), sep= "(?=_scaf)")
#rename weird col name
gc <- rename(gc,"PercGC" = "X..GCContent")

#average GC per scaffold per bin
gc_sum_perBin <- gc %>%
  group_by(Bin) %>% 
  summarise(value=mean(as.numeric(PercGC))) %>%
  ungroup()

#plot for all scaffolds to visualize
nxra_bin_vec <- unique(nxra_mags$V1)
gc_nxrA <- gc %>% filter(Bin %in% nxra_bin_vec)

#add mapping for highlight of specific genes
gc_nxrA$gene_color <- rep('grey', nrow(gc_nxrA))
nxra_scafs <- unique(nxra_scafs$V1)
gc_nxrA <- gc_nxrA %>%
  mutate(ID = str_replace_all(ID, ">", ""))
gc_nxrA <- gc_nxrA %>%
  mutate(gene_color = ifelse(ID %in% nxra_scafs, "red", gene_color))

#subset for just scaffold that nxrA is on
nxra_new <- unique(nxra_new$V1)
gc_nxrA_new <- gc_nxrA %>% filter(Bin %in% nxra_new)

## ggplot
dev.off()
plot <- gc_nxrA_new %>%
  ggplot(aes(x = ID, y = as.numeric(PercGC), group = 1)) + 
  geom_line(color = "grey", size = .3) +
  geom_point(size =.5, alpha = .5, color = gc_nxrA_new$gene_color) +
  labs(x = "Scaffold", y = "% GC") +
  guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.background = element_rect(color = "white"),
        legend.box.background = element_rect(fill = "transparent"),
        #panel.background = element_rect(fill = "transparent"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #panel.border = element_blank()) + #turn this off to get the outline back)
  scale_y_continuous(expand = c(0, 0)) + #turn this on to make it look aligned with ticks
  ggtitle("% GC nxrA MAGs") + #Change for top X grabbed
  facet_wrap(.~Bin, scales = "free") 
#scale_y_continuous(breaks = seq(0, 45, by=5))
#coord_flip()
plot

ggsave("output/gc_nxrA_new_phyla.png", plot, dpi = 500, width = 8)

############### Functional Redundancy of Metabolism #####################
library(devtools)
library(readr)
library(dplyr)
library(tidyverse)
library(vegan)
install_github("adsteen/funfunfun")

#read abundance input
abundance <- read_tsv("../Abundance_Barplot/Input/SFBMAGS-vs-Reads-CoverM.tsv")
#make first col rownames
abundance <- abundance %>%
  column_to_rownames(var = 'Genome')
#flip abundance for input needs of funfunfun
abundance <- t(abundance)
#drop unmapped
abundance <- abundance[ , -1]

#normalize the abundance following the example of funfunfun - dividing the values by the smallest value greater than 0
abundance_norm <- as.data.frame(round(abundance/min(abundance[abundance>0])))
#create sample effort based on example
sample_effort <- min(rowSums(abundance_norm))
#rarefy the abundance matrix based on example
abundance_norm_rare <- vegan::rrarefy(abundance_norm,sample_effort) # not clear about how this is calculated?
#sweep the data, divide the values by row wise sums
abundance_norm_rare_sweep <- sweep(abundance_norm_rare, 1, rowSums(abundance_norm_rare), '/')
#order the MAG columns alphabetically
abundance_norm_rare_sweep <- abundance_norm_rare_sweep[, order(colnames(abundance_norm_rare_sweep))]

#read trait input
traits <- read_tsv("input/METABOLIC_subset_fr.txt", col_names = TRUE)
traits <- traits[-1, ]
traits <- data.frame(traits, row.names = 1)
traits <- as.matrix(traits)
traits <- matrix(as.numeric(traits), nrow = nrow(traits), ncol = ncol(traits), dimnames = dimnames(traits))
#order the MAG rows alphabetically to match the abundance matrix
traits <- traits[order(rownames(traits)), ]
#calculate the functional redundancy metric
fr_sfb <- funfunfun::royalty_fr(abundance_norm_rare_sweep, traits, q = 0.5)
#change metadata columns
fr_sfb <- tidyr::separate(fr_sfb,sample,into = c("site","meta"), sep = '_filter')
fr_sfb <- fr_sfb %>%
  separate(site, into = c("month", "site"), sep = '_sed')
#remove blank
  

#plot
dev.off()
ggplot2::ggplot(fr_sfb,aes(x=trait,y=fr,color=site)) +
  geom_boxplot() +
  ylim(0,0.7) + #,color=sample
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

### functional redundancy example using funfunfun
abundance.matrix <- read.csv("https://raw.githubusercontent.com/adsteen/funfunfun/main/data/MAG_abundance_table.csv", row.names = 1 )
abundance.matrix.norm <- round(abundance.matrix/min(abundance.matrix[abundance.matrix>0])) 
sample.effort <- min(rowSums(abundance.matrix.norm))
abundance.matrix.norm.rare <- vegan::rrarefy(abundance.matrix.norm,sample.effort)
abundance.matrix.norm.rare.sweep <- sweep(abundance.matrix.norm.rare,1,rowSums(abundance.matrix.norm.rare),'/') #this is dividing the values by the row sum
#row_sum <- sum(df["Row2", ])


trait.matrix <- read.csv("https://raw.githubusercontent.com/adsteen/funfunfun/main/data/MAG_enzyme_gene_copies.csv", row.names = 1 )
#Notice that abundance matrix cols and trait rows are in same MAG order
fr <- funfunfun::royalty_fr(abundance.matrix, trait.matrix, q = 0.5)
fr <- tidyr::separate(fr,sample,into = c("site","size_fraction","depth"), sep = '_')
dev.off()
ggplot2::ggplot(fr,aes(x=trait,y=fr,color=depth))+geom_boxplot()+ylim(0,0.7)  

############### Functional Redundancy of Metabolism based on IMG data rather than METABOLIC #####################

#read trait input
traits_img <- readxl::read_excel("input/fr_genes_IMG_final.xlsx", sheet = "long_way")
#read in rare SFB MAGs
rare_MAGs <- read.delim(file = "../Abundance_Barplot/Output/MAGs_RareBiosphere.txt", header = FALSE)

#collapse so that presence >0 means path present
traits_img_sum <- traits_img %>%
  select(-`gene description`, -KO) %>%        # Drop non-numeric columns
  group_by(pathway) %>%                       # Group by pathway
  summarise(across(everything(),              # For each numeric column
                   ~ as.integer(any(. > 0)))) # 1 if any value > 0, else 0
#transpose
traits_img_sum <- t(traits_img_sum)
# Make first row the header
colnames(traits_img_sum) <- as.character(traits_img_sum[1, ])
# Get rid of first row
traits_img_sum <- traits_img_sum[-1, ]

#Make df a matrix and values num
traits_img_sum <- as.matrix(traits_img_sum)
traits_img_sum <- matrix(as.numeric(traits_img_sum), nrow = nrow(traits_img_sum), ncol = ncol(traits_img_sum), dimnames = dimnames(traits_img_sum))
#calculate the functional redundancy metric
fr_sfb <- funfunfun::royalty_fr(abundance_norm_rare_sweep, traits_img_sum, q = 0.5)

#remove any weight of abundance
#fr_sfb <- funfunfun::royalty_fr(abundance_norm_rare_sweep, traits_img_sum, q = 0.0)
#add extra weight to abundance
#fr_sfb <- funfunfun::royalty_fr(abundance_norm_rare_sweep, traits_img_sum, q = 1.0)

#change metadata columns
fr_sfb <- tidyr::separate(fr_sfb,sample,into = c("site","meta"), sep = '_filter')
fr_sfb <- fr_sfb %>%
  separate(site, into = c("month", "site"), sep = '_sed')

#plot
dev.off()
ggplot2::ggplot(fr_sfb,aes(x=trait,y=fr,color=site)) +
  geom_boxplot() +
  ylim(0,0.7) + #,color=sample
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

#view averages and sds
fr_sfb_stats <- fr_sfb %>%
  group_by(trait) %>%
  summarise(
    mean_fr = mean(fr) ,
    sd_fr = sd(fr)
  )

#view proportions of traits
# Number of rows
n_rows <- nrow(traits_img_sum)

# Proportion of 1s per column (trait)
trait_proportions <- colSums(traits_img_sum == 1) / n_rows * 100

# View the result
print(trait_proportions)

#view proportions of genes
#collapse so that presence >0 means path present
traits_img_genes <- traits_img %>%
  select(-pathway, -KO) %>%        # Drop non-numeric columns
  group_by(`gene description`) %>%                       # Group by pathway
  summarise(across(everything(),              # For each numeric column
                   ~ as.integer(any(. > 0)))) # 1 if any value > 0, else 0

#transpose
traits_img_genes <- t(traits_img_genes)
# Make first row the header
colnames(traits_img_genes) <- as.character(traits_img_genes[1, ])
# Get rid of first row
traits_img_genes <- traits_img_genes[-1, ]

#Make df a matrix and values num
traits_img_genes <- as.matrix(traits_img_genes)
traits_img_genes <- matrix(as.numeric(traits_img_genes), nrow = nrow(traits_img_genes), ncol = ncol(traits_img_genes), dimnames = dimnames(traits_img_genes))

# Number of rows
n_rows <- nrow(traits_img_genes)
# Proportion of 1s per column (trait)
gene_proportions <- colSums(traits_img_genes == 1) / n_rows * 100

# See proportion of pathways encoded by rare taxa
rare_ids <- rare_MAGs$V1
traits_img_rare <- traits_img_sum[rownames(traits_img_sum) %in% rare_ids, ]
rare_gene_perc <- colSums(traits_img_rare > 0, na.rm = TRUE) / nrow(traits_img_rare) * 100
rare_gene_perc

percent_encoded <- colMeans(traits_img_rare, na.rm = TRUE) * 100
#percent rare MAGs with nitrogen cycling genes
percent_rows_with_presence <- sum(rowSums(traits_img_rare[, 1:8] > 0, na.rm = TRUE) > 0) / nrow(traits_img_rare) * 100
#percent rare MAGs with nitrogen cycling genes
percent_rows_with_presence <- sum(rowSums(traits_img_rare[, 8:14] > 0, na.rm = TRUE) > 0) / nrow(traits_img_rare) * 100

#write rare MAGs traits table output
write_delim(as.data.frame(traits_img_rare), file = "output/traits_img_RareBiosphere.txt", quote = "none")
#### Unused

# #read trait input
# genes <- read_tsv("input/fr_genes_final.txt")
# #subset IMG data for genes of interest
# img_data_filt <- img_data %>%
#   filter(KO_Term %in% genes$KO)
# #map pathway name onto IMG data
# img_data_filt <- genes %>%
#   dplyr::select("KO", "Pathway2") %>%
#   right_join(img_data_filt, by = c("KO" = "Source"))

# # read in DRAM data
# DRAM <- readxl::read_excel("../../DRAM/metabolism_summary.xlsx", sheet = "Energy")
# # keep nitrogen and sulfur cycles
# DRAM <- DRAM %>%
#   filter(grepl("Nitrogen|Sulfur", header))
# # drop some modules because repeated or not of interest
# DRAM <- DRAM %>%
#   filter(!grepl("nitrite \\+ ammonia => nitrogen", module))
# DRAM <- DRAM %>%
#   filter(!grepl("Assimilatory nitrate reduction, nitrate => ammonia", module))
# 
# traits_dram <- read_tsv("input/fr_genes_DRAM.txt", col_names = FALSE)
# traits_dram <- as.data.frame(traits_dram)
# 
# # Combine row 1 and row 2, element-wise, separated by space
# combined_row <- paste(traits_dram[1, ], traits_dram[2, ], sep = " ")
# combined_row <- sub("\\[.*", "", combined_row)
# traits_dram[1, ] <- combined_row     # Overwrite row 1 with combined result
# # Make first row the header
# colnames(traits_dram) <- as.character(traits_dram[1, ])
# traits_dram <- traits_dram[-1, ]
# # Make first col row names
# traits_dram <- data.frame(traits_dram, row.names = 1)
# # Get rid of extra rows
# traits_dram <- traits_dram[-c(1:4), ]







