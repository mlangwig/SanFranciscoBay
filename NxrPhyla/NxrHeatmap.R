######################## Script adapted from J.D. Carlton, Baker Lab, UT Austin #######################################

library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(tidyverse)
library(magick)

#df <- read.csv("Input/input_heatmap.csv")
#head(df)

#df_num <- df %>% select(nxrA:mct)
#rownames(df_num) = df$MAG

######################### Obtain all annotations in MAGs from IMG data ####################################

# using the IMG tables of bins only annotations to get 1 big table

# Specify the directory where the files are located
directory <- "../../IMG/IMG_Tables_BinsOnly/"

# Get a list of all the files ending with tsv.BinsOnly
file_list <- list.files(path = directory, pattern = "tsv.BinsOnly$", full.names = TRUE)

# Read all the files and combine them into one data frame
img_data <- lapply(file_list, read.delim) %>%
  bind_rows()

# View the first few rows of the combined table
head(img_data)
tail(img_data)

# Read in KO IDs
ids <- read.delim(file = "Input/markergenes.txt", header = TRUE)
# Read in MAG list
mags <- read.delim(file = "Input/MAG_list.txt", header = FALSE)

#modify for parsing
img_data <- img_data %>%
  #mutate_all(na_if,"") %>% #add NA for absence of pfam
  mutate(Sample = Bin) %>%
  separate(Sample, into = c("Sample", NA), sep = "_SF") #%>% #get rid of bin name to just keep sample name
  #drop_na() #remove rows with NAs

# subset for just the cols you want
#img_data_sub <- img_data %>%
#  select(c("PFAM_ID", "Bin"))

# Make gene names more standardized
img_data <- img_data %>%
  mutate(KO_Term = gsub("KO:", "", KO_Term)) %>%
  mutate(PFAM_ID = gsub("pfam", "PF", PFAM_ID))

#get all IMG annotations for MAGs
img_data_large_mags <- img_data %>%
  filter(Bin %in% mags$V1)

# See all annotations for only nxr MAGs
#img_data_sub <- img_data %>%
#  filter(Bin %in% mags$V1)

# Filter IMG data for KOs
img_data_sub <- img_data %>%
  filter(PFAM_ID %in% ids$KO | KO_Term %in% ids$KO)

#Check counts
length(unique(ids$KO))
length(unique(img_data_sub$KO_Term))
# Both 140

# Filter for only nxr MAGs
img_data_sub <- img_data_sub %>%
  filter(Bin %in% mags$V1)

#Check counts
length(unique(mags$V1))
length(unique(img_data_sub$Bin))
# Both 39

# Add pathway metadata
img_data_sub <- ids %>%
  right_join(img_data_sub, by = c("KO" = "KO_Term"))

# Drop some pathways
img_data_sub <- img_data_sub %>%
  filter(!Cat_Broad %in% c("PPP", "Glycolysis", "WLP", "3HP"))

# K00374, K00855, and K01602 missing from MAGs

##################################################################

######################### Obtain all annotations in Cultured Refs from KEGG data ####################################

# Specify the directory where the files are located
directory2 <- "../../NarG_NxrA_phylogeny/Figure4/KEGG_output_culturedNOB/txt/"

# Get a list of all the files
file_list2 <- list.files(path = directory2, pattern = "^mod", full.names = TRUE)

# Read all the files and combine them into one data frame
refs <- lapply(file_list2, read.delim) %>%
  bind_rows()

# View the first few rows of the combined table
head(refs)
tail(refs)

#remove ko: from KEGG column data
refs <- refs %>%
  mutate(KEGG_ko = gsub("ko:", "", KEGG_ko)) %>%
  filter(!is.na(evalue))

# Add NCBI ID to MAGs
ref_mapping <- read.delim(file = "Input/RefGenomeIDs.txt", header = TRUE)
refs <- refs %>%
  mutate(NCBI_ID = query) %>%
  extract(NCBI_ID, into = "NCBI_ID", regex = "^([^_]+_[^_]+)", remove = TRUE)

refs <- ref_mapping %>%
  right_join(refs, by = c("Descriptive" = "NCBI_ID"))

# Filter refs for genes of interest

# collapse all KOs into one regex pattern like "K00370|K07306|K17050|..."
pattern <- paste0(ids$KO, collapse = "|")

refs_sub <- refs %>%
  filter(str_detect(KEGG_ko, pattern))

#Check counts
length(unique(refs_sub$KEGG_ko))
length(unique(ids$KO))
# References are missing 2 KOs: K00374 + K14470

######################### Create input matrix ####################################

###### MAGs
# get col of gene name
img_data_sub <- img_data_sub %>%
  separate(Name, into = c("Gene", "Gene_Desc"), sep = ";")

#get just info of interest for matrix
img_data_mags <- img_data_sub %>%
  dplyr::select(c("Gene", "Bin", "Order", "Category", "Cat_Broad"))

###### Refs
# do the same for the refs
met_data_refs <- refs_sub %>%
  dplyr::select(c("NCBI_ID", "KEGG_ko"))
# split the weird combo KO cols
met_data_refs <- met_data_refs %>%
  separate_rows(KEGG_ko, sep = ",")

# map metadata to ref KOs so I know which ones to drop
met_data_refs <- ids %>%
  right_join(met_data_refs, by = c("KO" = "KEGG_ko"))
# Drop some pathways
met_data_refs <- met_data_refs %>%
  filter(!Cat_Broad %in% c("PPP", "Glycolysis", "WLP", "3HP"))
# Drop KOs not of interest
met_data_refs <- met_data_refs %>%
  filter(!is.na(Cat_Broad))
# get col of gene name
met_data_refs <- met_data_refs %>%
  separate(Name, into = c("Gene", "Gene_Desc"), sep = ";")
# select only cols of interest
met_data_refs <- met_data_refs %>%
  dplyr::select(c("Gene", "NCBI_ID", "Order", "Category", "Cat_Broad"))
# change col name of ref table so can merge MAG and ref info
met_data_refs <- met_data_refs %>%
  rename("Bin" = "NCBI_ID")

# Combine ref and MAG data
mags_refs_metabolism <- rbind(img_data_mags, met_data_refs)

# Based on my phylogeny, rbcL in the SFB MAGs are all type IV aka not for C fixation. Clarifying that here:
mags_refs_metabolism <- mags_refs_metabolism %>%
  mutate(Gene = if_else(Gene == "rbcL" & grepl("Bin", Bin), "rbcL_typeIV", Gene))

# Create df of gene order
gene_order <- mags_refs_metabolism %>%
  select(Gene, Order) %>%
  distinct() %>%
  arrange(Order)

# Count occurrences per Gene-Bin
mags_refs_met_matrix <- mags_refs_metabolism %>%
  count(Bin, Gene) %>%
  pivot_wider(names_from = Gene, values_from = n, values_fill = 0)

# Set rownames to Bin
mags_refs_met_matrix <- mags_refs_met_matrix %>% 
  column_to_rownames("Bin")
# Reorder columns by gene_order
mags_refs_met_matrix <- mags_refs_met_matrix[, gene_order$Gene]

# Reorder rows by tax
mag_order <- read.delim(file = "Input/MAG_order.txt", header = TRUE)
# Reorder columns by gene_order
mag_order <- mag_order%>%
  select(-c("Tax"))
mags_refs_met_matrix <- mags_refs_met_matrix[mag_order$MAG,]

######################### Create multiple input matrices by pathway ####################################

category_info <- mags_refs_metabolism %>%
  select(Gene, Category, Cat_Broad) %>%
  distinct() %>%
  left_join(gene_order, by = "Gene") %>%
  arrange(Order)

#replace all gene presence values greater than 1 with 2
mags_refs_met_matrix_scale<-replace(mags_refs_met_matrix, mags_refs_met_matrix>1,2)

nitrogen <- as.matrix(t(mags_refs_met_matrix_scale[c(1:12)]))
other <- as.matrix(t(mags_refs_met_matrix_scale[c(13:17)]))
ETC <- as.matrix(t(mags_refs_met_matrix_scale[c(18:62)]))
#glyc <- as.matrix(t(mags_refs_met_matrix_scale[c(65:87)]))
#PPP <- as.matrix(t(mags_refs_met_matrix_scale[c(88:96)]))
TCA <- as.matrix(t(mags_refs_met_matrix_scale[c(63:75)]))
rTCA <- as.matrix(t(mags_refs_met_matrix_scale[c(76:84)]))
#WLP <- as.matrix(t(mags_refs_met_matrix_scale[c(120:127)]))
rGlyp <- as.matrix(t(mags_refs_met_matrix_scale[c(85:95)]))
CBB <- as.matrix(t(mags_refs_met_matrix_scale[c(96:99)]))
#THP <- as.matrix(t(mags_refs_met_matrix_scale[c(139)]))

#Heatmap for the full table
# pheatmap(df_scale,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 7, cellheight = 7,
#          gaps_row = c(2,9,10,11,12), gaps_col = c(12,17,62,85,94,107,112,117),
# fontsize = 7, color = colorRampPalette(c("white", "blue", "navy"))(50))

#heatmap by group to get distinct colors
nhm <- pheatmap::pheatmap(nitrogen,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 8, cellheight = 8, fontsize = 8,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#029E73"))(50), filename = "Output/nhm.png")
ohm <- pheatmap::pheatmap(other,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 8, cellheight = 8, fontsize = 8,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#f46a9b"))(50), filename = "Output/ohm.png")
etchm <- pheatmap::pheatmap(ETC,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 8, cellheight = 8, fontsize = 8,
                            show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                            color = colorRampPalette(c("white", "#7eb0d5"))(50), filename = "Output/etchm.png")
# ghm <- pheatmap::pheatmap(glyc,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
#                           show_rownames = 0, legend = FALSE, gaps_row = c(2,9,10,11,12),
#                           color = colorRampPalette(c("white", "purple"))(50))
# phm <- pheatmap::pheatmap(PPP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
#                           show_rownames = 0, legend = FALSE, gaps_row = c(2,9,10,11,12),
#                           color = colorRampPalette(c("white", "#FF4500"))(50))
thm <- pheatmap::pheatmap(TCA,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 8, cellheight = 8, fontsize = 8,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#989898"))(50), filename = "Output/thm.png")
rtchm <- pheatmap::pheatmap(rTCA,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 8, cellheight = 8, fontsize = 8,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#ffb55a"))(50), filename = "Output/rtchm.png")
# wlphm <- pheatmap::pheatmap(WLP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
#                           show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
#                           color = colorRampPalette(c("white", "#bd7ebe"))(50)) 
rglhm <- pheatmap::pheatmap(rGlyp,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 8, cellheight = 8, fontsize = 8,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#b2e061"))(50), filename = "Output/rglhm.png")
cbbhm <- pheatmap::pheatmap(CBB,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 8, cellheight = 8, fontsize = 8,
                            show_colnames = 1, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                            color = colorRampPalette(c("white", "#ffee65"))(50), filename = "Output/cbbhm.png")
# thphm <- pheatmap::pheatmap(THP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
#                             show_colnames = 1, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
#                             color = colorRampPalette(c("white", "#7fcdbb"))(50))

# test magick
imgs <- image_read(c("Output/nhm.png", "Output/ohm.png", "Output/etchm.png", "Output/thm.png",
                     "Output/rtchm.png", "Output/rglhm.png", "Output/cbbhm.png"))
stitched <- image_append(imgs, stack = TRUE)
image_write(stitched, "Output/combined_heatmap.png")

# #arrange heatmaps
# lm <- cbind(c(1,2,3,6,7,8,9,10,11))
# dev.off()
# g <- grid.arrange(grobs = list(nhm[[4]],
#                           ohm[[4]],
#                           etchm[[4]],
#                           # ghm[[4]],
#                           # phm[[4]],
#                           thm[[4]],
#                           rtchm[[4]],
#                           #wlphm[[4]],
#                           rglhm[[4]],
#                           cbbhm[[4]]),
#                           #thphm[[4]]),
#                           ncol = 1,
#              heights = c(1,1,4.2,1.7,1,1,1,1,1), # I try to arrange heights here so you can see the plot in plot zoom...but doesn't seem to translate well to ggsave
#              layout_matrix = lm)
# 
# #ggsave(filename = "Output/heatmap.png", g, dpi = 500, height = 12) #here I needed to establish height to make it non overlapping
# #ggsave(filename = "Output/heatmap.svg", g, dpi = 500, height = 12)
# 
# ggsave(filename = "Output/test_heatmap.png", g, dpi = 500, height = 12) #here I needed to establish height to make it non overlapping
# 
# 
# 
# 





