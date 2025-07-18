######################## Script adapted from J.D. Carlton, Baker Lab, UT Austin #######################################

library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)

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

# Filter IMG data for KOs
img_data_sub <- img_data %>%
  filter(PFAM_ID %in% ids$KO | KO_Term %in% ids$KO)

#Check counts
length(unique(ids$KO))
length(unique(img_data_sub$KO_Term))
# Both 139

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

# Drop PPP and Glycolysis
img_data_sub <- img_data_sub %>%
  filter(!Cat_Broad %in% c("PPP", "Glycolysis"))

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

##################################################################

#replace all gene presence values greater than 1 with 2
df_scale<-replace(df_num, df_num>1,2)

nitrogen <- as.matrix(t(df_scale[c(1:13)]))
other <- as.matrix(t(df_scale[c(14:18)]))
ETC <- as.matrix(t(df_scale[c(19:64)]))
glyc <- as.matrix(t(df_scale[c(65:87)]))
PPP <- as.matrix(t(df_scale[c(88:96)]))
TCA <- as.matrix(t(df_scale[c(97:109)]))
rTCA <- as.matrix(t(df_scale[c(110:119)]))
WLP <- as.matrix(t(df_scale[c(120:127)]))
rGlyp <- as.matrix(t(df_scale[c(128:134)]))
CBB <- as.matrix(t(df_scale[c(135:138)]))
THP <- as.matrix(t(df_scale[c(139)]))

#Heatmap for the full table
# pheatmap(df_scale,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 7, cellheight = 7,
#          gaps_row = c(2,9,10,11,12), gaps_col = c(12,17,62,85,94,107,112,117),
# fontsize = 7, color = colorRampPalette(c("white", "blue", "navy"))(50))

#heatmap by group to get distinct colors
nhm <- pheatmap::pheatmap(nitrogen,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#029E73"))(50))
ohm <- pheatmap::pheatmap(other,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#f46a9b"))(50))
etchm <- pheatmap::pheatmap(ETC,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                            show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                            color = colorRampPalette(c("white", "#7eb0d5"))(50))
# ghm <- pheatmap::pheatmap(glyc,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
#                           show_rownames = 0, legend = FALSE, gaps_row = c(2,9,10,11,12),
#                           color = colorRampPalette(c("white", "purple"))(50))
# phm <- pheatmap::pheatmap(PPP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
#                           show_rownames = 0, legend = FALSE, gaps_row = c(2,9,10,11,12),
#                           color = colorRampPalette(c("white", "#FF4500"))(50))
thm <- pheatmap::pheatmap(TCA,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#989898"))(50))
rtchm <- pheatmap::pheatmap(rTCA,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#ffb55a"))(50))
wlphm <- pheatmap::pheatmap(WLP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#bd7ebe"))(50)) 
rglhm <- pheatmap::pheatmap(rGlyp,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                          color = colorRampPalette(c("white", "#b2e061"))(50))
cbbhm <- pheatmap::pheatmap(CBB,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                            show_colnames = 0, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                            color = colorRampPalette(c("white", "#ffee65"))(50))
thphm <- pheatmap::pheatmap(THP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                            show_colnames = 1, legend = FALSE, gaps_col = c(1,3,10,11,13,14),
                            color = colorRampPalette(c("white", "#7fcdbb"))(50))
library(harrypotter)
pal <- hp(n = 8, house = "Gryffindor")
image(volcano, col = pal)

#arrange heatmaps
lm <- cbind(c(1,2,3,6,7,8,9,10,11))
dev.off()
g <- grid.arrange(grobs = list(nhm[[4]],
                          ohm[[4]],
                          etchm[[4]],
                          # ghm[[4]],
                          # phm[[4]],
                          thm[[4]],
                          rtchm[[4]],
                          wlphm[[4]],
                          rglhm[[4]],
                          cbbhm[[4]],
                          thphm[[4]]),
                          ncol = 1,
             heights = c(1,1,4.2,1.7,1,1,1,1,1), # I try to arrange heights here so you can see the plot in plot zoom...but doesn't seem to translate well to ggsave
             layout_matrix = lm)

ggsave(filename = "Output/heatmap.png", g, dpi = 500, height = 12) #here I needed to establish height to make it non overlapping
ggsave(filename = "Output/heatmap.svg", g, dpi = 500, height = 12)

