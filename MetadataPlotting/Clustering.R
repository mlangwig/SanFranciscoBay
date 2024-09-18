####################### Hierarchical Clustering Analysis ################################

library(vegan)
library(pvclust)
library(fpc)
library(ggplot2)
library(gridExtra)
library(dendextend)
library(dplyr)

##################### input #####################
genus <- read.table("../Abundance_Barplot/Input/SFBMAGS-vs-Reads-CoverM.tsv", header=TRUE, sep="\t", row.names = 1)
names <- read.table("input/site_name_change.txt", header = TRUE, sep="\t")

##################### input clean up #####################

#remove string from col names
names(genus) <- str_replace_all(names(genus), "_filter.METAGENOME.fastq.gz.Relative.Abundance....", "")
names(genus) <- str_replace_all(names(genus), "_sed_USGS", "")

#remove unmapped row
genus <- genus[-1, ]

#transpose
genust<-as.data.frame(t(genus))

#change the row names to be simplified
# Create a named vector where SiteComplex is the name and SiteSimple is the value
replacement_vector <- setNames(names$SiteSimple, names$SiteComplex)
# Replace the rownames of genust based on the replacement vector
rownames(genust) <- replacement_vector[rownames(genust)]



####################################### MAG abundance cluster #################################################

##################### Bray Curtis #####################

#generate cluster dendrogram
dist.mat<-vegdist(genust, method="bray")
#clust.res<-hclust(dist.mat)
dendG<-as.dendrogram(hclust(dist.mat))
clustG<-hclust(dist.mat, method = "ward.D")
dev.off()
#plot(clust.res, hang=-1, main="Bray Curtis SFB MAG Abundance")
plot(hclust(dist.mat), hang=-1, main="Bray Curtis SFB MAG Abundance" )
rect.hclust((clustG), k=3)

##################### Maximum #####################

# #different clustering method
# distanG<-dist(genust,method="maximum")
# #distanceP<-dist(pfamt,method="maximum")
# dendG<-as.dendrogram(hclust(distanG, method="ward.D"))
# #dendP<-as.dendrogram(hclust(distanceP, method="ward.D"))
# clustG<-hclust(distanG, method = "ward.D")
# #clustP<-hclust(distanceP,method = "ward.D")
# 
# #pdf("clustering.ward.maximum.pdf")
# #par(mfrow = c(2,1))
# dev.off()
# plot(hclust(distanG), hang=-1, main="Ward SFB MAG Abundance")
# rect.hclust((clustG), h=1.5)
# rect.hclust((clustG), k = 6)
# # plot(dendP, hang=-1, main=" Pfam" )
# # rect.hclust((clustP), k=4)


##################### testing MAG abundance clusters #####################

#testing with dist.mat, BC method
kbest.p<-3
cboot.hclust<-clusterboot(dist.mat,clustermethod=hclustCBI,
                          method="ward.D", k=kbest.p, B=1000, 
                          distances=TRUE)
groups<-cboot.hclust$result$partition

groups
mean<-cboot.hclust$bootmean
stab<-cboot.hclust$bootbrd
rec<-cboot.hclust$bootrecover
print(c( mean, stab, rec))

# #testing with distanG, Maximum method
# kbest.p<-2
# cboot.hclust<-clusterboot(distanG,clustermethod=hclustCBI,
#                           method="ward.D", k=kbest.p, B=1000, 
#                           distances=TRUE)
# groups<-cboot.hclust$result$partition
# 
# groups
# mean<-cboot.hclust$bootmean
# stab<-cboot.hclust$bootbrd
# rec<-cboot.hclust$bootrecover
# print(c( mean, stab,  rec))

##################### plotting data #####################

#par(mfrow = c(2,1), mar = c(5,2,1,0))
#dev.off()
dendG <- dendG %>%
  color_branches(k = 3) #%>%
  #set("branches_lwd", c(2,1,2)) #%>%
  #set("branches_lty", c(1,2,1))

dendG <- color_labels(dendG, k = 3)
# The same as:
# labels_colors(dend)  <- get_leaves_branches_col(dend)
dev.off()
par(mar = c(7, 4, 4, 2) + 3)
plot(dendG)

##################### plotting data with legend #####################

#Get the groups 
#First according to site 
s_type <- rep("Other",length(rownames(genust)))
is_x <- grepl("13_Jan", rownames(genust))
s_type[is_x] <- "Site A"
is_x <- grepl("13_July", rownames(genust))
s_type[is_x] <- "Site A"
is_x <- grepl("13_May", rownames(genust))
s_type[is_x] <- "Site A"
is_x <- grepl("13_Oct", rownames(genust))
s_type[is_x] <- "Site A"
is_x <- grepl("21_Jan", rownames(genust))
s_type[is_x] <- "Site B"
is_x <- grepl("21_July", rownames(genust))
s_type[is_x] <- "Site B"
is_x <- grepl("21_May", rownames(genust))
s_type[is_x] <- "Site B"
is_x <- grepl("21_Oct", rownames(genust))
s_type[is_x] <- "Site B"
is_x <- grepl("24_Jan", rownames(genust))
s_type[is_x] <- "Site C"
is_x <- grepl("24_July", rownames(genust))
s_type[is_x] <- "Site C"
is_x <- grepl("24_May", rownames(genust))
s_type[is_x] <- "Site C"
is_x <- grepl("24_Oct", rownames(genust))
s_type[is_x] <- "Site C"
is_x <- grepl("4_1_Jan", rownames(genust))
s_type[is_x] <- "Site D"
is_x <- grepl("4_1_July", rownames(genust))
s_type[is_x] <- "Site D"
is_x <- grepl("4_1_May", rownames(genust))
s_type[is_x] <- "Site D"
is_x <- grepl("4_1_Oct", rownames(genust))
s_type[is_x] <- "Site D"
is_x <- grepl("8_1_Jan", rownames(genust))
s_type[is_x] <- "Site E"
is_x <- grepl("8_1_July", rownames(genust))
s_type[is_x] <- "Site E"
is_x <- grepl("8_1_May", rownames(genust))
s_type[is_x] <- "Site E"
is_x <- grepl("8_1_Oct", rownames(genust))
s_type[is_x] <- "Site E"

s_type<-factor(s_type)

n_types <- length(unique(s_type))
cols_4 <- colorspace::rainbow_hcl(n_types, n=5, c = 120, l  = 60)
col_type <- cols_4[s_type]

#cutree 

k234<-cutree(dendG, k=2:4)
#color labels by site
labels_colors(dendG)<-col_type[order.dendrogram((dendG))]


#Get the groups 
#According to group detected by clusterboot 
cb_type <- rep("Other",length(rownames(genust)))
is_x <- grepl("13_Jan", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("13_July", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("13_May", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("13_Oct", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("21_Jan", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("21_July", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("21_May", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("21_Oct", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("24_Jan", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("24_July", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("24_May", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("24_Oct", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("4_1_Jan", rownames(genust))
cb_type[is_x] <- "Group 3"
is_x <- grepl("4_1_July", rownames(genust))
cb_type[is_x] <- "Group 3"
is_x <- grepl("4_1_May", rownames(genust))
cb_type[is_x] <- "Group 3"
is_x <- grepl("4_1_Oct", rownames(genust))
cb_type[is_x] <- "Group 3"
is_x <- grepl("8_1_Jan", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("8_1_July", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("8_1_May", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("8_1_Oct", rownames(genust))
cb_type[is_x] <- "Group 2"

cb_type<-factor(cb_type)
nb_types <- length(unique(cb_type))
colsnb_4 <- colorspace::diverge_hcl(nb_types, c = 100, l  = 85)


colnb_type <- colsnb_4[cb_type]

### plots
#par(mar = c(12,4,1,1))

# Specify the filename, width, and height of the PNG
png("output/dendrogram_sfb_abundance.png", width = 800, height = 400)

par(mar = c(7, 4, 4, 2) + 3)
plot(dendG, main="Bray Curtis SFB MAG Abundance")
colored_bars(cbind(k234[,3:1], colnb_type,col_type), dendG, rowLabels = c(paste0("k = ", 4:2), "Group" ,"Site"))

# Close the device to save the file
dev.off()

############################## Pfam-based clustering ####################################

############################## Creating the input ####################################

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

# subset for just the cols you want
img_data_sub <- img_data %>%
  select(c("PFAM_ID", "Bin"))

#modify for parsing
img_data_sub <- img_data_sub %>%
  mutate_all(na_if,"") %>% #add NA for absence of pfam
  separate(Bin, into = c("Sample", NA), sep = "_SF") %>% #get rid of bin name to just keep sample name
  drop_na() #remove rows with NAs

#pivot wider and sum pfams
img_data_sub_wide <- img_data_sub %>%
  count(Sample, PFAM_ID) %>%
  pivot_wider(names_from = PFAM_ID, values_from = n, values_fill = 0)

# Sample' is a non-numeric column, exclude it
row_sums <- rowSums(img_data_sub_wide[, sapply(img_data_sub_wide, is.numeric)])
# Print the row sums
print(row_sums)

pfam_normalized <- img_data_sub_wide %>%
  rowwise() %>%  # Apply the function row-wise
  mutate(total_pfam = sum(c_across(-Sample))) %>%  # Sum values across the row (ignoring Sample)
  mutate(across(-c(Sample, total_pfam), ~ . / total_pfam)) %>%  # Divide each value by the row sum
  ungroup() %>%  # Remove row-wise grouping
  select(-total_pfam)  # Optionally remove the total_pfam column
write_tsv(pfam_normalized, file = "output/pfam_normalized.tsv")

pfam_normalized <- pfam_normalized %>%
  column_to_rownames(var = "Sample")


############################## Bray Curtis clustering ####################################

#generate cluster dendrogram
dist.mat2<-vegdist(pfam_normalized, method="bray")
clust.res2<-hclust(dist.mat2)
dendG2<-as.dendrogram(hclust(dist.mat2))
clustG2<-hclust(dist.mat2, method = "ward.D")

dev.off()
plot(hclust(dist.mat2), hang=-1, main="Bray Curtis SFB MAG Pfams")
rect.hclust((clustG2), k=3)

#plot(clust.res2, hang=-1, main="Bray Curtis SFB Pfams")
#rect.hclust((clust.res2), k=3)

############################## Maximum clustering ####################################

#different clustering method
distanG<-dist(pfam_normalized,method="maximum")
#distanceP<-dist(pfamt,method="maximum")
dendG<-as.dendrogram(hclust(distanG, method="ward.D"))
#dendP<-as.dendrogram(hclust(distanceP, method="ward.D"))
clustG<-hclust(distanG, method = "ward.D")
#clustP<-hclust(distanceP,method = "ward.D")

#pdf("clustering.ward.maximum.pdf")
#par(mfrow = c(2,1))
dev.off()
plot(hclust(distanG), hang=-1, main="Ward SFB Pfams" )
rect.hclust((clustG), h=40000)
rect.hclust((clustG), k = 4)
# plot(dendP, hang=-1, main=" Pfam" )
# rect.hclust((clustP), k=4)


##################### Testing MAG pfam clusters #####################

#testing with dist.mat, BC method
kbest.p<-3
cboot.hclust<-clusterboot(dist.mat2,clustermethod=hclustCBI,
                          method="ward.D", k=kbest.p, B=1000, 
                          distances=TRUE)
groups<-cboot.hclust$result$partition

groups
mean<-cboot.hclust$bootmean
stab<-cboot.hclust$bootbrd
rec<-cboot.hclust$bootrecover
print(c( mean, stab, rec))

# #testing with distanG, Ward method
# kbest.p<-4
# cboot.hclust<-clusterboot(distanG,clustermethod=hclustCBI,
#                           method="ward.D", k=kbest.p, B=1000, 
#                           distances=TRUE)
# groups<-cboot.hclust$result$partition
# 
# groups
# mean<-cboot.hclust$bootmean
# stab<-cboot.hclust$bootbrd
# rec<-cboot.hclust$bootrecover
# print(c( mean, stab,  rec))

##################### plotting data #####################

#par(mfrow = c(2,1), mar = c(5,2,1,0))
#dev.off()
dendG2 <- dendG2 %>%
  color_branches(k = 3) #%>%
#set("branches_lwd", c(2,1,2)) #%>%
#set("branches_lty", c(1,2,1))

dendG2 <- color_labels(dendG2, k = 3)
# The same as:
# labels_colors(dend)  <- get_leaves_branches_col(dend)
dev.off()
par(mar = c(7, 4, 4, 2) + 3)
plot(dendG2)


##################### plotting data with legend #####################

#Get the groups 
#First according to site 
s_type <- rep("Other",length(rownames(pfam_normalized)))
is_x <- grepl("13_Jan", rownames(pfam_normalized))
s_type[is_x] <- "Site A"
is_x <- grepl("13_July", rownames(pfam_normalized))
s_type[is_x] <- "Site A"
is_x <- grepl("13_May", rownames(pfam_normalized))
s_type[is_x] <- "Site A"
is_x <- grepl("13_Oct", rownames(pfam_normalized))
s_type[is_x] <- "Site A"
is_x <- grepl("21_Jan", rownames(pfam_normalized))
s_type[is_x] <- "Site B"
is_x <- grepl("21_July", rownames(pfam_normalized))
s_type[is_x] <- "Site B"
is_x <- grepl("21_May", rownames(pfam_normalized))
s_type[is_x] <- "Site B"
is_x <- grepl("21_Oct", rownames(pfam_normalized))
s_type[is_x] <- "Site B"
is_x <- grepl("24_Jan", rownames(pfam_normalized))
s_type[is_x] <- "Site C"
is_x <- grepl("24_July", rownames(pfam_normalized))
s_type[is_x] <- "Site C"
is_x <- grepl("24_May", rownames(pfam_normalized))
s_type[is_x] <- "Site C"
is_x <- grepl("24_Oct", rownames(pfam_normalized))
s_type[is_x] <- "Site C"
is_x <- grepl("4_1_Jan", rownames(pfam_normalized))
s_type[is_x] <- "Site D"
is_x <- grepl("4_1_July", rownames(pfam_normalized))
s_type[is_x] <- "Site D"
is_x <- grepl("4_1_May", rownames(pfam_normalized))
s_type[is_x] <- "Site D"
is_x <- grepl("4_1_Oct", rownames(pfam_normalized))
s_type[is_x] <- "Site D"
is_x <- grepl("8_1_Jan", rownames(pfam_normalized))
s_type[is_x] <- "Site E"
is_x <- grepl("8_1_July", rownames(pfam_normalized))
s_type[is_x] <- "Site E"
is_x <- grepl("8_1_May", rownames(pfam_normalized))
s_type[is_x] <- "Site E"
is_x <- grepl("8_1_Oct", rownames(pfam_normalized))
s_type[is_x] <- "Site E"

s_type<-factor(s_type)

n_types <- length(unique(s_type))
cols_4 <- colorspace::rainbow_hcl(n_types, n=5, c = 120, l  = 60)
col_type <- cols_4[s_type]

#cutree 

k234<-cutree(dendG2, k=2:4)
#color labels by site
labels_colors(dendG2)<-col_type[order.dendrogram((dendG2))]


#Get the groups 
#According to group detected by clusterboot 
cb_type <- rep("Other",length(rownames(pfam_normalized)))
is_x <- grepl("13_Jan", rownames(pfam_normalized))
cb_type[is_x] <- "Group 1"
is_x <- grepl("13_July", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"
is_x <- grepl("13_May", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"
is_x <- grepl("13_Oct", rownames(pfam_normalized))
cb_type[is_x] <- "Group 3"
is_x <- grepl("21_Jan", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"
is_x <- grepl("21_July", rownames(pfam_normalized))
cb_type[is_x] <- "Group 3"
is_x <- grepl("21_May", rownames(pfam_normalized))
cb_type[is_x] <- "Group 1"
is_x <- grepl("21_Oct", rownames(pfam_normalized))
cb_type[is_x] <- "Group 3"
is_x <- grepl("24_Jan", rownames(pfam_normalized))
cb_type[is_x] <- "Group 3"
is_x <- grepl("24_July", rownames(pfam_normalized))
cb_type[is_x] <- "Group 3"
is_x <- grepl("24_May", rownames(pfam_normalized))
cb_type[is_x] <- "Group 3"
is_x <- grepl("24_Oct", rownames(pfam_normalized))
cb_type[is_x] <- "Group 1"
is_x <- grepl("4_1_Jan", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"
is_x <- grepl("4_1_July", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"
is_x <- grepl("4_1_May", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"
is_x <- grepl("4_1_Oct", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"
is_x <- grepl("8_1_Jan", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"
is_x <- grepl("8_1_July", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"
is_x <- grepl("8_1_May", rownames(pfam_normalized))
cb_type[is_x] <- "Group 3"
is_x <- grepl("8_1_Oct", rownames(pfam_normalized))
cb_type[is_x] <- "Group 2"

cb_type<-factor(cb_type)
nb_types <- length(unique(cb_type))
colsnb_4 <- colorspace::diverge_hcl(nb_types, c = 100, l  = 85)


colnb_type <- colsnb_4[cb_type]

### plots
#par(mar = c(12,4,1,1))

# Specify the filename, width, and height of the PNG
png("output/dendrogram_sfb_pfam.png", width = 800, height = 400)

par(mar = c(7, 4, 4, 2) + 3)
plot(dendG2, main="Bray Curtis SFB MAG Pfams")
colored_bars(cbind(k234[,3:1], colnb_type,col_type), dendG2, rowLabels = c(paste0("k = ", 4:2), "Group" ,"Site"))

# Close the device to save the file
dev.off()


