####################### Hierarchical Clustering Analysis ################################

library(vegan)
library(pvclust)
library(fpc)
library(ggplot2)
library(gridExtra)
library(dendextend)

##################### input #####################
genus <- read.table("../Abundance_Barplot/Input/SFBMAGS-vs-Reads-CoverM.tsv", header=TRUE, sep="\t", row.names = 1)

##################### input clean up #####################

#remove string from col names
names(genus) <- str_replace_all(names(genus), "_filter.METAGENOME.fastq.gz.Relative.Abundance....", "")
names(genus) <- str_replace_all(names(genus), "_sed_USGS", "")

#remove unmapped row
genus <- genus[-1, ]

#transpose
genust<-as.data.frame(t(genus))

##################### cluster #####################

#generate cluster dendrogram
dist.mat<-vegdist(genust, method="bray")
clust.res<-hclust(dist.mat)
dev.off()
plot(clust.res, hang=-1, main="Bray Curtis SFB MAG Abundance")
rect.hclust((clust.res), k=3)

#different clustering method
distanG<-dist(genust,method="maximum")
#distanceP<-dist(pfamt,method="maximum")
dendG<-as.dendrogram(hclust(distanG, method="ward.D"))
#dendP<-as.dendrogram(hclust(distanceP, method="ward.D"))
clustG<-hclust(distanG, method = "ward.D")
#clustP<-hclust(distanceP,method = "ward.D")

#pdf("clustering.ward.maximum.pdf")
#par(mfrow = c(2,1))
dev.off()
plot(hclust(distanG), hang=-1, main="Bray Curtis SFB MAG Abundance" )
rect.hclust((clustG), h=1.5)
rect.hclust((clustG), k = 6)
# plot(dendP, hang=-1, main=" Pfam" )
# rect.hclust((clustP), k=4)

##################### testing clusters #####################

#testing with distanG, Ward method
kbest.p<-2
cboot.hclust<-clusterboot(distanG,clustermethod=hclustCBI,
                          method="ward.D", k=kbest.p, B=1000, 
                          distances=TRUE)
groups<-cboot.hclust$result$partition

groups
mean<-cboot.hclust$bootmean
stab<-cboot.hclust$bootbrd
rec<-cboot.hclust$bootrecover
print(c( mean, stab,  rec))

#testing with dist.mat, BC method
kbest.p<-2
cboot.hclust<-clusterboot(dist.mat,clustermethod=hclustCBI,
                          method="ward.D", k=kbest.p, B=1000, 
                          distances=TRUE)
groups<-cboot.hclust$result$partition

groups
mean<-cboot.hclust$bootmean
stab<-cboot.hclust$bootbrd
rec<-cboot.hclust$bootrecover
print(c( mean, stab, rec))

##################### plotting data #####################

#par(mfrow = c(2,1), mar = c(5,2,1,0))
dendG <- dendG %>%
  color_branches(k = 2) %>%
  set("branches_lwd", c(2,1,2)) %>%
  set("branches_lty", c(1,2,1))

dendG <- color_labels(dendG, k = 2)
# The same as:
# labels_colors(dend)  <- get_leaves_branches_col(dend)
dev.off()
par(mar = c(7, 4, 4, 2) + 3)
plot(dendG)

##################### plotting data with legend #####################

#Get the groups 
#First according to site 
s_type <- rep("Other",length(rownames(genust)))
is_x <- grepl("SF_Jan12_13", rownames(genust))
s_type[is_x] <- "Site A"
is_x <- grepl("SF_Jul11_13", rownames(genust))
s_type[is_x] <- "Site A"
is_x <- grepl("SF_May12_13", rownames(genust))
s_type[is_x] <- "Site A"
is_x <- grepl("SF_Oct11_13", rownames(genust))
s_type[is_x] <- "Site A"
is_x <- grepl("SF_Jan12_21", rownames(genust))
s_type[is_x] <- "Site B"
is_x <- grepl("SF_Jul11_21", rownames(genust))
s_type[is_x] <- "Site B"
is_x <- grepl("SF_May12_21", rownames(genust))
s_type[is_x] <- "Site B"
is_x <- grepl("SF_Oct11_21", rownames(genust))
s_type[is_x] <- "Site B"
is_x <- grepl("SF_Jan12_24", rownames(genust))
s_type[is_x] <- "Site C"
is_x <- grepl("SF_Jul11_24", rownames(genust))
s_type[is_x] <- "Site C"
is_x <- grepl("SF_May12_24", rownames(genust))
s_type[is_x] <- "Site C"
is_x <- grepl("SF_Oct11_24", rownames(genust))
s_type[is_x] <- "Site C"
is_x <- grepl("SF_Jan12_4_1", rownames(genust))
s_type[is_x] <- "Site D"
is_x <- grepl("SF_Jul11_4_1", rownames(genust))
s_type[is_x] <- "Site D"
is_x <- grepl("SF_May12_4_1", rownames(genust))
s_type[is_x] <- "Site D"
is_x <- grepl("SF_Oct11_4_1", rownames(genust))
s_type[is_x] <- "Site D"
is_x <- grepl("SF_Jan12_8_1", rownames(genust))
s_type[is_x] <- "Site E"
is_x <- grepl("SF_Jul11_8_1", rownames(genust))
s_type[is_x] <- "Site E"
is_x <- grepl("SF_May12_8_1", rownames(genust))
s_type[is_x] <- "Site E"
is_x <- grepl("SF_Oct11_8_1", rownames(genust))
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
is_x <- grepl("SF_Jan12_13", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Jul11_13", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_May12_13", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Oct11_13", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Jan12_21", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Jul11_21", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_May12_21", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Oct11_21", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Jan12_24", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Jul11_24", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_May12_24", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Oct11_24", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Jan12_4_1", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Jul11_4_1", rownames(genust))
cb_type[is_x] <- "Group 1"
is_x <- grepl("SF_May12_4_1", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Oct11_4_1", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Jan12_8_1", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Jul11_8_1", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_May12_8_1", rownames(genust))
cb_type[is_x] <- "Group 2"
is_x <- grepl("SF_Oct11_8_1", rownames(genust))
cb_type[is_x] <- "Group 2"

cb_type<-factor(cb_type)
nb_types <- length(unique(cb_type))
colsnb_4 <- colorspace::diverge_hcl(nb_types, c = 100, l  = 85)


colnb_type <- colsnb_4[cb_type]

### plots
#par(mar = c(12,4,1,1))

dev.off()
par(mar = c(7, 4, 4, 2) + 3)
plot(dendG, main="SFB")
colored_bars(cbind(k234[,3:1], colnb_type,col_type), dendG, rowLabels = c(paste0("k = ", 4:2), "Group" ,"Site"))

############################## Metabolism-based clustering ####################################

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
img_data <- img_data %>%
  select(c("PFAM_ID", "Bin"))

#modify for parsing
img_data <- img_data %>%
  mutate_all(na_if,"") %>% #add NA for absence of pfam
  separate(Bin, into = c("Sample", NA), sep = "_SF") #get rid of bin name to just keep sample name















