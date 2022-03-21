library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)
library(stringr)

#read data from IMG tables (see readme for how this data was obtained)
data <- read.delim("../data/all_SFsamples_sulfur_binsOnly.header.tsv", sep = "\t", header = TRUE)

#read input from HMM search of custom sulfur database (see readme for how results were obtained)
data_hmm <- read.delim("../data/all_sulfur_hmmOutput_SFB.tsv")
#I forgot to remove the headers when concatenated the files so removing them here:
#This reduces data by 77. Makes sense because 1 is the header, and another is missing because no
#hits for soxC
data_hmm<-data_hmm[!grepl("protein", data_hmm$protein),]

#read in site mapping file
MAG_toSite<-read.delim("../data/MAGS_toSite.txt")

##NOTE TO SELF --> FIX 8_1_Jan_SF_Bin3 TO 8_1_Jan_SF_Bin3_1 IN HMM INPUT <-----------

#Also, note that this input has duplicate hits within a MAG aka it was sorted for best hits per 
#protein, NOT per MAG/genome. I'm keeping this because interested in the multi copy dsrCs in 
#some proteobacteria. We have a script to get best MAG hits though if needed.

#read in taxonomy data
tax<-read.delim("../../../SFB_tax_detailed.txt")

#length(unique(data$Product))

#################################  Sulfate reduction  ######################################

#Subset the data for just sulfate reduction genes
sulf_red <- data%>%filter(str_detect(KO_Term,"KO:K00958|KO:K00394|KO:K00395|KO:K11180|KO:K11181|KO:K07235|KO:K07236|KO:K07237|KO:K16885|KO:K16886|KO:K16887"))
#supplement above output with dsrC and dsrD because those are not included in IMG output
sulf_red_hmm <- data_hmm%>%filter(str_detect(accession, "dsrC_tusE|dsrE_protein.alignment"))

#separate bin and scaffold in hmm input
sulf_red_hmm<-separate(sulf_red_hmm,protein,into = c("Bin", "scaffold"), sep = "(?=_scaffold)")
#replace _scaffold to get just scaffold
sulf_red_hmm<- data.frame(lapply(sulf_red_hmm, function(x) {
  gsub("_scaffold", "scaffold", x)
}))

#generate column of 1s for value mapping
newcol<-rep(1,2143) #generate 65 1s
sulf_red$presence<-newcol

#generate same col for hmm data
newcol2<-rep(1,673) #generate 65 1s
sulf_red_hmm$presence<-newcol2

#sum the presence of marker genes by Bin name
sulf_red_sum <- sulf_red %>%
  group_by(Bin, KO_Term) %>% 
  summarise(presence=sum(presence))

sulf_red_hmm_sum <- sulf_red_hmm %>%
  group_by(Bin, accession) %>% 
  summarise(presence=sum(presence))

#concatenate sulf_red_sum and sulf_red_hmm_sum so can pivot wider
colnames(sulf_red_sum)[2]<-"accession"
sulf_red_KOhmm<-rbind(sulf_red_sum,sulf_red_hmm_sum)

#go from long to wide
sulf_red_KOhmm_wide <- sulf_red_KOhmm %>%
  dplyr::select(accession, presence, Bin) %>%
  pivot_wider(names_from = accession, values_from = presence, values_fill = 0)

#map taxonomy
sulf_red_KOhmm_wide <- tax %>%
  dplyr::select(Bin_name, GTDBtk) %>%
  right_join(sulf_red_KOhmm_wide, by = c("Bin_name" = "Bin"))

#rename to add gene name to KO
sulf_red_KOhmm_wide_renamed<-sulf_red_KOhmm_wide %>%
  rename( "sat_K00958"="KO:K00958",
          "aprA_K00394"="KO:K00394",
          "aprB_K00395"="KO:K00395",
          "dsrA_K11180"="KO:K11180",
          "dsrB_K11181"="KO:K11181",
          "dsrE_tusD_K07235"="KO:K07235",
          "dsrF_tusC_K07236"="KO:K07236",
          "dsrH_tusB_K07237"="KO:K07237",
          "qmoA_K16885"="KO:K16885",
          "qmoB_K16886"="KO:K16886",
          "qmoC_K16887"="KO:K16887"
  )

#grab complete sulfate reducers/oxidizers from all sites
complete_SRs<-sulf_red_KOhmm_wide_renamed%>%filter(sat_K00958>=1 & aprA_K00394>=1 & aprB_K00395>=1 & dsrA_K11180>=1 & dsrB_K11181>=1 & dsrC_tusE>=1)

#map whether Dsrs are reductive or oxidative according to Dsr phylogeny
dsrs_red<-read.delim("../data/reductive_dsrs.txt", header = FALSE)
dsrs_ox<-read.delim("../data/oxidative_dsrs.txt", header = FALSE)
#cat tables
dsr_types<-rbind(dsrs_red,dsrs_ox)
#rename cols
colnames(dsr_types)[2]<-"type"
#split V1 for Bin mapping
dsr_types<-separate(dsr_types,V1,into = c("Bin", "Scaffold"), sep = "(?=_scaffold)")
#drop column 2 of dsr_types
dsr_types<-dsr_types[,-2]
dsr_types<-distinct(dsr_types)
#map dsr type
complete_SRs <- dsr_types %>%
  dplyr::select(Bin, type) %>%
  right_join(complete_SRs, by = c("Bin" = "Bin_name"))

#dropping 13_July_SF_Bin18 and 8_1_Oct_SF_Bin28 because they didn't make it into my phylogeny
#and were not identified with the hmm searches
complete_SRs <- complete_SRs[-c(103, 104), ]

#get higher taxonomic rank/ditch taxonomy detail
complete_SRs<-separate(complete_SRs,GTDBtk,into = c("Class", "tax_detailed"), sep = "(?=;o)")


#how represented are complete SRs/SOs per site?
complete_SRs_tmp<-separate(complete_SRs,Bin_name,into = c("Site", "Bin"), sep = "(?=_SF)")
table(complete_SRs_tmp$Site)
#what taxa are complete SRs/SOs?
complete_SRs_tmp2<-separate(complete_SRs,GTDBtk,into = c("Taxa", "taxa_long"), sep = "(?=;o)")
table(complete_SRs_tmp2$Taxa)

##TOMORROW --> 
  #add a column of the site to complete_SRs
  #remake the long format one with site included in group_by
  #convert the oxidative values to negative or subset and make 2, one oxidative, one reductive
  #make the stacked barplot: y = sum_col, x = Site, fill = Class

#map site to complete_SRs
complete_SRs <- tax %>%
  dplyr::select(Bin_name, Site) %>%
  right_join(complete_SRs, by = c("Bin_name" = "Bin"))

#group by taxa and type and pivot long for plot 
#to show distribution of sulfate reducers vs sulfur oxidizers
complete_SRs_long <- complete_SRs %>%
  mutate(sum_col = 1) %>%
  group_by(Class, type, Site) %>%
  summarise(sum_col = sum(sum_col))

#convert oxidative values to negative?
complete_SRs_long_neg<-complete_SRs_long %>%
  mutate(sum_col = ifelse(type=="reductive", sum_col, -sum_col))

##make a stacked bar plot?
complete_SRs_long_neg$Site <- factor(complete_SRs_long_neg$Site, levels=rev(c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May", 
                                                                  "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                                                  "13_July", "13_Oct", "13_Jan", "13_May",
                                                                  "21_July", "21_Oct", "21_Jan", "21_May",
                                                                  "24_July", "24_Oct", "24_Jan", "24_May")))

my_colors <- pals::kelly(n=12)[2:13]

plot<-complete_SRs_long_neg %>%
  ggplot(aes(x = Site, y = as.numeric(sum_col), fill = Class, color = type)) + 
  #geom_line(aes(linetype = type, group = 2)) +
  ggtitle("Dissimilatory sulfate reduction and oxidation in SFB MAGs") +
  geom_bar(stat = "identity", alpha = .8) +
  scale_fill_manual(values = as.vector(my_colors)) +
  labs(x = "Site", y = "Number of MAGs") +
  coord_cartesian(ylim = c(-14, 4), expand = FALSE) +
  scale_y_continuous(breaks = seq(-14, 4, by = 2)) +
  geom_hline(yintercept=0) +
  #guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  coord_flip() 
  #scale_y_continuous(expand = c(0, 0))
plot

ggsave("Sulfur_OxRed_MAGs.png", plot, width = 10, height = 6, dpi = 500)

#Subset sites
SR_4_1<-sulf_red_KOhmm_wide_renamed[grepl("4_1", sulf_red_KOhmm_wide_renamed$Bin_name),]
SR_8_1<-sulf_red_KOhmm_wide_renamed[grepl("8_1", sulf_red_KOhmm_wide_renamed$Bin_name),]
SR_13<-sulf_red_KOhmm_wide_renamed[grepl("13_", sulf_red_KOhmm_wide_renamed$Bin_name),]
SR_21<-sulf_red_KOhmm_wide_renamed[grepl("21_", sulf_red_KOhmm_wide_renamed$Bin_name),]
SR_24<-sulf_red_KOhmm_wide_renamed[grepl("24_", sulf_red_KOhmm_wide_renamed$Bin_name),]

#4_1
##sum 4_1 all gene content columns
colSums(SR_4_1[,c(-1,-2)])

##get shorter taxonomy info
SR_4_1<-separate(SR_4_1,GTDBtk,into = c("Phyla", "trash"), sep = "(?=;c)")
SR_4_1<-subset(SR_4_1, select = -c(trash))

##4_1 by month
SR_4_1_July<-SR_4_1[grepl("July", SR_4_1$Bin_name),]
colSums(SR_4_1_July[,c(-1,-2)])
###count taxa in July 
table(SR_4_1_July$Phyla)
###that have at least 1 in all columns?
test <- SR_4_1%>%filter(SR_4_1$sat_K00958>=1 & SR_4_1$aprA_K00394>=1 & SR_4_1$aprB_K00395>=1 & dsrA_K11180>=1 & dsrB_K11181>=1 & dsrC_tusE>=1)

SR_4_1_Oct<-SR_4_1[grepl("Oct", SR_4_1$Bin_name),]
colSums(SR_4_1_Oct[,c(-1,-2)])

SR_4_1_Jan<-SR_4_1[grepl("Jan", SR_4_1$Bin_name),]
colSums(SR_4_1_Jan[,c(-1,-2)])

SR_4_1_May<-SR_4_1[grepl("May", SR_4_1$Bin_name),]
colSums(SR_4_1_May[,c(-1,-2)])


#8_1
##8.1 by month

#set row names as Sample Location
sulf_red_KOhmm_wide_renamed <- sulf_red_KOhmm_wide_renamed %>% remove_rownames %>% column_to_rownames(var = "Bin_name")

#set the proper order for sulfate reduction genes
sulf_red_KOhmm_wide_renamed_ordered<-select(sulf_red_KOhmm_wide_renamed,1,4,2,3,8,9,13,14,5,6,7,10,11,12)




#################################  Sulfur oxidation  ######################################
#Subset the data for just sulfur oxidation genes
sulfur_ox <- data%>%filter(str_detect(KO_Term,"KO:K17222|KO:K17223|KO:K17224|KO:K17225|KO:K17226|KO:K17227|KO:K17218|KO:K22622"))
#supplement above output with dsrC and dsrD because those are not included in IMG output

#generate column of 1s for value mapping
newcol<-rep(1,1434) #generate 65 1s
sulfur_ox$presence<-newcol
#sum the presence of marker genes by Bin name
sulfur_ox_sum <- sulfur_ox %>%
  group_by(Bin, KO_Term) %>% 
  summarise(presence=sum(presence))

#go from long to wide
sulfur_ox_sum_wide <- sulfur_ox_sum %>%
  dplyr::select(KO_Term, presence, Bin) %>%
  pivot_wider(names_from = KO_Term, values_from = presence, values_fill = 0)

#rename KOs
sulfur_ox_sum_wide_renamed<-sulfur_ox_sum_wide %>%
  rename( "soxA_K17222"="KO:K17222",
          "soxX_K17223"="KO:K17223",
          "soxB_K17224"="KO:K17224",
          "soxC_K17225"="KO:K17225",
          "soxY_K17226"="KO:K17226",
          "soxZ_K17227"="KO:K17227",
          "sqr_K17218"="KO:K17218"
  )

sulfur_ox_sum_wide_renamed <- tax %>%
  dplyr::select(Bin_name, GTDBtk) %>%
  right_join(sulfur_ox_sum_wide_renamed, by = c("Bin_name" = "Bin"))

#grab complete sulfate reducers/oxidizers from all sites
complete_SOs<-sulfur_ox_sum_wide_renamed%>%filter(soxA_K17222>=1 & soxX_K17223>=1 & soxB_K17224>=1 & soxC_K17225>=1 & soxY_K17226>=1 & soxZ_K17227>=1)

#map site
complete_SOs <- tax %>%
  dplyr::select(Bin_name, Site) %>%
  right_join(complete_SOs, by = c("Bin_name" = "Bin_name"))

#drop detailed taxonomy to only keep class
complete_SOs<-separate(complete_SOs,GTDBtk,into = c("Class", "tax_detailed"), sep = "(?=;o)")

#group by taxa and type and pivot long for plot 
#to show distribution of sulfate reducers vs sulfur oxidizers
complete_SOs_long <- complete_SOs %>%
  mutate(sum_col = 1) %>%
  group_by(Class, Site) %>%
  summarise(sum_col = sum(sum_col))


##make a stacked bar plot?
complete_SOs_long$Site <- factor(complete_SOs_long$Site, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May", 
                                                                              "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                                                              "13_July", "13_Oct", "13_Jan", "13_May",
                                                                              "21_July", "21_Oct", "21_Jan", "21_May",
                                                                              "24_July", "24_Oct", "24_Jan", "24_May"))

#harry potter color palette
install.packages("harrypotter")
library(harrypotter)
#my_colors <- pals::kelly(n=3)[11:13]

#plot
plot<-complete_SOs_long %>%
  ggplot(aes(x = Site, y = as.numeric(sum_col), fill = Class)) + 
  #geom_line(aes(linetype = type, group = 2)) +
  ggtitle("SOX System in SFB MAGs") +
  geom_bar(stat = "identity", alpha = .8) +
  scale_fill_hp(house = "Mischief", discrete = TRUE) +
  labs(x = "Site", y = "Number of MAGs") +
  #coord_cartesian(ylim = c(-14, 4), expand = FALSE) +
  #scale_y_continuous(breaks = seq(-14, 4, by = 2)) +
  #geom_hline(yintercept=0) +
  #guides(fill=guide_legend(override.aes = list(size=3))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.background = element_blank())
  #guides(color = guide_legend(override.aes = list(fill = NA))) +
  #coord_flip() 
#scale_y_continuous(expand = c(0, 0))
plot

ggsave("SOX_MAGs.png", plot, width = 9, dpi = 500)

############### which MAGs that are complete SRs/SOs using dsr also have complete SOX system? #######################

test <- complete_SRs %>%
  filter(Bin %in% complete_SOs$Bin_name)

test <- complete_SRs[complete_SRs$Bin %in% complete_SOs$Bin_name, ]

test2 <- dsr_types %>%
  mutate(ID = 1) %>%
  group_by(Bin) %>%
  summarise(ID = sum(ID)) %>%
  filter(ID == 2)







################### old
#grab only 4_1, K11180 and K11181 (dsrA and dsrB)
# SR_4_1<-sulf_red_sum_wide[, c(1, 8, 9)]
# SR_4_1<-SR_4_1[grepl("4_1", SR_4_1$Bin),]
# 
# #summarize by month
# SR_4_1_Jan<-SR_4_1[grepl("Jan", SR_4_1$Bin),]
# colSums(SR_4_1_Jan[, c("KO:K11180", "KO:K11181")])
# 
# SR_4_1_May<-SR_4_1[grepl("May", SR_4_1$Bin),]
# colSums(SR_4_1_May[, c("KO:K11180", "KO:K11181")])
# 
# SR_4_1_July<-SR_4_1[grepl("July", SR_4_1$Bin),]
# colSums(SR_4_1_July[, c("KO:K11180", "KO:K11181")])
# #map taxonomy to the data frame
# SR_4_1_July_tax <- tax %>%
#   dplyr::select(Bin_name, GTDBtk) %>%
#   right_join(SR_4_1_July, by = c("Bin_name" = "Bin"))
# 
# SR_4_1_Oct<-SR_4_1[grepl("Oct", SR_4_1$Bin),]
# colSums(SR_4_1_Oct[, c("KO:K11180", "KO:K11181")])

# grab only 4_1
# SO_4_1<-sulfur_ox_sum_wide[grepl("4_1", sulfur_ox_sum_wide$Bin),]
# 
# #summarize by month
# SO_4_1_Jan<-SO_4_1[grepl("Jan", SO_4_1$Bin),]
# colSums(SO_4_1_Jan[, c("KO:K17222", "KO:K17224", "KO:K17225")])
# 
# SO_4_1_May<-SO_4_1[grepl("May", SO_4_1$Bin),]
# colSums(SO_4_1_May[, c("KO:K17222", "KO:K17224", "KO:K17225")])
# 
# SO_4_1_July<-SO_4_1[grepl("July", SO_4_1$Bin),]
# colSums(SO_4_1_July[, c("KO:K17222", "KO:K17224", "KO:K17225")])
# #map taxonomy to the data frame
# SO_4_1_July_tax <- tax %>%
#   dplyr::select(Bin_name, GTDBtk) %>%
#   right_join(SO_4_1_July, by = c("Bin_name" = "Bin"))
# 
# SO_4_1_Oct<-SO_4_1[grepl("Oct", SO_4_1$Bin),]
# colSums(SO_4_1_Oct[, c("KO:K17222", "KO:K17224", "KO:K17225")])
# 
