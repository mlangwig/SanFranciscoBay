library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)

#read data from IMG tables (see readme for how this data was obtained)
data <- read.delim("../data/all_SFsamples_sulfur_binsOnly.header.tsv", sep = "\t", header = TRUE)

#read input from HMM search of custom sulfur database (see readme for how results were obtained)
data_hmm <- read.delim("../data/all_sulfur_hmmOutput_SFB.tsv")
#I forgot to remove the headers when concatenated the files so removing them here:
#This reduces data by 77. Makes sense because 1 is the header, and another is missing because no
#hits for soxC
data_hmm<-data_hmm[!grepl("protein", data_hmm$protein),]

#Also, note that this input has duplicate hits within a MAG aka it was sorted for best hits per 
#protein, NOT per MAG/genome. I'm keeping this because interested in the multi copy dsrCs in 
#some proteobacteria. We have a script to get best MAG hits though if needed.

#read in taxonomy data
tax<-read.delim("../../../SFB_tax_detailed.txt")

#length(unique(data$Product))

#Subset the data for just sulfate reduction genes
sulf_red <- data%>%filter(str_detect(KO_Term,"KO:K00958|KO:K00394|KO:K00395|KO:K11180|KO:K11181|KO:K07235|KO:K07236|KO:K07237|KO:K16885|KO:K16886|KO:K16887"))
#supplement above output with dsrC and dsrD because those are not included in IMG output
sulf_red_hmm <- data_hmm%>%filter(str_detect(accession, "dsrC_tusE|dsrE_protein.alignment"))


#generate column of 1s for value mapping
newcol<-rep(1,2143) #generate 65 1s
sulf_red$presence<-newcol

#sum the presence of marker genes by Bin name
sulf_red_sum <- sulf_red %>%
  group_by(Bin, KO_Term) %>% 
  summarise(presence=sum(presence))

#go from long to wide
sulf_red_sum_wide <- sulf_red_sum %>%
  dplyr::select(KO_Term, presence, Bin) %>%
  pivot_wider(names_from = KO_Term, values_from = presence, values_fill = 0)

#grab only 4_1, K11180 and K11181 (dsrA and dsrB)
SR_4_1<-sulf_red_sum_wide[, c(1, 8, 9)]
SR_4_1<-SR_4_1[grepl("4_1", SR_4_1$Bin),]

#summarize by month
SR_4_1_Jan<-SR_4_1[grepl("Jan", SR_4_1$Bin),]
colSums(SR_4_1_Jan[, c("KO:K11180", "KO:K11181")])

SR_4_1_May<-SR_4_1[grepl("May", SR_4_1$Bin),]
colSums(SR_4_1_May[, c("KO:K11180", "KO:K11181")])

SR_4_1_July<-SR_4_1[grepl("July", SR_4_1$Bin),]
colSums(SR_4_1_July[, c("KO:K11180", "KO:K11181")])
#map taxonomy to the data frame
SR_4_1_July_tax <- tax %>%
  dplyr::select(Bin_name, GTDBtk) %>%
  right_join(SR_4_1_July, by = c("Bin_name" = "Bin"))

SR_4_1_Oct<-SR_4_1[grepl("Oct", SR_4_1$Bin),]
colSums(SR_4_1_Oct[, c("KO:K11180", "KO:K11181")])


#set row names as Sample Location
sulf_red_sum_wide_names <- sulf_red_sum_wide %>% remove_rownames %>% column_to_rownames(var = "Bin")


#Subset the data for just sulfur oxidation genes
sulfur_ox <- data%>%filter(str_detect(KO_Term,"KO:K17222|KO:K17223|KO:K17224|KO:K17225|KO:K17226|KO:K17227|KO:K17218"))
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

#grab only 4_1
SO_4_1<-sulfur_ox_sum_wide[grepl("4_1", sulfur_ox_sum_wide$Bin),]

#summarize by month
SO_4_1_Jan<-SO_4_1[grepl("Jan", SO_4_1$Bin),]
colSums(SO_4_1_Jan[, c("KO:K17222", "KO:K17224", "KO:K17225")])

SO_4_1_May<-SO_4_1[grepl("May", SO_4_1$Bin),]
colSums(SO_4_1_May[, c("KO:K17222", "KO:K17224", "KO:K17225")])

SO_4_1_July<-SO_4_1[grepl("July", SO_4_1$Bin),]
colSums(SO_4_1_July[, c("KO:K17222", "KO:K17224", "KO:K17225")])
#map taxonomy to the data frame
SO_4_1_July_tax <- tax %>%
  dplyr::select(Bin_name, GTDBtk) %>%
  right_join(SO_4_1_July, by = c("Bin_name" = "Bin"))

SO_4_1_Oct<-SO_4_1[grepl("Oct", SO_4_1$Bin),]
colSums(SO_4_1_Oct[, c("KO:K17222", "KO:K17224", "KO:K17225")])
