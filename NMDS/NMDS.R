library(vegan)
library(ggplot2)
library(ggrepel)

######################### Input data relative abundance #####################################

#data<-read.delim2(file = "../../Abundance_Barplot/Input/SFBMAGS-vs-Reads-CoverM.tsv")
data<-read.delim2(file = "../../../Coverm/dereplicate/output_coverm.tsv")
str(data)
tax<-read.delim2(file = "../../Abundance_Barplot/Input/gtdbtk.r207.bacarc.summary.tsv")

#keep only relative abundance cols
data_filt <- data %>% 
  dplyr::select(contains("Relative.Abundance"))

#add back genome col
data_filt <- data_filt %>%
  dplyr::mutate(Genome = data[[1]]) %>%
  dplyr::select(Genome, everything())

data <- data_filt

#remove umapped row
data <- data[-1, ]

#remove extraneous strings from col names
colnames(data) <- gsub("_filter.METAGENOME.fastq.gz.Relative.Abundance....", "", colnames(data))

#add taxonomy
data <- tax %>%
  dplyr::select("user_genome", "classification") %>%
  right_join(data, by = c("user_genome" = "Genome"))

#clean up taxonomy
#phylum level
#data <- data %>% separate(classification, c("classification", NA), sep= ";c")
#data <- data %>% separate(classification, c(NA, "classification"), sep= ";")

#genus level
data <- data %>% separate(classification, c("classification", NA), sep= ";s")
data <- data %>% separate(classification, c(NA, "classification"), sep= ".*;")

#sum relative abundance by phylum
data <- data %>%
  mutate(across(starts_with("SF"), as.numeric)) %>%
  group_by(classification) %>%
  summarise(across(starts_with("SF"), sum))
str(data)

#transpose
data <- as.data.frame(t(data))
#make column names
data <- data %>%
  row_to_names(row_number = 1)
#fix numeric
data <- data %>%
  mutate(across(starts_with("p"), as.numeric)) 
str(data)

data <- data %>%
  mutate(across(starts_with("g"), as.numeric)) 
str(data)

data <- data %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))
str(data)

#order the data same as env data
mapping <- read.delim2("Input/mapping_samples.txt")
data <- data %>%
  rownames_to_column(var = "Site_Long")
#map the ordering variable
data <- mapping %>%
  right_join(data, by = "Site_Long")
data <- data %>%
  arrange(Order) %>%
  select(-c("Site_Short", "Order")) %>%
  column_to_rownames(var = "Site_Long")

######################### Geochemical data #####################################

# Tutorial from #http://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html#running-an-rda-in-r

#read in geochemical data
geo<-read_xlsx(path = "Input/Geochemical_Dataxlsx.xlsx")
#remove dist col
geo <- geo %>%
  select(-c('Dist'))

#order the geo data by order column
geo <- geo %>%
  arrange(Order)

#see if environmental variables are colinear
#subset data to assess for colinearity
env<-geo[7:24]

# We can visually look for correlations between variables:
heatmap(abs(cor(env)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("bottomright",
       #inset = c(-.5,0),
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

# Scale and center variables
env.z <- decostand(env, method = "standardize")

# Variables are now centered around a mean of 0
round(apply(env.z, 2, mean), 1)
# and scaled to have a standard deviation of 1
apply(env.z, 2, sd)
#remove correlated env data
env.z <- subset(env.z, select = -c(Al, Cl, Ntot, Ctot, porosity, Fe, Mn, Pb, Na, Cu, Temp))
#remove Mg just to see
#env.z <- subset(env.z, select = -c(Mg))

######################### Tutorial 1 - Example #####################################

# # Tutorial from #http://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html#running-an-rda-in-r
# 
# spe <- read.csv("~/Downloads/doubsspe.csv", row.names = 1)
# env <- read.csv("~/Downloads/doubsenv.csv", row.names = 1)
# 
# # Run the NMDS
# spe.nmds <- metaMDS(spe[, -8], distance = "bray", k = 2)
# 
# # Extract the results
# spe.nmds
# 
# # Assess the goodness of fit and draw a Shepard plot
# spe.nmds$stress
# stressplot(spe.nmds, main = "Shepard plot")
# 
# # Construct the biplot
# 
# plot(spe.nmds, 
#      type = "none", 
#      main = paste("NMDS/Bray - Stress=",round(spe.nmds$stress, 3)), 
#      xlab = c("NMDS1"), ylab = c("NMDS2"))
# points(scores(spe.nmds, display = "sites", choices = c(1, 2)),
#        pch = 21, col = "black", bg = "steelblue", cex = 1.2)
# text(scores(spe.nmds, display = "species", choices = c(1)), 
#      scores(spe.nmds, display = "species", choices = c(2)), 
#      labels = rownames(scores(spe.nmds, display = "species")), col = "red", cex = 0.8)

######################### Tutorial 1 - My Data #####################################

# Run the NMDS
set.seed(100)
#spe.nmds <- metaMDS(spe.rclr.mat, distance = "euclidean", trymax = 100)
spe.nmds <- metaMDS(data, distance = "bray", k = 2)
# Extract the results
spe.nmds

# Assess the goodness of fit and draw a Shepard plot
spe.nmds$stress
stressplot(spe.nmds, main = "Shepard plot")

# Construct the biplot

plot(spe.nmds, 
     type = "none", 
     main = paste("NMDS/Bray - Stress=",round(spe.nmds$stress, 3)), 
     xlab = c("NMDS1"), ylab = c("NMDS2"))
points(scores(spe.nmds, display = "sites", choices = c(1, 2)),
       pch = 21, col = "black", bg = "steelblue", cex = 1.2)
text(scores(spe.nmds, display = "species", choices = c(1)), 
     scores(spe.nmds, display = "species", choices = c(2)), 
     labels = rownames(scores(spe.nmds, display = "species")), col = "red", cex = 0.8)
text(scores(spe.nmds, display = "sites", choices = c(1))[, 1], 
     scores(spe.nmds, display = "sites", choices = c(2))[, 1], 
     labels = rownames(scores(spe.nmds, display = "sites")), 
     col = "black", cex = 0.8, pos = 4)  # Adjust position with 'pos' (4 = right)

######################### Tutorial 2 - My Data #####################################

#Tutuorial from https://jkzorz.github.io/2019/06/06/NMDS.html

#spe.rclr.mat <- as.matrix(spe.rclr)
#nmds code
set.seed(123)
# nmds = metaMDS(spe.rclr.mat, distance = "euclidean", trymax = 100)
#nmds = metaMDS(data, distance = "bray", k = 2)
spe.nmds
plot(spe.nmds)

stressplot(spe.nmds, main = "Shepard plot")

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(spe.nmds))
data.scores = as.data.frame(scores(spe.nmds)$sites)

#add metadata columns to data frame 
# data.scores$Sample = pc$Sample
# data.scores$Time = pc$Time
data.scores$Site_short = rownames(data.scores)
data.scores$Site <- rownames(data.scores)

#clean up the site names
data.scores <- data.scores %>%
  mutate(Site = str_replace(Site, "sed_USGS_", "")) %>%
  mutate(Site_short = str_replace(Site_short, "SF_.*_sed_USGS_", ""))

head(data.scores)

set.seed(150)
en = envfit(spe.nmds, env.z, permutations = 999, na.rm = TRUE)
en

plot(spe.nmds)
plot(en)

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

#assign colors to sites
my_colors = c("13" = "#CE9A28", "21" = "#28827A", "24" = "#3F78C1", "4_1" = "#4F508C", "8_1" = "#B56478")

gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Site_short), size = 3, alpha = 0.5) +  # Assign color only here
  geom_text_repel(aes(label = Site), size = 3, max.overlaps = Inf, colour = "black") + # Prevents text overlap
  scale_colour_manual(values = my_colors) + 
  #env data
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, linewidth = 1, alpha = 0.5, colour = "grey30") +
  geom_point(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2+0.04), 
            label = row.names(en_coord_cont), colour = "navy", fontface = "bold") + 
  # geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
  #           fontface = "bold", label = row.names(en_coord_cont)) + 
  #end env data
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")
        ) +
  labs(colour = "Site") +
  ggtitle(paste("NMDS/Bray - Stress =", round(spe.nmds$stress, 3))) +  # Add dynamic title
  expand_limits(x = c(min(data.scores$NMDS1) - 0.1, max(data.scores$NMDS1) + 0.1),
                y = c(min(data.scores$NMDS2) - 0.1, max(data.scores$NMDS2) + 0.1))
gg

ggsave("Output/nmds_env_sfb.svg", gg, dpi = 500)
ggsave("Output/nmds_env_sfb.png", gg, dpi = 500)

#ggsave("Output/nmds_env_sfb_genus.svg", gg, dpi = 500)
#ggsave("Output/nmds_env_sfb_genus.png", gg, dpi = 500, height = 7, width = 9)

#ggsave("Output/nmds_env_sfb_species.svg", gg, dpi = 500)
ggsave("Output/nmds_env_sfb_species.png", gg, dpi = 500)





