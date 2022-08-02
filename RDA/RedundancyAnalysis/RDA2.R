##Tutorial from:
#http://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html#running-an-rda-in-r
library(data.table)
library(tidyverse)
library(readxl)
library(vegan)

##### RDA

library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(pals)
library(RColorBrewer)

#read in abundance data
data<-read.delim2(file = "Input/input_AbsoluteAbundance.txt")

data$Site <- factor(data$Site, levels=c("4_1_July", "4_1_Oct", "4_1_Jan", "4_1_May",
                                                  "8_1_July", "8_1_Oct", "8_1_Jan", "8_1_May",
                                                  "13_July", "13_Oct", "13_Jan", "13_May",
                                                  "21_July", "21_Oct", "21_Jan", "21_May",
                                                  "24_July", "24_Oct", "24_Jan", "24_May"))

data<-data %>%
  mutate(newcol = 1)
data<-data %>% 
  rename(
    Presence = newcol
  )

##mutate so that you have table
data_summed <- data %>%
  group_by(Site, Phylum) %>% #grouping taxa together for the same site
  summarise(value=sum(Presence)) %>% #summing the group
  ungroup() %>% #don't group by Site anymore
  group_by(Site) %>% #now group by Phylum
  arrange(desc(Presence))


data<-read.delim2(file = "Input/input_table.txt")


#read in geochemical data
geo<-read_xlsx(path = "../../../Geochemical_Dataxlsx.xlsx")


#make the first column the row names
data<-data %>% remove_rownames %>% column_to_rownames(var="Site")


#check if there are 0s
sum(data == 0)
#calculate proportion of 0s in dataset
sum(data == 0)/(nrow(data) * ncol(data))
#large number of 0s so need Hellinger transformation
spe.hel <- decostand(data, method = "hellinger")

#see if environmental variables are colinear
#subset data to assess for colinearity
env<-geo[7:25]

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

#drop Al, Cl, Ntot, Ctot, porosity, Fe, Mn

# Scale and center variables
env.z <- decostand(env, method = "standardize")

# Variables are now centered around a mean of 0
round(apply(env.z, 2, mean), 1)
# and scaled to have a standard deviation of 1
apply(env.z, 2, sd)
#remove correlated env data
env.z <- subset(env.z, select = -c(Al, Cl, Ntot, Ctot, porosity, Fe, Mn))

# Model the effect of all environmental variables on fish
# community composition

spe.rda <- rda(spe.hel ~ ., data = env.z)
summary(spe.rda)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(spe.hel ~ 1, data = env.z), # lower model limit (simple!)
                      scope = formula(spe.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = TRUE) # change to TRUE to see the selection process!

# Check the new model with forward-selected variables
fwd.sel$call

# Write our new model
spe.rda.signif <- rda(spe.hel ~ Mg, data = env.z)
# check the adjusted R2 (corrected for the number of
# explanatory variables)
RsquareAdj(spe.rda.signif)
#test significance of RDA
anova.cca(spe.rda.signif, step = 1000)
anova.cca(spe.rda.signif, step = 1000, by = "term")
anova.cca(spe.rda.signif, step = 1000, by = "axis")


#PLOTTING
#dev.off()
# Type 1 scaling (Distances among objects reflect their similarities)
ordiplot(spe.rda.signif, scaling = 1, type = "text")

# Type 2 scaling (Angles between variables reflect their correlation)
ordiplot(spe.rda.signif, scaling = 2, type = "text")


