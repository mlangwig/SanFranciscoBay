##Tutorial from:
#http://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html#running-an-rda-in-r
library(data.table)
library(tidyverse)
library(readxl)
library(vegan)

##### RDA

#read in abundance data
data<-read_delim(file = "../../Abundance_Barplot/Input/input.tsv", header = FALSE)
data_log<-read.delim2(file = "../../Abundance_Barplot/Input/input_log.tsv", header = FALSE)
#read in geochemical data
geo<-read_xlsx(path = "../../../Geochemical_Dataxlsx.xlsx")

#transform so col names are row names
data_trans<-as.data.frame(t(data))
data_trans_log<-as.data.frame(t(data_log))

#remove taxonomy
data_trans<-data_trans[-1,]
data_trans_log<-data_trans_log[-1,]
#make columns col names
colnames(data_trans)<-data_trans[1,]
colnames(data_trans_log)<-data_trans_log[1,]
#remove duplicate
data_trans<-data_trans[-1,]
data_trans_log<-data_trans_log[-1,]
#make the first column the row names
data_trans<-data_trans %>% remove_rownames %>% column_to_rownames(var="Bin")
data_trans_log<-data_trans_log %>% remove_rownames %>% column_to_rownames(var="Bin")

#check if there are 0s
sum(data_trans == 0)

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


#convert to numeric
data_trans_num<-data_trans %>% mutate_if(is.character,as.numeric)
data_trans_log_num<-data_trans_log %>% mutate_if(is.character,as.numeric)

# Apply Hellinger transformation to correct for the double
# zero problem
#data.hel <- decostand(data_trans_num, method = "hellinger")

# Model the effect of all environmental variables on fish
# community composition

spe.rda <- rda(data_trans_num ~ ., data = env.z)
summary(spe.rda)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(data_trans_num ~ 1, data = env.z), # lower model limit (simple!)
                      scope = formula(spe.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = TRUE) # change to TRUE to see the selection process!

# Check the new model with forward-selected variables
fwd.sel$call

# Write our new model
spe.rda.signif <- rda(data_trans_num ~ Cu, data = env.z)
# check the adjusted R2 (corrected for the number of
# explanatory variables)
RsquareAdj(spe.rda.signif)
#test significance of RDA
anova.cca(spe.rda.signif, step = 1000)
anova.cca(spe.rda.signif, step = 1000, by = "term")
anova.cca(spe.rda.signif, step = 1000, by = "axis")


#PLOTTING
# Type 1 scaling
ordiplot(spe.rda.signif, scaling = 1, type = "text")
# Type 2 scaling
ordiplot(spe.rda.signif, scaling = 2, type = "text")




