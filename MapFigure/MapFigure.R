### SFB sampling site map

## Modified from original author, Dr. Mirna Vazquez Rosas Landa
## email: mirnavrl@austin.utexas.edu

#read input
BR <- readr::read_tsv("input.tsv")
raster_df <- readr::read_tsv("salinity.txt")

#summary
summary(BR)

#libraries
library(ggmap)
library(RColorBrewer)
library(pals)

#convert lat and long & salinity to numeric for plotting
BR$Latitude<-as.numeric(BR$Latitude)
BR$Longitude<-as.numeric(BR$Longitude)
BR$AvgSalinity<-as.numeric(BR$AvgSalinity)

#?
mp <- NULL 
nb.cols <- 17

#color scheme
#mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
mycolors <- pals::kelly(n=5)


#plotting
mapworld <- borders("world",colour = "gray80",fill="white") 
mp<-ggplot() + mapworld + ylim(-60,90) 
mp + 
  geom_point(data = BR, aes(x = Longitude, y = Latitude, color=factor(Site), 
                            size  = AvgSalinity)) +
  scale_fill_manual(values=c("blue","green","grey","yellow", "purple")) + #choose colors of sites
  theme_gray() +
  labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
  guides(color = guide_legend(title = "Sites")) + #change legend title
  coord_cartesian(xlim=c(-123,-121.5),ylim=c(37,38.5)) + #set limits of map
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12))


