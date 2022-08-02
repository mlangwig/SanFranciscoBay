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
library(sf)
library(dplyr)

# read in shape file
sf_map <- st_read("data/water_bodies_carto.shp") #read in shape file. Pay attention to CRS for coordinate system
sf_map <- sf_map %>% st_transform(crs = 4326) #make coordinate ref system lat / long

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

# map = shape file...extension .shp
# raster = salinity...extension .geotiff
# brick for average

#plotting
mapworld <- borders("world",colour = "gray80",fill="white") 
mp<-ggplot() + mapworld + ylim(-60,90) 
mp + 
  geom_point(data = BR, aes(x = Longitude, y = Latitude, color=factor(Site), 
                            size  = AvgSalinity)) +
  scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1")) + #choose colors of sites
  scale_size(limits = c(2, 29), breaks = c(2, 10, 20, 27, 28)) +
  theme_gray() +
  labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
  guides(color = guide_legend(title = "Site"), 
         size = guide_legend(title = "Average Salinity")) + #change legend title
  coord_cartesian(xlim=c(-123,-121.5),ylim=c(37,38.5)) + #set limits of map
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12))

### better map version
ggplot() + 
  geom_sf(data = sf_map) +
  geom_point(data = BR, aes(x = Longitude, y = Latitude, color=factor(Site), 
                            size  = AvgSalinity)) +
  scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1")) + #choose colors of sites
  scale_size(limits = c(2, 29), breaks = c(2, 10, 20, 27, 28)) +
  theme_gray() +
  labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
  guides(color = guide_legend(title = "Site"), 
         size = guide_legend(title = "Average Salinity")) + #change legend title
  #coord_cartesian(xlim=c(-123,-121.5),ylim=c(37,38.5)) + #set limits of map
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12))



