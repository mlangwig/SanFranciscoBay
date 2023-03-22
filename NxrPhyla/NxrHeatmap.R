######################## Script adapted from J.D. Carlton, Baker Lab, UT Austin #######################################

library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(dplyr)
df <- read.csv("Input/input_heatmap.csv")
head(df)

df_num <- df %>% select(nxrA:mct)
rownames(df_num) = df$MAG

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

