######################## Script adapted from J.D. Carlton, Baker Lab, UT Austin #######################################

library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)
library(gridExtra)
df <- read.csv("Input/input_heatmap.csv")
head(df)

df_num = data.frame(df[,3:136])
rownames(df_num) = df$MAG

#replace all gene presence values greater than 1 with 2
df_scale<-replace(df_num, df_num>1,2)

nitrogen <- as.matrix(t(df_scale[c(1:13)]))
other <- as.matrix(t(df_scale[c(14:18)]))
ETC <- as.matrix(t(df_scale[c(19:63)]))
glyc <- as.matrix(t(df_scale[c(64:86)]))
PPP <- as.matrix(t(df_scale[c(87:95)]))
TCA <- as.matrix(t(df_scale[c(96:108)]))
rTCA <- as.matrix(t(df_scale[c(109:118)]))
WLP <- as.matrix(t(df_scale[c(119:126)]))
rGlyp <- as.matrix(t(df_scale[c(127:133)]))
THP <- as.matrix(t(df_scale[c(134)]))

#Heatmap for the full table
# pheatmap(df_scale,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 7, cellheight = 7,
#          gaps_row = c(2,9,10,11,12), gaps_col = c(12,17,62,85,94,107,112,117),
# fontsize = 7, color = colorRampPalette(c("white", "blue", "navy"))(50))

#heatmap by group to get distinct colors
nhm <- pheatmap::pheatmap(nitrogen,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(2,9,10,11,12),
                          color = colorRampPalette(c("white", "#029E73"))(50))
ohm <- pheatmap::pheatmap(other,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(2,9,10,11,12),
                          color = colorRampPalette(c("white", "#f46a9b"))(50))
etchm <- pheatmap::pheatmap(ETC,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                            show_colnames = 0, legend = FALSE, gaps_col = c(2,9,10,11,12),
                            color = colorRampPalette(c("white", "#7eb0d5"))(50))
# ghm <- pheatmap::pheatmap(glyc,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
#                           show_rownames = 0, legend = FALSE, gaps_row = c(2,9,10,11,12),
#                           color = colorRampPalette(c("white", "purple"))(50))
# phm <- pheatmap::pheatmap(PPP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
#                           show_rownames = 0, legend = FALSE, gaps_row = c(2,9,10,11,12),
#                           color = colorRampPalette(c("white", "#FF4500"))(50))
thm <- pheatmap::pheatmap(TCA,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(2,9,10,11,12),
                          color = colorRampPalette(c("white", "#989898"))(50))
rtchm <- pheatmap::pheatmap(rTCA,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(2,9,10,11,12),
                          color = colorRampPalette(c("white", "#ffb55a"))(50))
wlphm <- pheatmap::pheatmap(WLP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(2,9,10,11,12),
                          color = colorRampPalette(c("white", "#bd7ebe"))(50)) 
rglhm <- pheatmap::pheatmap(rGlyp,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_colnames = 0, legend = FALSE, gaps_col = c(2,9,10,11,12),
                          color = colorRampPalette(c("white", "#b2e061"))(50))
thphm <- pheatmap::pheatmap(THP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                            show_colnames = 1, legend = FALSE, gaps_col = c(2,9,10,11,12),
                            color = colorRampPalette(c("white", "#ffee65"))(50))


#arrange heatmaps
lm <- rbind(c(1,2,3,6,7,8,9,10))
dev.off()
grid.arrange(grobs = list(nhm[[4]],
                          ohm[[4]],
                          etchm[[4]],
                          # ghm[[4]],
                          # phm[[4]],
                          thm[[4]],
                          rtchm[[4]],
                          wlphm[[4]],
                          rglhm[[4]],
                          thphm[[4]]),
             ncol = 1)#,
             #layout_matrix = lm)

g <- arrangeGrob(grobs = list(nhm[[4]],
                              ohm[[4]],
                              etchm[[4]],
                              # ghm[[4]],
                              # phm[[4]],
                              thm[[4]],
                              rtchm[[4]],
                              wlphm[[4]],
                              rglhm[[4]],
                              thphm[[4]]),
                 ncol = 1)
g

ggsave(filename = "Output/heatmap.png", g, dpi = 500, height = 10)
ggsave(filename = "Output/heatmap.svg", g, dpi = 500, height = 40)



