######################## Script adapted from JD Carlton, Baker lab #######################################

library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)
library(gridExtra)
df <- read.csv("Input/input_heatmap.csv")
head(df)

df_num = data.frame(df[,3:126])
rownames(df_num) = df$MAG

#replace all gene presence values greater than 1 with 2
df_scale<-replace(df_num, df_num>1,2)

nitrogen <- as.matrix(df_scale[c(1:12)])
other <- as.matrix(df_scale[c(13:17)])
ETC <- as.matrix(df_scale[c(18:62)])
glyc <- as.matrix(df_scale[c(63:85)])
PPP <- as.matrix(df_scale[c(86:94)])
TCA <- as.matrix(df_scale[c(95:107)])
rTCA <- as.matrix(df_scale[c(108:112)])
WLP <- as.matrix(df_scale[c(113:117)])
rGlyp <- as.matrix(df_scale[c(118:124)])

#Heatmap for the full table
#pheatmap(df_scale,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 7, cellheight = 7, 
#gaps_row = c(6,70), gaps_col = c(5,22,45,66,78,92,101,108),
#fontsize = 7, color = colorRampPalette(c("white", "blue", "navy"))(50))

nhm <- pheatmap::pheatmap(nitrogen,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_rownames = 0, legend = FALSE, gaps_row = c(2,9,10,11,12),
                          color = colorRampPalette(c("white", "#B22222"))(50))
ohm <- pheatmap::pheatmap(other,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_rownames = 0, legend = FALSE, gaps_row = c(6,70),
                          color = colorRampPalette(c("white", "blue"))(50))
etchm <- pheatmap::pheatmap(ETC,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                            show_rownames = 0, legend = FALSE, gaps_row = c(6,70),
                            color = colorRampPalette(c("white", "dark green"))(50))
ghm <- pheatmap::pheatmap(glyc,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_rownames = 0, legend = FALSE, gaps_row = c(6,70),
                          color = colorRampPalette(c("white", "purple"))(50))
phm <- pheatmap::pheatmap(PPP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_rownames = 0, legend = FALSE, gaps_row = c(6,70),
                          color = colorRampPalette(c("white", "#FF4500"))(50))
thm <- pheatmap::pheatmap(TCA,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_rownames = 0, legend = FALSE, gaps_row = c(6,70),
                          color = colorRampPalette(c("white", "goldenrod"))(50))
rtchm <- pheatmap::pheatmap(rTCA,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_rownames = 0, legend = FALSE, gaps_row = c(6,70),
                          color = colorRampPalette(c("white", "#6F4E37"))(50))
wlphm <- pheatmap::pheatmap(WLP,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_rownames = 0, legend = FALSE, gaps_row = c(6,70),
                          color = colorRampPalette(c("white", "magenta"))(50)) 
rglhm <- pheatmap::pheatmap(rGlyp,cluster_cols = F,cluster_rows = F,border_color = "black", cellwidth = 5, cellheight = 5, fontsize = 5,
                          show_rownames = 1, legend = FALSE, gaps_row = c(6,70),
                          color = colorRampPalette(c("white", "#474747"))(50))


#arrange heatmaps
lm <- rbind(c(1,2,3,4,5,6,7,8,9))
dev.off()
grid.arrange(grobs = list(nhm[[4]],
                          ohm[[4]],
                          etchm[[4]],
                          ghm[[4]],
                          phm[[4]],
                          thm[[4]],
                          rtchm[[4]],
                          wlphm[[4]],
                          rglhm[[4]]),
             layout_matrix = lm)

g <- arrangeGrob(grobs = list(ahm[[4]],
                              whm[[4]],
                              etchm[[4]],
                              ghm[[4]],
                              phm[[4]],
                              thm[[4]],
                              hhm[[4]],
                              dhm[[4]],
                              ehm[[4]]),
                 layout_matrix = lm)
ggsave(filename = "Output/heatmap.png", g, dpi = 500, width = 14, height = 10)
