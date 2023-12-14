
############################################################
# Purpose: 
#  1) compare number of metabolites in diverticulitis and controls
#  2) heatmap of highly prioritized metabolomic features from MAC
# Last update date: 9/25/2023
############################################################

# R --quiet
# 
# setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

library(readxl)
library(dplyr)
library(stringr)


#############################################################
# compared number of metablites in diverticulitis and control
#############################################################

abundances <- read.csv(file="./mbx/data/mac/abundances.csv")
abundances <- abundances[,-1] #56937   232
abundances[is.na(abundances)] <- 0
df1 <- abundances[rowSums(abundances[])>0,] #56914
count <- as.data.frame(colSums(df1 !=0))
names(count) <- "number"
count$sample <- rownames(count)

metadata <- read.csv(file="./mbx/data/mac/metadata_mbx.csv")
count_meta <- merge(count,metadata,by="sample")


group_by(count_meta, caco) %>%
  summarise(
    count = n(),
    mean = mean(number, na.rm = TRUE),
    sd = sd(number, na.rm = TRUE)
  )

  # Compute t-test
  res <- t.test(number ~ caco, data = count_meta, var.equal = TRUE)
  res

  # 	Two Sample t-test
  #
  # data:  number by caco
  # t = 0.81636, df = 230, p-value = 0.4151
  # alternative hypothesis: true difference in means is not equal to 0
  # 95 percent confidence interval:
  #  -800.653 1933.481
  # sample estimates:
  # mean in group Control  mean in group Divert
  #              51047.94              50481.53

#############################################################
# heatmap of highly priortized metabolomic features from MAC
#############################################################
abundances <- read.csv(file="./mbx/data/mac/abundances.csv")
abundances <- abundances[,-1] #56937   232
abundances[is.na(abundances)] <- 0

mbx <- as.data.frame(t(abundances)) # samples in rows

# Log transformation
LOG <- function(x) {
    y <- replace(x, x == 0, min(x[x>0]) / 2)
    return(log10(y))
}

mbx_log <- as.data.frame(apply(mbx, 2, LOG)) # 232 56937
mbx_log_scale <- as.data.frame(scale(mbx_log))
summary(mbx_log_scale[,1:10])
mbx_log_scale$sample <- rownames(mbx_log_scale)




# read in metadata
metadata <- read.csv(file="./mbx/data/mac/metadata_mbx.csv")
mbx_meta <- merge(mbx_log_scale,metadata,by="sample")
rownames(mbx_meta) <- mbx_meta$sample
mbx_meta <- mbx_meta[,-1]



# 11 232
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(seriation)
library(RColorBrewer)


#############################################################

# sort by V269 (Cer 18:1;O2/17:0)
mbx_meta <- mbx_meta[order(mbx_meta$caco,mbx_meta$V269),]

# select metabolites that I want to show in the heatmap
mbx_heatmap <- as.data.frame(t(mbx_meta[,c(269,14025,27879,13976,28549,1376,27175,27704,13981,13978,26971)]))
row.names(mbx_heatmap)[1] <- "Cer 18:1;O2/17:0"

# add m/z to names of unknowns 
# read in annotations
anno <- read.csv(file="./mbx/data/mac/annotations.csv")
anno <- anno[,-1]
t_anno <- as.data.frame(t(anno))
anno_subset <- as.data.frame(t(t_anno[,c(269,14025,27879,13976,28549,1376,27175,27704,13981,13978,26971)]))

mbx_heatmap <- cbind(anno_subset$MZ,mbx_heatmap)
rownames(mbx_heatmap) <- paste(rownames(mbx_heatmap),',m/z=',as.numeric(as.character(mbx_heatmap$`anno_subset$MZ`)),sep="")
mbx_heatmap <- mbx_heatmap[,-1]

# 11 232
#caco for annotation
caco_heatmap <- as.data.frame(mbx_meta$caco)
colnames(caco_heatmap) <- c('Status')
rownames(caco_heatmap) <- rownames(mbx_meta)

#anchor for annotation
# anchor_heatmap <- as.data.frame(mbx_meta$V269)
# colnames(anchor_heatmap) <- c('Cer 18:1;O2/17:0')
# rownames(anchor_heatmap) <- rownames(mbx_meta)
#
#
# annotation_heatmap <- mbx_meta[,colnames(mbx_meta) %in% c("caco","V269")]
# colnames(annotation_heatmap) <- c("Cer 18:1;O2/17:0","Status")


# summary(annotation_heatmap$'Cer 18:1;O2/17:0')
#
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 2.049   3.653   4.339   4.116   4.759   5.561


# ht = HeatmapAnnotation (df = annotation_heatmap,
#                         col = list (Status = c("Divert"="red","Control"="yellow"),
#                                     'Cer 18:1;O2/17:0' = colorRamp2(c(2, 3, 4, 5, 6), brewer.pal(n = 5, name = "Greens"))
#                                   ),
#                         annotation_name_gp= gpar(fontsize = 30),
#                         annotation_legend_param = list(title_gp = gpar(fontsize = 30),
#                                                labels_gp = gpar(fontsize = 20))
#                         )

ht = HeatmapAnnotation (df = caco_heatmap,
                        col = list(Status = c("Divert"="#4169E1","Control"="#E69F00")),
                        simple_anno_size = unit(3, "cm"),
                        annotation_name_gp= gpar(fontsize = 30),
                        annotation_legend_param = list(title_gp = gpar(fontsize = 30),
                                               labels_gp = gpar(fontsize = 20))
                        )


ht1_10species = Heatmap(mbx_heatmap, name = "Log10 z-score",top_annotation=ht,
                        col = colorRamp2(c(-6,-4,-2,0,2,4), rev(brewer.pal(n=6, name= "RdBu"))),
                        cluster_rows = FALSE, cluster_columns = FALSE, width = unit(80, "cm"),
                        heatmap_legend_param = list(color_bar = "continuous", at = c(-6,-4,-2,0,2,4), legend_direction = "horizontal",title_position="topcenter",
                        title_gp = gpar(fontsize = 30),labels_gp = gpar(fontsize = 50)),
                        row_names_gp = gpar(fontsize=50),
                        #column_names_gp  = gpar(fontsize = 8)
                        show_column_names = FALSE,
                      )

png("./mbx/results/mac/heatmap_mbx_by_cer_zscore.png", width = 1500, height = 300, units='mm', res = 300)
draw(ht1_10species, heatmap_legend_side = "bottom",annotation_legend_side = "left",
     #row_title = "Log10 abundances by Participant",
     #row_title_gp = gpar(fontsize = 26, fontface = "bold"),
     #column_title = "Cer 18:1;O2/17:0", column_title_side = "top",
     #column_title_gp = gpar(fontsize = 26, fontface = "bold")
   )
dev.off()


#############################################################
# heatmap of indolin
#############################################################
mbx_meta <- mbx_meta[order(mbx_meta$caco,mbx_meta$V180),]

# select metabolites that I want to show in the heatmap
mbx_heatmap <- as.data.frame(t(mbx_meta[,c(180,46426,14301,2947,4380,6944,39039,39953,3457,36334,37315)]))
row.names(mbx_heatmap)[1] <- "Indolin-2-one"

# 11 232
#caco for annotation
caco_heatmap <- as.data.frame(mbx_meta$caco)
colnames(caco_heatmap) <- c('Status')
rownames(caco_heatmap) <- rownames(mbx_meta)

ht = HeatmapAnnotation (df = caco_heatmap,
                        col = list(Status = c("Divert"="#4169E1","Control"="#E69F00")),
                        annotation_height=unit(6, "cm"),
                        annotation_name_gp= gpar(fontsize = 30),
                        annotation_legend_param = list(title_gp = gpar(fontsize = 30),
                                               labels_gp = gpar(fontsize = 20))
                        )


ht1_10species = Heatmap(mbx_heatmap, name = "Log10 abundances",top_annotation=ht,
                        col = colorRamp2(c(-6,-4,-2,0,2,4), rev(brewer.pal(n=6, name= "RdBu"))),
                        cluster_rows = FALSE, cluster_columns = FALSE, width = unit(80, "cm"),
                        heatmap_legend_param = list(color_bar = "continuous", at = c(-6,-4,-2,0,2,4), legend_direction = "horizontal",title_position="topcenter",
                        title_gp = gpar(fontsize = 30),labels_gp = gpar(fontsize = 20)),
                        row_names_gp = gpar(fontsize=50),
                        #column_names_gp  = gpar(fontsize = 8)
                        show_column_names = FALSE,
                      )

png("./mbx/results/mac/heatmap_mbx_by_indolin_zscore.png", width = 1200, height = 300, units='mm', res = 300)
draw(ht1_10species, heatmap_legend_side = "bottom",annotation_legend_side = "left",
   )
dev.off()

#############################################################
# heatmap of 1,7-Dimethyluric acid
#############################################################
29400,33,4513,3996,5126,46303,3604,4833

mbx_meta <- mbx_meta[order(mbx_meta$caco,mbx_meta$V29400),]

# select metabolites that I want to show in the heatmap
mbx_heatmap <- as.data.frame(t(mbx_meta[,c(29400,33,4513,3996,5126,46303,3604,4833)]))
row.names(mbx_heatmap)[1] <- "1,7-Dimethyluric acid"

# 11 232
#caco for annotation
caco_heatmap <- as.data.frame(mbx_meta$caco)
colnames(caco_heatmap) <- c('Status')
rownames(caco_heatmap) <- rownames(mbx_meta)

ht = HeatmapAnnotation (df = caco_heatmap,
                        col = list(Status = c("Divert"="#4169E1","Control"="#E69F00")),
                        annotation_height=unit(6, "cm"),
                        annotation_name_gp= gpar(fontsize = 30),
                        annotation_legend_param = list(title_gp = gpar(fontsize = 30),
                                               labels_gp = gpar(fontsize = 20))
                        )


ht1_10species = Heatmap(mbx_heatmap, name = "Log10 abundances",top_annotation=ht,
                        col = colorRamp2(c(-6,-4,-2,0,2,4), rev(brewer.pal(n=6, name= "RdBu"))),
                        cluster_rows = FALSE, cluster_columns = FALSE, width = unit(80, "cm"),
                        heatmap_legend_param = list(color_bar = "continuous", at = c(-6,-4,-2,0,2,4), legend_direction = "horizontal",title_position="topcenter",
                        title_gp = gpar(fontsize = 30),labels_gp = gpar(fontsize = 20)),
                        row_names_gp = gpar(fontsize=50),
                        #column_names_gp  = gpar(fontsize = 8)
                        show_column_names = FALSE,
                      )

png("./mbx/results/mac/heatmap_mbx_by_1,7-Dimethyluric_acid_zscore.png", width = 1200, height = 300, units='mm', res = 300)
draw(ht1_10species, heatmap_legend_side = "bottom",annotation_legend_side = "left",
   )
dev.off()
