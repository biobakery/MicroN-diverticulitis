#########################################################################################################
### Purpose: To correlate beta coefficients for species in association with diverticulitis and IBD 
### Last update date: 10/13/2023 
#########################################################################################################

#setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

library(readxl)
library(dplyr)
library(stringr)
library(lmerTest)
library(ggplot2)
library(tidyverse)


# read in all_results.tsv for diverticulitis
all_results <- read.table('./mgx/results/maaslin2/species/all_results.tsv', header=TRUE, sep='\t', check.names=TRUE, quote ="")
all_results_divert <- all_results[all_results$metadata=="caco",]


# read in results for IBD
all_results_ibd <- read.csv('./reference/HMP2_Nature/2018-01-01078D-s6/ST1.csv', header=TRUE, sep = ",", stringsAsFactors = FALSE)
all_results_ibd$feature <- all_results_ibd$Feature
all_results_ibd$feature <- sub(" ", "_", all_results_ibd$feature)
  
# merge results for diverticulitis and IBD
all <- merge(all_results_divert,all_results_ibd,by="feature")

# correlations between coefficients
library(Hmisc)
rcorr(all$coef,all$Coefficient..CD.)
rcorr(all$coef,all$Coefficient..UC.)

# scatter plots
library(ggplot2)
library(ggpubr)

colnames(all)[which(names(all) == "coef")] <- "Coefficient for diverticulitis"
colnames(all)[which(names(all) == "Coefficient..CD.")] <- "Coefficient for CD"
colnames(all)[which(names(all) == "Coefficient..UC.")] <- "Coefficient for UC"

ggp <- function(x,y) {
  plot <- ggplot(all, aes_string(x=x, y=y))+
    geom_point(alpha =2, shape = 21, size =10, stroke = 1, na.rm=TRUE,color='black',fill='#759DDB') +
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE,size=2)+
    theme_bw()+
    theme(legend.position = "bottom",
          legend.text = element_text(size=50), legend.title = element_text(size=50),
          axis.title=element_text(size=50,face='bold'),
          axis.text.x = element_text(size=50),
          axis.text.y = element_text(size=50),
          axis.ticks =element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_blank())+
    stat_cor(method = "spearman",size=12)
  filepath <- paste('./mgx/results/similarity_ibd/', x,'_',y,'.jpg', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 16, height = 15, dpi = 300)
  print(plot)
}# end

ggp("`Coefficient for diverticulitis`","`Coefficient for CD`")
ggp("`Coefficient for diverticulitis`","`Coefficient for UC`")
