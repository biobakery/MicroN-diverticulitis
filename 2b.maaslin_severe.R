##############################################################################################################################################
# Purpose: 
# 1) To identify microbial species and metabolites associated with diverticulitis severity using MaAsLin2
# 2) To show results in a heatmap
# Last update date: 10/11/2023 (formatting of heatmap)
#############################################################################################################################################

setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")


library(Maaslin2)
library(grid)
library(ggpubr)
library(gridExtra)
library(dplyr)

# read in taxonomy data
species_all <- read.table('./mgx/data/micron_species_filtered.tsv', header=TRUE, sep='\t')
rownames(species_all) <- species_all$UID
species_all <- species_all[,-1]


# read in metadata
meta <- read.table('./metadata/micron_metadata_severe.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID
meta <- meta[,-1]
#  Control Non-severe     Severe 
#118         98         19 

#meta <- meta[!meta$caco3 =="Control",]


meta_short <- meta[,colnames(meta) %in% c('caco3','age','race8905','bristolcat','abx')]





#############################################################################################################################################
                                                #### MaAsLin2 species #####
#############################################################################################################################################

# log-transform, default q=0.25
    Maaslin2(input_data       = species_all,
             input_metadata   = meta,
             output           = paste('./mgx/results/maaslin2/severe/species/', '/', sep=''),
             normalization    = 'TSS',
             standardize      = 'TRUE',
             transform        = 'LOG',
             analysis_method  = 'LM',
             min_abundance    = 0,
             min_prevalence   = 0,
             reference = c('bristolcat,Normal','abx,Non abx use','caco3,Non-severe'),
             cores            = 4)


   Maaslin2(input_data       = species_all,
            input_metadata   = meta_short,
            output           = paste('./mgx/results/maaslin2/severe/species_reduced/', '/', sep=''),
            normalization    = 'TSS',
            standardize      = 'TRUE',
            transform        = 'LOG',
            analysis_method  = 'LM',
            min_abundance    = 0,
            min_prevalence   = 0,
            reference = c('bristolcat,Normal','abx,Non abx use','caco3,Non-severe'),
            cores            = 4)

#############################################################################################################################################
                                        #### heatmap for species #####
#############################################################################################################################################
# combine results from maaslin for different models
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata,value))
  return(sig_result)
}

#sig_crude <- sig_results(dir='./mbx/results/maaslin2/individual/reduced/significant_results.tsv', exposure = "caco")
sig_adj   <- sig_results(dir='./mgx/results/maaslin2/severe/species/significant_results.tsv', exposure = "caco3")


all_results <- function(dir, exposure){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure)
  return(all_result)
}

sigall_adj   <- all_results(dir='./mgx/results/maaslin2/severe/species/all_results.tsv', exposure = "caco3")


sigfeature<-subset(sig_adj, select=c(feature))

join<-left_join(sigfeature, sigall_adj, by="feature")

join$name<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", join$feature)
join <- join %>% arrange(join$value, join$name)

bug_names<-join[1:(nrow(join)/2), (ncol(join)-1):ncol(join)]

# create heatmap for diet-taxonomy associations
level_x_order <- factor(join$value, level = c('Severe', 'Control'))
level_y_order <- factor(join$name, levels = bug_names$name)

join$stars <- cut(join$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))

p1<-ggplot(join, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = coef), color = "white",
            lwd = 1.5,
            linetype = 1) +
  coord_fixed(ratio = 0.5)+
  #scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  geom_text(aes(label=stars), color="black", size=20, show.legend = TRUE, vjust = 0.7) +
  xlab(NULL) +
  scale_y_discrete(limits=rev)+
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 50),
        legend.text = element_text(size = 30),
        legend.position = "right",
        plot.title = element_text(size=50),
        axis.title=element_text(size=50,face="bold"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=50, face="italic"),
        axis.text.x = element_text(size=50, angle = 90, vjust = 0.5, hjust=1, face="bold")) +
  labs(fill = "β cofficient")+
  guides(fill = guide_colorbar(barwidth = 2,
                               barheight = 30))
ggsave(file="./mgx/results/maaslin2/heatmap/heatmap_species_severe.jpg",plot=p1,width=18,height=18,dpi=300)





#############################################################################################################################################
                                                #### MaAsLin2 mbx #####
#############################################################################################################################################


mbx <- read.table('./mbx/data/mbx_maaslin_anno.tsv', header=TRUE, sep='\t')
rownames(mbx) <- mbx$sample
mbx <- mbx[,-1]
mbx <- mbx[order(rownames(mbx)),]

meta <- read.table('./mbx/data/mbx_maaslin_metadata_severe.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$sample
meta <- meta[,-1]
meta <- meta[order(rownames(meta)),]

meta_short <- meta[,colnames(meta) %in% c('caco3','age','race8905','bristolcat','abx')]


# log-transform, default q=0.25
    Maaslin2(input_data       = mbx,
             input_metadata   = meta,
             output           = paste('./mbx/results/maaslin2/individual/severe/full/', '/', sep=''),
             normalization    = 'TSS',
             standardize      = 'TRUE',
             transform        = 'LOG',
             analysis_method  = 'LM',
             min_abundance    = 0,
             min_prevalence   = 0,
             reference = c('bristolcat,Normal','abx,Non abx use','caco3,Non-severe'),
             cores            = 4)


     Maaslin2(input_data       = mbx,
              input_metadata   = meta_short,
              output           = paste('./mbx/results/maaslin2/individual/severe/reduced/', '/', sep=''),
              normalization    = 'TSS',
              standardize      = 'TRUE',
              transform        = 'LOG',
              analysis_method  = 'LM',
              min_abundance    = 0,
              min_prevalence   = 0,
              reference = c('bristolcat,Normal','abx,Non abx use','caco3,Non-severe'),
              cores            = 4)

#############################################################################################################################################
                                        #### heatmap for mbx #####
#############################################################################################################################################
# combine results from maaslin for different models
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata,value))
  return(sig_result)
}

#sig_crude <- sig_results(dir='./mbx/results/maaslin2/individual/reduced/significant_results.tsv', exposure = "caco")
sig_adj   <- sig_results(dir='./mbx/results/maaslin2/individual/severe/full/significant_results.tsv', exposure = "caco3")

#sig_crude$model <- "Reduced"

all_results <- function(dir, exposure){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure)
  return(all_result)
}

#sigall_crude <- all_results(dir='./mbx/results/maaslin2/individual/reduced/all_results.tsv', exposure = "caco")
sigall_adj   <- all_results(dir='./mbx/results/maaslin2/individual/severe/full/all_results.tsv', exposure = "caco3")


sigfeature<-subset(sig_adj, select=c(feature))

join<-left_join(sigfeature, sigall_adj, by="feature")

# get annotations
anno <- read.table('./mbx/data/mbx_maaslin_anno_name.tsv',header=TRUE,sep='\t')
anno$feature <- anno$mol_id

bind4_anno <- left_join(join,anno,by="feature")

bind4_anno <- bind4_anno %>% arrange(bind4_anno$value, bind4_anno$Metabolite.name)
#bind4_anno$name <- paste(bind4_anno$Metabolite.name, bind4_anno$mol_id)

bug_names<-bind4_anno[1:(nrow(bind4_anno)/2), (ncol(bind4_anno)-1):ncol(bind4_anno)]

# create heatmap for diet-taxonomy associations
level_x_order <- factor(bind4_anno$value, level = c('Severe', 'Control'))
level_y_order <- factor(bind4_anno$Metabolite.name, levels = bug_names$Metabolite.name)

bind4_anno$stars <- cut(bind4_anno$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))

p1<-ggplot(bind4_anno, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = coef), color = "white",
            lwd = 1.5,
            linetype = 1) +
  coord_fixed(ratio = 0.5)+
  #scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  geom_text(aes(label=stars), color="black", size=20, show.legend = TRUE,  vjust = 0.7) +
  xlab(NULL) +
  scale_y_discrete(limits=rev)+
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 50),
        legend.text = element_text(size = 30),
        legend.position = "right",
        plot.title = element_text(size=50),
        axis.title=element_text(size=50,face="bold"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=50),
        axis.text.x = element_text(size=50, angle = 90, vjust = 0.5, hjust=1, face="bold")) +
  labs(fill = "β cofficient")+
  guides(fill = guide_colorbar(barwidth = 2,
                               barheight = 30))
ggsave(file="./mbx/results/maaslin2/individual/heatmap/heatmap_mbx_severe.jpg",plot=p1,width=20,height=30,dpi=300)
