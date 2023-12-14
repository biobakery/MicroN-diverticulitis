##############################################################################################################################################
# Purpose: 
  # 1) To identify microbial species, pathways, enzymes, and metabolites associated with diverticulitis using MaAsLin2
  # 2) To show results in a heatmap
# Last update date: 10/11/2023 (formatting of heatmap)
#############################################################################################################################################
# salloc -p test --mem 64G -t 0-08:00
# module load R/4.2.2-fasrc01
# #install new packages
# export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
# R --quiet
# 
# rm(list=ls())
# getwd()
# setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
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

# read in pathway data
path <- read.table('./mgx/data/micron_path_filtered.tsv', header=TRUE, sep='\t')
rownames(path) <- path$UID
path <- path[,-1]

# read in ec data
ec <- read.table('./mgx/data/micron_ec_filtered.tsv', header=TRUE, sep='\t')
rownames(ec) <- ec$UID
ec <- ec[,-1]

# rename column names for ec data
colnames(ec) <- gsub("X","EC", colnames(ec))


# read in cazy dna data
# cazy <- read.table('./mgx/data/micron_cazy_filtered.tsv', header=TRUE, sep='\t')
# rownames(cazy) <- cazy$UID
# cazy <- cazy[,-1]

# read in metadata
meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID
meta <- meta[,-1]

meta_short <- meta[,colnames(meta) %in% c('caco','age','race8905','bristolcat','abx')]


#############################################################################################################################################
                                                #### MaAsLin2 species #####
#############################################################################################################################################

# log-transform, default q=0.25
    Maaslin2(input_data       = species_all,
             input_metadata   = meta,
             output           = paste('./mgx/results/maaslin2/species/', '/', sep=''),
             normalization    = 'TSS',
             standardize      = 'TRUE',
             transform        = 'LOG',
             analysis_method  = 'LM',
             min_abundance    = 0,
             min_prevalence   = 0,
             reference = c('bristolcat,Normal','abx,Non abx use'),
             cores            = 4)

   Maaslin2(input_data       = species_all,
            input_metadata   = meta_short,
            output           = paste('./mgx/results/maaslin2/species_reduced/', '/', sep=''),
            normalization    = 'TSS',
            standardize      = 'TRUE',
            transform        = 'LOG',
            analysis_method  = 'LM',
            min_abundance    = 0,
            min_prevalence   = 0,
            reference = c('bristolcat,Normal','abx,Non abx use'),
            cores            = 4)


 #############################################################################################################################################
                                                 #### MaAsLin2 path #####
 #############################################################################################################################################

 # log-transform, default q=0.25
     Maaslin2(input_data       = path,
              input_metadata   = meta,
              output           = paste('./mgx/results/maaslin2/path/', '/', sep=''),
              normalization    = 'TSS',
              standardize      = 'TRUE',
              transform        = 'LOG',
              analysis_method  = 'LM',
              min_abundance    = 0,
              min_prevalence   = 0,
             reference = c('bristolcat,Normal','abx,Non abx use'),
              cores            = 4)

    Maaslin2(input_data       = path,
             input_metadata   = meta_short,
             output           = paste('./mgx/results/maaslin2/path_reduced/', '/', sep=''),
             normalization    = 'TSS',
             standardize      = 'TRUE',
             transform        = 'LOG',
             analysis_method  = 'LM',
             min_abundance    = 0,
             min_prevalence   = 0,
            reference = c('bristolcat,Normal','abx,Non abx use'),
             cores            = 4)
    
    
# check parasutterella pathways
# path <- read.table('./mgx/data/micron_path_filtered_parasutterella.tsv', header=TRUE, sep='\t')
# rownames(path) <- path$UID
# path <- path[,-1]
#     Maaslin2(input_data       = path,
#              input_metadata   = meta,
#              output           = paste('./mgx/results/maaslin2/path_parasutterella/', '/', sep=''),
#              normalization    = 'TSS',
#              standardize      = 'TRUE',
#              transform        = 'LOG',
#              analysis_method  = 'LM',
#              min_abundance    = 0,
#              min_prevalence   = 0,
#              reference = c('bristolcat,Normal','abx,Non abx use'),
#              cores            = 4)
#     
#     
# # check parasutterella ecs
# ec <- read.table('./mgx/data/micron_ec_filtered_parasutterella.tsv', header=TRUE, sep='\t')
# rownames(ec) <- ec$UID
# ec<- ec[,-1]
# Maaslin2(input_data       = ec,
#          input_metadata   = meta,
#          output           = paste('./mgx/results/maaslin2/ec_parasutterella/', '/', sep=''),
#          normalization    = 'TSS',
#          standardize      = 'TRUE',
#          transform        = 'LOG',
#          analysis_method  = 'LM',
#          min_abundance    = 0,
#          min_prevalence   = 0,
#          reference = c('bristolcat,Normal','abx,Non abx use'),
         # cores            = 4)

#############################################################################################################################################
                                                #### MaAsLin2 ec #####
#############################################################################################################################################

# log-transform, default q=0.25
    Maaslin2(input_data       = ec,
             input_metadata   = meta,
             output           = paste('./mgx/results/maaslin2/ec/', '/', sep=''),
             normalization    = 'TSS',
             standardize      = 'TRUE',
             transform        = 'LOG',
             analysis_method  = 'LM',
             min_abundance    = 0,
             min_prevalence   = 0,
            reference = c('bristolcat,Normal','abx,Non abx use'),
             cores            = 4)

   Maaslin2(input_data       = ec,
            input_metadata   = meta_short,
            output           = paste('./mgx/results/maaslin2/ec_reduced/', '/', sep=''),
            normalization    = 'TSS',
            standardize      = 'TRUE',
            transform        = 'LOG',
            analysis_method  = 'LM',
            min_abundance    = 0,
            min_prevalence   = 0,
           reference = c('bristolcat,Normal','abx,Non abx use'),
            cores            = 4)

#############################################################################################################################################
                                                #### MaAsLin2 cazy dna #####
#############################################################################################################################################

# # log-transform, default q=0.25
#     Maaslin2(input_data       = cazy,
#              input_metadata   = meta,
#              output           = paste('./mgx/results/maaslin2/cazy/', '/', sep=''),
#              normalization    = 'TSS',
#              standardize      = 'TRUE',
#              transform        = 'LOG',
#              analysis_method  = 'LM',
#              min_abundance    = 0,
#              min_prevalence   = 0,
#             reference = c('bristolcat,Normal','abx,Non abx use'),
#              cores            = 4)
# 
#    Maaslin2(input_data       = cazy,
#             input_metadata   = meta_short,
#             output           = paste('./mgx/results/maaslin2/cazy_reduced/', '/', sep=''),
#             normalization    = 'TSS',
#             standardize      = 'TRUE',
#             transform        = 'LOG',
#             analysis_method  = 'LM',
#             min_abundance    = 0,
#             min_prevalence   = 0,
#            reference = c('bristolcat,Normal','abx,Non abx use'),
#             cores            = 4)


#############################################################################################################################################
                                              #### MaAsLin2 mbx #####
#############################################################################################################################################
mbx <- read.table('./mbx/data/mbx_maaslin_anno.tsv', header=TRUE, sep='\t')
rownames(mbx) <- mbx$sample
mbx <- mbx[,-1]
mbx <- mbx[order(rownames(mbx)),]

# mbx_class <- read.table('./mbx/data/class_maaslin.tsv', header=TRUE, sep='\t')
# rownames(mbx_class) <- mbx_class$sample
# mbx_class <- mbx_class[,-1]
# mbx_class <- mbx_class[order(rownames(mbx_class)),]
# 
# mbx_subclass <- read.table('./mbx/data/subclass_maaslin.tsv', header=TRUE, sep='\t')
# rownames(mbx_subclass) <- mbx_subclass$sample
# mbx_subclass <- mbx_subclass[,-1]
# mbx_subclass <- mbx_subclass[order(rownames(mbx_subclass)),]

meta <- read.table('./mbx/data/mbx_maaslin_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$sample
meta <- meta[,-1]
meta <- meta[order(rownames(meta)),]

meta_short <- meta[,colnames(meta) %in% c('caco','age','race8905','bristolcat','abx')]


# log-transform, default q=0.25
    Maaslin2(input_data       = mbx,
             input_metadata   = meta,
             output           = paste('./mbx/results/maaslin2/individual/full/', '/', sep=''),
             normalization    = 'TSS',
             standardize      = 'TRUE',
             transform        = 'LOG',
             analysis_method  = 'LM',
             min_abundance    = 0,
             min_prevalence   = 0,
             reference = c('bristolcat,Normal','abx,Non abx use'),
             cores            = 4)

   # Maaslin2(input_data       = mbx,
   #          input_metadata   = meta,
   #          output           = paste('./mbx/results/maaslin2/individual/full_noTSS/', '/', sep=''),
   #          normalization    = 'None',
   #          standardize      = 'TRUE',
   #          transform        = 'LOG',
   #          analysis_method  = 'LM',
   #          min_abundance    = 0,
   #          min_prevalence   = 0,
   #          reference = c('bristolcat,Normal','abx,Non abx use'),
   #          cores            = 4)

     Maaslin2(input_data       = mbx,
              input_metadata   = meta_short,
              output           = paste('./mbx/results/maaslin2/individual/reduced/', '/', sep=''),
              normalization    = 'TSS',
              standardize      = 'TRUE',
              transform        = 'LOG',
              analysis_method  = 'LM',
              min_abundance    = 0,
              min_prevalence   = 0,
              reference = c('bristolcat,Normal','abx,Non abx use'),
              cores            = 4)

    # Maaslin2(input_data       = mbx_class,
    #          input_metadata   = meta,
    #          output           = paste('./mbx/results/maaslin2/class/full/', '/', sep=''),
    #          normalization    = 'TSS',
    #          standardize      = 'TRUE',
    #          transform        = 'LOG',
    #          analysis_method  = 'LM',
    #          min_abundance    = 0,
    #          min_prevalence   = 0,
    #          reference = c('bristolcat,Normal','abx,Non abx use'),
    #          cores            = 4)
    # 
    # 
    #  Maaslin2(input_data       = mbx_class,
    #           input_metadata   = meta_short,
    #           output           = paste('./mbx/results/maaslin2/class/reduced/', '/', sep=''),
    #           normalization    = 'TSS',
    #           standardize      = 'TRUE',
    #           transform        = 'LOG',
    #           analysis_method  = 'LM',
    #           min_abundance    = 0,
    #           min_prevalence   = 0,
    #           reference = c('bristolcat,Normal','abx,Non abx use'),
    #           cores            = 4)
    # 
    # Maaslin2(input_data       = mbx_subclass,
    #          input_metadata   = meta,
    #          output           = paste('./mbx/results/maaslin2/subclass/full/', '/', sep=''),
    #          normalization    = 'TSS',
    #          standardize      = 'TRUE',
    #          transform        = 'LOG',
    #          analysis_method  = 'LM',
    #          min_abundance    = 0,
    #          min_prevalence   = 0,
    #          reference = c('bristolcat,Normal','abx,Non abx use'),
    #          cores            = 4)
    # 
    # 
    #  Maaslin2(input_data       = mbx_subclass,
    #           input_metadata   = meta_short,
    #           output           = paste('./mbx/results/maaslin2/subclass/reduced/', '/', sep=''),
    #           normalization    = 'TSS',
    #           standardize      = 'TRUE',
    #           transform        = 'LOG',
    #           analysis_method  = 'LM',
    #           min_abundance    = 0,
    #           min_prevalence   = 0,
    #           reference = c('bristolcat,Normal','abx,Non abx use'),
    #           cores            = 4)




#############################################################################################################################################
                                        #### heatmap for species (reduced + full)#####
#############################################################################################################################################
# combine results from maaslin for different models
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
  return(sig_result)
}

sig_crude <- sig_results(dir='./mgx/results/maaslin2/species_reduced/significant_results.tsv', exposure = "caco")
sig_adj   <- sig_results(dir='./mgx/results/maaslin2/species/significant_results.tsv', exposure = "caco")

sig_crude$model <- "Reduced"
sig_adj$model <- "Full"

all_results <- function(dir, exposure){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure, select =-c(value, N))
  return(all_result)
}

sigall_crude <- all_results(dir='./mgx/results/maaslin2/species_reduced/all_results.tsv', exposure = "caco")
sigall_adj   <- all_results(dir='./mgx/results/maaslin2/species/all_results.tsv', exposure = "caco")

join1<-full_join(sig_crude, sig_adj, by="feature")
sigfeature<-subset(join1, select=feature)

join1<-left_join(sigfeature, sigall_crude, by="feature")
join1$model <- "Reduced"
join2<-left_join(sigfeature, sigall_adj, by="feature")
join2$model <- "Full"

bind4<-rbind(join1,join2)

bind4$name<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", bind4$feature)
bind4 <- bind4 %>% arrange(bind4$model, bind4$name)
bug_names<-bind4[1:(nrow(bind4)/2), (ncol(bind4)-1):ncol(bind4)]

# create heatmap for diet-taxonomy associations
level_x_order <- factor(bind4$model, level = c('Reduced', 'Full'))
level_y_order <- factor(bind4$name, levels = bug_names$name)

bind4$stars <- cut(bind4$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))

# p1<-ggplot(bind4, aes(level_x_order, level_y_order)) +
#   geom_tile(aes(fill = coef), color = "white") +
#   scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
#   scale_size_continuous(range=c(2,30))+
#   geom_text(aes(label=stars), color="black", size=20, show.legend = TRUE) +
#   xlab("Model") +
#   scale_y_discrete(limits=rev)+
#   theme(panel.background = element_blank(),
#         legend.title = element_text(size = 50),
#         legend.text = element_text(size = 30),
#         legend.position = "bottom",
#         panel.spacing.y = unit(1.5, "lines"),
#         plot.title = element_text(size=50),
#         axis.title=element_text(size=50,face="bold"),
#         axis.title.y= element_blank(),
#         axis.text.y=element_text(size=50, face="bold.italic"),
#         axis.text.x = element_text(size=50, face="bold")) +
#   labs(fill = "β cofficient")+
#   guides(fill = guide_colorbar(barwidth = 30,
#                                barheight = 2))
# 
# ggsave(file="./mgx/results/maaslin2/heatmap/heatmap_species.jpg",plot=p1,width=25,height=25,dpi=300)


p1<-ggplot(bind4, aes(level_x_order, level_y_order)) +
  geom_tile(aes(fill = coef), color = "white",
            lwd = 1.5,
            linetype = 1) +
  coord_fixed(ratio = 0.5)+
  #scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  geom_text(aes(label=stars), color="black", size=20, show.legend = TRUE,vjust = 0.7) +
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
        axis.text.x = element_text(size=50, angle = 90, vjust = 0.5, hjust=1,face="bold")) +
  labs(fill = "β cofficient")+
  guides(fill = guide_colorbar(barwidth = 2,
                               barheight = 30))
ggsave(file="./mgx/results/maaslin2/heatmap/heatmap_species.jpg",plot=p1,width=25,height=25,dpi=300)

# #############################################################################################################################################
#                                         #### heatmap for species (full only)#####
# #############################################################################################################################################
# # combine results from maaslin for different models
# sig_results<- function(dir, exposure){
#   result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
#   sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
#   return(sig_result)
# }
# 
# sig_adj   <- sig_results(dir='./mgx/results/maaslin2/species/significant_results.tsv', exposure = "caco")
# 
# 
# all_results <- function(dir, exposure){
#   result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
#   all_result <- subset(result, metadata==exposure, select =-c(value, N))
#   return(all_result)
# }
# 
# sigall_adj   <- all_results(dir='./mgx/results/maaslin2/species/all_results.tsv', exposure = "caco")
# 
# sigfeature<-subset(sig_adj, select=feature)
# 
# join<-left_join(sigfeature, sigall_adj, by="feature")
# 
# bind4 <- join
# bind4$name<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", bind4$feature)
# bind4 <- bind4 %>% arrange(bind4$name)
# bug_names<-bind4[1:(nrow(bind4)/2), (ncol(bind4)-1):ncol(bind4)]
# 
# # create heatmap for diet-taxonomy associations
# level_y_order <- factor(bind4$name, levels = bug_names$name)
# 
# bind4$stars <- cut(bind4$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))
# p1<-ggplot(bind4, aes(level_y_order)) +
#   geom_tile(aes(fill = coef), color = "white") +
#   scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
#   geom_text(aes(label=stars), color="black", size=20, show.legend = TRUE) +
#   scale_y_discrete(limits=rev)+
#   theme(panel.background = element_blank(),
#         legend.title = element_text(size = 50),
#         legend.text = element_text(size = 30),
#         legend.position = "bottom",
#         plot.title = element_text(size=50),
#         axis.title=element_text(size=50,face="bold"),
#         axis.title.y= element_blank(),
#         axis.text.y=element_text(size=50, face="bold.italic"),
#         axis.text.x = element_text(size=50, face="bold")) +
#   labs(fill = "β cofficient")+
#   guides(fill = guide_colorbar(barwidth = 30,
#                                barheight = 2))
# ggsave(file="./mgx/results/maaslin2/heatmap/heatmap_species_full.jpg",plot=p1,width=25,height=25,dpi=300)


#############################################################################################################################################
                                        #### heatmap for pathways #####
#############################################################################################################################################
# combine results from maaslin for different models
# sig_results<- function(dir, exposure){
#   result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
#   sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
#   return(sig_result)
# }
# 
# sig_crude <- sig_results(dir='./mgx/results/maaslin2/path_reduced/significant_results.tsv', exposure = "caco")
# sig_adj   <- sig_results(dir='./mgx/results/maaslin2/path/significant_results.tsv', exposure = "caco")
# 
# sig_crude$model <- "Reduced"
# sig_adj$model <- "Full"
# 
# all_results <- function(dir, exposure){
#   result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
#   all_result <- subset(result, metadata==exposure, select =-c(value, N))
#   return(all_result)
# }
# 
# sigall_crude <- all_results(dir='./mgx/results/maaslin2/path_reduced/all_results.tsv', exposure = "caco")
# sigall_adj   <- all_results(dir='./mgx/results/maaslin2/path/all_results.tsv', exposure = "caco")
# 
# join1<-full_join(sig_crude, sig_adj, by="feature")
# sigfeature<-subset(join1, select=feature)
# 
# join1<-left_join(sigfeature, sigall_crude, by="feature")
# join1$model <- "Reduced"
# join2<-left_join(sigfeature, sigall_adj, by="feature")
# join2$model <- "Full"
# 
# bind4<-rbind(join1,join2)
# 
# bind4$name<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", bind4$feature)
# bind4 <- bind4 %>% arrange(bind4$model, bind4$name)
# bug_names<-bind4[1:(nrow(bind4)/2), (ncol(bind4)-1):ncol(bind4)]
# 
# # create heatmap for diet-taxonomy associations
# level_x_order <- factor(bind4$model, level = c('Reduced', 'Full'))
# level_y_order <- factor(bind4$name, levels = bug_names$name)
# 
# bind4$stars <- cut(bind4$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))
# p1<-ggplot(bind4, aes(level_x_order, level_y_order)) +
#   geom_tile(aes(fill = coef), color = "white") +
#   scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
#   geom_text(aes(label=stars), color="black", size=20, show.legend = TRUE) +
#   xlab("Model") +
#   scale_y_discrete(limits=rev)+
#   theme(panel.background = element_blank(),
#         legend.title = element_text(size = 50),
#         legend.text = element_text(size = 30),
#         legend.position = "bottom",
#         legend.justification = "bottom",
#         plot.title = element_text(size=50),
#         axis.title=element_text(size=50,face="bold"),
#         axis.title.y= element_blank(),
#         axis.text.y=element_text(size=50, face="bold.italic"),
#         axis.text.x = element_text(size=50, face="bold")) +
#   labs(fill = "β cofficient")+
#   guides(fill = guide_colorbar(barwidth = 30,
#                                barheight = 2))
# ggsave(file="./mgx/results/maaslin2/heatmap/heatmap_path.jpg",plot=p1,width=42,height=16,dpi=300)


#############################################################################################################################################
                                        #### heatmap for mbx #####
#############################################################################################################################################
# combine results from maaslin for different models
sig_results<- function(dir, exposure){
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
  return(sig_result)
}

sig_crude <- sig_results(dir='./mbx/results/maaslin2/individual/reduced/significant_results.tsv', exposure = "caco")
sig_adj   <- sig_results(dir='./mbx/results/maaslin2/individual/full/significant_results.tsv', exposure = "caco")

sig_crude$model <- "Reduced"
sig_adj$model <- "Full"

all_results <- function(dir, exposure){
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure, select =-c(value, N))
  return(all_result)
}

sigall_crude <- all_results(dir='./mbx/results/maaslin2/individual/reduced/all_results.tsv', exposure = "caco")
sigall_adj   <- all_results(dir='./mbx/results/maaslin2/individual/full/all_results.tsv', exposure = "caco")

join1<-full_join(sig_crude, sig_adj, by="feature")
sigfeature<-subset(join1, select=feature)

join1<-left_join(sigfeature, sigall_crude, by="feature")
join1$model <- "Reduced"
join2<-left_join(sigfeature, sigall_adj, by="feature")
join2$model <- "Full"

bind4<-rbind(join1,join2)

# remove duplicated metabolites (same HMDB IDs)
# Cer 18:1;O2/16:0 QI06417_HILIC_pos
# Cer 18:1;O2/16:0    QI18714_C8_pos


bind4 <- bind4[bind4$feature!="QI06417_HILIC_pos",]

# Two metabolites have different HMDB IDs and I will keep QI18892_HILIC_pos
# 1,7-Dimethyluric acid QI18892_HILIC_pos
# 1,7-Dimethyluric acid QI6992_C18_neg

# L-Urobilin    TF04_HILIC_pos
#L-Urobilin   QI49108_C18_neg

bind4 <- bind4[bind4$feature!="QI6992_C18_neg",]
bind4 <- bind4[bind4$feature!="QI49108_C18_neg",]

# get annotations
anno <- read.table('./mbx/data/mbx_maaslin_anno_name.tsv',header=TRUE,sep='\t')
anno$feature <- anno$mol_id

bind4_anno <- left_join(bind4,anno,by="feature")

bind4_anno <- bind4_anno %>% arrange(bind4_anno$model, bind4_anno$Metabolite.name)
#bind4_anno$name <- paste(bind4_anno$Metabolite.name, bind4_anno$mol_id)

bug_names<-bind4_anno[1:(nrow(bind4_anno)/2), (ncol(bind4_anno)-1):ncol(bind4_anno)]

# create heatmap for diet-taxonomy associations
level_x_order <- factor(bind4_anno$model, level = c('Reduced', 'Full'))
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
        axis.text.y=element_text(size=50),
        axis.text.x = element_text(size=50, angle = 90, vjust = 0.5, hjust=1, face="bold")) +
  labs(fill = "β cofficient")+
  guides(fill = guide_colorbar(barwidth = 2,
                               barheight = 30))
ggsave(file="./mbx/results/maaslin2/individual/heatmap/heatmap_mbx.jpg",plot=p1,width=20,height=45,dpi=300)




