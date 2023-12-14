##############################################################################################################################################
# R programs for Permanova
###                           permanova                            ####
### outcome:    species
### predictors: diverticulitis
### last update date: 10/6/2023 
#############################################################################################################################################
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

#setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")

library(scales)
library(tidyverse)
library(vegan)

#####################################
#######     Species           #######
#####################################

# read in taxonomy data
species_all <- read.table('./mgx/data/micron_species_filtered.tsv', header=TRUE, sep='\t')
rownames(species_all) <- species_all$UID
species_all <- species_all[,-1]

# read in metadata
meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID

species_all <- species_all[rownames(meta),]


# PERMANOVA using adonis2
caco_species <- as.data.frame(adonis2(species_all ~ caco, data = meta, permutations = 999, method = "bray"))
# R2=0.009265961, P=0.002



# x <- adonis(formula = taxa_bc ~ age, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# age <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ race8905, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# race <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ bmi19, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# bmi <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ abx, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# abx <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ bristolcat, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# bristol <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ smk19, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# smk <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ pmh, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# pmh <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ act17m, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# act <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ calor15n, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# calor <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ alco15n, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# alco <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ aofib15a, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# aofib <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ ahei2010_noETOH15, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# ahei <- x.data[!is.na(x.data$F.Model), ]
# 
# 
# # merge PERMANOVA results
# all_taxonomy <-rbind(caco,age,race,bmi,abx,bristol,smk,pmh,act,calor,alco,aofib,ahei)
# # label for variables
# component<- c('Diverticulitis diagnosis',
#               'Age',
#               'Race',
#               'Body mass index',
#               'Antibiotic use',
#               'Bristol stool scale',
#               'Smoking',
#               'Postmenopausal hormone use',
#               'Physical activity',
#               'Calorie intake',
#               'Alcohol intake',
#               'Dietary fiber',
#               'AHEI index')
# 
# all_taxonomy<-cbind(all_taxonomy, component)
# #all_taxonomy_diet<-all_taxonomy[1:5,]
# #all_taxonomy_diet$cat<-"Vitamin D status and components"
# #all_taxonomy_cov<-all_taxonomy[6:11,]
# #all_taxonomy_cov$cat<-"Covariates"
# #all_taxonomy_bio<-all_taxonomy[12:17,]
# #all_taxonomy_bio$cat<-"Biomarkers"
# #all_taxonomy_tax<-rbind(all_taxonomy_diet, all_taxonomy_cov, all_taxonomy_bio)
# #all_taxonomy_tax$cat<-as.factor(all_taxonomy_tax$cat)
# all_taxonomy_tax <- all_taxonomy
# # FDR adjustment for p-values
# all_taxonomy_tax$fdr <-
#   p.adjust(
#     all_taxonomy_tax$`Pr(>F)`,
#     method = "BH",
#     n = length(all_taxonomy_tax$`Pr(>F)`)
#   )
# # create figure to show PERMANOVA results
# all_taxonomy_tax$stars <- cut(all_taxonomy_tax$fdr, breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))
# level_y_order <- factor(all_taxonomy_tax$component, level = rev(all_taxonomy_tax$component))
# #color_code<-c("Vitamin D status and components"="#00AFBB",
#               #"Covariates"="#999999",
#               #"Biomarkers"="#E7B800")
# 
# plot <- ggplot(all_taxonomy_tax, aes(level_y_order, R2)) +
#   geom_bar(stat="identity",fill="cornflowerblue",color="black")+
#   scale_y_continuous(labels=scales::percent, limits = c(0, 0.02))+
#   geom_text(aes(label=stars), color="black", size=20) +
#   #(values = color_code)+
#   coord_flip()+
#   theme(panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          panel.background = element_blank(),
#          panel.border = element_rect(colour = "black", fill=NA, size=1),
#          legend.position = "bottom",
#         legend.text = element_text(size = 20),
#         legend.title = element_blank(),
#         axis.ticks =element_blank(),
#         plot.title = element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y= element_blank(),
#         axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
#         axis.text.y=element_text(size=20, face="bold"),
#         axis.text.x = element_text(size=20))
# print(plot)
# ggsave(filename='./mgx/results/permanova/permanova_species.jpg', plot=plot, width =15, height = 10, dpi = 300)
# write.csv(all_taxonomy_tax, file = "./mgx/results/permanova/permanova_species.csv")
# 
# 
# 
# ################################################################################################################
# 
# 
# 
# #####################################
# #######     path           #######
# #####################################
# 
# # read in pathway data
# path <- read.table('./mgx/data/micron_path_filtered.tsv', header=TRUE, sep='\t')
# rownames(path) <- path$UID
# path <- path[,-1]
# 
# 
# # read in metadata
# meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
# rownames(meta) <- meta$UID
# 
# path <- path[rownames(meta),]
# 
# # calculate bray-curtis
# taxa_bc <- vegdist(path, "bray")
# 
# 
# # PERMANOVA using adonis
# x <- adonis(formula = taxa_bc ~ caco, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# caco <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ age, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# age <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ race8905, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# race <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ bmi19, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# bmi <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ abx, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# abx <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ bristolcat, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# bristol <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ smk19, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# smk <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ pmh, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# pmh <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ act17m, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# act <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ calor15n, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# calor <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ alco15n, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# alco <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ aofib15a, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# aofib <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ ahei2010_noETOH15, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# ahei <- x.data[!is.na(x.data$F.Model), ]
# 
# 
# # merge PERMANOVA results
# all_taxonomy <-rbind(caco,age,race,bmi,abx,bristol,smk,pmh,act,calor,alco,aofib,ahei)
# # label for variables
# component<- c('Diverticulitis diagnosis',
#               'Age',
#               'Race',
#               'Body mass index',
#               'Antibiotic use',
#               'Bristol stool scale',
#               'Smoking',
#               'Postmenopausal hormone use',
#               'Physical activity',
#               'Calorie intake',
#               'Alcohol intake',
#               'Dietary fiber',
#               'AHEI index')
# 
# all_taxonomy<-cbind(all_taxonomy, component)
# #all_taxonomy_diet<-all_taxonomy[1:5,]
# #all_taxonomy_diet$cat<-"Vitamin D status and components"
# #all_taxonomy_cov<-all_taxonomy[6:11,]
# #all_taxonomy_cov$cat<-"Covariates"
# #all_taxonomy_bio<-all_taxonomy[12:17,]
# #all_taxonomy_bio$cat<-"Biomarkers"
# #all_taxonomy_tax<-rbind(all_taxonomy_diet, all_taxonomy_cov, all_taxonomy_bio)
# #all_taxonomy_tax$cat<-as.factor(all_taxonomy_tax$cat)
# all_taxonomy_tax <- all_taxonomy
# # FDR adjustment for p-values
# all_taxonomy_tax$fdr <-
#   p.adjust(
#     all_taxonomy_tax$`Pr(>F)`,
#     method = "BH",
#     n = length(all_taxonomy_tax$`Pr(>F)`)
#   )
# # create figure to show PERMANOVA results
# all_taxonomy_tax$stars <- cut(all_taxonomy_tax$fdr, breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))
# level_y_order <- factor(all_taxonomy_tax$component, level = rev(all_taxonomy_tax$component))
# #color_code<-c("Vitamin D status and components"="#00AFBB",
#               #"Covariates"="#999999",
#               #"Biomarkers"="#E7B800")
# 
# plot <- ggplot(all_taxonomy_tax, aes(level_y_order, R2)) +
#   geom_bar(stat="identity")+
#   scale_y_continuous(labels=scales::percent, limits = c(0, 0.03))+
#   geom_text(aes(label=stars), color="black", size=20) +
#   #(values = color_code)+
#   coord_flip()+
#   theme(panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          panel.background = element_blank(),
#          legend.position = "bottom",
#         legend.text = element_text(size = 20),
#         legend.title = element_blank(),
#         axis.ticks =element_blank(),
#         plot.title = element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y= element_blank(),
#         axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
#         axis.text.y=element_text(size=20, face="bold"),
#         axis.text.x = element_text(size=20))
# ggsave(filename='./mgx/results/permanova/permanova_path.jpg', plot=plot, width =15, height = 10, dpi = 300)
# write.csv(all_taxonomy_tax, file = "./mgx/results/permanova/permanova_path.csv")
# 
# 
# 
# #####################################
# #######     EC           #######
# #####################################
# 
# # read in pathway data
# ec <- read.table('./mgx/data/micron_ec_filtered.tsv', header=TRUE, sep='\t')
# rownames(ec) <- ec$UID
# ec <- ec[,-1]
# 
# 
# # read in metadata
# meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
# rownames(meta) <- meta$UID
# 
# ec <- ec[rownames(meta),]
# 
# # calculate bray-curtis
# taxa_bc <- vegdist(ec, "bray")
# 
# 
# # PERMANOVA using adonis
# x <- adonis(formula = taxa_bc ~ caco, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# caco <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ age, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# age <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ race8905, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# race <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ bmi19, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# bmi <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ abx, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# abx <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ bristolcat, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# bristol <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ smk19, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# smk <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ pmh, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# pmh <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ act17m, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# act <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ calor15n, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# calor <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ alco15n, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# alco <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ aofib15a, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# aofib <- x.data[!is.na(x.data$F.Model), ]
# 
# x <- adonis(formula = taxa_bc ~ ahei2010_noETOH15, data = meta, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# ahei <- x.data[!is.na(x.data$F.Model), ]
# 
# 
# # merge PERMANOVA results
# all_taxonomy <-rbind(caco,age,race,bmi,abx,bristol,smk,pmh,act,calor,alco,aofib,ahei)
# # label for variables
# component<- c('Diverticulitis diagnosis',
#               'Age',
#               'Race',
#               'Body mass index',
#               'Antibiotic use',
#               'Bristol stool scale',
#               'Smoking',
#               'Postmenopausal hormone use',
#               'Physical activity',
#               'Calorie intake',
#               'Alcohol intake',
#               'Dietary fiber',
#               'AHEI index')
# 
# all_taxonomy<-cbind(all_taxonomy, component)
# #all_taxonomy_diet<-all_taxonomy[1:5,]
# #all_taxonomy_diet$cat<-"Vitamin D status and components"
# #all_taxonomy_cov<-all_taxonomy[6:11,]
# #all_taxonomy_cov$cat<-"Covariates"
# #all_taxonomy_bio<-all_taxonomy[12:17,]
# #all_taxonomy_bio$cat<-"Biomarkers"
# #all_taxonomy_tax<-rbind(all_taxonomy_diet, all_taxonomy_cov, all_taxonomy_bio)
# #all_taxonomy_tax$cat<-as.factor(all_taxonomy_tax$cat)
# all_taxonomy_tax <- all_taxonomy
# # FDR adjustment for p-values
# all_taxonomy_tax$fdr <-
#   p.adjust(
#     all_taxonomy_tax$`Pr(>F)`,
#     method = "BH",
#     n = length(all_taxonomy_tax$`Pr(>F)`)
#   )
# # create figure to show PERMANOVA results
# all_taxonomy_tax$stars <- cut(all_taxonomy_tax$fdr, breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))
# level_y_order <- factor(all_taxonomy_tax$component, level = rev(all_taxonomy_tax$component))
# #color_code<-c("Vitamin D status and components"="#00AFBB",
#               #"Covariates"="#999999",
#               #"Biomarkers"="#E7B800")
# 
# plot <- ggplot(all_taxonomy_tax, aes(level_y_order, R2)) +
#   geom_bar(stat="identity")+
#   scale_y_continuous(labels=scales::percent, limits = c(0, 0.03))+
#   geom_text(aes(label=stars), color="black", size=20) +
#   #(values = color_code)+
#   coord_flip()+
#   theme(panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          panel.background = element_blank(),
#          panel.border = element_rect(colour = "black", fill=NA, size=2),
#          legend.position = "bottom",
#         legend.text = element_text(size = 20),
#         legend.title = element_blank(),
#         axis.ticks =element_blank(),
#         plot.title = element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y= element_blank(),
#         axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
#         axis.text.y=element_text(size=20, face="bold"),
#         axis.text.x = element_text(size=20))
# ggsave(filename='./mgx/results/permanova/permanova_ec.jpg', plot=plot, width =15, height = 10, dpi = 300)
# write.csv(all_taxonomy_tax, file = "./mgx/results/permanova/permanova_ec.csv")
