#########################################################################################################
### Purpose: Boxplot of selected species/pathways/ecs according to case/control status and severity
### Last update date: 10/13/2023 
#########################################################################################################

setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")


# salloc -p test --mem 64G -t 0-08:00
# module load R/4.2.2-fasrc01
# #install new packages
# export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
# R --quiet
# 
# rm(list=ls())
# getwd()
# setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")

library(grid)
library(ggplot2)
library(gridExtra)
library(dplyr)

# read in taxonomy data
species_all <- read.table('./mgx/data/micron_species_filtered.tsv', header=TRUE, sep='\t')
rownames(species_all) <- species_all$UID
species_all <- species_all[,-1]

# Log transformation
LOG <- function(x) {
    y <- replace(x, x == 0, min(x[x>0]) / 2)
    return(log2(y))
}

species_log <- apply(species_all, 2, LOG)

# Arc Sine Square Root Transformation
AST <- function(x) {
    y <- sign(x) * asin(sqrt(abs(x)))
    if(any(is.na(y))) {
        logging::logerror(
            paste0("AST transform is only valid for values between -1 and 1. ",
                   "Please select an appropriate normalization option or ",
                   "normalize your data prior to running."))
        stop()
    }
    return(y)
}
species_ast <- apply(species_all,2,AST)


# read in metadata
meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID
meta <- meta[,-1]

species_meta <- merge(species_ast,meta,by='row.names')
rownames(species_meta) <- species_meta[,1]
species_meta <- species_meta[,-1]

meta_severe <- read.table('./metadata/micron_metadata_severe.tsv', header=TRUE, sep='\t')
rownames(meta_severe) <- meta_severe$UID
meta_severe <- meta_severe[,-1]

species_meta <- merge(species_ast,meta,by='row.names')
rownames(species_meta) <- species_meta[,1]
species_meta <- species_meta[,-1]

severe <- as.data.frame(meta_severe[,colnames(meta_severe) %in% c('caco3')])
rownames(severe) <- rownames(meta_severe)
colnames(severe) <- 'caco3'
species_meta_2 <- merge(species_meta,severe,by='row.names')
rownames(species_meta_2) <- species_meta_2[,1]
species_meta_2 <- species_meta_2[,-1]



#############################################################################################################################################
                                        #### barplot for species according to caco#####
#############################################################################################################################################
scatter_plot <- function(y,ylab){
p1 <- ggplot(species_meta_2,aes_string(x='caco', y=y,fill='caco'))+
  geom_point(size=8, position=position_jitterdodge(),aes(fill=caco),alpha=1,shape=21,stroke=1,color='black',show.legend = FALSE)+
  geom_boxplot(alpha=0.5,size=1.5)+
  scale_fill_manual(values=c("#E69F00", "#4169E1"))+
  theme_classic() +
  xlab("Status") + ylab(ylab) +
  guides(legend.position=NULL)+
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 30),
      legend.position="none",
      axis.title=element_text(size=30,face="bold"),
      axis.text=element_text(size=30),
      axis.ticks=element_blank(),
      axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
      )
filepath <- paste('./mgx/results/maaslin2/boxplot/species/caco/', y, '.png', sep ="")
ggsave(filename=filepath, plot=p1, width = 10, height = 10, dpi = 300)
}

scatter_plot('Clostridium_citroniae',"Clostridium citroniae")
scatter_plot('Clostridium_lavalense',"Clostridium lavalense")
scatter_plot('Clostridium_clostridioforme',"Clostridium clostridioforme")
scatter_plot('Blautia_producta',"Blautia producta")
scatter_plot('Ruminococcus_gnavus',"Ruminococcus gnavus")
scatter_plot('Hungatella_hathewayi',"Hungatella hathewayi")
scatter_plot('Erysipelatoclostridium_ramosum',"Erysipelatoclostridium ramosum")
scatter_plot('Parasutterella_excrementihominis',"Parasutterella excrementihominis")
scatter_plot('Subdoligranulum_sp_APC924_74',"Subdoligranulum sp APC924 74")
scatter_plot('Eubacterium_eligens',"Eubacterium eligens")



#############################################################################################################################################
                                        #### barplot for species according to severity#####
#############################################################################################################################################
my_comparisons <- list( c("Control", "Non-severe"), c("Control", "Severe"), c("Non-severe", "Severe") )
scatter_plot <- function(y,ylab){
p1 <- ggplot(species_meta_2,aes_string(x='caco3', y=y,fill='caco3'))+
  geom_point(size=8, position=position_jitterdodge(),aes(fill=caco3),alpha=1,shape=21,stroke=1,color='black',show.legend = FALSE)+
  geom_boxplot(alpha=0.5,size=1.5)+
  scale_fill_manual(values=c("#E69F00", "#4169E1","darkslateblue"))+
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  xlab("Status") + ylab(ylab) +
  guides(legend.position=NULL)+
  stat_compare_means(comparisons = my_comparisons,
    aes(label = paste0("p = ", ..p.format..)), size=12)+
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 30),
      legend.position="none",
      axis.title=element_text(size=30,face="bold"),
      axis.text=element_text(size=30),
      axis.ticks=element_blank(),
      axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
      )
filepath <- paste('./mgx/results/maaslin2/boxplot/species/severity/', y, '.png', sep ="")
ggsave(filename=filepath, plot=p1, width = 10, height = 12, dpi = 300)
}
scatter_plot('Anaerotruncus_colihominis','Anaerotruncus colihominis')

#############################################################################################################################################
                                        #### barplot for pathways #####
#############################################################################################################################################
# read in path data
path <- read.table('./mgx/data/micron_path_filtered.tsv', header=TRUE, sep='\t')
rownames(path) <- path$UID
path <- path[,-1]

meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID
meta <- meta[,-1]

path_meta <- merge(path,meta,by='row.names')
rownames(path_meta) <- path_meta[,1]
path_meta <- path_meta[,-1]

scatter_plot <- function(y,ylab){
p1 <- ggplot(path_meta,aes_string(x='caco', y=y,fill='caco'))+
  geom_point(size=8, position=position_jitterdodge(),aes(fill=caco),alpha=1,shape=21,stroke=1,color='black',show.legend = FALSE)+
  geom_boxplot(alpha=0.5,size=1.5)+
  scale_fill_manual(values=c("#E69F00", "#4169E1"))+
  theme_classic() +
  xlab("Status") + ylab(ylab) +
  guides(legend.position=NULL)+
  #coord_cartesian(clip = "off")+
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 30),
      legend.position="none",
      axis.title=element_text(size=30,face="bold"),
      axis.text=element_text(size=30),
      axis.ticks=element_blank(),
      axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
      )
filepath <- paste('./mgx/results/maaslin2/boxplot/path/', y, '.png', sep ="")
ggsave(filename=filepath, plot=p1, width = 10, height = 12, dpi = 300)
}

scatter_plot('PWY.5677..succinate.fermentation.to.butanoate',"PWY-5677 Succinate fermentation to butanoate")
scatter_plot('PWY66.429..fatty.acid.biosynthesis.initiation..mitochondria.',"PWY66-429 Fatty acid biosynthesis initiation mitochondria")


#############################################################################################################################################
                                        #### density plot for pathways #####
#############################################################################################################################################
density_plot <- function(x,name){
p2 <- ggplot(data=path_meta, aes_string(x=x, group='caco', fill='caco')) +
    aes(y = after_stat(scaled))+
    geom_density(alpha=0.7,color=NA)+
    scale_fill_manual(values=c("#E69F00", "#4169E1"))+
    theme_classic() +
    guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 30),
        legend.position="none",
        axis.title=element_text(size=30,face="bold"),
        axis.text=element_text(size=30),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
        )
    filepath <- paste('./mgx/results/maaslin2/density_plot/path/',name,'_density.png',sep ="")
    ggsave(filename=filepath, plot=p2, width = 15, height = 12, dpi = 300)
}
density_plot('PWY66.429..fatty.acid.biosynthesis.initiation..mitochondria.',"PWY66-429")
density_plot('PWY.5677..succinate.fermentation.to.butanoate',"PWY-5677")
density_plot('PWY4FS.8..phosphatidylglycerol.biosynthesis.II..non.plastidic.',"PWY4FS-8")
density_plot('PWY4FS.7..phosphatidylglycerol.biosynthesis.I..plastidic.',"PWY4FS-7")



#############################################################################################################################################
                                        #### barplot for ecs #####
#############################################################################################################################################
# read in ec data
ec <- read.table('./mgx/data/micron_ec_filtered.tsv', header=TRUE, sep='\t')
rownames(ec) <- ec$UID
ec <- ec[,-1]

meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID
meta <- meta[,-1]

ec_meta <- merge(ec,meta,by='row.names')
rownames(ec_meta) <- ec_meta[,1]
ec_meta <- ec_meta[,-1]

scatter_plot <- function(y,ylab){
p1 <- ggplot(ec_meta,aes_string(x='caco', y=y,fill='caco'))+
  geom_point(size=8, position=position_jitterdodge(),aes(fill=caco),alpha=1,shape=21,stroke=1,color='black',show.legend = FALSE)+
  geom_boxplot(alpha=0.5,size=1.5)+
  scale_fill_manual(values=c("#E69F00", "#4169E1"))+
  theme_classic() +
  xlab("Status") + ylab(ylab) +
  guides(legend.position=NULL)+
  coord_cartesian(clip = "off")+
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 30),
      legend.position="none",
      axis.title=element_text(size=30,face="bold"),
      axis.text=element_text(size=30),
      axis.ticks=element_blank(),
      axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
      )
filepath <- paste('./mgx/results/maaslin2/boxplot/ec/', y, '.png', sep ="")
ggsave(filename=filepath, plot=p1, width = 10, height = 10, dpi = 300)
}

scatter_plot('X6.4.1.2',"Acetyl-CoA carboxylase")
scatter_plot('X6.3.4.14',"Biotin carboxylase")
scatter_plot('X1.4.1.16',"Diaminopimelate dehydrogenase")
scatter_plot('X1.3.1.31',"2-enoate reductase")

scatter_plot('X4.1.3.6',"citrate (pro-3S)-lyase")
scatter_plot('X6.1.1.17',"Glutamate-tRNA ligase")
scatter_plot('X2.8.3.10',"citrate CoA-transferase")

scatter_plot('X2.6.1.1',"aspartate transaminase")
scatter_plot('X2.6.1.37',"2-aminoethylphosphonate- \n pyruvate transaminase")

#############################################################################################################################################
                                        #### density plots for ecs #####
#############################################################################################################################################
density_plot <- function(x,name){
p2 <- ggplot(data=ec_meta, aes_string(x=x, group='caco', fill='caco')) +
    aes(y = after_stat(scaled))+
    geom_density(alpha=0.7,color=NA)+
    scale_fill_manual(values=c("#E69F00", "#4169E1"))+
    theme_classic() +
    #guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 30),
        #legend.position="none",
        axis.title=element_text(size=30,face="bold"),
        axis.text=element_text(size=30),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
        )
    filepath <- paste('./mgx/results/maaslin2/density_plot/ec/',name,'_density.png',sep ="")
    ggsave(filename=filepath, plot=p2, width = 15, height = 12, dpi = 300)
}
density_plot('X6.4.1.2',"Acetyl-CoA carboxylase")
density_plot('X6.3.4.14',"Biotin carboxylase")

all <- merge(ec,path,by='row.names')
rownames(all) <- all[,1]
all <- all[,-1]

all <- merge(species_ast, all, by='row.names')
rownames(all) <- all[,1]
all <- all[,-1]

library(Hmisc)
rcorr(all$`X6.4.1.2`,all$`PWY66.429..fatty.acid.biosynthesis.initiation..mitochondria.`) #0.33
rcorr(all$Parasutterella_excrementihominis,all$`PWY66.429..fatty.acid.biosynthesis.initiation..mitochondria.`) #0.03
rcorr(all$`PWY.5677..succinate.fermentation.to.butanoate`,all$`PWY66.429..fatty.acid.biosynthesis.initiation..mitochondria.`) #0.13
rcorr(all$Ruminococcus_gnavus,all$`PWY4FS.8..phosphatidylglycerol.biosynthesis.II..non.plastidic.`) #0.12
rcorr(all$Fusicatenibacter_saccharivorans,all$`PWY4FS.8..phosphatidylglycerol.biosynthesis.II..non.plastidic.`) #-0.21
rcorr(all$Escherichia_coli,all$`PWY4FS.8..phosphatidylglycerol.biosynthesis.II..non.plastidic.`) #0.26
