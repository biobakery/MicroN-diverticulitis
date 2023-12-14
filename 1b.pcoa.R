##############################################################################################################################################
# Purpose: PCA of Bray-Curtis dissmilarities 
# Last update date: 10/13/2023
#############################################################################################################################################

# salloc -p test --mem 64G -t 0-08:00
# module load gcc/9.3.0-fasrc01 R_packages/4.0.5-fasrc02
# 
# #install new packages
# export R_LIBS_USER=$HOME/apps/R_4.0.5:$R_LIBS_USER
# R --quiet
# 
 
# rm(list=ls())
# getwd()
# setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
# 

setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")


library(vegan)
library(ggplot2)
library(gridExtra)

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



# read in metadata
meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID


# metadata
meta_pcoa <- as.data.frame(meta$caco)
colnames(meta_pcoa) <- c('Status')
rownames(meta_pcoa) <- rownames(meta)





# bug subset
#bug_subset <- species_all[,colnames(species_all) %in% c('ruminococcus_bromii','eubacterium_rectale','prevotella_copri')]

# append bug_subset to metadata
#meta_pcoa <- merge(bug_subset, meta_pcoa, by = 'row.names')
#rownames(meta_pcoa) <- meta_pcoa[ ,1]
#meta_pcoa[ ,1] <- NULL

#######################################################################################
###########
# species #
###########
#### Note: I need to make sure that samples in rows
bugs_pcoa <- species_all
#samples in rows
bugs_pcoa <- capscale(bugs_pcoa ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )

write.csv(ord.bug.scores, file = "./mgx/data/pco.loading.csv")


# append metadata to scores
ord.scores.meta <- merge(ord.scores, meta_pcoa, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.scores.meta) <- ord.scores.meta[ ,1]
ord.scores.meta[ ,1] <- NULL

# number of samples
n_samples <- nrow(ord.scores.meta)

# PCoA colored by disease status

# Generate ordination plot objects and store all in list object
  # plot ordination
library(ggforce)
pcoa.plot.species <-
  ggplot( ord.scores.meta, aes(MDS1, MDS2) ) +
  geom_point(size = 5, alpha = 1, aes( fill =  Status),color='black',shape=21) +
  # geom_mark_ellipse(aes(color = Status,
  #                       label=Status),
  #                   expand = unit(0.5,"mm"),
  #                   label.buffer = unit(-5, 'mm'))+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 16, face = 'bold'),
        legend.text = element_text(size=12), legend.title = element_blank(),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'horizontal',
        plot.title = element_text( size = 20, face = 'bold')
        )+
        scale_fill_manual(values=c("#E69F00", "#4169E1"))+
        xlab(paste("PCo1 (",x_variance.bugs,"%)",sep = "")) + ylab(paste("PCo2 (",y_variance.bugs,"%)",sep = ""))
ggsave("./mgx/results/pcoa/pcoa_species.jpg", plot=pcoa.plot.species, width =7, height =6, dpi=300)

###########################################
# pcoa according to F prau
# append species to scores
ord.scores.species <- merge(ord.scores, species_all, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.scores.species) <- ord.scores.species[ ,1]
ord.scores.species[ ,1] <- NULL

library(viridis)
pcoa.plot.fprau <-
  ggplot( ord.scores.species, aes(MDS1, MDS2) ) +
  geom_point(size = 5, alpha = 1, aes( fill=  Faecalibacterium_prausnitzii),shape=21,color='black') +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 16, face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'horizontal',
        plot.title = element_text( size = 20, face = 'bold')
        )+
        scale_fill_gradientn(colours = viridis(10,direction=-1))+
        xlab(paste("PCo1 (",x_variance.bugs,"%)",sep = "")) + ylab(paste("PCo2 (",y_variance.bugs,"%)",sep = ""))
ggsave("./mgx/results/pcoa/pcoa_fprau.jpg", plot=pcoa.plot.fprau, width =7, height =6, dpi=300)


#############################################################################################################################################
########
# path #
########
#### Note: I need to make sure that samples in rows
bugs_pcoa <- path
#samples in rows
bugs_pcoa <- capscale(bugs_pcoa ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )

# append metadata to scores
ord.scores.meta <- merge(ord.scores, meta_pcoa, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.scores.meta) <- ord.scores.meta[ ,1]
ord.scores.meta[ ,1] <- NULL

# number of samples
n_samples <- nrow(ord.scores.meta)


# PCoA colored by disease status

pcoa.plot.path <-
  ggplot( ord.scores.meta, aes(MDS1, MDS2) ) +
  geom_point(size = 5, alpha = 1, aes( fill =  Status),color='black',shape=21) +
  theme_classic() +
  # geom_mark_ellipse(aes(fill= Status,
  #                       label=Status),
  #                   expand = unit(0.5,"mm"),
  #                   label.buffer = unit(-5, 'mm'))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 16, face = 'bold'),
        legend.text = element_text(size=12), legend.title = element_blank(),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'horizontal',
        plot.title = element_text( size = 20, face = 'bold')
        )+
        scale_fill_manual(values=c("#E69F00", "#4169E1"))+
        xlab(paste("PCo1 (",x_variance.bugs,"%)",sep = "")) + ylab(paste("PCo2 (",y_variance.bugs,"%)",sep = ""))
ggsave("./mgx/results/pcoa/pcoa_path.jpg", plot=pcoa.plot.path, width =7, height =6, dpi=300)





#############################################################################################################################################
########
# EC #
########
#### Note: I need to make sure that samples in rows
bugs_pcoa <- ec
#samples in rows
bugs_pcoa <- capscale(bugs_pcoa ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )

# append metadata to scores
ord.scores.meta <- merge(ord.scores, meta_pcoa, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.scores.meta) <- ord.scores.meta[ ,1]
ord.scores.meta[ ,1] <- NULL

# number of samples
n_samples <- nrow(ord.scores.meta)


# PCoA colored by disease status
pcoa.plot.ec <-
  ggplot( ord.scores.meta, aes(MDS1, MDS2) ) +
  geom_point(size = 5, alpha = 1, aes( fill =  Status),color='black',shape=21) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 16, face = 'bold'),
        legend.text = element_text(size=12), legend.title = element_blank(),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'horizontal',
        plot.title = element_text( size = 20, face = 'bold')
        )+
        scale_fill_manual(values=c("#E69F00", "#4169E1"))+
        xlab(paste("PCo1 (",x_variance.bugs,"%)",sep = "")) + ylab(paste("PCo2 (",y_variance.bugs,"%)",sep = ""))
ggsave("./mgx/results/pcoa/pcoa_ec.jpg", plot=pcoa.plot.ec, width =7, height =6, dpi=300)


#############################################################################################################################################
###########
# mbx #
###########
#### Note: I need to make sure that samples in rows
# read in mbx data
mbx <- read.csv(file="./mbx/data/mac/abundances.csv")
mbx <- mbx[,-1] #56937 232

mbx_t <- as.data.frame(t(mbx))
bugs_pcoa <- mbx_t
#samples in rows
#change NA to 0
bugs_pcoa <- dplyr::mutate_all(bugs_pcoa, ~replace(., is.na(.), 0))

# try sum normalization
bugs_pcoa_norm <- sweep(bugs_pcoa, 1, rowSums(bugs_pcoa), `/`)


###################### original values
bugs_pcoa <- capscale(bugs_pcoa ~ 1, distance = 'bray')
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )


# read in metadata from mbx
meta_mbx <- read.csv(file="./mbx/data/mac/metadata_mbx.csv")
rownames(meta_mbx) <- meta_mbx$sample
meta_mbx_pcoa <- as.data.frame(meta_mbx$caco)
colnames(meta_mbx_pcoa) <- c('Status')
rownames(meta_mbx_pcoa) <- rownames(meta_mbx)


# append metadata to scores
ord.scores.meta <- merge(ord.scores, meta_mbx_pcoa, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.scores.meta) <- ord.scores.meta[ ,1]
ord.scores.meta[ ,1] <- NULL

# number of samples
n_samples <- nrow(ord.scores.meta)


# PCoA colored by disease status
pcoa.plot.mbx <-
  ggplot( ord.scores.meta, aes(MDS1, MDS2) ) +
  geom_point(size = 5, alpha = 1, aes( fill =  Status),color='black',shape=21) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 16, face = 'bold'),
        legend.text = element_text(size=12), legend.title = element_blank(),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'horizontal',
        plot.title = element_text( size = 20, face = 'bold')
        )+
        scale_fill_manual(values=c("#E69F00", "#4169E1"))+
        xlab(paste("PCo1 (",x_variance.bugs,"%)",sep = "")) + ylab(paste("PCo2 (",y_variance.bugs,"%)",sep = ""))
ggsave("./mgx/results/pcoa/pcoa_mbx.jpg", plot=pcoa.plot.mbx, width =7, height =6, dpi=300)


###################### sum normalized
bugs_pcoa <- capscale(bugs_pcoa_norm ~ 1, distance = 'bray')
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )

# append metadata to scores
ord.scores.meta <- merge(ord.scores, meta_mbx_pcoa, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.scores.meta) <- ord.scores.meta[ ,1]
ord.scores.meta[ ,1] <- NULL

# number of samples
n_samples <- nrow(ord.scores.meta)


# PCoA colored by disease status
pcoa.plot.mbx <-
  ggplot( ord.scores.meta, aes(MDS1, MDS2) ) +
  geom_point(size = 5, alpha = 1, aes( fill =  Status) ,color='black',shape=21) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 16, face = 'bold'),
        legend.text = element_text(size=12), legend.title = element_blank(),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'horizontal',
        plot.title = element_text( size = 20, face = 'bold')
        )+
        scale_fill_manual(values=c("#E69F00", "#4169E1"))+
        xlab(paste("PCo1 (",x_variance.bugs,"%)",sep = "")) + ylab(paste("PCo2 (",y_variance.bugs,"%)",sep = ""))
ggsave("./mgx/results/pcoa/pcoa_mbx_norm.jpg", plot=pcoa.plot.mbx, width =7, height =6, dpi=300)



##############################################################################################################################################
# heatmap of 10 most abundant bugs accoring to PCoA
#############################################################################################################################################
#rank species data by mean abundance
mns <- colMeans(species_all, na.rm=TRUE)
order(-mns)
species_all_sort <- species_all[,order(-mns)]
colMeans(species_all_sort)


bugs_mds <- merge(x = species_all_sort, y = ord.scores, by = 'row.names', all = TRUE)
rownames(bugs_mds) <-bugs_mds[ ,1]
bugs_mds[ ,1] <- NULL

bugs_mds <- merge(x = bugs_mds, y = meta_pcoa, by = 'row.names', all = TRUE)
rownames(bugs_mds) <-bugs_mds[ ,1]
bugs_mds[ ,1] <- NULL


head(bugs_mds)
#sort by MDS1
bugs_mds <- bugs_mds[order(bugs_mds$Status,bugs_mds$MDS1),]

t_10species_sortmds1_heatmap <- as.data.frame(t(bugs_mds [, 1:20]))
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(seriation)
library(RColorBrewer)

#only MDS1 for annotation
caco_heatmap <- as.data.frame(bugs_mds$Status)
colnames(caco_heatmap) <- c('Status')
rownames(caco_heatmap) <- rownames(bugs_mds)


mds1_heatmap <- as.data.frame(bugs_mds$MDS1)
colnames(mds1_heatmap) <- c('MDS1')
rownames(mds1_heatmap) <- rownames(bugs_mds)


mds1_heatmap <- bugs_mds[,colnames(bugs_mds) %in% c("Status","MDS1")]


summary(mds1_heatmap$MDS1)

Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
-1.3432 -0.5358 -0.1315  0.0000  0.4227  1.8527



ht = HeatmapAnnotation (df = mds1_heatmap,
                        col = list (Status = c("Divert"="red","Control"="yellow"),
                                    MDS1 = colorRamp2(c(-2, -1, 0, 1, 2), brewer.pal(n = 5, name = "Greens"))
                                  ),
                        annotation_name_gp= gpar(fontsize = 30),
                        annotation_legend_param = list(title_gp = gpar(fontsize = 30),
                                               labels_gp = gpar(fontsize = 20))
                        )


ht1_10species = Heatmap(t_10species_sortmds1_heatmap, name = "Relative Abundance",top_annotation=ht,
                        col = colorRamp2(c(0, 0.1, 0.2, 0.3, 0.4), rev(brewer.pal(n=5, name= "RdBu"))),
                        cluster_rows = FALSE, cluster_columns = FALSE, width = unit(80, "cm"),
                        heatmap_legend_param = list(color_bar = "continuous", at = c(0, 0.1, 0.2, 0.3, 0.4), legend_direction = "horizontal",title_position="topcenter",
                        title_gp = gpar(fontsize = 30),labels_gp = gpar(fontsize = 20)),
                        row_names_gp = gpar(fontface = "italic",fontsize=50),
                        column_names_gp  = gpar(fontsize = 8)
                      )

png("./mgx/results/pcoa/heatmap_10species_byMDS1.png", width = 1300, height = 600, units='mm', res = 300)
draw(ht1_10species, heatmap_legend_side = "bottom",annotation_legend_side = "left",
     row_title = "Relative Abundance by Participant",
     row_title_gp = gpar(fontsize = 26, fontface = "bold"),
     column_title = "PCo1", column_title_side = "top",
     column_title_gp = gpar(fontsize = 26, fontface = "bold")
   )
dev.off()
