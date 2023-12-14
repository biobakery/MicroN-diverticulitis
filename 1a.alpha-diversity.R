
#########################################################################################################
### Purpose: To compare alpha diversity indexes between diverticulitis and controls
### Last update date: 10/13/2023 
#########################################################################################################

# module load R/4.2.2-fasrc01
# #install new packages
# export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
# R --quiet
# 
# 
# 
# 
# rm(list=ls())
# getwd()
#setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")


library(ggplot2)
library(vegan)
library(ggpubr)


# read in taxonomy data
species_all <- read.table('./mgx/data/micron_species_filtered.tsv', header=TRUE, sep='\t')
rownames(species_all) <- species_all$UID
species_all <- species_all[,-1]

# read in metadata
meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID


# metadata
meta_pcoa <- as.data.frame(meta$caco)
colnames(meta_pcoa) <- c('Status')
rownames(meta_pcoa) <- rownames(meta)


shannon <- as.data.frame(diversity(species_all,"shannon"))
colnames(shannon) <- "shannon"

richness <- as.data.frame(specnumber(species_all))
colnames(richness) <- "richness"

simp <- as.data.frame(diversity(species_all, "simpson"))
colnames(simp) <- "simp"

alpha <- cbind(shannon,richness,simp)

alpha$evenness <- alpha$shannon/log(alpha$richness)

alpha.meta <- merge(alpha,meta_pcoa,by="row.names")
rownames(alpha.meta) <- alpha.meta[,1]
alpha.meta[,1] <- NULL

alpha.meta$mean=mean(alpha.meta$richness)
alpha.meta$se=sd(alpha.meta$richness) / sqrt(length(alpha.meta$richness))
alpha.meta$LCL=alpha.meta$mean-1.96*alpha.meta$se
alpha.meta$UCL=alpha.meta$mean+1.96*alpha.meta$se



# Split the values by group
shannon.spl <- split(alpha.meta$shannon, alpha.meta$Status)
# Kolmogorov-Smironv test
shannon.pv <- ks.test(shannon.spl$Divert, shannon.spl$Control)$p.value
print(shannon.pv) #0.3982005

simp.spl <- split(alpha.meta$simp, alpha.meta$Status)
simp.pv <- ks.test(simp.spl$Divert, simp.spl$Control)$p.value
print(simp.pv) #0.3255683

inv.spl <- split(alpha.meta$inv, alpha.meta$Status)
inv.pv <- ks.test(inv.spl$Divert, inv.spl$Control)$p.value
print(inv.pv) #0.3255683

rich.spl <- split(alpha.meta$richness, alpha.meta$Status)
rich.pv <- ks.test(rich.spl$Divert, rich.spl$Control)$p.value
print(rich.pv) #0.03

even.spl <- split(alpha.meta$evenness, alpha.meta$Status)
even.pv <- ks.test(even.spl$Divert, even.spl$Control)$p.value
print(even.pv) #0.49


res <- wilcox.test(inv ~ Status, data = alpha.meta,
                   exact = FALSE)
res

###############################################################
######################  richness       #######################
##############################################################
p1 <- ggplot(alpha.meta,aes(x = Status, y = richness,fill=Status))+
  #geom_crossbar(stat = "summary",fun.data = "mean_se", fun.args = list(mult = 1.96), width = 0.2,alpha=1) +
  geom_jitter(size=6, position=position_jitter(0.2),aes(color=Status),alpha=0.75)+
  geom_boxplot(alpha=0.75)+
  scale_fill_manual(values=c("#E69F00", "#4169E1"))+
  scale_color_manual(values=c("#E69F00", "#4169E1"))+
theme_classic() +
theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
      plot.title = element_text(size = 30),
      legend.title = element_text(size=30),
      legend.text = element_text(size=30),
      axis.title=element_text(size=30,face="bold"),
      axis.text=element_text(size=30),
      axis.ticks=element_blank(),
      axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
      ) +
xlab("Status") +
ylab("Species richness")+
stat_compare_means(
  aes(label = paste0("p = ", ..p.format..)), size=12)
print(p1)

ggsave('./mgx/results/vis/richness.jpg', plot=p1, width = 10, height = 10, dpi = 300)


###############################################################
######################  shannon      #######################
##############################################################
p2 <- ggplot(alpha.meta,aes(x = Status, y = shannon,fill=Status))+
  geom_jitter(size=6, position=position_jitter(0.2),aes(color=Status),alpha=0.75)+
  geom_boxplot(alpha=0.75)+
  scale_fill_manual(values=c("#E69F00", "#4169E1"))+
  scale_color_manual(values=c("#E69F00", "#4169E1"))+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 30),
        legend.title = element_text(size=30),
        legend.text = element_text(size=30),
        axis.title=element_text(size=30,face="bold"),
        axis.text=element_text(size=30),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
  ) +
  xlab("Status") +
  ylab("Shannon index")+
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)), size=12)
print(p2)

ggsave('./mgx/results/vis/shannon.jpg', plot=p2, width = 10, height = 10, dpi = 300)


###############################################################
######################  simpson      #######################
##############################################################
p3 <- ggplot(alpha.meta,aes(x = Status, y = simp,fill=Status))+
  geom_jitter(size=3, position=position_jitter(0.2),aes(color=Status),alpha=0.75)+
  geom_boxplot(alpha=0.75)+
  scale_fill_manual(values=c("#E69F00", "#4169E1"))+
  scale_color_manual(values=c("#E69F00", "#4169E1"))+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 30),
        legend.title = element_text(size=30),
        legend.text = element_text(size=30),
        axis.title=element_text(size=30,face="bold"),
        axis.text=element_text(size=30),
        axis.ticks=element_blank(),
        axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
  ) +
  xlab("Status") +
  ylab("Simpson index")+
  stat_compare_means(
    aes(label = paste0("p = ", ..p.format..)), size=12)
print(p3)

ggsave('./mgx/results/vis/simp.jpg', plot=p3, width = 10, height = 10, dpi = 300)

########################################################
library(ggpubr)
p <- ggboxplot(alpha.meta, x = "Status", y = "richness",
               color = "Status", palette = "jco",
               add = "jitter")
#  Add p-value
p + stat_compare_means()
