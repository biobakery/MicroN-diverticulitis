################################################################################
## Purpose: Scatter plots of correlations between species and metabolites in divert and control ###
## Last update date: 10/11/2023 (formatting)
################################################################################

# salloc -p test --mem 64G -t 0-08:00
# module load R/4.2.2-fasrc01
# #install new packages
# export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
# R --quiet
# 
# setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

library(ggplot2)
library(ggpubr)

################################################################################
###                           species (cohort id)                           ####
################################################################################

taxa <- read.table('./mbx/data/halla/raw/all/taxa.tsv', header=TRUE, sep="\t")
rownames(taxa) <- taxa$Feature
taxa_t <- as.data.frame(t(taxa))
taxa_t <- taxa_t[-1,]

# change to numeric
for(i in 1:ncol(taxa_t))
{
  taxa_t[ , i] <- as.numeric(as.character(taxa_t[, i]))
}

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
taxa_t_ast <- apply(taxa_t,2,AST)


################################################################################
###                          metabolites (cohort id)                        ####
################################################################################

mbx <- read.table('./mbx/data/halla/raw/all/mbx.tsv', header=TRUE, sep="\t")
mbx_t <- as.data.frame(t(mbx))
colnames(mbx_t) <- mbx_t[1,]
mbx_t <- mbx_t[-1,]

# change to numeric
for(i in 1:ncol(mbx_t))
{
  mbx_t[ , i] <- as.numeric(as.character(mbx_t[, i]))
}

# Log transformation
LOG <- function(x) {
    y <- replace(x, x == 0, min(x[x>0]) / 2)
    return(log2(y))
}

mbx_t_log <- apply(mbx_t, 2, LOG)


################################################################################
###                         merge species with metabolites                  ####
################################################################################

mbx_taxa <- merge(taxa_t_ast,mbx_t_log,by='row.names')
rownames(mbx_taxa) <- mbx_taxa[,1]
mbx_taxa <- mbx_taxa[,-1]

rownames(mbx_taxa) <- gsub("X", "", rownames(mbx_taxa))


################################################################################
###                           metadata                       ####
################################################################################
meta <- read.table('./mbx/data/mbx_maaslin_metadata.tsv', header=TRUE, sep='\t')
meta <- meta[,colnames(meta) %in% c("sample","caco")]
colnames(meta)[1] <- 'SAMPLEID'

# change to cohort id
# idkey for mbx_id and cohort id
id_sampleid_mbx <- read.csv(file="./idkey/id_sampleid_mbx.csv",row.names=1)
meta<- merge(meta,id_sampleid_mbx,by='SAMPLEID')
rownames(meta) <- meta$id


################################################################################
###                           merge                        ####
################################################################################
all <- merge(mbx_taxa,meta,by='row.names')
rownames(all) <- all[,1]
all <- all[,-1]

################################################################################
###                      scatter plot                        ####
################################################################################
ggp_by_group <- function(x,y) {
  plot <- ggplot(all, aes_string(x=x, y=y))+
    geom_point(alpha =2, shape = 21, size =10, stroke = 1, na.rm=TRUE,color='black',aes(fill=caco)) +
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE,size=5,aes(color=caco))+
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
    scale_fill_manual(values = c("#E69F00", "#4169E1"), labels=c("Control", "Divert"))+
    scale_color_manual(values = c("#E69F00", "#4169E1"), labels=c("Control", "Divert"))
    #stat_cor(method = "spearman",size=12)
  filepath <- paste('./mbx/results/mgx_mbx_scatter/by_caco/', x,'_',y,'.png', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 15, height = 16, dpi = 300)
  print(plot)
}# end

# ggp_by_group('`Deoxycholic acid`','Clostridium_clostridioforme')
# ggp_by_group('`Deoxycholic acid`','Ruminococcus_gnavus')
# ggp_by_group('`Sebacic acid`','Ruminococcus_gnavus')
# ggp_by_group('DPA','Ruminococcus_gnavus')
# ggp_by_group('`Undecanedioic acid`','Clostridium_clostridioforme')
# ggp_by_group('`Undecanedioic acid`','Firmicutes_bacterium_CAG_129')

ggp_by_group('`alpha-Muricholic acid`','Bilophila_wadsworthia')
ggp_by_group('`Cholic acid`','Bilophila_wadsworthia')

ggp_by_group('`Phe-Tyr`','Bacteroides_fragilis')
ggp_by_group('`Thr-Phe`','Bacteroides_thetaiotaomicron')
ggp_by_group('`Ser-Ile`','Ruminococcaceae_bacterium_D16')

ggp_by_group('`CE 16:0`','Firmicutes_bacterium_CAG_56')
ggp_by_group('`TG 50:2`','Firmicutes_bacterium_CAG_56')

# ggp_by_group('`Cer 18:1;O2/17:0`','Bacteroides_thetaiotaomicron')
# ggp_by_group('`Cer 18:1;O2/16:0`','Bacteroides_thetaiotaomicron')
# ggp_by_group('`Cer 18:1;O2/17:0`','Bacteroides_ovatus')
# ggp_by_group('`Cer 18:1;O2/16:0`','Bacteroides_ovatus')
