################################################################################
# Purpose: microbe-metabolite correlations in all samples                        
# Last update date: 10/3/2023
################################################################################



#setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

library(ggplot2)
library(ggpubr)

################################################################################
###             association file generated using halla                      ####
################################################################################

r <- read.table('./mbx/data/halla/raw/all/synthetic_output/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")

cer <- r[r$Y_features == "Cer 18:1;O2/17:0",]
rank <- order(cer$association)
cer_sort <- cer[rank,]

parasutterella <- r[r$X_features == "Parasutterella_excrementihominis",]
rank <- order(parasutterella$association)
parasutterella_sort <- parasutterella[rank,]

dimethy <- r[r$Y_features == "1,7-Dimethyluric acid",]
rank <- order(dimethy$association)
dimethy_sort <- dimethy[rank,]

rgnavus <- r[r$X_features == "Ruminococcus_gnavus",]
rank <- order(rgnavus$association)
rgnavus_sort <- rgnavus[rank,]

dimethy <- r[r$Y_features == "1,7-Dimethyluric acid",]
rank <- order(dimethy$association)
dimethy_sort <- dimethy[rank,]


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
  filepath <- paste('./mbx/results/mgx_mbx_scatter/all/', x,'_',y,'.jpg', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 16, height = 15, dpi = 300)
  print(plot)
}# end

ggp('`Cer 18:1;O2/17:0`','Subdoligranulum_sp')
ggp('`1,7-Dimethyluric acid`','Blautia_sp_CAG_257')
ggp('`1,7-Dimethyluric acid`','Ruminococcus_gnavus')
ggp('`1,7-Dimethyluric acid`','Flavonifractor_plautii')
ggp('`1,7-Dimethyluric acid`','Parasutterella_excrementihominis')
ggp('`3-Methyladipic acid or Pimelic acid`','Parasutterella_excrementihominis')

ggp('`Deoxycholic acid`','Blautia_sp_CAG_257')
ggp('Tyramine','Ruminococcus_gnavus')
ggp('`alpha-Muricholic acid`','Parasutterella_excrementihominis')
ggp('`1,7-Dimethyluric acid`','Oscillibacter_sp_57_20')
ggp('Histamine','Oscillibacter_sp_57_20')
ggp('`N-Acetylhistamine`','Oscillibacter_sp_57_20')

ggp('`Cer 18:1;O2/17:0`','Bacteroides_thetaiotaomicron')


#rcorr(all$`Cer 18:1;O2/16:0`,all$Bacteroides_thetaiotaomicron)
#rcorr(all$`Cer 18:1;O2/16:0`,all$Bacteroides_ovatus)
