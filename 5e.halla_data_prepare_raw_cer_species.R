################################################################################
# Purpose: prepare data of species and prioritized ceramide metabolites for HALLA                        
# Last update date: 9/25/2023
################################################################################

#setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")


################################################################################
###                           species                              ####
################################################################################
# read in species
species_all <- read.table('./mgx/data/micron_species_filtered.tsv', header=TRUE, sep='\t')
# 189 species 235 samples
# UID

################################################################################
###                           all metabolites               ####
################################################################################
abundances <- read.csv(file="./mbx/data/mac/abundances.csv") # 56937 metabolites 232 samples
abundances <- abundances[,-1] #56937   232
abundances[is.na(abundances)] <- 0

mbx <- as.data.frame(t(abundances)) # samples in rows

# Log transformation
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log10(y))
}

mbx_log <- as.data.frame(apply(mbx, 2, LOG)) # 232 56937

# keep metabolites that were prioritized 
mbx_subset <- as.data.frame(t(mbx_log[,c(269,14025,27879,13976,28549,1376,27175,27704,13981,13978,26971)]))
row.names(mbx_subset)[1] <- "Cer 18:1;O2/17:0"



################################################################################
###                           idkey file                         ####
################################################################################
# idkey for mgx_id and cohort_id
id_UID <- read.csv(file="./idkey/id_UID.csv",row.names=1)
dup <- id_UID[duplicated(id_UID$id),]

# idkey for mbx_id and cohort id
id_sampleid_mbx <- read.csv(file="./idkey/id_sampleid_mbx.csv",row.names=1)

# combined id mapping file
id_mgx_mbx <- merge(id_UID, id_sampleid_mbx, by="id")

################################################################################
###                           case control matching file         ####
################################################################################
caco_matrix <- read.csv(file="./idkey/caco_matrix.csv")
caco_data <- caco_matrix[caco_matrix$caco==2,]
caco_data$id <- substr(caco_data$id,1,6) #keep first 6 digits
caco_data$matchid <- substr(caco_data$matchid,1,6) #keep first 6 digits


################################################################################
###                           metadata                       ####
################################################################################
meta <- read.table('./mbx/data/mbx_maaslin_metadata.tsv', header=TRUE, sep='\t')
meta <- meta[,colnames(meta) %in% c("sample","caco")]
colnames(meta)[1] <- 'SAMPLEID'

# change to cohort id
meta<- merge(meta,id_sampleid_mbx,by='SAMPLEID')
rownames(meta) <- meta$id
divert_id <- meta[meta$caco=="Divert",]$id
control_id <- meta[meta$caco=="Control",]$id

################################################################################
###                           id match                         ####
################################################################################
# change mbx id to match mgx file
species_all <- merge(species_all,id_UID,by="UID")
species_all <- species_all[,!colnames(species_all) %in% c("UID","SAMPLEID")]

mbx <- as.data.frame(t(mbx_subset))
mbx$SAMPLEID <- rownames(mbx)


mbx <- merge(id_sampleid_mbx,mbx,by="SAMPLEID")
mbx <- mbx[,!colnames(mbx) %in% "SAMPLEID"]

inta <- intersect(species_all$id,mbx$id)
species_for_halla <- species_all[species_all$id %in% inta,]


rownames(species_for_halla) <- species_for_halla$id
species_for_halla <- species_for_halla[,!colnames(species_for_halla) %in% "id"]
species_for_halla <- species_for_halla[order(rownames(species_for_halla)),]

taxa <- as.data.frame(t(species_for_halla))
taxa$Feature <- rownames(taxa)
library(dplyr)
taxa <- taxa %>%
  select(Feature, everything())

write.table(taxa, file = "./mbx/data/halla/raw/species_ceramide/taxa.tsv", row.names=FALSE, sep="\t")


mbx_for_halla <- mbx[mbx$id %in% inta,] #226 samples
rownames(mbx_for_halla) <- mbx_for_halla$id
mbx_for_halla <- mbx_for_halla[order(rownames(mbx_for_halla)),]
mbx_for_halla <- mbx_for_halla[,!colnames(mbx_for_halla) %in% "id"]


mbx_for_halla <- as.data.frame(t(mbx_for_halla))

mbx_for_halla$Feature <- rownames(mbx_for_halla)
library(dplyr)
mbx_for_halla<- mbx_for_halla %>%
  select(Feature, everything())

write.table(mbx_for_halla, file = "./mbx/data/halla/raw/species_ceramide/mbx_cer.tsv", row.names=FALSE, sep="\t")




######use halla on hutlab11 server
wma@hutlab11.rc.fas.harvard.edu
cd /net/fs-huttenhower01/srv/export/huttenhower_lab/share/tools/hutlab/src/modules
source /net/fs-huttenhower01/srv/export/huttenhower_lab/share/tools/hutlab/src/hutlabrc.sh
hutlab load centos7/python3/biobakery_workflows/3.0.0-beta-devel-dependsUpdate


cd /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mbx/data/halla/raw
halla -x ./species_ceramide/mbx_cer.tsv -y ./species_ceramide/taxa.tsv -m spearman -o ./species_ceramide/synthetic_output/ 
  
 
cd /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mbx/data/halla/example
halla -x ./X_dataset.txt -y ./Y_dataset.txt -m spearman -o ./synthetic_output/ 
  
cd /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mbx/data/halla/yiqing
halla -x ./predict_carbs.txt -y ./predict_taxa.txt -m spearman -o ./synthetic_output/ 
  


