################################################################################
# Purpose: prepare data of ECs and prioritized ceramide metabolites for HALLA                        
# Last update date: 9/25/2023
# Note: I tested using orginal mbx data (non-log transformed) and had exactly the same halla plot
################################################################################

#setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")


################################################################################
###                           species                              ####
################################################################################
# read in ECs
ec_all <- read.table('./mgx/data/micron_ec_filtered.tsv', header=TRUE, sep='\t')
# 1745 ecs 235 samples
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
mbx_subset <- as.data.frame(t(mbx[,c(269,14025,27879,13976,28549,1376,27175,27704,13981,13978,26971)]))
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
ec_all <- merge(ec_all,id_UID,by="UID")
ec_all <- ec_all[,!colnames(ec_all) %in% c("UID","SAMPLEID")]


mbx <- as.data.frame(t(mbx_subset))
mbx$SAMPLEID <- rownames(mbx)


mbx <- merge(id_sampleid_mbx,mbx,by="SAMPLEID")
mbx <- mbx[,!colnames(mbx) %in% "SAMPLEID"]

inta <- intersect(ec_all$id,mbx$id)
ec_for_halla <- ec_all[ec_all$id %in% inta,]


rownames(ec_for_halla) <- ec_for_halla$id
ec_for_halla <- ec_for_halla[,!colnames(ec_for_halla) %in% "id"]
ec_for_halla <- ec_for_halla[order(rownames(ec_for_halla)),]

# remove X in column names
colnames(ec_for_halla) <- gsub("X","", colnames(ec_for_halla))

# change ec labels to full names
# bring in the dna_ec key
dnaec_key <-  read.table( '/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/microbiome/mlvs_fiber_crp/data/channing_output/outcome/ecdna/dna_lookup.txt', sep = "\t", header = T, row.names = 1)
# compatible names
dnaec_key$correct <- make.names(dnaec_key$ec_name)
dnaec_key <- as.data.frame(t (dnaec_key))
dnaec_key[] <-lapply(dnaec_key, as.character)

# duplicate the one with formatted name as well
dnaec_key2 <- dnaec_key

#finally keep only the R compatible name
dnaec_key <- dnaec_key[rownames(dnaec_key) %in% 'correct', ]

# dnaec_key has more columns, so just keep cols in common
common_col <- intersect(colnames(dnaec_key),colnames(ec_for_halla))

#remove_col <- ec_for_halla[,!colnames(ec_for_halla) %in% common_col]
# Wenjie: may need to update the dnaec_key file

ec_for_halla <- ec_for_halla[,colnames(ec_for_halla) %in% common_col]
dnaec_key <- dnaec_key[,colnames(dnaec_key) %in% common_col]

ec_for_halla_name <- rbind(dnaec_key,ec_for_halla)
colnames(ec_for_halla_name) <- ec_for_halla_name[1,]
ec_for_halla_name <- ec_for_halla_name[-1,] #226 1463

taxa <- as.data.frame(t(ec_for_halla_name))
taxa$Feature <- rownames(taxa)
library(dplyr)
taxa <- taxa %>%
  select(Feature, everything())

write.table(taxa, file = "./mbx/data/halla/raw/ec_ceramide/ec.tsv", row.names=FALSE, sep="\t")


mbx_for_halla <- mbx[mbx$id %in% inta,] #226 samples
rownames(mbx_for_halla) <- mbx_for_halla$id
mbx_for_halla <- mbx_for_halla[order(rownames(mbx_for_halla)),]
mbx_for_halla <- mbx_for_halla[,!colnames(mbx_for_halla) %in% "id"]


mbx_for_halla <- as.data.frame(t(mbx_for_halla))

mbx_for_halla$Feature <- rownames(mbx_for_halla)
library(dplyr)
mbx_for_halla<- mbx_for_halla %>%
  select(Feature, everything())

write.table(mbx_for_halla, file = "./mbx/data/halla/raw/ec_ceramide/mbx_cer.tsv", row.names=FALSE, sep="\t")




######use halla on hutlab11 server
wma@hutlab11.rc.fas.harvard.edu
cd /net/fs-huttenhower01/srv/export/huttenhower_lab/share/tools/hutlab/src/modules
source /net/fs-huttenhower01/srv/export/huttenhower_lab/share/tools/hutlab/src/hutlabrc.sh
hutlab load centos7/python3/biobakery_workflows/3.0.0-beta-devel-dependsUpdate

cd /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mbx/data/halla/raw
halla -x ./ec_ceramide/ec.tsv -y ./ec_ceramide/mbx_cer.tsv -m spearman -o ./ec_ceramide/synthetic_output
halla -x ./ec_ceramide/mbx_cer.tsv -y ./ec_ceramide/ec100.tsv -m spearman -o ./ec_ceramide/synthetic_output_ec100


#hallagram -i ./all/synthetic_output/ --cbar_label 'Pairwise Spearman' --x_dataset_label 'Species' --y_dataset_label 'Metabolites' 
