################################################################################
# Purpose: prepare data for HALLA (raw species and metabolites)                      ####
# Last update date: 9/25/2023
################################################################################


setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")


################################################################################
###                           species                              ####
################################################################################
# read in species
species_all <- read.table('./mgx/data/micron_species_filtered.tsv', header=TRUE, sep='\t')
# 189 species 235 samples
# UID

################################################################################
###                           annotated metabolites               ####
################################################################################
mbx_anno_df <- read.csv(file="./mbx/data/mbx_anno_df.csv",row.names=1) # 552 metabolites 232 samples
mbx_anno_df <- mbx_anno_df[,7:dim(mbx_anno_df)[2]]

# 552 metabolites 232 samples
# Metablite

t_mbx_anno_df <- as.data.frame(t(mbx_anno_df))
t_mbx_anno_df_1 <- t_mbx_anno_df[-1,]
# 232 552


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

mbx <- as.data.frame(t(mbx_anno_df))
mbx$SAMPLEID <- rownames(mbx)

name <- as.data.frame(mbx[1,])
mbx_noname <- mbx[-1,]
mbx_noname <- merge(id_sampleid_mbx,mbx_noname,by="SAMPLEID")
mbx_noname <- mbx_noname[,!colnames(mbx_noname) %in% "SAMPLEID"]

inta <- intersect(species_all$id,mbx_noname$id)
species_for_halla <- species_all[species_all$id %in% inta,]


rownames(species_for_halla) <- species_for_halla$id
species_for_halla <- species_for_halla[,!colnames(species_for_halla) %in% "id"]
species_for_halla <- species_for_halla[order(rownames(species_for_halla)),]

taxa <- as.data.frame(t(species_for_halla))
taxa$Feature <- rownames(taxa)
library(dplyr)
taxa <- taxa %>%
  select(Feature, everything())

write.table(taxa, file = "./mbx/data/halla/raw/all/taxa.tsv", row.names=FALSE, sep="\t")


mbx_for_halla <- mbx_noname[mbx_noname$id %in% inta,] #226 samples
rownames(mbx_for_halla) <- mbx_for_halla$id
mbx_for_halla <- mbx_for_halla[order(rownames(mbx_for_halla)),]
mbx_for_halla <- mbx_for_halla[,!colnames(mbx_for_halla) %in% "id"]

name <- name[,!colnames(name) %in% "SAMPLEID"]

mbx_for_halla <- rbind(name,mbx_for_halla)
mbx_for_halla <- as.data.frame(t(mbx_for_halla))

names(mbx_for_halla)[names(mbx_for_halla) == "Metablite"] <- "Feature"
write.table(mbx_for_halla, file = "./mbx/data/halla/raw/all/mbx.tsv", row.names=FALSE, sep="\t")

################################################################################
###                           divert                         ####
################################################################################
taxa_divert <- taxa[,colnames(taxa) %in% divert_id]
taxa_divert$Feature <- rownames(taxa_divert)
taxa_divert <- taxa_divert %>%
  select(Feature, everything())

mbx_for_halla_divert <-mbx_for_halla[,colnames(mbx_for_halla) %in% divert_id]
test <- as.data.frame(mbx_for_halla[,colnames(mbx_for_halla) %in% c("Feature")])
colnames(test) <- "Feature"
mbx_for_halla_divert <- cbind(test,mbx_for_halla_divert)


################################################################################
###                          control                         ####
################################################################################
taxa_control <- taxa[,colnames(taxa) %in% control_id]
t_taxa_control <- as.data.frame(t(taxa_control))
t_taxa_control$id <- rownames(t_taxa_control)
t_taxa_control <- merge(t_taxa_control,caco_data,by="id")
rownames(t_taxa_control) <- t_taxa_control$matchid
t_taxa_control <- t_taxa_control[,!colnames(t_taxa_control) %in% c("id","caco","matchid")]
t_taxa_control <- t_taxa_control[order(rownames(t_taxa_control)),]

taxa_control_2 <- as.data.frame(t(t_taxa_control))

taxa_control_2$Feature <- rownames(taxa_control_2)
taxa_control_2 <- taxa_control_2 %>%
  select(Feature, everything())

#######################################
mbx_for_halla_control <-mbx_for_halla[,colnames(mbx_for_halla) %in% control_id]
t_mbx_control <- as.data.frame(t(mbx_for_halla_control))
t_mbx_control$id <- rownames(t_mbx_control)
t_mbx_control <- merge(t_mbx_control,caco_data,by="id")
rownames(t_mbx_control) <- t_mbx_control$matchid
t_mbx_control <- t_mbx_control[,!colnames(t_mbx_control) %in% c("id","caco","matchid")]
t_mbx_control <- t_mbx_control[order(rownames(t_mbx_control)),]

mbx_control_2 <- as.data.frame(t(t_mbx_control))
test <- as.data.frame(mbx_for_halla[,colnames(mbx_for_halla) %in% c("Feature")])
colnames(test) <- "Feature"
mbx_for_halla_control <- cbind(test,mbx_control_2)


################################################################################
###                           matched data                         ####
################################################################################
inta <- intersect(colnames(taxa_divert),colnames(taxa_control_2))
taxa_divert_match <- taxa_divert[,colnames(taxa_divert) %in% inta]
taxa_control_match <- taxa_control_2[,colnames(taxa_control_2) %in% inta]

table(colnames(taxa_divert_match) == colnames(taxa_control_match))
table(rownames(taxa_divert_match) == rownames(taxa_control_match))
write.table(taxa_divert_match, file = "./mbx/data/halla/raw/divert/taxa_divert.tsv", row.names=FALSE, sep="\t")
write.table(taxa_control_match, file = "./mbx/data/halla/raw/control/taxa_control.tsv", row.names=FALSE, sep="\t")

inta2 <- intersect(colnames(mbx_for_halla_divert),colnames(mbx_for_halla_control))
mbx_divert_match <- mbx_for_halla_divert[,colnames(mbx_for_halla_divert) %in% inta2]
mbx_control_match <- mbx_for_halla_control[,colnames(mbx_for_halla_control) %in% inta2]
table(colnames(mbx_divert_match) == colnames(mbx_control_match))
table(rownames(mbx_divert_match) == rownames(mbx_control_match))
write.table(mbx_divert_match, file = "./mbx/data/halla/raw/divert/mbx_divert.tsv", row.names=FALSE, sep="\t")
write.table(mbx_control_match, file = "./mbx/data/halla/raw/control/mbx_control.tsv", row.names=FALSE, sep="\t")




# take difference
taxa_diff <- taxa_divert_match[,2:dim(taxa_divert_match)[2]]-taxa_control_match[,2:dim(taxa_control_match)[2]]
taxa_diff <- cbind(taxa_divert_match[,1],taxa_diff)
names(taxa_diff)[1] <- 'Feature'
write.table(taxa_diff, file = "./mbx/data/halla/raw/diff/taxa_diff.tsv", row.names=FALSE, sep="\t")

i <- c(2:dim(mbx_divert_match)[2])                                  # Specify columns you want to change
mbx_divert_match[ , i] <- apply(mbx_divert_match[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
mbx_control_match[ , i] <- apply(mbx_control_match[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))

mbx_diff <- mbx_divert_match[,2:dim(mbx_divert_match)[2]]-mbx_control_match[,2:dim(mbx_control_match)[2]]
mbx_diff <- cbind(mbx_divert_match[,1],mbx_diff)
names(mbx_diff)[1] <- 'Feature'
write.table(mbx_diff, file = "./mbx/data/halla/raw/diff/mbx_diff.tsv", row.names=FALSE, sep="\t")


######use halla on hutlab11 server
wma@hutlab11.rc.fas.harvard.edu
cd /net/fs-huttenhower01/srv/export/huttenhower_lab/share/tools/hutlab/src/modules
source /net/fs-huttenhower01/srv/export/huttenhower_lab/share/tools/hutlab/src/hutlabrc.sh
hutlab load centos7/python3/biobakery_workflows/3.0.0-beta-devel-dependsUpdate

cd /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mbx/data/halla/raw
halla -x ./divert/taxa_divert.tsv   -y ./divert/mbx_divert.tsv   -m spearman -o ./divert/synthetic_output
halla -x ./control/taxa_control.tsv -y ./control/mbx_control.tsv -m spearman -o ./control/synthetic_output
halla -x ./all/taxa.tsv -y ./all/mbx.tsv -m spearman -o ./all/synthetic_output
