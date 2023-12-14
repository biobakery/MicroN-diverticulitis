
##############################################################################################################################################
# Purpose: To prepare data for metabolomics analysis
# Last update date: 10/13/2023 (update metadata using the finalized data from Channing)
#############################################################################################################################################

# salloc -p test --mem 64G -t 0-08:00
# module load R/4.2.2-fasrc01
# #install new packages
# export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
# R --quiet
# 
# setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

library(readxl)
library(dplyr)
library(stringr)

################################################################################
#ID mapping
idkey <- read.csv(file="./idkey/lab20212.diverticulitis.manifest.mapping.csv")
idkey <- idkey[!idkey$PARTICIPANTID=="",c("PARTICIPANTID","SAMPLEID")]
idkey$id <- substr(idkey$PARTICIPANTID,1,6) #keep first 6 digits
idkey <- idkey[,colnames(idkey) %in% c("id","SAMPLEID")]
write.csv(idkey,file="./idkey/id_sampleid_mbx.csv", na="")

################################################################################

HILIC_pos <- read_excel('./mbx/data/original_data_lock/22_0720_MICRO-N_diverticulitis_HILIC-pos.xlsx',skip=6)
HILIC_neg <- read_excel('./mbx/data/original_data_lock/22_0720_MICRO-N_diverticulitis_HILIC-neg.xlsx',skip=6)
C8_pos <- read_excel('./mbx/data/original_data_lock/22_0720_MICRO-N_diverticulitis_C8-pos.xlsx',skip=6)
C8_neg <- read_excel('./mbx/data/original_data_lock/22_0916_MICRO-N_diverticulitis_C18-neg.xlsx',skip=6)

mbx_files <- c('HILIC_pos', 'HILIC_neg', 'C8_pos', 'C8_neg')


for (i in 1:length(mbx_files))
{
  df <- NULL
  df <- get (mbx_files[i])
  df <- df[apply(is.na(df), 1, sum) < dim(df)[2],] # number of NAs in the row less than number of columns
  df <- df[,apply(is.na(df), 2, sum) < dim(df)[1]]
  #df <- dplyr::mutate_all(df, ~replace(., is.na(.), 0))
  my.name <- mbx_files[i]
  assign(my.name, df)
}


  mbx_full_df <- full_join(HILIC_pos, HILIC_neg)
  mbx_full_df <- full_join(mbx_full_df, C8_pos) #29283   274
  mbx_full_df <- full_join(mbx_full_df, C8_neg) #56948   274


#remove QC samples
test <- mbx_full_df[,grep("V000", colnames(mbx_full_df))] #56948 237
mbx_full_df <- cbind(mbx_full_df[,1:7],test) #56948 244

mbx_full_df$mol_id <- paste(mbx_full_df$Compound_ID, mbx_full_df$Method, sep = '_')
mbx_full_df <- mbx_full_df[!mbx_full_df$mol_id == '_',]
mbx_full_df$mol_id <- str_replace(mbx_full_df$mol_id, '-', '_')

mbx_full_df <- dplyr::mutate_all(mbx_full_df, ~replace(., is.na(.), 0))

mbx_full_df <- mbx_full_df[!mbx_full_df$HMDB_ID =="internal standard",] #56937 245
rownames(mbx_full_df) <- mbx_full_df$mol_id
mbx_full_df <- dplyr::select(mbx_full_df, -mol_id)


write.table(mbx_full_df, file = "./mbx/data/mbx_full_df.tsv", row.names=TRUE, sep="\t")


################################################################################
# prepare input files for MACARRoN

mbx_full_df <- read.csv(file="./mbx/data/mbx_full_df.tsv", header=TRUE, sep="\t") #56937 244
mbx_mol_info <- mbx_full_df[,1:7] #56937 7
mbx_full_df <- mbx_full_df[,8:dim(mbx_full_df)[2]] #56937 237
#existing_sampleid <- colnames(mbx_full_df)
#write.csv(existing_sampleid, file="./mbx/data/existing_sampleid.csv",row.names=F)


############
# 1-abundances
############
mbx_full_df_export <- mbx_full_df
one <- as.data.frame(t(mbx_full_df_export))
one$SAMPLEID <- rownames(one)

#average duplicates
idkey <- read.csv(file="./idkey/id_sampleid_mbx.csv")
idkey <- idkey[,-1]

two <- merge(idkey,one,by="SAMPLEID")
rownames(two) <- two$SAMPLEID

dup <- two[duplicated(two$id),]

# SAMPLEID     id
# 63  V000407912 620785
# 78  V000420784 689113
# 80  V000423257 640855
# 85  V000427089 748749
# 168 V000540410 568046

# take average of the dups
df_620785 <- two[two$id=="620785",]
df_620785 <- df_620785[,!colnames(df_620785) %in% c("id","SAMPLEID")]
df_620785 <- as.data.frame(t(as.data.frame(colMeans(df_620785))))
rownames(df_620785) <- "V000407912"


df_689113 <- two[two$id=="689113",]
df_689113 <- df_689113[,!colnames(df_689113) %in% c("id","SAMPLEID")]
df_689113 <- as.data.frame(t(as.data.frame(colMeans(df_689113))))
rownames(df_689113) <- "V000420784"


df_640855 <- two[two$id=="640855",]
df_640855 <- df_640855[,!colnames(df_640855) %in% c("id","SAMPLEID")]
df_640855 <- as.data.frame(t(as.data.frame(colMeans(df_640855))))
rownames(df_640855) <- "V000423257"

df_748749 <- two[two$id=="748749",]
df_748749 <- df_748749[,!colnames(df_748749) %in% c("id","SAMPLEID")]
df_748749 <- as.data.frame(t(as.data.frame(colMeans(df_748749))))
rownames(df_748749) <- "V000427089"

df_568046 <- two[two$id=="568046",]
df_568046 <- df_568046[,!colnames(df_568046) %in% c("id","SAMPLEID")]
df_568046 <- as.data.frame(t(as.data.frame(colMeans(df_568046))))
rownames(df_568046) <- "V000540410"

three <- rbind(df_620785,df_689113,df_640855,df_748749,df_568046)

two <- two[!two$id=="620785",]
two <- two[!two$id=="689113",]
two <- two[!two$id=="640855",]
two <- two[!two$id=="748749",]
two <- two[!two$id=="568046",]

two <- two[,!colnames(two) %in% c("id","SAMPLEID")]
four <- rbind(two,three) # 232 56937
mbx_full_df_export <- as.data.frame(t(four))


abudances <- mbx_full_df_export
rownames(abudances) <- NULL
addcol <- as.data.frame(rownames(abudances))
colnames(addcol) <- "Feature"
abudances <- cbind(addcol,abudances)
abudances[abudances== 0] <- NA
write.csv(abudances, file="./mbx/data/mac/abundances.csv",row.names=F)



############
# 2-annotations
############
# Feature(1,2...) HMDB_ID Metabolite mz RT Method
annotations <- mbx_mol_info[,c("HMDB_ID","Metablite","MZ","RT","Method")]
names(annotations)[names(annotations) == "Metablite"] <- "Metabolite"
rownames(annotations) <- NULL
addcol <- as.data.frame(rownames(annotations))
colnames(addcol) <- "Feature"
annotations <- cbind(addcol,annotations)
annotations[annotations == 0] <- NA
annotations[is.na(annotations)] <- ""
write.csv(annotations, file="./mbx/data/mac/annotations.csv",row.names=F)


############
# 3-metadata
############
metadata <- read.table('./metadata/micron.divert.metadata.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
#categorical variables
metadata <- metadata %>%
  mutate(caco = case_when(
  caco==1 ~ 'Divert',
  caco==2 ~ 'Control'))
metadata$caco<- factor(metadata$caco,
                                levels = c('Control', 'Divert'))

metadata <- metadata %>%
  mutate(pmh = case_when(
    pmh==1 ~ 'Nonuser',
    pmh==2 ~ 'Hormone user'))
metadata$pmh<- factor(metadata$pmh,
                       levels = c('Nonuser', 'Hormone user'))


metadata <- metadata %>%
  mutate(smk19 = case_when(
    smk19=='no' ~ 'Non smoker',
    smk19=='yes' ~ 'Smoker'))
metadata$smk19<- factor(metadata$smk19,
                      levels = c('Non smoker', 'Smoker'))

metadata <- metadata %>%
  mutate(abx = case_when(
    abx1m==0 ~ 'Non abx use',
    abx1m==1 ~ 'abx use'))
metadata$abx<- factor(metadata$abx,
                      levels = c('Non abx use', 'abx use'))

metadata <- metadata %>%
  mutate(bristolcat = case_when(
    bristolcat==1 ~ 'Hard',
    bristolcat==2 ~ 'Normal',
    bristolcat==3 ~ 'Soft',
    TRUE ~ 'Normal'))
metadata$bristolcat<- factor(metadata$bristolcat,
                      levels = c('Normal', 'Hard', 'Soft'))


metadata$race8905<- factor(metadata$race8905,
                      levels = c('white', 'black'))


metadata <- merge(idkey,metadata,by="id")

#metadata <- metadata[,!colnames(metadata) %in% c("id","abx1m","qq","kitret_date","incident","retmo19")]
metadata <- metadata[,colnames(metadata) %in% c("SAMPLEID","caco","race8905","pmh","smk19","act17m","calor15n","alco15n","aofib15a","ahei2010_noETOH15","bristolcat","abx","age","bmi19")]

metadata <- metadata[metadata$SAMPLEID %in% colnames(mbx_full_df_export),]
names(metadata)[names(metadata) == "SAMPLEID"] <- "sample"
write.csv(metadata, file="./mbx/data/mac/metadata_mbx.csv",row.names=F)

############
# 4-taxonomy
############
library(Macarron)
?Macarron
# MACARRoN can only be installed in R 4.2 and newer, so I use local R
# Import annotation table as a dataframe
annotations_df <- read.csv(file="./mbx/data/mac/annotations.csv", row.names = 1)

#source("/n/holystore01/LABS/huttenhower_lab/Users/abhosle/Macarron_tool/decorateID.R")
source("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/code/decorateID.R")

# Create taxonomy file
met_taxonomy <- decorateID(annotations_df)
write.csv(met_taxonomy, file="./mbx/data/mac/taxonomy.csv")



################################################################################
# Annotated metabolites for MaAsLin2 #
################################################################################
##annotated metabolites
abundances_df = read.csv(file = "./mbx/data/mac/abundances.csv", row.names = 1) # setting features as rownames
abundances_df <- dplyr::mutate_all(abundances_df, ~replace(., is.na(.), 0))
mbx_anno <- cbind(mbx_mol_info,abundances_df)
mbx_anno <- mbx_anno[!mbx_anno$Metablite == 0,]
write.csv(mbx_anno,file="./mbx/data/mbx_anno_df.csv")


#names of the 552 annotated metablites
#dup <- mbx_anno[duplicated(mbx_anno$Metablite),]
mbx_maaslin_anno_name <- as.data.frame(mbx_anno$Metablite)
rownames(mbx_maaslin_anno_name) <- rownames(mbx_anno)
colnames(mbx_maaslin_anno_name) <- "Metabolite name"
mbx_maaslin_anno_name$mol_id <- rownames(mbx_maaslin_anno_name)
write.table(mbx_maaslin_anno_name,file="./mbx/data/mbx_maaslin_anno_name.tsv", row.names=FALSE, sep="\t")

mbx_maaslin_anno <- mbx_anno[,8:dim(mbx_anno)[2]]
mbx_maaslin_anno <- as.data.frame(t(mbx_maaslin_anno))
mbx_maaslin_anno$sample <- rownames(mbx_maaslin_anno)
mbx_maaslin_anno <- mbx_maaslin_anno %>%
  select(sample, everything())

write.table(mbx_maaslin_anno,file="./mbx/data/mbx_maaslin_anno.tsv", row.names=FALSE, sep="\t")
# 232 samples, 552 metabolites


############
# metadata #
############
metadata <- read.table('./metadata/micron.divert.metadata.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
idkey = read.csv(file = "./idkey/id_sampleid_mbx.csv")
idkey <- idkey[,-1]
metadata <- merge(idkey,metadata,by="id")

#categorical variables
metadata <- metadata %>%
  mutate(caco = case_when(
  caco==1 ~ 'Divert',
  caco==2 ~ 'Control'))
metadata$caco<- factor(metadata$caco,
                                levels = c('Control', 'Divert'))

metadata <- metadata %>%
  mutate(pmh = case_when(
    pmh==1 ~ 'Nonuser',
    pmh==2 ~ 'Hormone user'))
metadata$pmh<- factor(metadata$pmh,
                       levels = c('Nonuser', 'Hormone user'))


metadata <- metadata %>%
  mutate(smk19 = case_when(
    smk19=='no' ~ 'Non smoker',
    smk19=='yes' ~ 'Smoker'))
metadata$smk19<- factor(metadata$smk19,
                      levels = c('Non smoker', 'Smoker'))

metadata <- metadata %>%
  mutate(abx = case_when(
    abx1m==0 ~ 'Non abx use',
    abx1m==1 ~ 'abx use'))
metadata$abx<- factor(metadata$abx,
                      levels = c('Non abx use', 'abx use'))

metadata <- metadata %>%
  mutate(bristolcat = case_when(
    bristolcat==1 ~ 'Hard',
    bristolcat==2 ~ 'Normal',
    bristolcat==3 ~ 'Soft',
    TRUE ~ 'Normal'))
metadata$bristolcat<- factor(metadata$bristolcat,
                      levels = c('Normal', 'Hard', 'Soft'))


metadata$race8905<- factor(metadata$race8905,
                      levels = c('white', 'black'))

dup <- metadata[duplicated(metadata$id),]

#metadata <- metadata[,!colnames(metadata) %in% c("id","abx1m","qq","kitret_date","incident","retmo19")]
metadata <- metadata[,colnames(metadata) %in% c("SAMPLEID","caco","race8905","pmh","smk19","act17m","calor15n","alco15n","aofib15a","ahei2010_noETOH15","bristolcat","abx","age","bmi19")]

metadata <- metadata[metadata$SAMPLEID%in%colnames(mbx_full_df_export),]
names(metadata)[names(metadata) == "SAMPLEID"] <- "sample"
write.table(metadata, file="./mbx/data/mbx_maaslin_metadata.tsv",row.names=F,sep="\t")


############
# metadata with multiple, single, control#
############
# metadata <- read.table('./metadata/micron.divert.metadata.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
# idkey = read.csv(file = "./idkey/id_sampleid_mbx.csv")
# idkey <- idkey[,-1]
# metadata <- merge(idkey,metadata,by="id")
# 
# #categorical variables
# metadata <- metadata %>%
#   mutate(caco3 = case_when(
#   dvamult2yr==1 ~ 'Single',
#   dvamult2yr==2 ~ 'Multiple',
#   TRUE ~ 'Control'))
# metadata$caco3<- factor(metadata$caco3,
#                                 levels = c('Control', 'Single', 'Multiple'))
# 
# metadata <- metadata %>%
#   mutate(pmh = case_when(
#     pmh==1 ~ 'Nonuser',
#     pmh==2 ~ 'Hormone user'))
# metadata$pmh<- factor(metadata$pmh,
#                        levels = c('Nonuser', 'Hormone user'))
# 
# 
# metadata <- metadata %>%
#   mutate(smk19 = case_when(
#     smk19==1 ~ 'Non smoker',
#     smk19==2 ~ 'Smoker'))
# metadata$smk19<- factor(metadata$smk19,
#                       levels = c('Non smoker', 'Smoker'))
# 
# metadata <- metadata %>%
#   mutate(abx = case_when(
#     abx1m==0 ~ 'Non abx use',
#     abx1m==1 ~ 'abx use'))
# metadata$abx<- factor(metadata$abx,
#                       levels = c('Non abx use', 'abx use'))
# 
# metadata <- metadata %>%
#   mutate(bristolcat = case_when(
#     bristolcat==1 ~ 'Hard',
#     bristolcat==2 ~ 'Normal',
#     bristolcat==3 ~ 'Soft',
#     TRUE ~ 'Normal'))
# metadata$bristolcat<- factor(metadata$bristolcat,
#                       levels = c('Normal', 'Hard', 'Soft'))
# 
# 
# metadata$race8905<- factor(metadata$race8905,
#                       levels = c('white', 'black'))
# 
# dup <- metadata[duplicated(metadata$id),]
# 
# #metadata <- metadata[,!colnames(metadata) %in% c("id","abx1m","qq","kitret_date","incident","retmo19")]
# metadata <- metadata[,colnames(metadata) %in% c("SAMPLEID","caco3","race8905","pmh","smk19","act17m","calor15n","alco15n","aofib15a","ahei2010_noETOH15","bristolcat","abx","age","bmi19")]
# 
# metadata <- metadata[metadata$SAMPLEID%in%colnames(mbx_full_df_export),]
# names(metadata)[names(metadata) == "SAMPLEID"] <- "sample"
# write.table(metadata, file="./mbx/data/mbx_maaslin_metadata_multiple.tsv",row.names=F,sep="\t")



############
# metadata with severe, non-severe, control#
############
metadata <- read.table('./metadata/micron.divert.metadata.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
idkey = read.csv(file = "./idkey/id_sampleid_mbx.csv")
idkey <- idkey[,-1]
metadata <- merge(idkey,metadata,by="id")

#categorical variables
metadata <- metadata %>%
  mutate(caco3 = case_when(
  severe==1 ~ 'Severe',
  severe==0 ~ 'Non-severe',
  TRUE ~ 'Control'))
metadata$caco3<- factor(metadata$caco3,
                                levels = c('Control', 'Non-severe', 'Severe'))

metadata <- metadata %>%
  mutate(pmh = case_when(
    pmh==1 ~ 'Nonuser',
    pmh==2 ~ 'Hormone user'))
metadata$pmh<- factor(metadata$pmh,
                       levels = c('Nonuser', 'Hormone user'))


metadata <- metadata %>%
  mutate(smk19 = case_when(
    smk19=='no' ~ 'Non smoker',
    smk19=='yes' ~ 'Smoker'))
metadata$smk19<- factor(metadata$smk19,
                      levels = c('Non smoker', 'Smoker'))

metadata <- metadata %>%
  mutate(abx = case_when(
    abx1m==0 ~ 'Non abx use',
    abx1m==1 ~ 'abx use'))
metadata$abx<- factor(metadata$abx,
                      levels = c('Non abx use', 'abx use'))

metadata <- metadata %>%
  mutate(bristolcat = case_when(
    bristolcat==1 ~ 'Hard',
    bristolcat==2 ~ 'Normal',
    bristolcat==3 ~ 'Soft',
    TRUE ~ 'Normal'))
metadata$bristolcat<- factor(metadata$bristolcat,
                      levels = c('Normal', 'Hard', 'Soft'))


metadata$race8905<- factor(metadata$race8905,
                      levels = c('white', 'black'))

dup <- metadata[duplicated(metadata$id),]

#metadata <- metadata[,!colnames(metadata) %in% c("id","abx1m","qq","kitret_date","incident","retmo19")]
metadata <- metadata[,colnames(metadata) %in% c("SAMPLEID","caco3","race8905","pmh","smk19","act17m","calor15n","alco15n","aofib15a","ahei2010_noETOH15","bristolcat","abx","age","bmi19")]

metadata <- metadata[metadata$SAMPLEID%in%colnames(mbx_full_df_export),]
names(metadata)[names(metadata) == "SAMPLEID"] <- "sample"
write.table(metadata, file="./mbx/data/mbx_maaslin_metadata_severe.tsv",row.names=F,sep="\t")











# 
# 
# 
# 
# ####create_annoted_dataframe.R
#   info <- mbx_mol_info
#   dat_df <- mbx_full_df
#   annotated_mols <- row.names(mbx_mol_info[!mbx_mol_info$Metablite == 0,])
#   annotated_mols <- annotated_mols[2:length(annotated_mols)]
#   resolve_duplicates <- function(dat_df, annotated_mols, info, rows_or_cols){
#     duplicates <- info[!info$Metablite == '' & duplicated(info$Metablite),]$Metablite
#     if(rows_or_cols == 'rows'){
#       for(mol in duplicates){
#         duplicated_mols <- row.names(info[info$Metablite %in% mol,])
#         annotated_mols <- annotated_mols[!annotated_mols %in% duplicated_mols]
#         annotated_mols <- c(annotated_mols, names(which.max(apply(dat_df[duplicated_mols,],1,median))))
#       }
#     }
#     if(rows_or_cols == 'cols'){
#       for(mol in duplicates){
#         duplicated_mols <- row.names(info[info$Metablite %in% mol,])
#         annotated_mols <- annotated_mols[!annotated_mols %in% duplicated_mols]
#         annotated_mols <- c(annotated_mols, names(which.max(apply(dat_df[,duplicated_mols],2,median))))
#       }
#     }
#     return(annotated_mols)
#   }
# 
#   if(annotated_mols %in% row.names(dat_df)){
#     annotated_mols <- intersect(annotated_mols, row.names(dat_df))
#     annotated_mols <- resolve_duplicates(dat_df, annotated_mols, info, 'rows')
#     dat_df <- dat_df[annotated_mols,]
#     row.names(dat_df) <- make.unique(info[annotated_mols,]$Metablite)
#   }else{
#     annotated_mols <- intersect(annotated_mols, colnames(dat_df))
#     annotated_mols <- resolve_duplicates(dat_df, annotated_mols, info, 'cols')
#     dat_df <- dat_df[,annotated_mols]
#     colnames(dat_df) <- make.unique(info[annotated_mols,]$Metablite)
#   }
#   return(dat_df)
#   }
