
##############################################################################################################################################
# Purpose: To prepare metadata for multivariate analysis
# Last update date: 10/12/2023 (update using the finalized data from Channing)
#############################################################################################################################################
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")


############
# metadata for caco#
############
metadata <- read.table('./metadata/micron.divert.metadata.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
idkey = read.csv(file = "./idkey/id_UID_235.csv")
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



metadata <- metadata[,colnames(metadata) %in% c("UID","caco","race8905","pmh","smk19","act17m","calor15n","alco15n","aofib15a","ahei2010_noETOH15","bristolcat","abx","age","bmi19")]

write.table(metadata, file="./metadata/micron_metadata.tsv",row.names=F,sep="\t")


############
# metadata for severity#
############
metadata <- read.table('./metadata/micron.divert.metadata.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
idkey = read.csv(file = "./idkey/id_UID_235.csv")
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



metadata <- metadata[,colnames(metadata) %in% c("UID","caco3","race8905","pmh","smk19","act17m","calor15n","alco15n","aofib15a","ahei2010_noETOH15","bristolcat","abx","age","bmi19")]

write.table(metadata, file="./metadata/micron_metadata_severe.tsv",row.names=F,sep="\t")





############
# food data#
############
metadata <- read.table('./metadata/micron.divert.food.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
idkey = read.csv(file = "./idkey/id_UID_235.csv")
idkey <- idkey[,-1]
metadata <- merge(idkey,metadata,by="id")

metadata <- metadata[,colnames(metadata) %in% c("UID","procm15","rmeat15","ormeat15","othfish15","wine15","beer15", "liq15","tea15", "coff15" , "fruj15","yelveg15","toma15","lfveg15","othveg15","rgrain15", "snack15","sugbev15","lowbev15","pizza15","eggs15","marg15","lowdai15","dess15","wgrain15","fries15","nuts15","crsoup15", "saldre15","fish15","poul15", "butter15","hidai15","fruit15" ,"f115v","f215v")]
write.table(metadata, file="./metadata/micron_food.tsv",row.names=F,sep="\t")
