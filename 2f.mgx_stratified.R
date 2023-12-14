#########################################################################################################
### Purpose: stratified pathway
### Last update date: 10/13/2023 
#########################################################################################################


setwd("C:/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")


#setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Microbiome data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(tidyr)


######MicroN TAXA#####
idkey <- read.csv(file="./idkey/lab20200.diverticulitis.manifest.csv")
idkey$id <- substr(idkey$PARTICIPANTID,1,6) #keep first 6 digits
idkey <- idkey[order(idkey$MANUFACTURERBARCODE),]
idkey <- idkey[idkey$MANUFACTURERBARCODE != 360521094,]
idkey <- idkey[,colnames(idkey) %in% c("id","MANUFACTURERBARCODE") ]
idkey$UID <- idkey$MANUFACTURERBARCODE
idkey <- idkey[,!colnames(idkey) %in% c("MANUFACTURERBARCODE")]
write.csv(idkey,file="./idkey/id_UID.csv", na="")

dup <- idkey[duplicated(idkey$id),]
id       UID
18  849435 311586345
16  486227 311603716
204 530623 311670140
112 789590 354217942
111 818734 360522188
# 5 duplicated cohort ids


#######################
#########pathways######
#######################

path1 <- read.csv(file="./mgx/data/pathabundance_relab.tsv", header=TRUE, sep="\t")
#nrow=23359

#remove X0360521094 since it does not have any taxa
path2 <- path1[, !colnames(path1) %in% c("X0360521094_Abundance") ]

#merge with id
names(path2)[names(path2) == "X..Pathway"] <- "Pathway"
rownames(path2) <- path2$Pathway
path2$Pathway <- NULL
path_t <- as.data.frame(t(as.matrix(path2))) # transpose
path_t$UID <- gsub("X0","", rownames(path_t))
path_t$UID <- gsub("_Abundance","",path_t$UID)
path_t <- merge(idkey,path_t,by='UID') #n=240 23361

rownames(path_t) <- path_t$UID


dup <- path_t[duplicated(path_t$id),]
UID     id
13  311586345 849435
21  311603716 486227
32  311670140 530623
120 354217942 789590
190 360522188 818734

# take average of the dups
df_849435 <- path_t[path_t$id=="849435",]
df_849435 <- df_849435[,!colnames(df_849435) %in% c("id","UID")]
df_849435 <- as.data.frame(t(as.data.frame(colMeans(df_849435))))
rownames(df_849435) <- "311586345"

df_486227 <- path_t[path_t$id=="486227",]
df_486227 <- df_486227[,!colnames(df_486227) %in% c("id","UID")]
df_486227 <- as.data.frame(t(as.data.frame(colMeans(df_486227))))
rownames(df_486227) <- "311603716"

df_530623 <- path_t[path_t$id=="530623",]
df_530623 <- df_530623[,!colnames(df_530623) %in% c("id","UID")]
df_530623 <- as.data.frame(t(as.data.frame(colMeans(df_530623))))
rownames(df_530623) <- "311670140"

df_789590 <- path_t[path_t$id=="789590",]
df_789590 <- df_789590[,!colnames(df_789590) %in% c("id","UID")]
df_789590 <- as.data.frame(t(as.data.frame(colMeans(df_789590))))
rownames(df_789590) <- "354217942"

df_818734 <- path_t[path_t$id=="818734",]
df_818734 <- df_818734[,!colnames(df_818734) %in% c("id","UID")]
df_818734 <- as.data.frame(t(as.data.frame(colMeans(df_818734))))
rownames(df_818734) <- "360522188"

path_1 <- rbind(df_849435,df_486227,df_530623,df_789590,df_818734)
path_1$UID <- rownames(path_1)
path_1 <- path_1 %>%
  select(UID, everything())

path_t <- path_t[!path_t$id=="849435",]
path_t <- path_t[!path_t$id=="486227",]
path_t <- path_t[!path_t$id=="530623",]
path_t <- path_t[!path_t$id=="789590",]
path_t <- path_t[!path_t$id=="818734",]

path_t <- path_t [, !colnames(path_t) %in% c("id")]

path_final <- rbind(path_t,path_1) #235 23360 (23359 paths)

write.table(path_final,"./mgx/data/micron_path_stratified.tsv", row.names = FALSE,sep="\t")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Metadata ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
meta <- read.csv(file="./metadata/micron_metadata.tsv", header=TRUE, sep="\t")
meta <- meta[,colnames(meta) %in% c("UID","caco")]

humann_input <- merge(meta,path_final,by="UID")

humann_input <- as.data.frame(t(humann_input)) #23361 235

write.table(humann_input,"./mgx/data/humann_bar_input.tsv", row.names = TRUE, sep="\t",col.names=FALSE)


##humann3
source /net/fs-huttenhower01/srv/export/huttenhower_lab/share/tools/hutlab/src/hutlabrc.sh
hutlab load centos7/python3/humann3/3.5-devel

grep 'PWY66-429' /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/data/humann_bar_input.tsv | less -S

humann_barplot --input /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/data/humann_bar_input.tsv --focal-metadata caco --last-metadata caco --output /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/results/humann_barplot/PWY66-429.png --focal-feature PWY66-429 --sort sum metadata


grep 'PWY-5677' /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/data/humann_bar_input.tsv | less -S

humann_barplot --input /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/data/humann_bar_input.tsv --focal-metadata caco --last-metadata caco --output /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/results/humann_barplot/PWY-5677.png --focal-feature PWY-5677 --sort sum metadata

humann_barplot --input /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/data/humann_bar_input.tsv --focal-metadata caco --last-metadata caco --output /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/results/humann_barplot/PWY66-391.png --focal-feature PWY66-391 --sort sum metadata


grep 'PWY4FS-8' /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/data/humann_bar_input.tsv | less -S

humann_barplot --input /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/data/humann_bar_input.tsv --focal-metadata caco --last-metadata caco --output /n/holystore01/LABS/huttenhower_lab/Users/wma/divert/mgx/results/humann_barplot/PWY4FS-8.png --focal-feature PWY4FS-8 --sort sum metadata
