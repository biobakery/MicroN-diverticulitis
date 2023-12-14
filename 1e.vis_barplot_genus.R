#########################################################################################################
### Purpose: To visualize top 20 genus in diverticulitis and controls ###
### Last update date: 10/13/2023 
#########################################################################################################


setwd("C:/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Microbiome data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)


#######################
####### Genus #########
#######################
taxa <- read.csv(file="./mgx/data/micron_species_filtered.tsv", header=TRUE, sep="\t")
rownames(taxa) <- taxa$UID
taxa <- taxa[,-1] # 235 189
species <- as.data.frame(t(taxa))

taxaname <- read.csv(file="./mgx/data/micron_taxaname_filtered.csv")
taxaname <- taxaname[,3:8]

species_category <- cbind(species,taxaname)

species_genus  <- species_category[,!colnames(species_category) %in% c("Family","Order","Class","Phylum","Domain")]

genus <- aggregate(species_genus[,1:235], list(species_genus$Genus), FUN=sum)

rownames(genus) <- genus$Group.1
genus <- genus[,-1]

genus_t <- as.data.frame(t(genus)) # 235 71

# sort genus by mean abundance
mns <- colMeans(genus_t, na.rm=TRUE)
order(-mns)
genus_sort <- genus_t[,order(-mns)]

genus_sort <- genus_sort[order(-genus_sort$Bacteroides),]

# sum normalize - transform to relative abundance
genus_sort <- sweep(genus_sort, 1, rowSums(genus_sort), `/`)

# read in metadata
meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID


# merge with caco
taxa <- merge(genus_sort,meta[,colnames(meta) %in% c("caco","UID")],by="row.names")
rownames(taxa) <- taxa[,1]
taxa <- taxa[,-1]
taxa <- taxa[,!colnames(taxa) %in% c('UID')]


###################################################################
####### boxplot of genus Bacteroides and Faecalibacterium #########
###################################################################
box_plot <- function(y,ylab){
p1 <- ggplot(taxa,aes_string(x='caco', y=y,fill='caco'))+
  geom_point(size=8, position=position_jitterdodge(),aes(fill=caco),alpha=1,shape=21,stroke=1,color='black',show.legend = FALSE)+
  geom_boxplot(alpha=0.5,size=1.5)+
  scale_fill_manual(values=c("#E69F00", "#4169E1"))+
  theme_classic() +
  #xlab("Status") + ylab(ylab) +
  guides(legend.position=NULL)+
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 30),
      legend.position="none",
      axis.title=element_text(size=30,face="bold"),
      axis.text=element_text(size=30),
      axis.ticks=element_blank(),
      axis.line = element_line(colour = 'black', size=0.75, linetype='solid')
      )
filepath <- paste('./mgx/results/vis/genus_boxplot/', y, '.png', sep ="")
ggsave(filename=filepath, plot=p1, width = 10, height = 10, dpi = 300)
}

box_plot('Bacteroides',"Bacteroides")
box_plot('Faecalibacterium',"Faecalibacterium")




############################################################
####### boxplot of 20 genera in divert and control #########
############################################################

taxa_divert <- taxa[taxa$caco=="Divert",]
taxa_control <- taxa[taxa$caco=="Control",]

taxa_divert <- taxa_divert[,!colnames(taxa_divert) %in% c('caco')]
taxa_control <- taxa_control[,!colnames(taxa_control) %in% c('caco')]

taxa_divert  <- taxa_divert [order(-taxa_divert$Bacteroides),]
taxa_control  <- taxa_control[order(-taxa_control$Bacteroides),]


taxa_divert <- as.data.frame(t(taxa_divert))
taxa_control <- as.data.frame(t(taxa_control))

# keep 20 genus, and group genus beyond 20 to "other"
if (nrow(taxa_divert)>20)
{
  other <- as.data.frame(t(colSums(taxa_divert[21:nrow(taxa_divert),])))
  rownames(other) <- "Other"
  taxa_divert <- rbind(taxa_divert[1:20,], other)
}

if (nrow(taxa_control)>20)
{
  other <- as.data.frame(t(colSums(taxa_control[21:nrow(taxa_control),])))
  rownames(other) <- "Other"
  taxa_control <- rbind(taxa_control[1:20,], other)
}

taxa_divert=taxa_divert[order(nrow(taxa_divert):1),]
taxa_control=taxa_control[order(nrow(taxa_control):1),]

##########
#barplot#
#########
bug_colors <- list(
  # Bacteroidetes (red/yellow)
  # Firmicutes (blue/purple)
  "Bacteroides"             = "red",
  "Faecalibacterium"        = "dark blue",
  "Blautia"                 = "lightskyblue",
  "Alistipes"               = "pink",
  "Firmicutes_unclassified" = "mediumslateblue",
  "Bifidobacterium"         = "green",
  "Ruminococcus"            = "steelblue",
  "Parabacteroides"         = "deep pink",
  "Eubacterium"             = "dodgerblue",
  "Clostridium"             = "cornflowerblue",
  "Roseburia"               = "blue",
  "Lachnospiraceae_unclassified"="darkslateblue",
  "Anaerostipes"            = "navy",
  "Akkermansia"             = "yellowgreen",
  "Ruminococcaceae_unclassified" = "blue3",
  "Streptococcus"           = "deepskyblue",
  "Fusicatenibacter"        = "royalblue",
  "Prevotella"              = "palevioletred",
  "Dorea"                   = "steelblue1",
  "Collinsella"             = "dodgerblue1",
  "Other"                   = "#aaaaaa")


bugcolor <- scale_fill_manual(values = unlist(bug_colors))

##########################################################
#start of barplot function
barplot <- function (data,figname){
stratified <- data
stratified$genus <- row.names(stratified)
stratified_names <- c(rownames(stratified))
stratified$genus <- factor(stratified$genus, ordered = TRUE, levels = stratified_names)

# melt
stratified_melt <- melt(stratified, id.vars = 'genus')
dim(stratified_melt) # 2457 3


bar<-ggplot(data=stratified_melt, aes(x=variable, y=value, fill=genus)) +
  geom_bar(stat="identity")+
  theme_classic() +
  theme(plot.title = element_text(size = 30),
        legend.title = element_text(size=30),
        legend.text = element_text(size=30),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=30,face="bold"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=30),
        axis.ticks=element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
        ) +
  guides(fill = guide_legend(ncol = 1,keyheight = 0.55)) +
  scale_y_continuous(expand = c(0, 0)) +
  bugcolor +
  xlab("Sample") +
  ylab("Relative abundance")
print(bar)

filepath <- paste('./mgx/results/vis/', figname,sep = "")
ggsave(filename=filepath, plot=bar, width = 20, height = 15, dpi = 300)
}# end of stacked barplot function

barplot(taxa_divert,"genus_divert_barplot.jpg")
barplot(taxa_control,"genus_control_barplot.jpg")
