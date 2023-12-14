################################################################################
###                           combined permanova                            ####
### outcome:    fecal metabolomics (known)
### predictors: microbiome (10 PCs) and metadata
### last update date: 10/13/2023 (rank factors by R2)
################################################################################

#salloc -p test --mem 64G -t 0-08:00
#module load R/4.2.2-fasrc01
#install new packages
#export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
#R --quiet

#setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

library(scales)
library(tidyverse)
library(vegan)

################################################################################
###                           species                              ####
################################################################################
# read in species
species_all <- read.table('./mgx/data/micron_species_filtered.tsv', header=TRUE, sep='\t')
# 189 species 235 samples

################################################################################
###                           annotated metabolites               ####
################################################################################
mbx_anno_df <- read.csv(file="./mbx/data/mbx_anno_df.csv",row.names=1) # 552 metabolites 232 samples
mbx_anno_df <- mbx_anno_df[,8:dim(mbx_anno_df)[2]]


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
###                           metadata                       ####
################################################################################
meta <- read.table('./mbx/data/mbx_maaslin_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$sample
# 232 samples


################################################################################
###                           id matched files                        ####
################################################################################
# change mbx id to match mgx file
species_all <- merge(id_mgx_mbx,species_all,by="UID")
species_all <- species_all[,!colnames(species_all) %in% c("UID","id")]

mbx <- as.data.frame(t(mbx_anno_df))
mbx$SAMPLEID <- rownames(mbx)

inta <- intersect(species_all$SAMPLEID,mbx$SAMPLEID)
species <- species_all[species_all$SAMPLEID %in% inta,]

rownames(species) <- species$SAMPLEID

species <- species[,-1]
species <- species[order(rownames(species)),]

mbx <- mbx[mbx$SAMPLEID %in% inta,]
mbx <- mbx[,!colnames(mbx) %in% "SAMPLEID"]
mbx <- mbx[order(rownames(mbx)),]

# replace 0
IMPUTE <- function(x) {
    y <- replace(x, x == 0, min(x[x>0]) / 2)
    return(y)
}

# Log transformation
LOG <- function(x) {
    y <- replace(x, x == 0, min(x[x>0]) / 2)
    return(log2(y))
}

mbx_log <- apply(mbx, 2, LOG)
mbx_imp <- apply(mbx, 2, IMPUTE)


meta <- meta[rownames(meta) %in% inta,]
meta <- meta[order(rownames(meta)),]


################################################################################
###                           get 10 PCs                        ####
################################################################################

# get 15 PCos
data_bray <- vegdist(species, method = 'bray')
cmdscale(data_bray, k = 15, eig = TRUE)$GOF # gives idea of total variance explained by first 15 eigenvectors
data_bray_pcoa <- cmdscale(data_bray, k = (nrow(species) - 1), eig = TRUE) # ordination

# extract % var explained
pcos <- as.data.frame(data_bray_pcoa$points[, 1:15]) %>%
     rename_with(~ paste0("PC", str_remove(.x, pattern = "V")))


################################################################################
###                           PERMANOVA                    ####
################################################################################
meta_pcos <- merge(meta,pcos,by="row.names")
rownames(meta_pcos) <- meta_pcos[,1]
meta_pcos <- meta_pcos[,-1]

table(rownames(mbx) == rownames(meta_pcos))


# PERMANOVA using adonis
# x <- adonis(formula = mbx_bray ~ caco, data = meta_pcos, permutations = 999)
# x.data <- as.data.frame(x$aov.tab)
# caco <- x.data[!is.na(x.data$F.Model), ]


# all
all <- as.data.frame(adonis2(
     mbx ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15
     + caco + age + race8905 + bmi19 + abx + bristolcat + smk19 + pmh + act17m + alco15n + aofib15a + calor15n,
     data = meta_pcos, permutations = 2000, method = "bray"
))%>%
     filter(!is.na(F)) %>%
     summarize(total = sum(R2)) %>%
     mutate(var = "all") %>%
     column_to_rownames("var") %>%
     rename(R2 = total)
# this is equivalent to use by=NULL

# microbiome 15 pc
pc <- as.data.frame(adonis2(
     mbx ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15,
     data = meta_pcos, permutations = 999, method = "bray"
))%>%
     filter(!is.na(F)) %>%
     summarize(total = sum(R2)) %>%
     mutate(var = "pc") %>%
     column_to_rownames("var") %>%
     rename(R2 = total)

# divert
caco <- as.data.frame(adonis2(mbx ~ caco, data = meta_pcos, permutations = 999, method = "bray"))

# age
age <- as.data.frame(adonis2(mbx ~ age, data = meta_pcos, permutations = 999, method = "bray"))

# race
race <- as.data.frame(adonis2(mbx ~ race8905, data = meta_pcos, permutations = 999, method = "bray"))

# bmi
bmi <- as.data.frame(adonis2(mbx ~ bmi19, data = meta_pcos, permutations = 999, method = "bray"))

# abx
abx <- as.data.frame(adonis2(mbx ~ abx, data = meta_pcos, permutations = 999, method = "bray"))

# Bristol
bristol <- as.data.frame(adonis2(mbx ~ bristolcat, data = meta_pcos, permutations = 999, method = "bray"))

# smk
smk <- as.data.frame(adonis2(mbx ~ smk19, data = meta_pcos, permutations = 999, method = "bray"))

# pmh
pmh <- as.data.frame(adonis2(mbx ~ pmh, data = meta_pcos, permutations = 999, method = "bray"))

# act
act <- as.data.frame(adonis2(mbx ~ act17m, data = meta_pcos, permutations = 999, method = "bray"))

# alc
alc <- as.data.frame(adonis2(mbx ~ alco15n, data = meta_pcos, permutations = 999, method = "bray"))

# aofib
aofib <- as.data.frame(adonis2(mbx ~ aofib15a, data = meta_pcos, permutations = 999, method = "bray"))

# calorie
calor <- as.data.frame(adonis2(mbx ~ calor15n, data = meta_pcos, permutations = 999, method = "bray"))


# extract adonis results
adonis_all <- bind_rows(all, pc, caco, age, race, bmi, abx, bristol, smk, pmh, act, alc, aofib, calor) %>%
     rownames_to_column("var") %>%
     filter(!str_detect(string = var, pattern = "Total|Residual"))

# label for variables
component<- c('Combined model',
             'Microbiome',
             'Diverticulitis diagnosis',
             'Age',
             'Race',
             'Body mass index',
             'Antibiotic use',
             'Bristol stool scale',
             'Smoking',
             'Postmenopausal hormone use',
             'Physical activity',
             'Alcohol intake',
             'Dietary fiber',
             'Calorie intake')

adonis_all <- cbind(adonis_all, component)


# FDR adjustment for p-values
adonis_all$fdr <-
 p.adjust(
   adonis_all$`Pr(>F)`,
   method = "BH",
   n = length(adonis_all$`Pr(>F)`)
 )
# create figure to show PERMANOVA results
# combined model and pc do not have p value
adonis_all_star <- adonis_all

adonis_all_star$stars <- cut(adonis_all_star$fdr, breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))

adonis_all_star$stars[is.na(adonis_all_star$stars)] <- '***'

# rank factors by R2
adonis_all_star <- adonis_all_star[order(adonis_all_star$R2,decreasing=TRUE),]


level_y_order <- factor(adonis_all_star$component, level = rev(adonis_all_star$component))

plot <- ggplot(adonis_all_star, aes(level_y_order, R2)) +
 geom_bar(stat="identity",fill="cornflowerblue",color="black")+
 scale_y_continuous(labels=scales::percent, limits = c(0, 0.35))+
 geom_text(aes(label=stars), color="black", size=10, hjust = -0.5) +
 coord_flip()+
 theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
       panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
       legend.text = element_text(size = 20),
       legend.title = element_blank(),
       axis.ticks =element_blank(),
       plot.title = element_blank(),
       axis.title=element_blank(),
       axis.line = element_line(colour = 'black', size=0.75, linetype='solid'),
       axis.text.y=element_text(size=20, face="bold"),
       axis.text.x = element_text(size=20)
     )
ggsave(filename='./mbx/results/permanova/permanova_combine.jpg', plot=plot, width =15, height = 10, dpi = 300)
write.csv(adonis_all_star, file = "./mbx/results/permanova/permanova_combine.csv")
