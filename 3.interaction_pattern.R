
#########################################################################################################
### Purpose: To explore interactions between diet and microbiome composition on diverticulitis ###
### Last update date: 10/13/2023 (update using the finalized metadata)
#########################################################################################################


# salloc -p test --mem 64G -t 0-08:00
# module load R/4.2.2-fasrc01
# #install new packages
# export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
# R --quiet
# 
# 
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

# 
# rm(list=ls())
# getwd()
# setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")


library(vegan)
library(aod)
library(dplyr)
library(interactionR)
library(ggplot2)


########################################################################################
#######################       read in dataset ########################
########################################################################################
# read in taxonomy data
species_all <- read.table('./mgx/data/micron_species_filtered.tsv', header=TRUE, sep='\t')
rownames(species_all) <- species_all$UID
species_all <- species_all[,-1]

# read in metadata
meta <- read.table('./metadata/micron_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$UID

#### Note: I need to make sure that samples in rows
bugs_pcoa <- species_all
#samples in rows
bugs_pcoa <- capscale(bugs_pcoa ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )

ord.bug.scores <- ord.bug.scores[order(ord.bug.scores$MDS1),]


# append metadata to scores
ord.scores.meta <- merge(ord.scores, meta, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.scores.meta) <- ord.scores.meta[ ,1]
ord.scores.meta[ ,1] <- NULL

# number of samples
n_samples <- nrow(ord.scores.meta)

# association dataset
association <- merge(ord.scores,meta,by='row.names')
rownames(association) <- association[,1]
association <- association[,-1]

association <- merge(association,species_all,by='row.names')


association <- association %>%
  mutate(caco = case_when(
    caco=='Divert' ~ 1,
    caco=='Control' ~ 0))
association$caco<- factor(association$caco,
                          levels = c(0, 1))

rownames(association) <- association[,1]
association$UID <- rownames(association)


##############################################################
######### food groups and prudent/western dietary pattern ######
##############################################################
food <- read.table('./metadata/micron_food.tsv', header=TRUE, sep='\t')
rownames(food) <- food$UID


#dup1 <- association[duplicated(rownames(association)),]
#dup2 <- food[duplicated(rownames(food)),]


##########################
######### merge data######
##########################
association_food <- merge(association,food,by="UID")

# MDS2
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-2.01107 -0.47062  0.05567  0.00000  0.55535  1.28761

# f115v
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-1.82218 -0.49711  0.07473  0.18193  0.79275  3.02685

# Bacteroides_uniformis
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00000 0.02232 0.06286 0.08962 0.13132 0.48374



association_food <- association_food %>%
  mutate(mds2m = case_when(
  MDS2<=  0.05567 ~ 0,
  MDS2> 0.05567 ~ 1))

association_food <- association_food %>%
  mutate(bunif2c = case_when(
  Bacteroides_uniformis<=  0.06286 ~ 0,
  Bacteroides_uniformis> 0.06286 ~ 1))

association_food <- association_food %>%
  mutate(f115vm = case_when(
  f115v<0.07473 ~ 0,
  f115v>=0.07473 ~ 1))

association_food<- association_food %>%
  mutate(aofibm = case_when(
  aofib15a<22.82 ~ 0,
  aofib15a>=22.82 ~ 1))

# create joint categories
association_food <- association_food %>%
  mutate(mds2_f115v = case_when(
  mds2m==1 & f115vm==0 ~ 0,
  mds2m==1 & f115vm==1 ~ 1,
  mds2m==0 & f115vm==0 ~ 2,
  mds2m==0 & f115vm==1 ~ 3
))
association_food$mds2_f115v<- factor(association_food$mds2_f115v,
                                levels = c(0, 1, 2, 3))

association_food <- association_food %>%
  mutate(pcopri2c = case_when(
    Prevotella_copri>0 ~ 1,
    Prevotella_copri==0 ~ 0))

#####################################################################################

##########################
######### MDS1*f115v######
##########################
int_mds1_f115v <- glm(caco ~ MDS1 *f115v  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_mds1_f115v <- as.data.frame(summary(int_mds1_f115v)$coefficients)
print(coef_mds1_f115v)

##########################
######### MDS1*f215v######
##########################
int_mds1_f215v <- glm(caco ~ MDS1 *f215v  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_mds1_f215v <- as.data.frame(summary(int_mds1_f215v)$coefficients)
print(coef_mds1_f215v)


##########################
######### MDS2*f115v######
##########################
int_mds2_f115v <- glm(caco ~ MDS2 *f115v  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_mds2_f115v <- as.data.frame(summary(int_mds2_f115v)$coefficients)
print(coef_mds2_f115v)
#0.049 (consistent in crude models)

##########################
######### MDS2*f215v######
##########################
int_mds2_f215v <- glm(caco ~ MDS2 *f215v  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_mds2_f215v <- as.data.frame(summary(int_mds2_f215v)$coefficients)
print(coef_mds2_f215v)
#0.45

#####check pco loadings
# F prau and B uniformis have top loadings for both MDS1 and MDS2
##########################
library(viridis)
pcoa.plot.bunif <-
  ggplot(association_food, aes(MDS1, MDS2) ) +
  geom_point(size = 5, alpha = 1, aes( fill=  Bacteroides_uniformis),shape=21,color='black') +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 16, face = 'bold'),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'horizontal',
        plot.title = element_text( size = 20, face = 'bold')
        )+
        scale_fill_gradientn(colours = viridis(10,direction=-1))+
        xlab(paste("PCo1 (10%)",sep = "")) + ylab(paste("PCo2 (7%)",sep = ""))
ggsave("./mgx/results/pcoa/pcoa_bunif.jpg", plot=pcoa.plot.bunif, width =7, height =6, dpi=300)

library(Hmisc)
rcorr(association_food$MDS2,association_food$Bacteroides_uniformis) # r=0.71
rcorr(association_food$MDS1,association_food$Bacteroides_uniformis) # r=0.47




##########################
#Faecalibacterium_prausnitzii*f115v
##########################
int_fprau_f115v <- glm(caco ~ Faecalibacterium_prausnitzii*f115v  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_fprau_f115v <- as.data.frame(summary(int_fprau_f115v)$coefficients)
print(coef_fprau_f115v)
#0.78

##########################
#Prevotella_copri*f115v
##########################
int_pcopri_f115v <- glm(caco ~ Prevotella_copri*f115v  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_pcopri_f115v <- as.data.frame(summary(int_pcopri_f115v)$coefficients)
print(coef_pcopri_f115v)
#0.11

##########################
#Bacteroides_uniformis*f115v
##########################
int_bunif_f115v <- glm(caco ~ Bacteroides_uniformis*f115v  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_bunif_f115v <- as.data.frame(summary(int_bunif_f115v)$coefficients)
print(coef_bunif_f115v)
#0.03 (consistent in crude models)


##########################
#Bacteroides_uniformis*aofib15a
int_bunif_aofib <- glm(caco ~ Bacteroides_uniformis*aofib15a  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_bunif_aofib <- as.data.frame(summary(int_bunif_aofib)$coefficients)
print(coef_bunif_aofib)
#0.10

##########################
#MDS2*aofib15a
##########################
int_mds2_aofib <- glm(caco ~ MDS2*aofib15a  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_mds2_aofib <- as.data.frame(summary(int_mds2_aofib)$coefficients)
print(coef_mds2_aofib)
#0.05 (consistent in crude models)

##########################
#pcopri*aofib15a
##########################
int_pcopri_aofib <- glm(caco ~ Prevotella_copri*aofib15a  + age + race8905 + pmh + smk19 + act17m + calor15n +alco15n +abx + bristolcat +bmi19 , data = association_food, family = "binomial")
coef_pcopri_aofib <- as.data.frame(summary(int_pcopri_aofib)$coefficients)
print(coef_pcopri_aofib)
#0.20 (consistent in crude models)


##########################
######### MDS2*f115 subgroups######
##########################
mds2_high <- glm(caco ~ f115v +age + race8905 + pmh + act17m + calor15n +alco15n +abx + bristolcat +bmi19, data = association_food[association_food$mds2m==1,], family = "binomial")
summary(mds2_high)
exp(cbind(OR = coef(mds2_high), confint(mds2_high)))
#f115v            0.56384635 3.057522e-01  1.0005403

mds2_low <- glm(caco ~ f115v +age + pmh + act17m + calor15n +alco15n +abx + bristolcat +bmi19, data = association_food[association_food$mds2m==0,], family = "binomial")
summary(mds2_low)
exp(cbind(OR = coef(mds2_low), confint(mds2_low)))
#f115v            0.7707855 0.4452062426   1.297847

mds2_f115v_joint <- glm(caco ~ mds2_f115v +age + race8905 + pmh + act17m + calor15n +alco15n +abx + bristolcat +bmi19, data = association_food, family = "binomial")
summary(mds2_f115v_joint)
exp(cbind(OR = coef(mds2_f115v_joint), confint(mds2_f115v_joint)))
# mds2_f115v1      0.6378836 0.283404805  1.418874
# mds2_f115v2      0.7670696 0.354651189  1.647206
# mds2_f115v3      0.5018575 0.207789085  1.189332

##########################
######### Bacteroides_uniformis*f115v subgroups######
##########################
bunif_high <- glm(caco ~ f115v +age + race8905 + pmh + act17m + calor15n +alco15n +abx + bristolcat +bmi19, data = association_food[association_food$bunif2c==1,], family = "binomial")
summary(bunif_high)
exp(cbind(OR = coef(bunif_high), confint(bunif_high)))
#0.6422827 0.368135297   1.084958


bunif_low <- glm(caco ~ f115v +age + pmh + act17m + calor15n +alco15n +abx + bristolcat +bmi19, data = association_food[association_food$bunif2c==0,], family = "binomial")
summary(bunif_low)
exp(cbind(OR = coef(bunif_low), confint(bunif_low)))
#f115v            0.8578810 4.916693e-01  1.471513

##########################
#MDS2*aofib15a subgroups
##########################
mds2_high <- glm(caco ~ aofib15a +age + race8905 + pmh + act17m + calor15n +alco15n +abx + bristolcat +bmi19, data = association_food[association_food$mds2m==1,], family = "binomial")
summary(mds2_high)
exp(cbind(OR = coef(mds2_high)*5, confint(mds2_high)*5))
#aofib15a         0.706092869 4.898572e-01 9.982169e-01

mds2_low <- glm(caco ~ aofib15a +age + pmh + act17m + calor15n +alco15n +abx + bristolcat +bmi19, data = association_food[association_food$mds2m==0,], family = "binomial")
summary(mds2_low)
exp(cbind(OR = coef(mds2_low)*5, confint(mds2_low)*5))
#aofib15a         0.88327528 6.407932e-01 1.197863e+00


##########################
#pcopri*aofib15a subgroups
##########################
pcopri_high <- glm(caco ~ aofib15a +age + pmh + act17m + calor15n +alco15n +abx +bmi19, data = association_food[association_food$pcopri2c==1,], family = "binomial")
summary(pcopri_high)
exp(cbind(OR = coef(pcopri_high)*5, confint(pcopri_high)*5))
#aofib15a       9.485942e-01 3.597780e-01 2.310326e+00

pcopri_low <- glm(caco ~ aofib15a +age + race8905+ pmh + act17m + calor15n +alco15n +abx + bristolcat +bmi19, data = association_food[association_food$pcopri2c==0,], family = "binomial")
summary(pcopri_low)
exp(cbind(OR = coef(pcopri_low)*5, confint(pcopri_low)*5))
#aofib15a         0.78896175 6.181327e-01 9.960789e-01
