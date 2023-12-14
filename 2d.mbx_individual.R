#########################################################################################################
### Purpose: Boxplot of individual metabolites according to case/control status and severity
### Last update date: 10/13/2023 
#########################################################################################################


#setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

library(readxl)
library(dplyr)
library(stringr)
library(lmerTest)
library(ggplot2)
library(ggpubr)

################################################################################
###                           Data preparation                              ####
################################################################################
# read in annotated metabolites
mbx_anno_df <- read.csv(file="./mbx/data/mbx_anno_df.csv") #552 239
rownames(mbx_anno_df) <- mbx_anno_df[,1]
mbx_anno_df <- mbx_anno_df[,-1]

mbx_anno_df <- mbx_anno_df[,8:dim(mbx_anno_df)[2]] #552 232

# transform data
t_mbx_anno_df <- as.data.frame(t(mbx_anno_df)) #232 samples
#t_mbx_anno_df <- log2(t_mbx_anno_df)
# change -Inf to 0
#t_mbx_anno_df[t_mbx_anno_df== "-Inf"] <- 0

# TSS
mbx_norm <-t(apply(t_mbx_anno_df, 1, function(x) x/sum(x)))

# Log transformation
LOG <- function(x) {
    y <- replace(x, x == 0, min(x[x>0]) / 2)
    return(log2(y))
}

mbx_norm_log <- apply(t_mbx_anno_df, 2, LOG)

# read in metadata
meta <- read.table('./mbx/data/mbx_maaslin_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$sample
meta <- meta[,-1] # 232 13

# merge metabolites with metadata
mbx <- merge(mbx_norm_log,meta,by="row.names")
rownames(mbx) <- mbx[,1]
mbx <- mbx[,-1]

meta_severe <- read.table('./mbx/data/mbx_maaslin_metadata_severe.tsv', header=TRUE, sep='\t')
rownames(meta_severe) <- meta_severe[,1]
meta_severe <- meta_severe[,-1]
severe <- as.data.frame(meta_severe[,colnames(meta_severe) %in% c('caco3')])
rownames(severe) <- rownames(meta_severe)
colnames(severe) <- 'caco3'

mbx2 <- merge(mbx,severe,by='row.names')
rownames(mbx2) <- mbx2[,1]
mbx2 <- mbx2[,-1]


#############################################################################################################################################
                                        #### barplot for metabolites according to caco #####
#############################################################################################################################################
boxplot_plot <- function(y,ylab){
p1 <- ggplot(mbx2,aes_string(x='caco', y=y,fill='caco'))+
  geom_point(size=8, position=position_jitterdodge(),aes(fill=caco),alpha=1,shape=21,stroke=1,color='black',show.legend = FALSE)+
  geom_boxplot(alpha=0.5,size=1.5)+
  scale_fill_manual(values=c("#E69F00", "#4169E1"))+
  theme_classic() +
  xlab("Status") + ylab(ylab) +
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
filepath <- paste('./mbx/results/maaslin2/individual/boxplot/', ylab, '.png', sep ="")
ggsave(filename=filepath, plot=p1, width = 10, height = 10, dpi = 300)
}

boxplot_plot('QI44345_HILIC_pos',"3-Methylhistidine")
boxplot_plot('QI06254_HILIC_pos',"Cer 18:1;O2/17:0")
boxplot_plot('TF16_HILIC_neg',"3-Methyladipic acid or Pimelic acid")
boxplot_plot('QI32557_C18_neg',"Stearic acid")
boxplot_plot('QI18892_HILIC_pos',"1,7-Dimethyluric acid")
boxplot_plot('QI18714_C8_pos',"Cer 18:1;O2/16:0")
boxplot_plot('QI55567_C18_neg',"PGF2beta")
boxplot_plot('TF120_HILIC_neg',"N-Acetylglutamic acid")
boxplot_plot('QI28729_C18_neg',"DPA")
boxplot_plot('QI24126_C18_neg',"9-cis-Retinoic acid")
boxplot_plot('QI07890_C8_pos',"LPC 18:3")
boxplot_plot('QI06417_HILIC_pos',"Cer 18:1;O2/16:0")



#############################################################################################################################################
                                        #### barplot for metabolites according to severity #####
#############################################################################################################################################

my_comparisons <- list( c("Control", "Non-severe"), c("Control", "Severe"), c("Non-severe", "Severe") )
boxplot_plot <- function(y,ylab){
p1 <- ggplot(mbx2,aes_string(x='caco3', y=y,fill='caco3'))+
  geom_point(size=8, position=position_jitterdodge(),aes(fill=caco3),alpha=1,shape=21,stroke=1,color='black',show.legend = FALSE)+
  geom_boxplot(alpha=0.5,size=1.5)+
  scale_fill_manual(values=c("#E69F00", "#4169E1","darkslateblue"))+
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  xlab("Status") + ylab(ylab) +
  guides(legend.position=NULL)+
  stat_compare_means(comparisons = my_comparisons,
    aes(label = paste0("p = ", ..p.format..)), size=12)+
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
filepath <- paste('./mbx/results/maaslin2/individual/severe/boxplot/', y, '.png', sep ="")
ggsave(filename=filepath, plot=p1, width = 10, height = 12, dpi = 300)
}
boxplot_plot('TF2_HILIC_neg','1-Methyluric acid')
boxplot_plot('QI22732_C8_pos','Cer 18:1;O2/2:0')
boxplot_plot('QI06254_HILIC_pos',"Cer 18:1;O2/17:0")
boxplot_plot('QI30733_HILIC_pos',"5-Hydroxy-tryptophan")

















#############
##regression
#############
#fit <-  lm(formula = all[,1] ~ caco +race8905 + pmh + smk19 + act17m + calor15n + alco15n + aofib15a + ahei2010_noETOH15 + bristolcat + abx + bmi19, data=all)
for(i in 1:552){
  fit <-  lm(formula = mbx[,i] ~ caco + age + race8905 + pmh + smk19 + act17m + calor15n + alco15n + aofib15a + ahei2010_noETOH15 + bristolcat + abx + bmi19, data=mbx)
  out <- summary(fit)
  x <- out$coefficients[row.names(out$coefficients)=='cacoDivert',]
  if (i == 1){
    coef <- as.data.frame(t(x))
  } else {
    coef_c = as.data.frame(t(x))
    coef <- rbind(coef,coef_c)
  }
}

coef <- coef

rownames(coef) <- colnames(mbx[,1:552])
coef <- coef[order(coef$`Pr(>|t|)`),]

coef$fdr <-
  p.adjust(
    coef$`Pr(>|t|)`,
    method = "BH",
    n = length(coef$`Pr(>|t|)`)
  )

  p.adjust(
    0.00262506251380625,
    method = "BH")

# sig <- coef[coef$`Pr(>|t|)`<0.05,]
# sig_class <- sig_class[order(-sig_class$Estimate),]
#
# sig_class$UCI <- sig_class$Estimate+1.96*sig_class$`Std. Error`
# sig_class$LCI <- sig_class$Estimate-1.96*sig_class$`Std. Error`
#
# sig_class$Class <- rownames(sig_class)




################################################################################
# QI44345_HILIC_pos 3-Methylhistidine
# QI18892_HILIC_pos 1,7-Dimethyluric acid

# QI04432_HILIC_pos Dihydrocholesterol
# QI25330_HILIC_pos Kynurenic acid
# TF48_HILIC_neg Butyric acid
# QI07781_HILIC_pos 2-Oxindole-3-acetate
# QI06105_HILIC_pos 5-Hydroxytryptophol

boxplot <- function(metabolite,name){
plot <- ggplot(mbx, aes(x=caco, y=metabolite, fill=caco)) +
    geom_boxplot()
filepath <- paste('./mbx/results/', name, sep ="" )
ggsave(filename=filepath, plot=plot, width = 10, height = 10, dpi = 300)
}

boxplot(mbx$QI04432_HILIC_pos,"mbx_dihydrocholesterol.jpg")
boxplot(mbx$QI25330_HILIC_pos,"mbx_kynurenic.jpg")
boxplot(mbx$TF48_HILIC_neg,"mbx_butyric.jpg")
boxplot(mbx$QI07781_HILIC_pos,"mbx_acetate.jpg")
boxplot(mbx$QI06105_HILIC_pos,"mbx_hydroxytryptophol.jpg")


boxplot(mbx$QI44345_HILIC_pos,"mbx_3-Methylhistidine.jpg")



fit <-  lm(formula = mbx$QI04432_HILIC_pos ~ caco + age + race8905, data=mbx)
summmary(fit)

res <- wilcox.test(QI06105_HILIC_pos~ caco, data = mbx,
                   exact = FALSE)
res
