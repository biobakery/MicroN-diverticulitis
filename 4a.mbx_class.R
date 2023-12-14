

#########################################################################################################
# Purpose: compare mean levels of metabolites by class in diverticulitis and controls ###
# Last update date: 9/25/2023
#########################################################################################################


# salloc -p test --mem 64G -t 0-08:00
# module load R/4.2.2-fasrc01
# #install new packages
# export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
# R --quiet
# 
# 
# setwd("/n/holystore01/LABS/huttenhower_lab/Users/wma/divert/")
setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

library(readxl)
library(dplyr)
library(stringr)
library(lmerTest)
library(ggplot2)
library(data.table)

################################################################################
###                           Data preparation                              ####
################################################################################
# read in annotated metabolites
mbx_anno_df <- read.csv(file="./mbx/data/mbx_anno_df.csv") #552 239
mbx_anno_df <- mbx_anno_df[,-1]

# remove those without HMDB_ID since they are not assigned to class

mbx_anno_df <- mbx_anno_df[!is.na(mbx_anno_df$HMDB_ID),] #529 239

# read in taxonomy
taxonomy_df = read.csv(file = "./mbx/data/mac/taxonomy.csv", row.names = 1)

mbx_anno_df <- merge(mbx_anno_df,taxonomy_df,by="HMDB_ID") # 527 241


#name_class <- mbx_anno_df[,colnames(mbx_anno_df) %in% c('Metablite','Sub_Class','Class')]
#write.table(name_class,file="./mbx/data/name_class.tsv", row.names=FALSE, sep="\t")



# remove metabolite with missing class
mbx_anno_df <- mbx_anno_df[!mbx_anno_df$Class=='',] #523 241

mbx_anno_df <- mbx_anno_df[,8:dim(mbx_anno_df)[2]] #523 234

# transform data
t_mbx_anno_df <- as.data.frame(t(mbx_anno_df)) #234 (including class and subclass) 523

# sum according to Class
sum_class <- aggregate(mbx_anno_df[,1:232], list(mbx_anno_df$Class), FUN=sum)
rownames(sum_class) <- sum_class[,1]
sum_class <- sum_class[,-1]

sum_subclass <- aggregate(mbx_anno_df[,1:232], list(mbx_anno_df$Sub_Class), FUN=sum)
rownames(sum_subclass) <- sum_subclass[,1]
sum_subclass <- sum_subclass[-1,-1]

t_sum_class <- as.data.frame(t(sum_class))
t_sum_subclass <- as.data.frame(t(sum_subclass))



######################################################
# count number of metabolites in each class/subclass #
######################################################
number_class <- mbx_anno_df[,1:232] %>% group_by(mbx_anno_df$Class) %>%tally()

number_subclass <- mbx_anno_df[,1:232] %>% group_by(mbx_anno_df$Sub_Class) %>%tally()


############################################
# save for MaAsLin2
class_maaslin <- t_sum_class
class_maaslin$sample <- rownames(class_maaslin)
class_maaslin <- class_maaslin %>%
  select(sample, everything())
write.table(class_maaslin,file="./mbx/data/class_maaslin.tsv", row.names=FALSE, sep="\t")

subclass_maaslin <- t_sum_subclass
subclass_maaslin$sample <- rownames(subclass_maaslin)
subclass_maaslin <- subclass_maaslin %>%
  select(sample, everything())
write.table(subclass_maaslin,file="./mbx/data/subclass_maaslin.tsv", row.names=FALSE, sep="\t")

###########################################

# do log transformation on the sum of metabolites
log_sum_class <- log2(t_sum_class)
log_sum_subclass <- log2(t_sum_subclass)

# change -Inf to 0
log_sum_class[log_sum_class== "-Inf"] <- 0
log_sum_subclass[log_sum_subclass== "-Inf"] <- 0

# read in metadata
meta <- read.table('./mbx/data/mbx_maaslin_metadata.tsv', header=TRUE, sep='\t')
rownames(meta) <- meta$sample
meta <- meta[,-1] # 232 13

# merge metabolites with metadata
all_class <- merge(log_sum_class,meta,by="row.names")
all_subclass <- merge(log_sum_subclass,meta,by="row.names")

rownames(all_class) <- all_class[,1]
all_class <- all_class[,-1]

rownames(all_subclass) <- all_subclass[,1]
all_subclass <- all_subclass[,-1]

##############################################################################
##class## 46
##linear regression - metabolite classes significantly associated with caco
##############################################################################
#fit <-  lm(formula = all[,1] ~ caco +race8905 + pmh + smk19 + act17m + calor15n + alco15n + aofib15a + ahei2010_noETOH15 + bristolcat + abx + age + bmi19, data=all)
for(i in 1:46){
  fit <-  lm(formula = all_class[,i] ~ caco, data=all_class) #+ age + race8905 + bristolcat + abx
  out <- summary(fit)
  x <- out$coefficients[row.names(out$coefficients)=='cacoDivert',]
  if (i == 1){
    coef <- as.data.frame(t(x))
  } else {
    coef_c = as.data.frame(t(x))
    coef <- rbind(coef,coef_c)
  }
}

coef_class <- coef

rownames(coef_class) <- colnames(all_class[,1:46])

sig_class <- coef_class[coef_class$`Pr(>|t|)`<0.05,]
sig_class <- sig_class[order(-sig_class$Estimate),]

sig_class$UCI <- sig_class$Estimate+1.96*sig_class$`Std. Error`
sig_class$LCI <- sig_class$Estimate-1.96*sig_class$`Std. Error`

sig_class$Class <- rownames(sig_class)

###########################################################################################
###Subclass## 84
##linear regression - metabolite subclasses significantly associated with caco
##############################################################################


for(i in 1:84){
  fit <-  lm(formula = all_subclass[,i] ~ caco, data=all_subclass)
  out <- summary(fit)
  x <- out$coefficients[row.names(out$coefficients)=='cacoDivert',]
  if (i == 1){
    coef <- as.data.frame(t(x))
  } else {
    coef_c = as.data.frame(t(x))
    coef <- rbind(coef,coef_c)
  }
}

coef_subclass <- coef

rownames(coef_subclass) <- colnames(all_subclass[,1:84])

sig_subclass <- coef_subclass[coef_subclass$`Pr(>|t|)`<0.05,]
sig_subclass <- sig_subclass[order(-sig_subclass$Estimate),]

sig_subclass$UCI <- sig_subclass$Estimate+1.96*sig_subclass$`Std. Error`
sig_subclass$LCI <- sig_subclass$Estimate-1.96*sig_subclass$`Std. Error`

sig_subclass$Sub_class <- rownames(sig_subclass)


################################################################################
###                          barplot means and 95% CIs by subclass           ####
################################################################################
subclass_sig_data <- all_subclass[,colnames(all_subclass) %in% rownames(sig_subclass)]
a <- match(rownames(sig_subclass), colnames(subclass_sig_data))
subclass_sig_data <- subclass_sig_data[,a]

#attach number of metabolites in each subclass
t_number_subclass <- as.data.frame(t(number_subclass))
colnames(t_number_subclass) <- t_number_subclass[1,]
t_number_subclass <- t_number_subclass[-1,]
t_number_subclass <- t_number_subclass[,colnames(t_number_subclass) %in% rownames(sig_subclass)]


#attach p value
pvalue <- as.data.frame(t(sig_subclass$`Pr(>|t|)`))
rownames(pvalue) <- 'pvalue'
colnames(pvalue) <- rownames(sig_subclass)


test <- rbind(subclass_sig_data,t_number_subclass,pvalue)


t_test <- as.data.frame(t(test))
t_test$stars <- cut(as.numeric(as.character(t_test$pvalue)), breaks=c(-Inf, 0.001, 0.01, 0.05), label=c("***", "**", "*"))


rownames(t_test) <- paste(rownames(t_test),as.numeric(as.character(t_test$n)),sep="/")
rownames(t_test) <- paste(rownames(t_test),t_test$stars,sep="")

subclass_sig_data <- as.data.frame(t(t_test)) # 235 12
subclass_sig_data <- subclass_sig_data[1:232,]

for(i in 1:ncol(subclass_sig_data)) 
{ 
  subclass_sig_data[ , i] <- as.numeric(as.character(subclass_sig_data[, i])) 
}


subclass_sig_data <- cbind(subclass_sig_data,all_subclass$caco)
colnames(subclass_sig_data)[13] <- 'caco'



################################################################################
###                          barplot means and 95% CIs by subclass (raw data)####
################################################################################
data <- subclass_sig_data
data$subject <- rownames(data)
data_long <- melt(setDT(data), id.vars = c("subject","caco"))
level_y_order <- factor(data_long$variable, level = rev(unique(data_long$variable)))

p1 <- ggplot(data_long,aes(x = value, y = level_y_order,fill=caco))+
  geom_point(size=3,aes(color=caco),position = position_dodge(width = 0.5), alpha=0.75)+
  geom_boxplot(alpha=0.75)+
  xlab("Mean (95% CI)") +
  ylab("Metabolite subclass") +
  scale_fill_manual(values=c("#E69F00", "#4169E1"))+
  scale_color_manual(values=c("#E69F00", "#4169E1"))+
theme_classic() +
theme(plot.title = element_text(size = 30),
      legend.title = element_text(size=30),
      legend.text = element_text(size=30),
      axis.title=element_text(size=50,face="bold"),
      axis.text=element_text(size=50),
      axis.ticks=element_blank(),
      axis.line = element_line(colour = 'black', linewidth=0.75, linetype='solid')
      )
print(p1)
ggsave('./mbx/results/mbx_subclass_boxplot.jpg', plot=p1, width = 30, height = 12, dpi = 300)


################################################################################
mean_subclass <- aggregate(subclass_sig_data[,1:12], list(subclass_sig_data$caco), FUN=mean)
rownames(mean_subclass) <- mean_subclass[,1]
mean_subclass <- mean_subclass[,-1]
t_mean_subclass <- as.data.frame(t(mean_subclass))
t_mean_subclass$metabolite <- rownames(t_mean_subclass)
# change to long format

long_mean <- melt(setDT(t_mean_subclass), id.vars = c("metabolite"), variable.name = "caco")
colnames(long_mean)[3] <- 'mean'


st.err <- function(x) {
    sd(x)/sqrt(length(x))
     }
se_subclass <- aggregate(subclass_sig_data[,1:12], list(subclass_sig_data$caco), st.err)
rownames(se_subclass) <- se_subclass[,1]
se_subclass <- se_subclass[,-1]
t_se_subclass <- as.data.frame(t(se_subclass))
t_se_subclass$metabolite <- rownames(t_se_subclass)
long_se <- melt(setDT(t_se_subclass), id.vars = c("metabolite"), variable.name = "caco")
colnames(long_se)[3] <- 'se'

long_mean_se <- cbind(long_mean,long_se[,3])


long_mean_se$UCI <- long_mean_se$mean+1.96*long_mean_se$se
long_mean_se$LCI <- long_mean_se$mean-1.96*long_mean_se$se

#write.table(long_mean_se,file="./mbx/data/subclass_mean_se.tsv", row.names=FALSE, sep="\t")
#long_mean_se <- read.table('./mbx/data/subclass_mean_se.tsv', header=TRUE, sep='\t')

level_y_order <- factor(long_mean_se$metabolite, level = rev(unique(long_mean_se$metabolite)))

plot <- ggplot(long_mean_se, aes(x=mean,y=level_y_order,fill = caco),group=caco) +
  geom_point(aes(x=mean, y=level_y_order,color=factor(caco)), size = 10, position = position_dodge(width = 0.5))+
  geom_errorbar(aes(xmin=LCI, xmax=UCI), width=.5, position=position_dodge(width = 0.5),size=1)+
  xlab("Mean (95% CI)") +
  ylab("Metabolite subclass") +
  scale_color_manual(values=c("#E69F00", "#4169E1"))+
  theme(legend.position="right",
        legend.text = element_text(size=20), legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        axis.title = element_text(size=50,face = 'bold'),
        axis.text = element_text(size=50,face = 'bold'),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
  )
print(plot)
filepath <- paste('./mbx/results/', 'mbx_subclass_mean.jpg', sep ="" )
ggsave(filename=filepath, plot=plot, width = 30, height = 12, dpi = 300)




################################################################################
###                          Visualization of beta coefficient by class    ####
################################################################################
level_y_order <- factor(sig_class$Class, level = rev(unique(sig_class$Class)))

plot <- ggplot(sig_class, aes(x=Estimate,y=level_y_order,fill = Estimate > 0)) +
  geom_bar(stat = "identity", width=0.5, position = position_dodge())+
  geom_errorbar(aes(xmin=LCI, xmax=UCI), width=0.1, position=position_dodge(width = 0.5),size=0.75)+
  geom_vline(xintercept=0, size=1)+
  xlab("Regression coefficient") +
  ylab("Metabolite class") +
  scale_fill_manual(values = c("TRUE" = "firebrick",
                               "FALSE" = "royalblue")) +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        axis.title = element_text(size=50,face = 'bold'),
        axis.text = element_text(size=50,face = 'bold'),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
  )
print(plot)
filepath <- paste('./mbx/results/', 'mbx_class.jpg', sep ="" )
ggsave(filename=filepath, plot=plot, width = 30, height = 10, dpi = 300)

################################################################################
###                          Visualization of beta coefficient by subclass ####
################################################################################
level_y_order <- factor(sig_subclass$Sub_class, level = rev(unique(sig_subclass$Sub_class)))

plot <- ggplot(sig_subclass, aes(x=Estimate,y=level_y_order,fill = Estimate > 0)) +
  geom_bar(stat = "identity", width=0.5, position = position_dodge())+
  geom_errorbar(aes(xmin=LCI, xmax=UCI), width=0.1, position=position_dodge(width = 0.5),size=0.75)+
  geom_vline(xintercept=0, size=1)+
  xlab("Regression coefficient") +
  ylab("Metabolite subclass") +
  scale_fill_manual(values = c("TRUE" = "firebrick",
                               "FALSE" = "royalblue")) +
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        axis.title = element_text(size=50,face = 'bold'),
        axis.text = element_text(size=50,face = 'bold'),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
  )
print(plot)
filepath <- paste('./mbx/results/', 'mbx_subclass.jpg', sep ="" )
ggsave(filename=filepath, plot=plot, width = 30, height = 12, dpi = 300)


################################################################################
###                          barplot means and 95% CIs by subclass           ####
################################################################################
subclass_sig_data <- all_subclass[,colnames(all_subclass) %in% rownames(sig_subclass)]
a <- match(rownames(sig_subclass), colnames(subclass_sig_data))
subclass_sig_data <- subclass_sig_data[,a]

#attach number of metabolites in each subclass
t_number_subclass <- as.data.frame(t(number_subclass))
colnames(t_number_subclass) <- t_number_subclass[1,]
t_number_subclass <- t_number_subclass[-1,]
t_number_subclass <- t_number_subclass[,colnames(t_number_subclass) %in% rownames(sig_subclass)]

test <- rbind(subclass_sig_data,t_number_subclass)
t_test <- as.data.frame(t(test))
rownames(t_test) <- paste(rownames(t_test),t_test$n,sep="/")



subclass_sig_data <- as.data.frame(t(t_test))



subclass_sig_data <- cbind(subclass_sig_data,all_subclass$caco)
colnames(subclass_sig_data)[13] <- 'caco'



mean_subclass <- aggregate(subclass_sig_data[,1:12], list(subclass_sig_data$caco), FUN=mean)
rownames(mean_subclass) <- mean_subclass[,1]
mean_subclass <- mean_subclass[,-1]
t_mean_subclass <- as.data.frame(t(mean_subclass))
t_mean_subclass$metabolite <- rownames(t_mean_subclass)
# change to long format

long_mean <- melt(setDT(t_mean_subclass), id.vars = c("metabolite"), variable.name = "caco")
colnames(long_mean)[3] <- 'mean'


st.err <- function(x) {
    sd(x)/sqrt(length(x))
     }
se_subclass <- aggregate(subclass_sig_data[,1:12], list(subclass_sig_data$caco), st.err)
rownames(se_subclass) <- se_subclass[,1]
se_subclass <- se_subclass[,-1]
t_se_subclass <- as.data.frame(t(se_subclass))
t_se_subclass$metabolite <- rownames(t_se_subclass)
long_se <- melt(setDT(t_se_subclass), id.vars = c("metabolite"), variable.name = "caco")
colnames(long_se)[3] <- 'se'

long_mean_se <- cbind(long_mean,long_se[,3])


long_mean_se$UCI <- long_mean_se$mean+1.96*long_mean_se$se
long_mean_se$LCI <- long_mean_se$mean-1.96*long_mean_se$se

#write.table(long_mean_se,file="./mbx/data/subclass_mean_se.tsv", row.names=FALSE, sep="\t")
#long_mean_se <- read.table('./mbx/data/subclass_mean_se.tsv', header=TRUE, sep='\t')

level_y_order <- factor(long_mean_se$metabolite, level = rev(unique(long_mean_se$metabolite)))

plot <- ggplot(long_mean_se, aes(x=mean,y=level_y_order,fill = caco),group=caco) +
  geom_point(aes(x=mean, y=level_y_order,color=factor(caco)), size = 10, position = position_dodge(width = 0.5))+
  geom_errorbar(aes(xmin=LCI, xmax=UCI), width=.5, position=position_dodge(width = 0.5),size=1)+
  xlab("Mean (95% CI)") +
  ylab("Metabolite subclass") +
  scale_color_manual(values=c("#E69F00", "#4169E1"))+
  theme(legend.position="right",
        legend.text = element_text(size=20), legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(1, "lines"),
        axis.title = element_text(size=50,face = 'bold'),
        axis.text = element_text(size=50,face = 'bold'),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
  )
print(plot)
filepath <- paste('./mbx/results/', 'mbx_subclass_mean.jpg', sep ="" )
ggsave(filename=filepath, plot=plot, width = 30, height = 12, dpi = 300)
