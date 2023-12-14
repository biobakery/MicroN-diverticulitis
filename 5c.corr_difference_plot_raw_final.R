
##################################################################################################
# Purpose: Visualize differences in correlations between species and metabolites between diverticulitis and controls
# Last update date: 10/16/2023
# Note: I first got correlations between each metabolite and bug pair in divert and controls using HALLA and then use this script to visualize differences in the correlations
##################################################################################################


setwd("/Users/wm897/Dropbox (Partners HealthCare)/MGH/AC/AC/diverticulitis/NHS2/MicroN/")

library(readxl)
library(dplyr)
library(stringr)
library(lmerTest)
library(ggplot2)
library(tidyverse)


# read in all_associations.txt got from running HALLA
r_divert <- read.table('./mbx/data/halla/raw/divert/synthetic_output/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
r_control <- read.table('./mbx/data/halla/raw/control/synthetic_output/all_associations.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")

r_diff <- r_divert %>% right_join(r_control, by=c("X_features","Y_features"))
r_diff$diff <- r_diff$association.x-r_diff$association.y

##################################################################################################
#######reduce dimensionality by removing bugs or metabolites that do not have any big values (>0.45)
##################################################################################################
test <- r_diff
test$abs_diff <- abs(test$diff)


test$X_features_factor <- as.factor(test$X_features)
test$Y_features_factor <- as.factor(test$Y_features)


rows <- c()
row_names <- c()
for (i in 1:length(levels(test$X_features_factor))){
  include = 0
  name = NULL
  for (j in 1:length(levels(test$Y_features_factor))){
    if (test[(i-1)*505 + j,]$abs_diff > 0.45) {
      include <- 1
      name <- test[(i-1)*505 + j,]$X_features
    }
  }
  print(name)
  if (include == 1) {
    rows <- append(rows, i)
    row_names <- append(row_names, name)
  }
}


cols = c()
col_names =c()
for (j in 1:length(levels(test$Y_features_factor))){
  include = 0
  name = NULL
  for (i in 1:length(levels(test$X_features_factor))){
    if (test[(i-1)*505 + j,]$abs_diff > 0.45) {
      include = 1
      name <- test[(i-1)*505 + j,]$Y_features
    }
  }
  if (include == 1) {
    cols = append(cols, j)
    col_names <- append(col_names, name)
  }
}


test_reduce <- test[(test$X_features %in% row_names),]
test_reduce2 <- test_reduce[(test_reduce$Y_features %in% col_names),]

test_reduce2$stars <- cut(test_reduce2$abs_diff, breaks=c(0,0.1,0.3,0.5,1), label=c("","*","**","***"))

# check the extremes
#check <- test_reduce2[order(test_reduce2$abs_diff),]
#write.table(check,"./mbx/data/halla/raw/association_diff_filter_sort.tsv", row.names = FALSE,sep="\t")

#test_reduce2 <- read.table('./mbx/data/halla/raw/association_diff_filter_sort.tsv', header=TRUE, sep='\t')

test_reduce2 <- test_reduce2 %>%
  mutate(group = case_when(
    association.x > 0 & association.y > 0 ~ '+,+',
    association.x < 0 & association.y < 0 ~ '-,-',
    association.x > 0 & association.y < 0 ~ '+,-',
    association.x < 0 & association.y > 0 ~ '-,+'))
test_reduce2$group<- factor(test_reduce2$group,
                        levels = c('+,+', '-,-','+,-','-,+'))

################################################################################
###class###
################################################################################
class <- read.table('./mbx/data/name_class.tsv',header=TRUE, sep='\t')
colnames(class)[1] <- 'Y_features'
class <- class[!duplicated(class$Y_features), ]

# remove metabolites that do not have Class
class <- class[!class$Class=="",]

test <- test_reduce2
test2 <- merge(test,class,by="Y_features")



########################################################################################################
################################################################################
###clustering by diff in correlations###
################################################################################

a <- test2[,colnames(test2) %in% c("X_features","Y_features","diff")]

b <- a |> 
  pivot_wider(names_from = "Y_features", values_from = "diff") |>
  as.data.frame() |>
  column_to_rownames("X_features") |>
  as.matrix()

c <- b |> dist() |> hclust()
order_species <- c$labels[c$order]


d <- t(b)

e <- d |> dist() |> hclust()
order_mbx <- e$labels[e$order]


test2$X_features_factor<-factor(test2$X_features)
test2$Y_features_factor<-factor(test2$Y_features)
all  <- test2 
all$X_features_factor <- factor(all$X_features_factor, levels=order_species)
all$Y_features_factor <- factor(all$Y_features_factor, levels=order_mbx)
all$abs_diff <- abs(all$diff)




####add annotation by class#####
mbx <- class[class$Y_features %in% order_mbx,]
mbx <- mbx %>% arrange(factor(mbx$Y_features, levels = order_mbx))
mbx$Y_features_factor <- factor(mbx$Y_features, levels=order_mbx)

h1 <- ggplot(mbx)+
  geom_bar(mapping = aes(x=1, y = Y_features_factor, fill = Class), 
           stat = "identity")+
  theme_classic(base_size = 50) +
  theme(
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 50),
        legend.text = element_text(size = 50),
        panel.spacing.y = unit(1.5, "lines"),
        axis.text.y = element_text(size=50),
        axis.text.x = element_blank(),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank()
  )+
  guides(colour = guide_legend(override.aes = list(size=30))) 
print(h1)

h1_legend <- as_ggplot(get_legend(h1))
print(h1_legend)

h2 <- all |>
  as_tibble() |>
  ggplot(aes(X_features_factor, Y_features_factor)) +
  geom_point(aes( color = factor(group), size = abs_diff))+
  scale_size_continuous(range=c(2,30))+
  labs(size = "Absolute difference", color = "Group") +
  theme_classic(base_size = 50) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=50),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.position = "top",
        panel.spacing.y = unit(1.5, "lines"),
        axis.text.y = element_blank(),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', linewidth=0.75, linetype='solid')
  )+
  guides(colour = guide_legend(override.aes = list(size=30))) 
print(h2)

h2_legend <- as_ggplot(get_legend(h2))
print(h2_legend)

library(cowplot)


h1 <- h1 + theme(legend.position = "none")
h2 <- h2 + theme(legend.position = "none")

filepath <- paste('./mbx/data/halla/raw/', 'h1.jpg', sep ="" )
ggsave(filename=filepath, plot=h1, width = 30, height = 48, dpi = 300)

filepath <- paste('./mbx/data/halla/raw/', 'h2.jpg', sep ="" )
ggsave(filename=filepath, plot=h2, width = 30, height = 48, dpi = 300)


# combine the two plots
plot_combine <- plot_grid(h1, h2, align = "h",axis = "l", rel_widths = c(1, 4))
print(plot_combine)

h3 <- plot_grid(plot_combine, h1_legend, align = "h",axis = "l",rel_widths = c(6, 1),rel_heights = c(4,1))
print(h3)

h4 <- plot_grid(h2_legend, h3, align = "v",axis = "b", ncol= 1,rel_heights = c(1,12))
print(h4)


filepath <- paste('./mbx/data/halla/raw/', 'r_diff_plot_gt_0.45_raw_by_diff.jpg', sep ="" )
ggsave(filename=filepath, plot=h4, width = 78, height = 60, dpi = 300,limitsize = FALSE)






########################################################################################################
################################################################################
###clustering by group###
################################################################################

a <- test_reduce2[,colnames(test_reduce2) %in% c("X_features","Y_features","group")]

a <- a %>%
  mutate(group = case_when(
    group == '+,+' ~ 1,
    group == '-,-' ~ 2,
    group == '+,-' ~ 3,
    group == '-,+' ~ 4))
a$group<- factor(a$group,levels = c('1', '2','3','4'))

b <- a |> 
  pivot_wider(names_from = "Y_features", values_from = "group") |>
  as.data.frame() |>
  column_to_rownames("X_features") |>
  as.matrix()

c <- b |> dist() |> hclust()
order_species <- c$labels[c$order]


d <- t(b)

e <- d |> dist() |> hclust()
order_mbx <- e$labels[e$order]


test_reduce2$X_features_factor<-factor(test_reduce2$X_features)
test_reduce2$Y_features_factor<-factor(test_reduce2$Y_features)
all  <- test_reduce2 
all$X_features_factor <- factor(all$X_features_factor, levels=order_species)
all$Y_features_factor <- factor(all$Y_features_factor, levels=order_mbx)
all$abs_diff <- abs(all$diff)


plot <- all |>
  as_tibble() |>
  ggplot(aes(X_features_factor, Y_features_factor)) +
  geom_point(aes( color = factor(group), size = abs_diff))+
  scale_size_continuous(range=c(2,20))+
  labs(size = "Absolute difference", color = "Group") +
  theme_classic(base_size = 30) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=30),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30),
        legend.position = "top",
        panel.spacing.y = unit(1.5, "lines"),
        axis.text = element_text(size=30),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )+
  guides(colour = guide_legend(override.aes = list(size=20))) 
print(plot)

filepath <- paste('./mbx/data/halla/raw/', 'r_diff_plot_gt_0.45_raw_by_group.jpg', sep ="" )
ggsave(filename=filepath, plot=plot, width = 30, height = 48, dpi = 300)





