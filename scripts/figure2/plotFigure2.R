#written by Noah Friedman

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(plyr)
library(data.table); setDTthreads(6)
library(stringr)
library(ggrepel)

#adjust this as needed
plottingDataPath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure2/FIGURE2_PLOTTING_FILES/plotDataFiles/'
plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure2/FIGURE2_PLOTTING_FILES/figurePdfs/'

plot_figure_2a <- function(df){
  p <- ggplot(df, aes(x=burdenType, fill=mutType, y=frac))+
    geom_bar(stat = "identity", position="fill")+
    xlab('mutation burden')+
    ylab('fraction of drivers')+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    scale_fill_viridis_d()+
    scale_x_discrete(limits = rev(levels(df$burdenType)))+
    xlab('tumor type')
  return(p)
}

figure2aDataFrame <- read.table(paste(plottingDataPath, 'figure_2a.tsv', sep=''), sep='\t', header=TRUE)
plt2a <- plot_figure_2a(figure2aDataFrame)
saveFilePath = paste(plottingFilePath, 'figure2a.pdf')
ggsave(saveFilePath, plot=plt2a,  width = 3, height = 4)

#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2B

plot_figure_2b <- function(df){
  p <- ggplot()+
    geom_smooth(data = df[df$hypermutationStatus == 'hypermutated',], aes(x=TMB, y=fracDriverRelated, colour='cancer type related genes'), method='loess', span=1)+
    geom_smooth(data = df[df$hypermutationStatus == 'hypermutated',], aes(x=TMB, y=fracDriverUnrelated, colour='cancer type unrelated genes'), method='loess', span=1)+
    ylab('N putative drivers/N all mutations')+
    scale_color_manual(values=c('Purple', 'Orange'))+
    ylim(0,1)+
    theme_classic()+
    ggtitle('All genes')
  return(p)
}

figure2bDataFrame <- read.table(paste(plottingDataPath, 'figure_2b.tsv', sep=''), sep='\t', header=TRUE)
plt2b <- plot_figure_2b(figure2bDataFrame)
saveFilePath = paste(plottingFilePath, 'figure2b.pdf')
ggsave(saveFilePath, plot=plt2b,  width = 5, height = 4)

#
#######
#############
#####################
##################################
#####################
#############
########
#
#FIGURE 2C

plot_figure_2c <- function(df){
  p <- ggplot(df, aes(x=perCase_c1, y=perCase_c2, color=pathway))+
    geom_point(alpha=0.5)+
    geom_point(data=df[df$qVal < .05,], pch=21, colour='black')+ #outlines
    geom_text_repel(data=df[(df$perCase_c1 > .2) | (df$perCase_c2 > .2), ], aes(label=Allele), size=2)+
    xlim(0, .75)+
    ylim(0, .75)+
    theme_classic()+
    geom_segment(aes(x=0,y=0,xend=0.5, yend=0.5), colour='black', linetype='dashed')+
    scale_color_manual(values = c('gray', 'blue', 'orange'))+
    ylab('Fraction of Colorectal cases mutated')+
    xlab('Fraction of Endometrial cases mutated')+
    ggtitle('MSI Indels')
  return(p)
}

figure2cDataFrame <- read.table(paste(plottingDataPath, 'figure_2c.tsv', sep=''), sep='\t', header=TRUE)
plt2c <- plot_figure_2c(figure2cDataFrame[(figure2cDataFrame$c1_cancerType == 'Endometrial Cancer') & (figure2cDataFrame$c2_cancerType == 'Colorectal Cancer'), ])
saveFilePath = paste(plottingFilePath, 'figure2c.pdf')
ggsave(saveFilePath, plot=plt2c,  width = 5, height = 5)

#Correlation scores
cor(df[df$GeneType == 'oncogene',]$perCaseColorectal, df[df$GeneType == 'oncogene',]$perCaseEndometrial, method='pearson')
cor(df[df$GeneType == 'tsg',]$perCaseColorectal, df[df$GeneType == 'tsg',]$perCaseEndometrial, method='pearson')

#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2D

make_figure_2d <- function(df){
  p <- ggplot(df, aes(x=factor(label, levels=c('B2M', 'PTEN', 'TP53', 'ARID1A', 'APC',
                                               'related_tsg', 'other_tsg', 'PIK3CA', 'TERT', 'related_oncogene', 'other_oncogene', '1 or 2 VUS', '1 or 2 silent')),
                      y=isTrans))+
    stat_summary()+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    xlab('Gene Type')+
    ylab('Fraction phasable mutations in trans')
  return(p)
}

figure2dDataFrame <- read.table(paste(plottingDataPath, 'figure_2d.tsv', sep=''), sep='\t', header=TRUE)
plt2d <- make_figure_2d(figure2dDataFrame)
saveFilePath = paste(plottingFilePath, 'figure2d.pdf')
ggsave(saveFilePath, plot=plt2d,  width = 5, height = 5)

#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2E

plot_figure_2e <- function(df){
  
  fun.1 <- function(x) x^2 
  
  p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))+
    geom_point(data=df, aes(x=nAllele/nGene, y=nDoubleHit/nBiallelicLoss, size=nDoubleHit, colour=geneType))+
    geom_text_repel(data=df[df$nDoubleHit >= 4,], aes(x=nAllele/nGene, y=nDoubleHit/nBiallelicLoss, label=allele), force=5)+
    xlab('Fraction of all drivers in gene\ncaused by double hit allele')+
    ylab('Fraction of all biallelic inactivation of gene\nattributable to double hit')+
    ylim(0,1)+
    xlim(0,1)+
    stat_function(fun = fun.1, linetype='dashed')+
    theme_classic()
  return(p)
}

figure2eDataFrame <- read.table(paste(plottingDataPath, 'figure_2e.tsv', sep=''), sep='\t', header=TRUE)
plt2e <- plot_figure_2e(figure2eDataFrame)
saveFilePath = paste(plottingFilePath, 'figure2e.pdf')
ggsave(saveFilePath, plot=plt2e,  width = 5, height = 5)



#
#####
#################
##############################
#################
######
#
#FINAL PLOT

#First pad plots

#2a
RightPadRatio <- 1
LowerPadRatio <- 0
padded2a <- plot_grid(plot_grid(plt2a, ggplot(), ncol=2, rel_widths = c(1, RightPadRatio)),
                      ggplot(), nrow=2, rel_heights = c(1, LowerPadRatio))
#2b
f1RightPadRatio <- 0
f1LowerPadRatio <- .25
padded2b <- plot_grid(plot_grid(plt2b, ggplot(), ncol=2, rel_widths = c(1, f1RightPadRatio)),
                      ggplot(), nrow=2, rel_heights = c(1, f1LowerPadRatio))

#2c
f1RightPadRatio <- 0
f1LowerPadRatio <- 0
padded2c <- plot_grid(plot_grid(plt2c, ggplot(), ncol=2, rel_widths = c(1, f1RightPadRatio)),
                      ggplot(), nrow=2, rel_heights = c(1, f1LowerPadRatio))

#2d
f1RightPadRatio <- 0
f1LowerPadRatio <- 0
padded2d <- plot_grid(plot_grid(plt2d, ggplot(), ncol=2, rel_widths = c(1, f1RightPadRatio)),
                      ggplot(), nrow=2, rel_heights = c(1, f1LowerPadRatio))

#2e
f1RightPadRatio <- 0
f1LowerPadRatio <- 0
padded2e <- plot_grid(plot_grid(plt2e, ggplot(), ncol=2, rel_widths = c(1, f1RightPadRatio)),
                      ggplot(), nrow=2, rel_heights = c(1, f1LowerPadRatio))


combinedPdfPlot <- plot_grid(padded2a, padded2b, 
                             padded2c, padded2d, padded2e,
                             ggplot(), ggplot(), ggplot(), ggplot(), #include extra plots for padding purposes
                             nrow=3, ncol=3)
saveFilePath = paste(plottingFilePath, 'figure2.pdf')
ggsave(saveFilePath,
       plot=combinedPdfPlot,  width = 15, height = 15, units = c("in"), limitsize = FALSE)



