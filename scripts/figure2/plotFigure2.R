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
library(ggpubr)

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
    xlab('tumor type')+
    ggtitle('2a.')
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

#Todo add a per/Mb value here

plot_figure_2b <- function(df){
  pTmb <- ggplot(df, aes(x=TMB))+
    geom_histogram()+emptyTheme+ylab('n cases')
  p <- ggplot()+
    geom_smooth(data = df, aes(x=TMB, y=(1e6*nRelatedDriver)/relatedGeneSize, colour='cancer type related genes'), method='loess', span=1)+
    geom_smooth(data = df, aes(x=TMB, y=(1e6*nUnrelatedDrivers)/(1.3e6-relatedGeneSize), colour='cancer type unrelated genes'), method='loess', span=1)+
    ylab('Mutations per MB of gene content')+
    scale_color_manual(values=c('Purple', 'Orange'))+
    theme_classic()+
    xlab('case-TMB')+
    ggtitle('2b.')
  leg <- get_legend(p)
  p <- p + theme(legend.position = 'none')
  alignedPlot <- plot_grid(p, pTmb, nrow=2, rel_heights=c(1,.25))
  alignedPlotWithLegend <- plot_grid(alignedPlot, leg, ncol=2, rel_widths = c(1,.5))
  return(alignedPlotWithLegend)
}

figure2bDataFrame <- read.table(paste(plottingDataPath, 'figure_2b.tsv', sep=''), sep='\t', header=TRUE)
plt2b <- plot_figure_2b(figure2bDataFrame)
saveFilePath = paste(plottingFilePath, 'figure2b.pdf')
ggsave(saveFilePath, plot=plt2b,  width = 4.5, height = 4)

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
  
  p1 <- ggplot(df, aes(x=mutBurdenPathway, y=nTruncTsg/nTotal, color=cancerType))+
    geom_boxplot()+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    geom_text_repel(data=figure2cDataFrame[figure2cDataFrame$nTruncTsg/figure2cDataFrame$nTotal > .175,], aes(label=gene))+
    ggtitle('TSGs')+
    xlab('mutation burden & pathway')+
    ylab('Fraction cases with truncating mutations')+
    ylim(0,.75)
  
  leg <- get_legend(p1)
  p1 <- p1 + theme(legend.position = 'none')
  
  p2 <- ggplot(df, aes(x=mutBurdenPathway, y=nTruncOncogene/nTotal, color=cancerType))+
    geom_boxplot()+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    geom_text_repel(data=figure2cDataFrame[figure2cDataFrame$nTruncOncogene/figure2cDataFrame$nTotal > .2,], aes(label=gene))+
    ggtitle('Oncogenes')+
    xlab('mutation burden & pathway')+
    ylab('Fraction cases with truncating mutations')+
    ylim(0, .75)+
    theme(legend.position = 'none')
  
  p3 <- ggplot(df, aes(x=mutBurdenPathway, y=nVus/nTotal, color=cancerType))+
    geom_boxplot()+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    geom_text_repel(data=figure2cDataFrame[figure2cDataFrame$nVus/figure2cDataFrame$nTotal > .2,], aes(label=gene))+
    ggtitle('VUS mutations')+
    xlab('mutation burden & pathway')+
    ylab('Fraction cases with VUS-missense mutations')+
    ylim(0,.75)+ 
    theme(legend.position = 'none')
  
  alignedPlot <- plot_grid(p1, p2, p3, leg, ncol=4, rel_widths = c(1,1,1,.3))
  plotWithTitle <- plot_grid( ggplot()+ggtitle('2c.'), alignedPlot, nrow=2, rel_heights = c(.2,1))
  return(plotWithTitle)
}

figure2cDataFrame <- read.table(paste(plottingDataPath, 'figure_2c.tsv', sep=''), sep='\t', header=TRUE)
plt2c <- plot_figure_2c(figure2cDataFrame)
saveFilePath = paste(plottingFilePath, 'figure2c.pdf')
ggsave(saveFilePath, plot=plt2c,  width = 15, height = 5)



#
######
###############
########################
##################################
########################
###############
######
#

#FIGURE 2D Summarizing composite mutations in hypermutated tumors

plot_figure_2d <- function(df){
  ggplot()+
    stat_summary(data=df[df$geneType == 'oncogene',], aes(x='1', y=isComposite))+
    stat_summary(data=df[df$geneType == 'tsg_missense',], aes(x='2', y=isComposite))+
    stat_summary(data=df[df$geneType == 'tsg_truncating',], aes(x='3', y=isComposite))+
    stat_summary(data=df[df$hypermutationInduced == 'Almost certain',], aes(x='4', y=isComposite))+
    stat_summary(data=df[df$hypermutationInduced == 'Unlikely',], aes(x='5', y=isComposite))+
    
    stat_summary(data=df[(df$related == 'related') & (df$geneType != 'oncogene'),], aes(x='6', y=isComposite))+
    stat_summary(data=df[(df$related == 'related') & (df$geneType == 'oncogene'),], aes(x='7', y=isComposite))+
    stat_summary(data=df[(df$related == 'not-related') & (df$geneType != 'oncogene'),], aes(x='8', y=isComposite))+
    stat_summary(data=df[(df$related == 'not-related') & (df$geneType == 'oncogene'),], aes(x='9', y=isComposite))+
    theme_classic()+
    scale_x_discrete(labels=c('1' = "oncogene", "2" = "TSG missense",
                              "3" = "TSG truncating", "4" = 'hypermutation induced', 
                              '5' = 'not hypermutation induced', '6' = 'related tsg', '7' = 'related oncogene',
                              '8' = 'unrelated tsg', '9' = 'unrelated oncogene'))+
    
    theme(axis.text.x = element_text(angle=90))+
    coord_cartesian(ylim= c(0,.5))+
    ylab('n times driver observed in composite/\nn times driver observed overall')+
    xlab('mutation type')+
    ggtitle('2d.')
}

figure2dDataFrame <- read.table(paste(plottingDataPath, 'figure_2d.tsv', sep=''), sep='\t', header=TRUE)
plt2d <- plot_figure_2d(figure2dDataFrame)
saveFilePath = paste(plottingFilePath, 'figure2d.pdf')
ggsave(saveFilePath, plot=plt2d,  width = 4, height = 4.5)

#
######
###############
########################
##################################
########################
###############
######
#

#FIGURE 2E Phasing

make_figure_2e <- function(df){
  p <- ggplot(df, aes(x=factor(label, levels=c('B2M', 'PTEN', 'TP53', 'ARID1A', 'APC',
                                               'related_tsg', 'other_tsg', 'PIK3CA', 'TERT', 'related_oncogene', 'other_oncogene', '1 or 2 VUS', '1 or 2 silent')),
                      y=isTrans))+
    stat_summary()+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    xlab('Gene Type')+
    ylab('Fraction phasable mutations in trans')+
    ggtitle('2e.')
  return(p)
}

figure2eDataFrame <- read.table(paste(plottingDataPath, 'figure_2e.tsv', sep=''), sep='\t', header=TRUE)
plt2e <- make_figure_2e(figure2eDataFrame)
saveFilePath = paste(plottingFilePath, 'figure2e.pdf')
ggsave(saveFilePath, plot=plt2e,  width = 5, height = 5)



#
#####
#################
##############################
#################
######
#

