#written by Noah Friedman
#the code to plot all the figures for figure 1
library(ggplot2)
library(grid)
require(cowplot)
library(dplyr)
library(data.table); setDTthreads(6)
library(gdata)
library(ggpubr)

#Pan plot set up info

#adjust this as needed
plottingDataPath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/plotDataFiles/'
plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/figurePdfs/'

#
###
######
###########
#######
###
#
#PLOT FIGURE 1A


plot_figure_1a <- function(df){
  bottomScaleFactor <- 10
  p <- ggplot(df, aes(x=reorder(label, fracHypermutatedOrdering)))+
    geom_bar(aes(y=fracHypermutated), stat='identity')+
    
    theme(axis.text.x = element_text(angle=90))+
    theme(axis.ticks.x = element_blank())+
    scale_fill_manual(values=c('#858585', 'light gray'))+
    theme_classic()+
    xlab('Cancer Type')+
    ylab('Fraction of cases hypermutated')+
    coord_flip()+
    ggtitle('1a.')
  return(p)
}

figure1aDataFrame <- read.table(paste(plottingDataPath, 'figure_1a.tsv', sep=''), sep='\t', header=TRUE)
plt1a <- plot_figure_1a(figure1aDataFrame)
saveFilePath = paste(plottingFilePath, 'figure1a.pdf')
ggsave(saveFilePath, plot=plt1a,  width = 6, height = 4)

#
###
######
###########
#######
###
#
#PLOT FIGURE 1B

plot_figure_1b <- function(df){
  p <- ggplot(df, aes(x=1, fill=signature, y=nHyperHigh/sum(df$nHyperHigh)))+
    geom_bar(stat='identity')+
    theme_classic()+
    scale_fill_manual(values=c("#FF0000","#267574",
                               'gray', "#ADFF2F","#00dd5d", '#ffb347', '#2A52BE', "#FFF600"))+
    xlab('')+
    theme(axis.text.x = element_blank())+
    ylab('fraction of all hypermutated cases')+
    ggtitle('1b.')
  return(p)
}

figure1bDataFrame <- read.table(paste(plottingDataPath, 'figure_1b.tsv', sep=''), sep='\t', header=TRUE)
plt1b <- plot_figure_1b(figure1bDataFrame)
saveFilePath = paste(plottingFilePath, 'figure1b.pdf')
ggsave(saveFilePath,
       plot=plt1b,  width = 3, height = 5)

#
###
######
###########
#######
###
#
#PLOT FIGURE 1C

plot_figure_1c <- function(df){
  plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nOncMuts))+
    geom_violin(bw = 1, aes(fill=factor(cohort,
                                  levels = c('normal_Colorectal', 'hyper_Colorectal', 'normal_Endometrial', 'hyper_Endometrial',
                                  'normal_Glioma', 'hyper_Glioma', 'normal_Other', 'hyper_Other'))))+
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    ylab('N Driver Mutations')+
    xlab('Cancer Type')+
    coord_cartesian(ylim=c(0,50))+
    scale_fill_manual(values =c('orange', '#b36200', 'lavender', '#301934', '#add8e6', 'blue', 'gray', '#333333'), name='Cohort')+
    ggtitle('1c.')
  
  return(plt)
}

figure1cDataFrame <- read.table(paste(plottingDataPath, 'figure_1c.tsv', sep=''), sep='\t', header=TRUE)
plt1c <- plot_figure_1c(figure1cDataFrame)
saveFilePath = paste(plottingFilePath, 'figure1c.pdf')
ggsave(saveFilePath,
       plot=plt1c,  width = 6, height = 4, units = c("in"), limitsize = FALSE)


#P values: (change to $nHotspots if desired)
pEndometrial <- t.test(figure1cDataFrame[figure1cDataFrame$cohort == 'hyper_Endometrial',]$nOncMuts,
                       figure1cDataFrame[figure1cDataFrame$cohort == 'normal_Endometrial',]$nOncMuts)$p.value
pColorectal <- t.test(figure1cDataFrame[figure1cDataFrame$cohort == 'hyper_Colorectal',]$nOncMuts,
                      figure1cDataFrame[figure1cDataFrame$cohort == 'normal_Colorectal',]$nOncMuts)$p.value
pGlioma <- t.test(figure1cDataFrame[figure1cDataFrame$cohort == 'hyper_Glioma',]$nOncMuts,
                  figure1cDataFrame[figure1cDataFrame$cohort == 'normal_Glioma',]$nOncMuts)$p.value
pOther <- t.test(figure1cDataFrame[figure1cDataFrame$cohort == 'hyper_Other',]$nOncMuts,
                 figure1cDataFrame[figure1cDataFrame$cohort == 'normal_Other',]$nOncMuts)$p.value
print(paste('p values: ', 'endometrial:', pEndometrial, 'colorectal:', pColorectal, 'glioma:', pGlioma, 'other:', pOther))

#
####
#########
##############
#########
####
#

#plot figure 1d

plot_figure_1d <- function(df){
  pTmb <- ggplot(df, aes(x=expectedOncogenicSNP))+
    geom_histogram()+emptyTheme+ylab('n cases')
  p <- ggplot(df)+
    stat_summary_bin(bins=10, aes(x=expectedOncogenicSNP, y = obsOncogenicSNP, colour = 'All SNV drivers'))+
    stat_summary_bin(bins=10, aes(x=expectedTruncatingOncogene, y = obsStopGainOncogene, colour = 'Stop-gain oncogene'))+
    stat_summary_bin(bins=10, aes(x=expectedTruncatingTSG, y = obsStopGainTSG, colour = 'Stop-gain TSG'))+
    stat_summary_bin(bins=10, aes(x=expectedHotspot, y = obsHotspot, colour = 'Hotspot drivers'))+
    geom_segment(aes(x=30, y=30, yend=0, xend=0), linetype='dashed')+
    xlab('Expected')+
    ylab('Observed')+
    theme_classic()+
    ggtitle('1d.')+
    coord_fixed()
  leg <- get_legend(p)
  p <- p + theme(legend.position = 'none')
  alignedPlot <- plot_grid(p, pTmb, nrow=2, rel_heights=c(1,.25))
  alignedPlotWithLegend <- plot_grid(alignedPlot, leg, ncol=2, rel_widths = c(1,.5))
  return(alignedPlotWithLegend)
}

figure1dDataFrame <- read.table(paste(plottingDataPath, 'figure_1d.tsv', sep=''), sep = '\t', header=TRUE)
plt1d <- plot_figure_1d(figure1dDataFrame)
saveFilePath = paste(plottingFilePath, 'figure1d.pdf')
ggsave(saveFilePath,
       plot=plt1d,  width = 3, height = 4, units = c("in"), limitsize = FALSE)

#
####
#########
##############
#########
####
#


#note this omits some outliers
plot_figure_1e <- function(df){
  pTmb <- ggplot(df, aes(x=nIndels/30))+
    geom_histogram()+emptyTheme+ylab('n cases')+xlab('Indels/MB')+
    xlim(0,50)
    #scale_x_log10()
  p <- ggplot()+
    
    stat_summary_bin(data=df[(df$dominantSignature == 'MMR'),], bins=20, aes(x=nIndels/30, y=(OncogeneObs/.782), colour='Oncogene_obs'), alpha=0.5)+
    stat_summary_bin(data=df[(df$dominantSignature == 'MMR'),], bins=20, aes(x=nIndels/30, y=(TSGObs/.501), colour='TSG_obs'), alpha=0.5)+
    geom_segment(data=df[(df$dominantSignature == 'MMR'),], aes(x=50, y=50, yend=0, xend=0), linetype='dashed')+
    geom_histogram(data=df[(df$dominantSignature == 'MMR'),], aes(x=nIndels/30, y=-..count..))+
    coord_fixed()+
    ylab('Indels/MB')+
    theme_classic()+
    xlab('Indels/MB neutral genes')+
    xlim(0,50)+
    ggtitle('1e')+
    scale_y_continuous(breaks = c(-50, 0, 30, 60), labels=c("n cases=50", "0", "30", "60"))
  leg <- get_legend(p)
  p <- p + theme(legend.position = 'none')
  alignedPlot <- p
  alignedPlotWithLegend <- plot_grid(alignedPlot, leg, ncol=2, rel_widths = c(1,.5))
  return(alignedPlotWithLegend)
}

#plot figure 1e
figure1eDataFrame <- read.table(paste(plottingDataPath, 'figure_1e.tsv', sep=''), sep = '\t', header=TRUE)
plt1e <- plot_figure_1e(figure1eDataFrame[figure1eDataFrame$dominantSignature == 'MMR',])
saveFilePath = paste(plottingFilePath, 'figure1e.pdf')
ggsave(saveFilePath,
       plot=plt1e,  width = 5, height = 4, units = c("in"), limitsize = FALSE)

#
###
######
#############
##################
##########################
#################
############
#######
###

