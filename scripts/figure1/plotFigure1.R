#written by Noah Friedman
#the code to plot all the figures for figure 1
library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(stringr)
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
    coord_flip()
  return(p)
}

figure1aDataFrame <- read.table(paste(plottingDataPath, 'figure_1a.tsv', sep=''), sep='\t', header=TRUE)
plt1a <- plot_figure_1a(figure1aDataFrame)
saveFilePath = paste(plottingFilePath, 'figure1a.pdf')
ggsave(saveFilePath, plot=p,  width = 6, height = 4)

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
    ylab('fraction of all hypermutated cases')
  return(p)
}

figure1bDataFrame <- read.table(paste(plottingDataPath, 'figure_1b.tsv', sep=''), sep='\t', header=TRUE)
plt1b <- plot_figure_1b(figure1bDataFrame)
saveFilePath = paste(plottingFilePath, 'figure1b.pdf')
ggsave(saveFilePath,
       plot=p,  width = 3, height = 5)

#
###
######
###########
#######
###
#
#PLOT FIGURE 1C

plot_figure_1c <- function(df){
  #plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nHotspots))+
  plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nOncMuts))+
    #geom_boxplot(fatten = NULL, outlier.shape=NA)+
    geom_violin(bw = 1, aes(fill=factor(cohort,
                                  levels = c('normal_Colorectal', 'hyper_Colorectal', 'normal_Endometrial', 'hyper_Endometrial',
                                  'normal_Glioma', 'hyper_Glioma', 'normal_Other', 'hyper_Other'))))+
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    #scale_colour_manual(values =  c('black', "#267574", 'gray', '#ADFF2F', "#9acd32", '#2A52BE'), name="Dominant\nSignature")+
    ylab('N Driver Mutations')+
    #ylab('N hotspot mutations')+
    xlab('Cancer Type')+
    coord_cartesian(ylim=c(0,50))+
    #geom_jitter(aes(colour=factor(cohort,
    #                              levels = c('normal_Colorectal', 'hyper_Colorectal', 'normal_Endometrial', 'hyper_Endometrial',
    #                                         'normal_Glioma', 'hyper_Glioma', 'normal_Other', 'hyper_Other'))),
    #            shape=16, position=position_jitter(0.1), alpha=0.75)+
    
    scale_fill_manual(values =c('orange', '#b36200', 'lavender', '#301934', '#add8e6', 'blue', 'gray', '#333333'), name='Cohort')
  
  return(plt)
}

figure1cDataFrame <- read.table(paste(plottingDataPath, 'figure_1c.tsv', sep=''), sep='\t', header=TRUE)
plt1c <- plot_figure_1c(figure1cDataFrame)
saveFilePath = paste(plottingFilePath, 'figure1c.pdf')
ggsave(saveFilePath,
       plot=plt,  width = 6, height = 4, units = c("in"), limitsize = FALSE)


#P values: (change to $nHotspots if desired)
pEndometrial <- t.test(figure1dDataFrame[figure1dDataFrame$cohort == 'hyper_Endometrial',]$nOncMuts,
       figure1dDataFrame[figure1dDataFrame$cohort == 'normal_Endometrial',]$nOncMuts)$p.value
pColorectal <- t.test(figure1dDataFrame[figure1dDataFrame$cohort == 'hyper_Colorectal',]$nOncMuts,
                       figure1dDataFrame[figure1dDataFrame$cohort == 'normal_Colorectal',]$nOncMuts)$p.value
pGlioma <- t.test(figure1dDataFrame[figure1dDataFrame$cohort == 'hyper_Glioma',]$nOncMuts,
                       figure1dDataFrame[figure1dDataFrame$cohort == 'normal_Glioma',]$nOncMuts)$p.value
pOther <- t.test(figure1dDataFrame[figure1dDataFrame$cohort == 'hyper_Other',]$nOncMuts,
                       figure1dDataFrame[figure1dDataFrame$cohort == 'normal_Other',]$nOncMuts)$p.value
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
  p <- ggplot(df, aes(x=nmut))+
    geom_smooth(aes(y = expectedOncogenicSNP, colour = 'Expected'))+
    geom_smooth(aes(y = obsOncogenicSNP, colour = 'Observed'))+
    scale_colour_manual(values=c('gray', 'black'))+
    ylab('Putative SNP drivers')+
    xlab('nmut in IMPACT-341 genes')+
    theme_classic()
  return(p)
}

figure1dDataFrame <- read.table(paste(plottingDataPath, 'figure_1d.tsv', sep=''), sep = '\t', header=TRUE)
plt1d <- plot_figure_1d(figure1dDataFrame)
saveFilePath = paste(plottingFilePath, 'figure1d.pdf')
ggsave(saveFilePath,
       plot=p,  width = 3, height = 4, units = c("in"), limitsize = FALSE)

#
####
#########
##############
#########
####
#

plot_figure_1e <- function(df){
  p <- ggplot()+
    geom_smooth(data= df[df$dominantSignature == 'Signature.MMR',], aes(x=tmb, y=OncogeneObs - OncogeneExp, linetype='Oncogene', colour='MMR'), span=1)+
    geom_smooth(data= df[df$dominantSignature == 'Signature.MMR',], aes(x=tmb, y=TSGObs - TSGExp, linetype='TSG', colour='MMR'), span=1)+
    
    geom_smooth(data= df[df$dominantSignature == 'Signature.11',], aes(x=tmb, y=OncogeneObs - OncogeneExp, linetype='Oncogene', colour='TMZ'), span=1)+
    geom_smooth(data= df[df$dominantSignature == 'Signature.11',], aes(x=tmb, y=TSGObs - TSGExp, linetype='TSG', colour='TMZ'), span=1)+
    
    geom_smooth(data= df[df$dominantSignature == 'Signature.10',], aes(x=tmb, y=OncogeneObs - OncogeneExp, linetype='Oncogene', colour='POLE'), span=1)+
    geom_smooth(data= df[df$dominantSignature == 'Signature.10',], aes(x=tmb, y=TSGObs - TSGExp, linetype='TSG', colour='POLE'), span=1)+
    
    geom_smooth(data= df[df$dominantSignature == 'Signature.APOBEC',], aes(x=tmb, y=OncogeneObs - OncogeneExp, linetype='Oncogene', colour='APOBEC'), span=1)+
    geom_smooth(data= df[df$dominantSignature == 'Signature.APOBEC',], aes(x=tmb, y=TSGObs - TSGExp, linetype='TSG', colour='APOBEC'), span=1)+
    
    geom_smooth(data= df[df$dominantSignature == 'Signature.7',], aes(x=tmb, y=OncogeneObs - OncogeneExp, linetype='Oncogene', colour='UV'), span=1)+
    geom_smooth(data= df[df$dominantSignature == 'Signature.7',], aes(x=tmb, y=TSGObs - TSGExp, linetype='TSG', colour='UV'), span=1)+
    
    geom_smooth(data= df[df$dominantSignature == 'Signature.SMOKING',], aes(x=tmb, y=OncogeneObs - OncogeneExp, linetype='Oncogene', colour='SMOKING'), span=1)+
    geom_smooth(data= df[df$dominantSignature == 'Signature.SMOKING',], aes(x=tmb, y=TSGObs - TSGExp, linetype='TSG', colour='SMOKING'), span=1)+
    
    scale_linetype_manual(values=c("dotted", "solid"))+
    xlim(0,200)+
    #scale_color_manual(values=c('#CB9D06', '#4682b4'))+
    ylab('Difference between observed & expected indels')+
    theme_classic()+
    scale_color_manual(values=c("#FF0000","#267574",
                                "#ADFF2F", '#ffb347', '#2A52BE', "#FFF600"))
    ggtitle('All Exomes')
  return(p)
}

plot_figure_1e_msi_only <- function(df){
  p <- ggplot(df[df$dominantSignature == 'Signature.MMR',], aes(x=tmb))+
    geom_smooth(aes(y=OncogeneExp, linetype='expected', colour='Oncogene'), span=1)+
    geom_smooth(aes(y=OncogeneObs, linetype='observed', colour='Oncogene'), span=1)+
    geom_smooth(aes(y=TSGExp, linetype='expected', colour='TSG'), span=1)+
    geom_smooth(aes(y=TSGObs, linetype='observed', colour='TSG'), span=1)+
    scale_linetype_manual(values=c("dotted", "solid"))+
    scale_color_manual(values=c('#CB9D06', '#4682b4'))+
    ylab('N indels')+
    theme_classic()+
    ggtitle('MSI only')+
    xlim(0,200)
  return(p)
}

#plot figure 1e
figure1eDataFrame <- read.table(paste(plottingDataPath, 'figure_1e.tsv', sep=''), sep = '\t', header=TRUE)
plt1e <- plot_figure_1e(figure1eDataFrame)
saveFilePath = paste(plottingFilePath, 'figure1e.pdf')
ggsave(saveFilePath,
       plot=p,  width = 3, height = 4, units = c("in"), limitsize = FALSE)

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

#FINAL PLOT

#First pad plots

#1a
RightPadRatio <- 0
LowerPadRatio <- .5
padded1a <- plot_grid(ggplot()+ggtitle('a.'), plot_grid(plt1a, ggplot(), ncol=2, rel_widths = c(1, RightPadRatio)),
                      ggplot(), nrow=3, rel_heights = c(.2, 1, LowerPadRatio))
#1b
f1RightPadRatio <- 1.5
f1LowerPadRatio <- 0
padded1b <- plot_grid(ggplot()+ggtitle('b.'), plot_grid(plt1b, ggplot(), ncol=2, rel_widths = c(1, f1RightPadRatio)),
                      ggplot(), nrow=3, rel_heights = c(.2, 1, LowerPadRatio))

#1c
f1RightPadRatio <- 0
f1LowerPadRatio <- .2
padded1c <- plot_grid(ggplot()+ggtitle('c.'), plot_grid(plt1c, ggplot(), ncol=2, rel_widths = c(1, f1RightPadRatio)),
                      ggplot(), nrow=3, rel_heights = c(.2, 1, LowerPadRatio))

#1d
f1RightPadRatio <- .5
f1LowerPadRatio <- 0
padded1d <- plot_grid(ggplot()+ggtitle('d.'), plot_grid(plt1d, ggplot(), ncol=2, rel_widths = c(1, f1RightPadRatio)),
                      ggplot(), nrow=3, rel_heights = c(.2, 1, LowerPadRatio))

#1e
f1RightPadRatio <- .5
f1LowerPadRatio <- 0
padded1e <- plot_grid(ggplot()+ggtitle('e.'), plot_grid(plt1e, ggplot(), ncol=2, rel_widths = c(1, f1RightPadRatio)),
                      ggplot(), nrow=3, rel_heights = c(.2, 1, LowerPadRatio))


combinedPdfPlot <- plot_grid(padded1a, padded1b, 
                             padded1c, padded1d, padded1e,
                             ggplot(), ggplot(), ggplot(), ggplot(), #include extra plots for padding purposes
                             nrow=3, ncol=3)
saveFilePath = paste(plottingFilePath, 'figure1.pdf')
ggsave(saveFilePath,
       plot=combinedPdfPlot,  width = 20, height = 20, units = c("in"), limitsize = FALSE)



