#written by Noah Friedman

plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/'
library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(plyr)
library(data.table); setDTthreads(6)
library(stringr)
library(ggrepel)

emptyTheme <- theme(axis.line = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())


df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/figure2aGeneTypes.tsv', sep='\t', header=TRUE)

p <- ggplot(df, aes(x=burdenType, fill=mutType, y=frac))+
  geom_bar(stat = "identity", position="fill")+
  xlab('mutation burden')+
  ylab('fraction of drivers')+
  theme(axis.text.x = element_text(angle=90))+
  scale_fill_viridis_d()+
  scale_x_discrete(limits = rev(levels(df$burdenType)))+
  xlab('tumor type')+
  emptyTheme

ggsave('~/Desktop/plot.pdf', plot=p,  width = 3, height = 4, units = c("in"))

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
    emptyTheme+
    ggtitle('All genes')
  return(p)
}

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/relatedUnrelated.tsv', sep='\t', header=TRUE)

p <- plot_figure_2b(df)
ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 4, units = c("in"))

p <- ggplot()+
  geom_smooth(data = df[df$cancerTypeAdj == 'Glioma',], aes(x=TMB, y=fracDriverRelated, colour='Related genes: Glioma'), method='loess', span=1)+
  geom_smooth(data = df[df$cancerTypeAdj == 'Glioma',], aes(x=TMB, y=fracDriverUnrelated, colour='Unrelated genes: Glioma'), method='loess', span=1)+
  geom_smooth(data = df[df$cancerTypeAdj == 'Other',], aes(x=TMB, y=fracDriverRelated, colour='Related genes: Other'), method='loess', span=1)+
  geom_smooth(data = df[df$cancerTypeAdj == 'Other',], aes(x=TMB, y=fracDriverUnrelated, colour='Unrelated genes: Other'), method='loess', span=1)+
  scale_color_manual(values=c('Purple', '#b19cd9', '#fed8b1', 'Orange'))+
  ylab('N putative drivers/N all mutations')+
  emptyTheme+
  ylim(0,1)

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 4, units = c("in"))

#
#######
#############
#####################
##################################
#####################
#############
########
#



#DEPRECATED
make_dnds_plot <- function(dndsData, title){
  
  
  emptyTheme <- theme(axis.line = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
  
  capThresh <- 1e-9
  signifThresh = .0001
  plotThreshHyper <- .001
  plotThreshNormal <- .000001
  minCoord = .9
  
  maxCoord = 2-log10(capThresh)
  
  ggplot()+
    geom_text_repel(data = dndsData[((dndsData$qglobal_cv.Normal > signifThresh) & (dndsData$qglobal_cv.Hyper < plotThreshHyper)),],
                    aes(x=1- log10(qglobal_cv.Normal + capThresh), y=1-log10(qglobal_cv.Hypermutated+ capThresh),
                        label=gene_name, colour=cancerType))+
    geom_point(data = dndsData,
               aes(x=1- log10(qglobal_cv.Normal + capThresh), y=1-log10(qglobal_cv.Hypermutated + capThresh),
                   colour=cancerType))+
    
    scale_x_continuous(breaks=c(1,2,3,5,1-log10(capThresh)), labels=c('1', '2', '3', '5', '>10'))+
    scale_y_continuous(breaks=c(1,2,3,5,1-log10(capThresh)), labels=c('1', '2', '3', '5', '>10'))+
    
    #BELLS and whistles for the plot
    geom_segment(aes(x=minCoord, xend=maxCoord, y=1- log10(.01), yend= 1- log10(.01)), colour='black', linetype=2)+
    geom_segment(aes(x=1- log10(.01), xend=1- log10(.01), y=minCoord, yend= maxCoord), colour='black', linetype=2)+
    
    geom_segment(aes(x=minCoord, xend=maxCoord, y=1- log10(.1), yend= 1- log10(.1)), colour='black')+
    geom_segment(aes(x=1- log10(.1), xend=1- log10(.1), y=minCoord, yend= maxCoord), colour='black')+
    ylab('1 minus log(q) value in hypermutators')+
    xlab('1 minus log(q) value in non-hypermutators')+
    emptyTheme+
    labs(caption='runDnDsCv.R\nmake_dnds_figure_2c.ipynb')+
    ggtitle(title)
}

#current file lives @ '/Users/friedman/Desktop/WORK/dataForLocalPlotting/dndscvSummary.tsv'
figure2aDataFrame <- read.table(paste(plottingFilePath, 'figure2aDNDSSummary.tsv', sep='') , sep = '\t', header=TRUE)
p <- make_dnds_plot(figure2aDataFrame, 'DNDS-CV Comparing Hypermutated and\n Non-Hypermutated Cancers')
ggsave(paste(plottingFilePath, 'figure2a.pdf', sep=''), plot=p,  width = 6, height = 6)

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

df <- read.table('/Users/friedman/Desktop/hypermutationProjectFinal/tables/table2.tsv', sep='\t', header=TRUE)

p <- ggplot(df, aes(x=TMB))+
        geom_smooth(aes(y=NEUTRAL_TRUNCATING_RATE, color='Neutral'), span=1, se=FALSE)+
        geom_smooth(aes(y=ESSENTIAL_TRUNCATING_RATE, color='Essential'), span=1, se=FALSE)+
        geom_smooth(aes(y=TSG_TRUNCATING_RATE, color='TSG'), span=1, se=FALSE)+
        geom_smooth(aes(y=ONCOGENE_TRUNCATING_RATE,  color='Oncogene'), span=1, se=FALSE)+
  
        stat_summary_bin(aes(y=NEUTRAL_TRUNCATING_RATE, color='Neutral'), bins=8, alpha=0.5)+
        stat_summary_bin(aes(y=ESSENTIAL_TRUNCATING_RATE, color='Essential'), bins=8, alpha=0.5)+
        stat_summary_bin(aes(y=TSG_TRUNCATING_RATE, color='TSG'), bins=8, alpha=0.5)+
        stat_summary_bin(aes(y=ONCOGENE_TRUNCATING_RATE,  color='Oncogene'), bins=8, alpha=0.5)+      
  
        emptyTheme+
        scale_color_manual(values=c('Green', 'Gray', 'Red', 'Blue'))+
        ylab('Rate of truncating mutation')+
        xlab('TMB')+
        xlim(0,200)+
        ylim(0,40)

ggsave('~/Desktop/plot.pdf', plot=p,  width = 3, height = 4)



#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2C

plot_figure_2c <- function(df){
  p <- ggplot(df, aes(x=perCase_c1, y=perCase_c2, color=pathway))+
    geom_point(alpha=0.5)+
    geom_point(data=df[df$qVal < .05,], pch=21, colour='black')+ #outlines
    geom_text_repel(data=df[(df$perCase_c1 > .2) | (df$perCase_c2 > .2), ], aes(label=Allele), size=2)+
    xlim(0, .75)+
    ylim(0, .75)+
    emptyTheme+
    geom_segment(aes(x=0,y=0,xend=0.5, yend=0.5), colour='black', linetype='dashed')+
    scale_color_manual(values = c('gray', 'blue', 'orange'))+
    ylab('Fraction of Colorectal cases mutated')+
    xlab('Fraction of Endometrial cases mutated')+
    ggtitle('MSI Indels')
  return(p)
}

#Correlation scores
cor(df[df$GeneType == 'oncogene',]$perCaseColorectal, df[df$GeneType == 'oncogene',]$perCaseEndometrial, method='pearson')
cor(df[df$GeneType == 'tsg',]$perCaseColorectal, df[df$GeneType == 'tsg',]$perCaseEndometrial, method='pearson')

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/craigStylePlot.tsv', sep='\t', header=TRUE)
p <- plot_figure_2c(df[(df$c1_cancerType == 'Endometrial Cancer') & (df$c2_cancerType == 'Colorectal Cancer'), ])

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5)
#
######
###############
########################
##################################
########################
###############
######
#
#FIGURE 2D: phasing

make_figure_2d <- function(df){
  p <- ggplot(df, aes(x=factor(label, levels=c('B2M', 'PTEN', 'TP53', 'ARID1A', 'APC',
                                               'related_tsg', 'other_tsg', 'PIK3CA', 'TERT', 'related_oncogene', 'other_oncogene', '1 or 2 VUS', '1 or 2 silent')),
                      y=isTrans))+
    stat_summary()+
    emptyTheme+
    theme(axis.text.x = element_text(angle=90))+
    xlab('Gene Type')+
    ylab('Fraction phasable mutations in trans')
  return(p)
}

df <- read.table('/Users/friedman/Desktop/hypermutationProjectFinal/files/infoFiles/phasingSummary.tsv', sep='\t', header=TRUE)

p <- make_figure_2d(df)
ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5)

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
    emptyTheme
  return(p)
}


df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/doubleHitPlot.tsv', sep='\t', header=TRUE)
p <- plot_figure_2e(df)
ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5)


#
#
#DEPRECATED?

make_dnds_weak_drivers_plot <- function(dndsData, title){
  plotThresh <- .1
  minCoord = .9
  
  dndsData <- dndsData[(dndsData$qglobal_cv.Normal > .01) & (dndsData$qglobal_cv.Normal > .01),]
  ggplot()+
    geom_text_repel(data = dndsData[(dndsData$qglobal_cv.Normal < plotThresh) | (dndsData$qglobal_cv.Hypermutated < plotThresh),],
                    aes(x=1- log10(qglobal_cv.Normal), y=1-log10(qglobal_cv.Hypermutated), label=gene_name, colour=cancerType))+
    geom_point(data = dndsData[(dndsData$qglobal_cv.Normal >= plotThresh) & (dndsData$qglobal_cv.Hypermutated >= plotThresh),],
               aes(x=1- log10(qglobal_cv.Normal), y=1-log10(qglobal_cv.Hypermutated), colour=cancerType))+
    
    geom_point(aes(x=1 - mean(log(dndsData[dndsData$qglobal_cv.Normal > 0,]$qglobal_cv.Normal)),
                   y=1 - mean(log(dndsData[dndsData$qglobal_cv.Hypermutated > 0,]$qglobal_cv.Hypermutated))), colour='orange', size=3)+
    #geom_text_repel(aes(x=1 - mean(log(dndsData[dndsData$qglobal_cv.Normal > 0,]$qglobal_cv.Normal)),
    #                    y=1 - mean(log(dndsData[dndsData$qglobal_cv.Hypermutated > 0,]$qglobal_cv.Hypermutated))), label='mean q_val', colour='orange')+
    
    xlim(minCoord,3)+
    ylim(minCoord,3)+
    geom_segment(aes(x=minCoord, xend=3, y=1- log10(.01), yend= 1- log10(.01)), colour='black', linetype=2)+
    geom_segment(aes(x=1- log10(.01), xend=1- log10(.01), y=minCoord, yend= 3), colour='black', linetype=2)+
    geom_segment(aes(x=minCoord, xend=3, y=1- log10(.1), yend= 1- log10(.1)), colour='black')+
    geom_segment(aes(x=1- log10(.1), xend=1- log10(.1), y=minCoord, yend= 3), colour='black')+
    ylab('1 minus q value in hypermutators')+
    xlab('1 minus q value in non-hypermutators')+
    ggtitle(title)
}




#CURRENTLY CLONALITY PLOT
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/clonalSummary.tsv', sep='\t', header=TRUE)
p <- ggplot()+
  
  stat_summary(data=df[(df$driverType == 'driver') & (df$mutationType == 'Endogenous') &
                         ((df$Hugo_Symbol == 'POLE')),], aes(x=1, y=clonal), colour='#228B22')+
  stat_summary(data=df[(df$driverType == 'driver') & (df$mutationType == 'Endogenous') &
                         ((df$Hugo_Symbol == 'MLH1') | (df$Hugo_Symbol == 'MSH2') | (df$Hugo_Symbol == 'MSH6')| (df$Hugo_Symbol == 'PMS2')),],
               aes(x=2, y=clonal), colour='#228B22')+
  
  stat_summary(data=df[(df$driverType == 'driver') & (df$mutationType == 'Endogenous'),], aes(x=3, y=clonal), colour='#228B22')+
  stat_summary(data=df[(df$driverType == 'driver') & (df$mutationType == 'Endogenous') & (df$geneType == 'TSG'),], aes(x=4, y=clonal), colour='#228B22')+
  stat_summary(data=df[(df$driverType == 'driver') & (df$mutationType == 'Endogenous') & (df$geneType == 'Oncogene'),], aes(x=5, y=clonal), colour='#228B22')+
  stat_summary(data=df[(df$driverType == 'VUS') & (df$mutationType == 'Endogenous'),], aes(x=6, y=clonal), colour='#228B22')+
  
  stat_summary(data=df[(df$driverType == 'driver') & (df$mutationType == 'TMZ-glioma'),], aes(x=7, y=clonal), colour='purple')+
  stat_summary(data=df[(df$driverType == 'driver') & (df$mutationType == 'TMZ-glioma') & (df$geneType == 'TSG'),], aes(x=8, y=clonal), colour='purple')+
  stat_summary(data=df[(df$driverType == 'driver') & (df$mutationType == 'TMZ-glioma') & (df$geneType == 'Oncogene'),], aes(x=9, y=clonal), colour='purple')+
  stat_summary(data=df[(df$driverType == 'VUS') & (df$mutationType == 'TMZ-glioma'),], aes(x=10, y=clonal), colour='purple')+
  
  ylim(0,1)+
  ylab('Fraction of mutations clonal')+
  scale_x_continuous(breaks=c(1,2,3,4,5,
                              6,7,8,9,10), labels=c( 'POLE', 'MMR genes', 'All Drivers', 'TSG Drivers', 'Oncogene Drivers', 'All VUS',
                                                     'All Drivers', 'TSG Drivers', 'Oncogene Drivers', 'All VUS'))+
  theme(axis.text.x = element_text(angle=90))+
  xlab('MMR/POLE tumors                TMZ Glioma')+
  emptyTheme+
  labs(caption='plotFigure2.R\nmutation_clonality_analyses_4a.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 4, height = 4)









