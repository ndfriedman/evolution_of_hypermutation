#written by Noah Friedman

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggpubr)
library(ggrepel)

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

plottingDataPath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure2/FIGURE2_PLOTTING_FILES/plotDataFiles/'
plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure2/FIGURE2_PLOTTING_FILES/figurePdfs/'

#
###
######
#############
##################
##########################
#################
############
#######
##
#Figure S2 (i) indel and stop gain prevalence by signature type

plot_figure_s2_i <- function(df){
  p <- ggplot(df, aes(x=dominantSignatureAdj, fill=variable, y = frac))+
    geom_bar(stat='identity')+
    theme(axis.text.x = element_text(angle=90))+
    scale_fill_viridis_d()+
    emptyTheme+
    xlab('signature aetiology')+
    ylab('fraction of drivers')
  return(p)
}

figureS2iDf <- read.table(paste(plottingDataPath, 'figureS2_i.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_figure_s2_i(figureS2iDf)
saveFilePath = paste(plottingFilePath, 'figureS2_i.pdf')
ggsave(saveFilePath, plot=p,  width = 4, height = 3)


#
######
###########
###################
############################
####################################
###########################################
####################################
############################
###################
###########
######
#
#Figure S2 (ii) Fraction drivers in related genes

plot_figure_s2_ii <- function(df){
  p <- ggplot(df[(df$nUnrelatedExpected > .01),],
            aes(x=nmut_IM341, colour=dominantSignature))+
  #geom_smooth(aes(y=(nRelatedObserved - nRelatedExpected)/(nRelatedObserved - nRelatedExpected + nUnrelatedObserved - nUnrelatedExpected),  linetype='related'), method='loess', span=1, se=FALSE)+
  geom_smooth(aes(y=(nRelatedObserved)/(nRelatedObserved + nUnrelatedObserved),  linetype='unrelated'), method='loess', span=1, se=FALSE)+
  scale_x_log10()+
  emptyTheme+
  ylab('fraction drivers in related genes')+
  xlab('mutations in IM-341 genes')+
  ylim(0,1)+
  scale_colour_manual(values = c("#FF0000","#267574",
                                 'gray', "#ADFF2F", '#2A52BE', "#FFF600"))
  return(p)
}

figureS2iiDf <- read.table(paste(plottingDataPath, 'figureS2_ii.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_figure_s2_ii(figureS2iiDf)
saveFilePath = paste(plottingFilePath, 'figureS2_ii.pdf')
ggsave(saveFilePath, plot=p,  width = 4, height = 3.5, units = c("in"))


#
######
###########
###################
############################
####################################
###########################################
####################################
############################
###################
###########
######
#
#Figure S2 (iii) DNDS significance values across 'related' and 'unrelated genes'

#df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/plotGeneFractions.tsv', sep = '\t', header=TRUE)

plot_figure_s2_iii <- function(df){
  p <- ggplot(df[(df$CANCER_TYPE == 'Glioma') | (df$CANCER_TYPE == 'Colorectal Cancer') 
                 | (df$CANCER_TYPE == 'Endometrial Cancer'), ],
              aes(x=cancerTypeGene, y=dndsIsSignificantHyper, fill=CANCER_TYPE))+
    stat_summary(geom='bar')+
    ylim(0, 1)+
    theme(axis.text.x = element_text(angle=90))+
    ylab('Fraction genes signficant by DNDS')+
    emptyTheme+
    ggtitle('Hypermutated')
  return(p)
}

figureS2iiiDf <- read.table(paste(plottingDataPath, 'figureS2_iii.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_figure_s2_iii(figureS2iiiDf)
saveFilePath = paste(plottingFilePath, 'figureS2_iii.pdf')
ggsave(saveFilePath, plot=p,  width = 4, height = 4, units = c("in"))

#
######
###########
###################
############################
####################################
###########################################
####################################
############################
###################
###########
######
#
#Figure S2 (iv) Gene mutation enrichment in hypermutated tumors

plot_figure_s2_iv <- function(df){
  plot_supplementary_enrichment_figure <- function(df, sigType, ratioThreshold=100, fracHyperThreshold=.25){
    df <- df[df$signature == sigType,]
    hyperToNormalRatio = median(df$ratio)
    p <- ggplot()+
      geom_point(data=df, aes(x=fracHyper, y=ratio), colour='gray')+
      geom_segment(aes(x=0,xend=.5, y=hyperToNormalRatio, yend=hyperToNormalRatio, colour='median enrichment'), linetype='dashed')+
      geom_text_repel(data=df[(df$ratio > ratioThreshold) & (df$fracHyper > fracHyperThreshold),], 
                      aes(x=fracHyper, y=ratio, label=gene), force=10)+
      ylab('enrichment of putative drivers\nin hypermutated tumors')+
      xlab('percent of hypermutated cases with\nputative driver mutation')+
      scale_y_log10(limits=c(.1, 100000))+
      xlim(0,1)+
      emptyTheme+
      ggtitle(sigType)
    return(p)
  }
  
  alignedPlot <- plot_grid(plot_supplementary_enrichment_figure(df, 'POLE', ratioThreshold = 75),
                           plot_supplementary_enrichment_figure(df, 'MSI', ratioThreshold = 75),
                           plot_supplementary_enrichment_figure(df, 'APOBEC', ratioThreshold = 75)
                           , ncol=3)
  return(alignedPlot)
}

figureS2ivDf <- read.table(paste(plottingDataPath, 'figureS2_iv.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_iv(figureS2ivDf) #FYI trying to display the variable p will cause an R viewport error (it saves fine though)
saveFilePath = paste(plottingFilePath, 'figureS2_iv.pdf')
ggsave(saveFilePath, plot=p,  width = 15, height = 4, units = c("in"))

#
######
################
##########################
#####################################
##########################
################
#######
#
#Figure S2 (v) mutation enrichment and pole pentanucleotides



#
######
################
##########################
#####################################
##########################
################
#######
#
#Figure S2 (v) cumulative hotspots explained by related and unrelated mutations

#
######
################
##########################
#####################################
##########################
################
#######
#
#Figure S2 (vi)



#
######
################
##########################
#####################################
##########################
################
#######
#
#Figure S2 (vii) Figure 2c style plot for all msi cancer types with substantial number of samples

#NOTE: the 'plot_comp' function is used for both this figure and figure S2_viii
plot_comp <- function(df, title='MSI', thresh=.2){
  p <- ggplot(df, aes(x=perCase_c1, y=perCase_c2))+#, color=pathway))+
    geom_point(alpha=0.5, colour='#D3D3D3')+
    geom_point(data=df[df$qVal < .05,], pch=21, colour='black')+ #outlines
    geom_text_repel(data=df[(df$perCase_c1 > thresh) | (df$perCase_c2 > thresh), ], aes(label=Allele), size=2)+
    xlim(0, .75)+
    ylim(0, .75)+
    emptyTheme+
    geom_segment(aes(x=0,y=0,xend=0.5, yend=0.5), colour='black', linetype='dashed')+
    scale_color_manual(values = c('gray', 'blue', 'orange'))+
    xlab(paste(df[1,]$c1_cancerType, 'n cases:', df[1,]$total_C1))+
    ylab(paste(df[1,]$c2_cancerType, 'n cases:', df[1,]$total_C2))+
    ggtitle(title)
  return(p)
}

plot_figure_s2_vii <- function(df){
  l <- list()
  for(comp in unique(df$comp)){
    dfHere <- df[df$comp == comp,]
    l[[comp]] <- plot_comp(dfHere)
  }
  p <- plot_grid(plotlist=l)
  return(p)
}

figureS2viiDf <- read.table(paste(plottingDataPath, 'figureS2_vii.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_vii(figureS2viiDf) #FYI trying to display the variable p will cause an R viewport error (it saves fine though)
saveFilePath = paste(plottingFilePath, 'figureS2_vii.pdf')
ggsave(saveFilePath, plot=p,  width = 15, height = 10, units = c("in"))


#
######
################
##########################
#####################################
##########################
################
#######
#
#Figure S2 (viii)

#NOTE: the 'plot_comp' function used for this figure is defined in S2_vii

plot_figure_s2_viii <- function(df){
  p <- plot_comp(df, title='POLE', thresh=.3)
  return(p)
}

figureS2viiiDf <- read.table(paste(plottingDataPath, 'figureS2_viii.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_viii(figureS2viiiDf) #FYI trying to display the variable p will cause an R viewport error (it saves fine though)
saveFilePath = paste(plottingFilePath, 'figureS2_viii.pdf')
ggsave(saveFilePath, plot=p,  width = 4, height = 4, units = c("in"))



#
######
################
##########################
#####################################
##########################
################
#######
#
#Figure S2 (ix)

plot_figure_s2_ix <- function(df){
  hyperPlot <- ggplot(df, aes(x=cancerType))+
    geom_bar(aes(y=1, fill='LOH'), stat='identity')+
    geom_bar(aes(y=nHyperComposite/totalHyper, fill='composite'), stat='identity')+
    theme(axis.text.x = element_text(angle=90))+
    emptyTheme+
    ylab('')+
    ggtitle('hypermutated cases')

  normalPlot <- ggplot(df, aes(x=cancerType))+
   geom_bar(aes(y=1, fill='LOH'), stat='identity')+
    geom_bar(aes(y=nNormalComposite/totalNormal, fill='composite'), stat='identity')+
    theme(axis.text.x = element_text(angle=90))+
    emptyTheme+
    ggtitle('non-hypermutated cases')+
    ylab('fraction')

  combinedPlot <- plot_grid(normalPlot, hyperPlot, ncol=2)
  return(combinedPlot)
}

figureS2ixDf <- read.table(paste(plottingDataPath, 'figureS2_ix.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_ix(figureS2ixDf) #FYI trying to display the variable p will cause an R viewport error (it saves fine though)
saveFilePath = paste(plottingFilePath, 'figureS2_ix.pdf')
ggsave(saveFilePath, plot=p,  width = 6, height = 4, units = c("in"))

#
####
###########
######################
#############
#########
####


#
####
###########
######################
#############
#########
####


#
####
###########
######################
#############
#########
####
#TEMP
emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())


df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/genesMutatedExomeCurves.tsv', sep = '\t', header=TRUE)
p <- ggplot(df, aes(x=n, y=frac, group=cohort, color=cohort))+
  geom_path()+
  scale_color_manual(values=c('#D3D3D3', 'red', '#C0C0C0', 'pink', '#808080', 'orange'))+
  emptyTheme+
  xlab('<----less frequently mutated genes   more frequently mutated genes ---->')+
  ylab('Fraction of mutations')+
  labs(caption='tcga_gene_mutation_figure_5_unkonwn.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 7, height = 3, units = c("in"))

#
######
###########
###################
############################
####################################
###########################################
####################################
############################
###################
###########
######
#


#
######
###########
###################
############################
####################################
###########################################
####################################
############################
###################
###########
######
#

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/relatedGeneFracs.tsv', sep = '\t', header=TRUE)
p <- ggplot(df)+
  stat_summary(aes(x=cohortThresh, y=fracRelated, fill=hypermutantStats), geom = 'bar')+
  stat_summary(aes(x=cohortThresh, y=fracRelated, group=cancerType, color=cancerType), geom='point')+
  emptyTheme+
  scale_x_discrete(
                 labels=c("1%", "", "2%", "", "4%", "", '10%', ''
                ))+
  theme(axis.ticks.x = element_blank())+
  ylim(0,1)+
  ylab('Fraction drivers in related genes')+
  xlab('Threshold of mutated cases for a gene to be considered related')+
  scale_fill_manual(values=c('black', 'gray'))+
  ggtitle('Related genes and hypermutation')+
  labs(caption='plotSupplementaryFiguresFig2.R\nmake_obs_vs_expected_by_gene_type_figure2b.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 6, height = 5, units = c("in"))


#
######
###########
###################
############################
####################################
###########################################
####################################
############################
###################
###########
######
#

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/tsgGeneFracs.tsv', sep = '\t', header=TRUE)
p <- ggplot(df, aes(x=tmb, y=fracTsg, group=signatureCohort, color=signatureCohort))+
  geom_smooth(se=FALSE)+
  scale_x_log10()+
  ylim(0,1)+
  ylab('Fraction of drivers in TSGs')+
  emptyTheme+
  labs(caption='plotSupplementaryFiguresFig2.R\nmake_obs_vs_expected_by_gene_type_figure2b.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 4, height = 4, units = c("in"))


#
##########
#####################
################################
############################################
################################
#####################
##########
#

#DEPRECATED FIGURE 1

#
#



#DEPRECATED?

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



#df <- read.table('~/Desktop/WORK/dataForLocalPlotting/relatedUnrelated.tsv', sep='\t', header=TRUE)

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

#DNDS supplemental figure
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/dndsSupplement.tsv', sep = '\t', header=TRUE)

p <- ggplot(df[(df$CANCER_TYPE == 'Glioma') | (df$CANCER_TYPE == 'Colorectal Cancer') 
               | (df$CANCER_TYPE == 'Endometrial Cancer'), ],
            aes(x=cancerTypeGene, y=dndsIsSignificantHyper, fill=CANCER_TYPE))+
  stat_summary(geom='bar')+
  ylim(0, 1)+
  theme(axis.text.x = element_text(angle=90))+
  ylab('Fraction genes signficant by DNDS')+
  emptyTheme+
  ggtitle('Hypermutated')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5, units = c("in"))

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/dndsAndGeneMutations.tsv', sep = '\t', header=TRUE)
p <- ggplot(df, aes(x=scoreAdj, y=frac, group=cohort, color=hypermutantStatus))+
  geom_path()+
  geom_segment(aes(x=0.05, xend=0.05, y=0, yend=1), colour='gray', alpha=0.5, linetype='dashed')+
  scale_x_log10(limits =c(1e-10, 1), breaks=c(1e-10, 1e-5, 1e-3, 5e-2))+
  emptyTheme+
  theme(axis.text.x=element_text(angle=90))+
  xlab('DNDS q value')+
  ylab('Fraction of all drivers')+
  labs(caption = 'plotSupplementaryFiguresFig2.R\nmake_dnds_figures_2a_2c.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 3, units = c("in"))

#
#
#
#
#
#
#
#
#
#
#
#
#
###MUTATION ATTRIBUTION

plotPie <- function(df, type, includeLegend = FALSE){
  p <- ggplot(df[df$geneType == type,], aes(x= '', fill = hypermutationInduced))+
    geom_bar(width=1)+
    coord_polar("y", start=0)+
    emptyTheme+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())+
    ggtitle(type)+
    scale_fill_manual(values=c('black', '#A9A9A9', '#C0C0C0'))
  if(includeLegend == FALSE){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}


dfAttribution <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/mutAttributionData.tsv', sep='\t', header=TRUE)
alignedP <- plot_grid(
  
  get_legend(plotPie(dfAttribution, 'tsgTruncating', includeLegend = TRUE)),
  plotPie(dfAttribution, 'oncogene'),
  plotPie(dfAttribution, 'tsgMissense'),
  plotPie(dfAttribution, 'tsgTruncating'),
  nrow=4)

ggsave('~/Desktop/plot.pdf', plot=alignedP,  width = 2, height = 4, units = c("in"))




p <- ggplot(df[(df$cancerType == 'Colorectal Cancer'),], aes(x=n, y=frac,
                                                             color=cohort, group=cohort))+
  geom_path()+
  emptyTheme+
  xlab('nth gene')
scale_x_log10()

p <- ggplot(df, aes(x=n, y=frac, color=cohort, group=cohort))+
  geom_smooth(se=FALSE)+
  emptyTheme+
  ylab('fraction of all drivers')+
  xlab('nth most commonly mutated gene')+
  scale_color_manual(values=c('red', '#ffcccb', 'blue', '#add8e6'))+
  scale_x_log10(breaks=c(1,10,50,100))+
  ggtitle('Diversity of drivers')+
  labs(caption = 'plotSupplementaryFiguresFig2.R\nmake_obs_vs_expected_by_gene_type_figure2b.ipynb')






