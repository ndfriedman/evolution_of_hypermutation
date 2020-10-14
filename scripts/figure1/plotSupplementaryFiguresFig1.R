#written by Noah Friedman
#scripts for plotting supplementary figures for figure 1

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(stringr)
library(gdata)
library(ggpubr)

emptyTheme <- theme(axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

#adjust this as needed
plottingDataPath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/plotDataFiles/'
plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/figurePdfs/'

# Figure S1(a) Distributions of TMB by cancer type: refer to code in scripts/utilityScripts/plotAndDefineHypermutationThresholds.R

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
# Figure S1 (b) Dominant Signatures in Cohort by Cancer Type

plot_figure_s1_b <- function(df){
  p <- ggplot(df, aes(x=reorder(cancerType, -orderingVal), y= frac))+
    
    geom_bar(aes(fill=
                   factor(signature, levels=c('APOBEC', 'MMR',
                                              'other', 'POLE', 'POLE&MMR', 'SMOKING', 'TMZ', 'UV')))
             , stat='identity')+
    theme_classic()+
    geom_text(aes(label=nCasesHyperHighType, y=1.05), angle=60)+
    theme(axis.text.x = element_text(angle=90))+
    theme(axis.ticks.x = element_blank())+
    xlab('Cancer type')+
    ylab('Fraction of cases')+
    scale_fill_manual(values=c("#FF0000","#267574",
                               'gray', "#ADFF2F", '#ffb347', '#2A52BE', "#FFF600"))+
    guides(fill=guide_legend(title="Signature"))+
    ggtitle('S1(b.)')
  return(p)
}

figureS1bDf <- read.table(paste(plottingDataPath, 'figure_1b.tsv', sep=''), sep='\t', header=TRUE)
pltS1b <- plot_figure_s1_b(figureS1bDf)
saveFilePath = paste(plottingFilePath, 'figureS1_b.pdf')
ggsave(saveFilePath,
       plot=pltS1b,  width = 5, height = 5)



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
#Figure S1 (c) Ratio of passengers to driver by TMB

plot_figure_s1_c <- function(df){
  ggplot(df, aes(x=nMut, y=nDriver/(nMut - nDriver), colour=dominantSignature))+
    geom_smooth(se = FALSE)+
    theme_classic()+
    xlab('mutations in IMPACT-341 genes')+
    ylab('Ratio of drivers/VUS')+
    scale_colour_manual(values=c("#FF0000",
                                 'black', "#267574", 'gray', "#ADFF2F",
                                 '#ffb347', '#2A52BE', "#FFF600"))+
    ggtitle('S1(c.)')
}

figureS1cDf <- read.table(paste(plottingDataPath, 'figureS1_c.tsv', sep=''), sep='\t', header=TRUE)
pltS1c <- plot_figure_s1_c(figureS1cDf)
saveFilePath = paste(plottingFilePath, 'figureS1_c.pdf')
ggsave(saveFilePath,
       plot=pltS1c,  width = 5, height = 3)
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
#Figure S1 (z1) observed vs expected by tmb tier
#TODO assign this to a figure
plot_figure_s1_z1 <- function(df){
  
  plot_fig <- function(df, title){
    p <- ggplot(df, aes(x=nmut))+
      geom_smooth(aes(y = expectedOncogenicSNP, colour = 'Expected'))+
      geom_smooth(aes(y = obsOncogenicSNP, colour = 'Observed'))+
      scale_colour_manual(values=c('gray', 'black'))+
      ylab('Putative SNP drivers')+
      xlab('n SNV in IMPACT-341 genes')+
      theme_classic()+
      ggtitle(title)+
      coord_cartesian(xlim=c(0,300), ylim=c(0,50))
    return(p)
  }
  
  l <- list()
  for(sig in unique(df$dominantSignature)){
    dfHere <- df[df$dominantSignature == sig,]
    l[[sig]] <- plot_fig(dfHere, sig)
  }
  gridP <- plot_grid(plotlist=l)
  finalP <- plot_grid(ggplot()+ggtitle('S1(z1.)'), gridP, nrow=2, rel_heights = c(.1,1))
  return(finalP)
}


#this figure is made using the Figure 1d data fram
figureS1z1Df <- read.table(paste(plottingDataPath, 'figure_1d.tsv', sep=''), sep='\t', header=TRUE)
pltS1z1 <- plot_figure_s1_z1(figureS1z1Df[figureS1z1Df$dominantSignature != 'other',])
saveFilePath = paste(plottingFilePath, 'figureS1_z1.pdf')
ggsave(saveFilePath,
       plot=pltS1z1,  width = 10, height = 10)

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

plot_figure_s1_z2 <- function(df){
  
  plot_fig <- function(df, sig){
    p <- ggplot(df, aes(x=nIndels/30))+
      geom_smooth(aes(y=OncogeneObs, linetype='Observed',  color='oncogene'), span=1)+
      geom_smooth(aes(y=OncogeneExp, linetype='Expected',  color='oncogene'), span=1)+
      geom_smooth(aes(y=TSGObs, linetype='Observed',  color='TSG'), span=1)+
      geom_smooth(aes(y=TSGExp, linetype='Expected', color='TSG'), span=1)+
      scale_linetype_manual(values=c("dotted", "solid"))+
      theme_classic()+
      ylab('N indels ')+
      xlab('N indels/MB')+
      coord_cartesian(xlim = c(0,40), ylim=c(0,40))+
      ggtitle(sig)
  }
  
  l <- list()
  for(sig in unique(df$dominantSignature)){
    dfHere <- df[df$dominantSignature == sig,]
    l[[sig]] <- plot_fig(dfHere, sig)
  }
  gridP <- plot_grid(plotlist=l)
  finalP <- plot_grid(ggplot()+ggtitle('S1(z2.)'), gridP, nrow=2, rel_heights = c(.1,1))
  
  ggtitle('All Exomes')
  return(finalP)
}

#uses the figure 1e data frame
figureS1z2Df <- read.table(paste(plottingDataPath, 'figure_1e.tsv', sep=''), sep='\t', header=TRUE)
#REMOVE one weird apobec case
pltS1z2 <- plot_figure_s1_z2(figureS1z2Df[figureS1z2Df$Tumor_Sample_Barcode != 'TCGA-D8-A27V-01A-12D-A17D-09',])
saveFilePath = paste(plottingFilePath, 'figureS1_z2.pdf')
ggsave(saveFilePath,
       plot=pltS1z2,  width = 10, height = 10)

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
#Figure S1 (d) Distributions of observed and expected mutations in TCGA

plot_figure_s1_d <- function(df){
  p1 <- ggplot(df)+
    geom_density(aes(x=impactVUSSNPExpected, linetype='expected'))+
    geom_density(aes(x=vusObserved, linetype='observed'))+
    emptyTheme+
    scale_x_log10()+
    xlab('n mutations')+
    geom_text(aes(x=30, y=1.75), label='ns')+
    ggtitle('VUS mutations\n in IMPACT genes')+
    scale_linetype_manual(values = c("dashed", "solid"))+
    theme(axis.text.x = element_text(angle=90))
  p2 <- ggplot(df)+
    geom_density(aes(x=oncogenicSNPExpected, linetype='expected'))+
    geom_density(aes(x=oncogenicObserved, linetype='observed'))+
    emptyTheme+
    scale_x_log10()+
    xlab('n mutations')+
    geom_text(aes(x=3, y=1.5), label='*')+
    ggtitle('Drivers')+
    scale_linetype_manual(values = c("dashed", "solid"))+
    theme(axis.text.x = element_text(angle=90))
  alignedPlot <- plot_grid(p1, p2, ncol=2)
  plotWithTitle <- plot_grid(ggplot()+ggtitle('S1(d)'), alignedPlot, rel_heights = c(.1,1), nrow=2)
  return(plotWithTitle)
}

figureS1dDf <- read.table(paste(plottingDataPath, 'figureS1_d.tsv', sep=''), sep='\t', header=TRUE)
pltS1d <- plot_figure_s1_d(figureS1dDf)
saveFilePath = paste(plottingFilePath, 'figureS1_d.pdf')
ggsave(saveFilePath,
       plot=pltS1d,  width = 7, height = 3)

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
#Figure S1 (e) Mutation rates in essential, neutral and TSGs

plot_figure_s1_e <- function(df){
  p <- ggplot(df)+
    geom_violin(bw=1, aes(x='TSG', y=mutRateTSG*1e6))+
    geom_violin(bw = 1, aes(x='neutral', y=mutRateNotEssential*1e6))+
    geom_violin(bw = 1, aes(x='essential', y=mutRateEssential*1e6))+
    theme_classic()+
    stat_summary(aes(x='TSG', y=mutRateTSG*1e6))+
    stat_summary(aes(x='neutral', y=mutRateNotEssential*1e6))+
    stat_summary(aes(x='essential', y=mutRateEssential*1e6))+
    ylab('rate of truncating mutation/MB')+
    coord_cartesian(ylim=c(0,10))+
    xlab('gene type')+
    ggtitle('S1(e)')
  return(p)
}

figureS1eDf <- read.table(paste(plottingDataPath, 'figureS1_e.tsv', sep=''), sep = '\t', header=TRUE)
pltS1e <- plot_figure_s1_e(figureS1eDf)
saveFilePath = paste(plottingFilePath, 'figureS1_e.pdf')
ggsave(saveFilePath,
       plot=pltS1e,  width = 3, height = 4, units = c("in"), limitsize = FALSE)

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
#Figure S1 (f) indel and stop gain prevalence by signature type
plot_figure_s1_f <- function(df){
  pTrunc <- ggplot(df, aes(x=reorder(dominantSignature, isTruncating), y=isTruncating))+
    stat_summary(geom='bar', fill='#838383')+
    stat_summary(geom='errorbar', colour='gray')+
    theme(axis.text.x = element_text(angle=90))+
    ylab('fraction stop gain mutations')+
    xlab('dominant signature')+
    coord_cartesian(ylim=c(0,.3))+
    emptyTheme

  pIndel <- ggplot(df, aes(x=reorder(dominantSignature, isIndel), y=isIndel))+
    stat_summary(geom='bar', fill='#838383')+
    stat_summary(geom='errorbar', colour='gray')+
    theme(axis.text.x = element_text(angle=90))+
    ylab('fraction INDELs')+
    xlab('dominant signature')+
    coord_cartesian(ylim=c(0,.3))+
    emptyTheme
  alignedPlot <- plot_grid(pTrunc, pIndel, ncol=2)
  plotWithTitle <- plot_grid(ggplot()+ggtitle('S1(f)'), alignedPlot, rel_heights = c(.1,1), nrow=2)
  return(plotWithTitle)
}

figureS1fDf <- read.table(paste(plottingDataPath, 'figureS1_f.tsv', sep=''), sep='\t', header=TRUE)
pltS1f <- plot_figure_s1_f(figureS1fDf)
saveFilePath = paste(plottingFilePath, 'figureS1_f.pdf')
ggsave(saveFilePath,
       plot=pltS1f,  width = 5, height = 3)
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
#
#Figure S1 (g) propensity of signatures to cause oncogenic, truncating and hotspot mutations
plot_figure_s1_g <- function(df, mutType, legend=FALSE){
  
  barColorPalette = c(
    "#00DFFF", "#FF0000", "#FF1493", 
    "#267574","#D3D3D3","#ADFF2F",
    "#2A52BE","#FFA500", "purple", "#FFF600"
  )
  p <- ggplot(df[df$mutType == mutType,], aes(x=reorder(Signature_Name, frac), y=frac, fill=colorName))+
    geom_bar(stat='identity')+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    xlab('Signature Name')+
    emptyTheme+
    ylim(0, .11)+
    ylab(paste('chance of a nonsynonymous\nmutation being', mutType))+
    scale_fill_manual(values=barColorPalette)
  if(legend == FALSE){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}

figureS1gDf <- read.table(paste(plottingDataPath, 'figureS1_g.tsv', sep=''), sep='\t', header=TRUE)
pltS1g <- plot_grid(plot_figure_s1_g(figureS1gDf, 'truncating'), 
                   plot_figure_s1_g(figureS1gDf, 'oncogenic'),
                   plot_figure_s1_g(figureS1gDf, 'hotspot'),
                   get_legend(plot_figure_s1_g(figureS1gDf, 'hotspot', legend = TRUE)),
                   ncol=4)

pltS1g <- plot_grid(ggplot()+ggtitle('S1(g)'), pltS1g, rel_heights = c(.1,1), nrow=2)
saveFilePath = paste(plottingFilePath, 'figureS1_g.pdf')
ggsave(saveFilePath, plot=pltS1g,  width = 17, height = 4, units = c("in"))


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
#



#Figure S1 (h) nucleosome positioning and mutation rate

make_figure_s1_h <- function(df){
  dyadPlot <- ggplot(df, aes(x=closestNucleosomeDistance, y=(..count..)/sum(..count..)))+
    stat_bin(geom='point', bins=2000)+
    xlim(-1000, 1000)+
    ylab('Fraction of all mutations')+
    xlab('Distance from dyad center')+
    emptyTheme+
    ggtitle('Distance from nucleosome dyad\nand mutation rate')

  barPlot <- ggplot()+
    stat_summary(data = df[df$driver == 'True',], aes(x='_driver', y=isCloseToNucleosome, fill='putative driver'), geom='bar')+
    stat_summary(data = df[df$driver == 'False',], aes(x='_non-driver', y=isCloseToNucleosome, fill='VUS'), geom='bar')+
    stat_summary(data = df[df$truncatingType == 'truncatingTSG',], aes(x='_truncating TSG', y=isCloseToNucleosome, fill='putative driver'), geom='bar')+
    stat_summary(data = df[df$truncatingType == 'truncatingOncogene',], aes(x='truncating oncogene', y=isCloseToNucleosome, fill='VUS'), geom='bar')+
    ylab('Fraction mutations proximal to dyad')+
    theme(axis.text.x = element_text(angle=90))+
    ggtitle('Nucleosome impact on mutations')+
    scale_fill_manual(values=c('dark gray', 'light gray'))+
    emptyTheme+
    xlab('mutationType')
  alignedPlot <- plot_grid(dyadPlot, barPlot, ncol=2)
  plt <- plot_grid(ggplot()+ggtitle('S1(h)'), alignedPlot, rel_heights = c(.1,1), nrow=2)
  return(plt)
}

#NOTE THIS csv can stress memory, be careful
figureS1hDf <- read.table(paste(plottingDataPath, 'figureS1_h.csv', sep=''), sep=',', header=TRUE)
pltS1h <- make_figure_s1_h(figureS1hDf)
saveFilePath = paste(plottingFilePath, 'figureS1_h.pdf')
ggsave(saveFilePath,
       plot=pltS1h,  width = 8, height = 4)
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
#Figure S1 (i) expression and mutation rate

make_figure_s1_i <- function(df){
  pltBySignature <- ggplot(df)+
    stat_summary(aes(x=Signature, y=mutRateExpressed*1e6, colour='expressed'))+
    stat_summary(aes(x=Signature, y=mutRateNotExpressed*1e6, colour='not expressed'))+
    emptyTheme+
    theme(axis.text.x = element_text(angle=90))+
    ylab('mutations per megabase')+
    xlab('dominant signature')
  
  pltByMutationType <- ggplot(df)+
    stat_summary(aes(x='VUS', y=mutRateVUSExpressed*1e6, colour='expressed'))+
    stat_summary(aes(x='VUS', y=mutRateVUSNotExpressed*1e6, colour='not expressed'))+
    stat_summary(aes(x='Truncating', y=mutRateTruncatingExpressed*1e6, colour='expressed'))+
    stat_summary(aes(x='Truncating', y=mutRateTruncatingNotExpressed*1e6, colour='not expressed'))+
    stat_summary(aes(x='Driver', y=mutRateDriverExpressed*1e6, colour='expressed'))+
    stat_summary(aes(x='Driver', y=mutRateDriverNotExpressed*1e6, colour='not expressed'))+
    emptyTheme+
    theme(axis.text.x = element_text(angle=90))+
    xlab('mutation type')+
    ylab('mutations per megabase')
  
  alignedPlot <- plot_grid(pltBySignature, pltByMutationType, ncol=2)
  plt <- plot_grid(ggplot()+ggtitle('S1(i)'), alignedPlot, rel_heights = c(.1,1), nrow=2)
  return(plt)
}

figureS1iDf <- read.table(paste(plottingDataPath, 'figureS1_i.tsv', sep=''), sep='\t', header=TRUE)
pltS1i <- make_figure_s1_i(figureS1iDf)
saveFilePath = paste(plottingFilePath, 'figureS1_i.pdf')
ggsave(saveFilePath,
       plot=pltS1i,  width = 8, height = 4)


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
#Figure S1 (j) replication timing and mutation rate

make_figure_s1_j <- function(df){
  ggplot()+
    stat_summary(data=df[df$driverType == 'driver',], aes(x='__driver', y=replicationTime, color='driver'))+
    stat_summary(data=df[df$driverType == 'VUS',], aes(x='__VUS', y=replicationTime, color='passenger'))+
    stat_summary(data=df[df$truncatingType == 'truncatingOncogene',], aes(x='truncating_oncogene', y=replicationTime, color='passenger'))+
    stat_summary(data=df[df$truncatingType == 'truncatingTSG',], aes(x='_truncating_TSG', y=replicationTime, color='driver'))+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    ylab('Replication Time:\n<-early                                    late->')+
    xlab('mutation type')+
    scale_color_viridis_d(option='C')+
    ggtitle('S1(j)')
}

figureS1jDf <- read.table(paste(plottingDataPath, 'figureS1_j.tsv', sep=''), sep='\t', header=TRUE)
pS1j <- make_figure_s1_j(figureS1jDf)
saveFilePath = paste(plottingFilePath, 'figureS1_j.pdf')
ggsave(saveFilePath,
       plot=pS1j,  width = 4, height = 4)






































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
#
#Figure S1

plot_figure_s1_g <- function(df){
  p3 <- 
    ggplot(df[df$vusObserved > 100,])+
    geom_boxplot(aes(x=1, y=truncatingEssentialExpected, colour='expected'))+ 
    geom_boxplot(aes(x=2, y=truncatingEssentialObserved, colour='observed'))+
    scale_y_log10()+
    emptyTheme+
    ggtitle('Truncating mutations\nin essential genes')+
    theme(axis.text.x = element_text(angle=90))
}
#We use the same figureS1f tsv for figure S1g
figureS1gDf <- read.table(paste(plottingDataPath, 'figureS1_f.tsv', sep=''), sep='\t', header=TRUE)


#
####
##########
####################
#############################
########################################
####################################################
#DEPRECATED FILES HERE FOR NOW, will be deleted

#WORK on Monday Dec 2

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/essentialGenesMuts.tsv', ,sep = '\t', header=TRUE)

ggplot(df, aes(x= nmut))+
  
  #geom_smooth(aes(y=nEssentialDoubleNormed, colour='essentialDoubleRate'), span=25)+
  #geom_smooth(aes(y=nNeutralDoubleNormed, colour='neutralDoubleRate'), span=25)+
  #geom_smooth(aes(y=nTSGDoubleNormed, colour='tsgDoubleRate'), span=25)
  
  geom_point(aes(y=nEssentialDoubleNormed, colour='essentialDoubleRate'))+
  geom_point(aes(y=nNeutralDoubleNormed, colour='neutralDoubleRate'))+
  geom_point(aes(y=nTSGDoubleNormed, colour='tsgDoubleRate'))

p <- ggplot()+
  geom_smooth(data = df[(df$signature == 'MMR') & (df$isOutlier == 'False'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_mmr'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_mmr'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_mmr'), method = 'lm', se=FALSE)+
  
  stat_summary_bin(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_mmr'), alpha=0.35)+
  stat_summary_bin(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_mmr'), alpha=0.35)+
  stat_summary_bin(data = df[(df$signature == 'MMR')& (df$isOutlier == 'False'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_mmr'), alpha=0.35)+
  
  geom_smooth(data = df[(df$signature == 'POLE') & (df$isOutlier == 'False'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_pole'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_pole'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_pole'), method = 'lm', se=FALSE)+
  
  stat_summary_bin(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_pole'), alpha=0.35, bins=5)+
  stat_summary_bin(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_pole'), alpha=0.35, bins=5)+
  stat_summary_bin(data = df[(df$signature == 'POLE')& (df$isOutlier == 'False'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_pole'), alpha=0.35, bins=5)+
  
  xlab('TMB')+
  ylab('Difference between observed and\nexpected n truncating muts')+
  scale_color_manual(values=c('#013220', '#037D50', '#FF6347',
                              '#FF8C00', '#82CFFD', '#0D4F8B'))+
  labs(colour = 'Signature and\ngene type')+
  emptyTheme+
  scale_x_log10()

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5, units = c("in"))


#
##
######
###############
########################
################
######
##
#


#TEMP
##
#######
#################
###########################
#################
########
##

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/observedVsExpectedByType.tsv', sep='\t', header=TRUE)
s=1
p <- ggplot(df[(df$nmut > 5),], aes(x=nmut))+
  geom_smooth(aes(y = nmutNormal/dif, colour='Non-hypermutated driver genes'), method='loess', span=s)+
  geom_smooth(aes(y = nmutNormalAndHyperStrong/dif, colour='Non-hypermutated driver genes\nand strong hypermutant drivers'), method='loess', span=s)+
  geom_smooth(aes(y = nmutNormalAndHyperAll/dif, colour='Non-hypermutated driver genes\nand strong hypermutant drivers\nand weak hypermutant drivers'), method='loess', span=s)+
  scale_colour_manual(values=c('black', 'red', 'orange'))+
  ylim(0,1)+
  emptyTheme+
  xlab('Nmut in IM-341')+
  ylab('Fraction of divergence between\nobserved and expected\nexplained by:')+
  labs(caption='plotFigure1.R\nmake_obs_vs_expected_by_gene_type.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5, units = c("in"))

#
###
########
###############
########
###
#




##
####
#######
####TEMP DISPLAY AREA
p1 <- ggplot(dfObsExp[dfObsExp$isMsi == 'True',], aes(x=nmut))+
  geom_smooth(aes(y = expectedOncogenicAll, colour = 'Expected'))+
  geom_smooth(aes(y = obsOncogenic, colour = 'Observed'))+
  scale_colour_manual(values=c('black', 'gray'))+
  ylab('All putative driver mutations')+
  xlab('nmut in IMPACT-341 genes')+
  emptyTheme

p2 <- ggplot(dfObsExp, aes(x=nmut))+
  geom_smooth(aes(y = expectedOncogenicSNP, colour = 'Expected'))+
  geom_smooth(aes(y = obsOncogenicSNP, colour = 'Observed'))+
  scale_colour_manual(values=c('black', 'gray'))+
  ylab('Putative SNP drivers')+
  xlab('nmut in IMPACT-341 genes')+
  emptyTheme

p3 <- ggplot(dfObsExp[dfObsExp$isMsi == 'True',], aes(x=nmut))+
  geom_smooth(aes(y = expectedOncogenicIndel, colour = 'Expected'))+
  geom_smooth(aes(y = obsOncogenicINDEL, colour = 'Observed'))+
  scale_colour_manual(values=c('black', 'gray'))+
  ylab('Putative INDEL drivers')+
  ggtitle('[MSI ONLY]')+
  xlab('nmut in IMPACT-341 genes')+
  emptyTheme

p4 <- ggplot(dfObsExp, aes(x=nmut))+
  geom_smooth(aes(y = expectedHotspot, colour = 'Expected'))+
  geom_smooth(aes(y = obsHotspot, colour = 'Observed'))+
  scale_colour_manual(values=c('black', 'gray'))+
  ylab('Putative hotspot drivers')+
  xlab('nmut in IMPACT-341 genes')+
  emptyTheme

p5 <- ggplot(dfObsExp, aes(x=nmut))+
  geom_smooth(aes(y = expectedTruncatingTSG, colour = 'Expected Truncating TSG'))+
  geom_smooth(aes(y = expectedTruncatingOncogene, colour = 'Expected Truncating Oncogene'))+
  geom_smooth(aes(y = obsStopGainTSG, colour = 'Observed TSG stop gain'))+
  geom_smooth(aes(y = obsStopGainOncogene, colour = 'Observed oncogene stop gain'))+
  ylab('Stop-gain mutations')+
  xlab('nmut in IMPACT-341 genes')+
  emptyTheme

plt <- plot_grid(p1, p2, p3, p4, p5, ncol=5, rel_widths = c(.75, .75, .75, .75, 1))
ggsave('~/Desktop/plot.pdf',
       plot=plt,  width = 20, height = 4, units = c("in"))

p1 <- ggplot(dfObsExp, aes(x=nmut, y=obsOncogenic - expectedOncogenicAll, group=dominantSignature, color=dominantSignature))+
  geom_smooth(method='lm')+
  ylab('observed minus expected\n ALL PUTATIVE DRIVERS')+
  emptyTheme



p3 <- ggplot(dfObsExp, aes(x=nmut, y=obsHotspot - expectedHotspot, group=dominantSignature, color=dominantSignature))+
  geom_smooth(method='lm')+
  ylab('observed minus expected\n ALL PUTATIVE HOTSPOT DRIVERS')+
  emptyTheme

p4 <- ggplot(dfObsExp, aes(x=nmut, y=obsStopGainTSG - expectedTruncatingTSG, group=dominantSignature, color=dominantSignature))+
  geom_smooth(method='lm')+
  ylab('observed minus expected\n ALL PUTATIVE STOP GAIN DRIVERS')+
  emptyTheme

plt <- plot_grid(p1, p2, p3, p4, ncol=4)
ggsave('~/Desktop/plot.pdf',
       plot=plt,  width = 15, height = 4, units = c("in"))

#
#######
##################
###############################
#########
#################
#

