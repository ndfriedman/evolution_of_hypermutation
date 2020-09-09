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

emptyTheme <- theme(axis.line = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())

#SUPPLEMENTARY FIGURES


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
#Fraction mutations in panel that are drivers

plot_figure_s1_ii <- function(df){
  ggplot(df, aes(x=nMut, y=fracDriver, colour=dominantSignature))+
    geom_smooth(se = FALSE)+
    emptyTheme+
    xlab('mutations in IM-341 genes')+
    ylab('fraction tmb from drivers')+
    scale_colour_manual(values=c("#FF0000",
      'black', "#267574", 'gray', "#ADFF2F",
      '#ffb347', '#2A52BE', "#FFF600"))
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/fractionDrivers.tsv', sep='\t', header=TRUE)
p <- plot_figure_s1_ii(df)
ggsave('~/Desktop/plot.pdf',
       plot=p,  width = 5, height = 3)


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
#SIGNATURES BY CANCER TYPE SUPPLEMENTARY FIGURE

plot_signatures_figure_sup <- function(df){
  p <- ggplot(df, aes(x=reorder(cancerType, -orderingVal), y= frac))+
    
    geom_bar(aes(fill=
                   factor(signature, levels=c('APOBEC', 'MMR',
                                              'other', 'POLE', 'POLE&MMR', 'SMOKING', 'TMZ', 'UV')))
             , stat='identity')+
    geom_text(aes(label=nCasesHyperHighType, y=1.05), angle=60)+
    theme(axis.text.x = element_text(angle=90))+
    theme(axis.ticks.x = element_blank())+
    xlab('Cancer type')+
    ylab('Fraction of cases')+
    scale_fill_manual(values=c("#FF0000","#267574",
                               'gray', "#ADFF2F","#00dd5d", '#ffb347', '#2A52BE', "#FFF600"))+
    guides(fill=guide_legend(title="Signature"))+
    emptyTheme
  return(p)
}

df <- read.table('/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/figure1bSignatureSummary.tsv',
                                sep='\t', header=TRUE)

pSup <- plot_signatures_figure_sup(df)
ggsave('~/Desktop/plot.pdf',
       plot=pSup,  width = 5, height = 5)

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
#MUTAITON TYPES BY SIGNATURE
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/mutationTypeSummary.tsv',
                 sep='\t', header=TRUE)

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
ggsave('~/Desktop/plot.pdf',
       plot=alignedPlot,  width = 5, height = 3)


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

#Nucleosomes supplemental figure

df = read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/nucleosomeDyadOnfo.tsv', sep=',', header=TRUE)
dyadPlot <- ggplot(df, aes(x=closestNucleosomeDistance, y=(..count..)/sum(..count..)))+
         stat_bin(geom='point', bins=2000)+
         xlim(-1000, 1000)+
         ylab('Fraction of all mutations')+
         xlab('Distance from dyad center')+
         emptyTheme+
         ggtitle('Distance from nucleosome dyad\nand mutation rate')+
         labs(caption='plotSupplementaryFiguresFig1.R\ngenerate_figure1_supplementary_figures.ipynb')

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
ggsave('~/Desktop/plot.pdf', plot=alignedPlot,  width = 8, height = 4, units = c("in"))


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

#ONCOGENICITY OF SIGNATURES
plot_oncogenicity_figure <- function(df, mutType, legend=FALSE){

  
  barColorPalette = c(
    "#00DFFF", "#FF0000", "#FF1493", 
    "#267574","#D3D3D3","#ADFF2F",
    "#2A52BE","#FFA500", "#FFF600"
  )
  p <- ggplot(df[df$mutType == mutType,], aes(x=reorder(Signature_Name, frac), y=frac, fill=colorName))+
    geom_bar(stat='identity')+
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

df <- read.table('/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/sig_oncogenicity_figure1s?.tsv', ,sep = '\t', header=TRUE)

fullP <- plot_grid(plot_oncogenicity_figure(df, 'truncating'), 
                   plot_oncogenicity_figure(df, 'oncogenic'),
                   plot_oncogenicity_figure(df, 'hotspot'),
                   get_legend(plot_oncogenicity_figure(df, 'hotspot', legend = TRUE)),
                   ncol=4)

ggsave('~/Desktop/plot.pdf', plot=fullP,  width = 17, height = 4, units = c("in"))

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


#PIE CHART INSET
sigsInset <- ggplot(dfSigs, aes(x=1, y=nCases, fill=
                                  factor(signature, levels=c('APOBEC', 'MMR',
                                                             'other', 'POLE', 'POLE&MMR', 'SMOKING', 'TMZ', 'UV'))))+
  geom_bar(stat='identity')+
  coord_polar("y", start=0)+
  guides(fill=guide_legend(title="Signature"))+
  scale_fill_manual(values=c("#FF0000","#267574",
                             'gray', "#ADFF2F","#00dd5d", '#ffb347', '#2A52BE', "#FFF600"))+
  emptyTheme+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  xlab('')+
  ylab('')

ggsave(paste(plottingFilePath, 'figure1c_pieChart.pdf', sep=''),
       plot=sigsInset,  width = 5, height = 5)

