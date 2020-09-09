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


#2A supplements

#########TYPE OF MUTATION BY SIGNATURE TYPE

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/fig2aSupplement.tsv', sep='\t', header=TRUE)

p <- ggplot(df, aes(x=dominantSignatureAdj, fill=variable, y = frac))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle=90))+
  scale_fill_viridis_d()+
  emptyTheme+
  xlab('signature aetiology')+
  ylab('fraction of drivers')+
  labs(caption='make_type_of_genes_mutated_figure2a.ipynb')

ggsave('~/Desktop/plot.pdf', plot=p,  width = 4, height = 3, units = c("in"))

###############################
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
#S1


df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/plotGeneFractions.tsv', sep = '\t', header=TRUE)

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


ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 3, units = c("in"))

#
####
###########
######################
#############
#########
####
#

#Related unrelated obs Expected supplement
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/obsExpectedRelatedGenesSup.tsv', sep = '\t', header=TRUE)
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

ggsave('~/Desktop/plot.pdf', plot=p,  width = 4, height = 3.5, units = c("in"))



#
######
################
##########################
#####################################
##########################
################
#######
#
#2C (i)

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/craigStylePlotPOLE.tsv', sep = '\t', header=TRUE)

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


#2C (i)

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/allComps.tsv', sep = '\t', header=TRUE)

l <- list()
for(comp in unique(df$comp)){
  dfHere <- df[df$comp == comp,]
  l[[comp]] <- plot_comp(dfHere)
}

p <- plot_grid(plotlist=l)

ggsave('~/Desktop/plot.pdf', plot=p,  width = 15, height = 10, units = c("in"))


#2C (? test)

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/apobecCompTest.tsv', sep = '\t', header=TRUE)
p <- plot_comp(df, title='APOBEC', thresh=.3)

ggsave('~/Desktop/plot.pdf', plot=p,  width = 4, height = 4, units = c("in"))


#2C (percentages of genes hyper vs non hyper)
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

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/percentComps.tsv', sep = '\t', header=TRUE)

alignedPlot <- plot_grid(plot_supplementary_enrichment_figure(df, 'POLE', ratioThreshold = 75),
                         plot_supplementary_enrichment_figure(df, 'MSI', ratioThreshold = 75),
                         plot_supplementary_enrichment_figure(df, 'APOBEC', ratioThreshold = 75)
                         , ncol=3)

ggsave('~/Desktop/plot.pdf', plot=alignedPlot,  width = 15, height = 4, units = c("in"))

#
####
###########
######################
#############
#########
####

#Inactivation method 

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/inactivationMethod.tsv', sep = '\t', header=TRUE)
hp <- ggplot(df, aes(x=cancerType))+
  geom_bar(aes(y=1, fill='LOH'), stat='identity')+
  geom_bar(aes(y=nHyperComposite/totalHyper, fill='composite'), stat='identity')+
  theme(axis.text.x = element_text(angle=90))+
  emptyTheme+
  ylab('')+
  ggtitle('hypermutated cases')

np <- ggplot(df, aes(x=cancerType))+
  geom_bar(aes(y=1, fill='LOH'), stat='identity')+
  geom_bar(aes(y=nNormalComposite/totalNormal, fill='composite'), stat='identity')+
  theme(axis.text.x = element_text(angle=90))+
  emptyTheme+
  ggtitle('non-hypermutated cases')+
  ylab('fraction')

combinedPlot <- plot_grid(np, hp, ncol=2)

ggsave('~/Desktop/plot.pdf', plot=combinedPlot,  width = 6, height = 4, units = c("in"))


#
####
###########
######################
#############
#########
####

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




