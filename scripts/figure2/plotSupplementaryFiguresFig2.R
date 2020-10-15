#written by Noah Friedman

library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(ggpubr)
library(ggrepel)

plottingDataPath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure2/FIGURE2_PLOTTING_FILES/plotDataFiles/'
plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure2/FIGURE2_PLOTTING_FILES/figurePdfs/'

emptyTheme <- theme(axis.line = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())
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
#Figure S2 (a) indel and stop gain prevalence by signature type

plot_figure_s2_a <- function(df){
  p <- ggplot(df, aes(x=dominantSignatureAdj, fill=variable, y = frac))+
    geom_bar(stat='identity')+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    scale_fill_viridis_d()+
    xlab('signature aetiology')+
    ylab('fraction of drivers')+
    ggtitle('S2(a.)')
  return(p)
}

figureS2aDf <- read.table(paste(plottingDataPath, 'figureS2_a.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_figure_s2_a(figureS2aDf)
saveFilePath = paste(plottingFilePath, 'figureS2_a.pdf')
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
#Figure S2 (b) Attribution of mutations to signatures.
plot_figure_s2_b <- function(df){
  
  plot_pie_chart <- function(df, title){
    ggplot(df, aes(x=1, fill=hypermutationInduced))+
      geom_bar()+
      coord_polar("y", start=0)+
      emptyTheme+
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank())+
      xlab('')+ylab('')+
      scale_fill_manual(values=c('black', 'gray', '#d3d3d3'))+
      ggtitle(title)
  }
  
  p1 <- plot_pie_chart(df[df$geneType == 'VUS',], 'VUS')
  p2 <- plot_pie_chart(df[df$geneType == 'Oncogene',], 'Oncogene')
  p3 <- plot_pie_chart(df[df$geneType == 'TSG_missense',], 'TSG missense')
  p4 <- plot_pie_chart(df[df$geneType == 'TSG_truncating',], 'TSG truncating')
  leg <- get_legend(p1)
  p1 <- p1 + theme(legend.position = 'none')
  p2 <- p2 + theme(legend.position = 'none')
  p3 <- p3 + theme(legend.position = 'none')
  p4 <- p4 + theme(legend.position = 'none')
  alignedPlot <- plot_grid(p1, p2, p3, p4)
  legPlot <- plot_grid(alignedPlot, leg, rel_widths = c(1,.5))
  finalPlot <- plot_grid(ggplot()+ggtitle('S2(b.)'), legPlot, nrow=2, rel_heights = c(.1,1))
  return(finalPlot)
}

figureS2bDf <- read.table(paste(plottingDataPath, 'figureS2_b.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_figure_s2_b(figureS2bDf)
saveFilePath = paste(plottingFilePath, 'figureS2_b.pdf')
ggsave(saveFilePath, plot=p,  width = 6, height = 6)

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
#Figure S2 (c) All possible drivers

plot_figure_s2_c <- function(df){
  pAll <- ggplot(df, aes(x=1, fill=mutationType, y=number))+
            geom_bar(aes(y=number), stat='identity')+
            emptyTheme+
            ylab('n possible mutations')+
            ggtitle('non-\nsynonymous')+
            xlab('')+
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
  pDrivers <- ggplot(df[df$mutationType != 'VUS',], aes(x=1, fill=mutationType, y=number))+
    geom_bar(aes(y=number), stat='identity')+
    emptyTheme+
    ggtitle('driver')+
    xlab('')+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    ylab('n possible mutations')+
    theme(legend.position = 'none')
  
  leg <- get_legend(pAll)
  pAll <- pAll + theme(legend.position = 'none')
  alignedPlot <- plot_grid(pAll, pDrivers, leg, ncol=3, rel_widths = c(1,1,.75))
  finalPlot <-  plot_grid(ggplot()+ggtitle('S2(c.)'), alignedPlot, nrow=2, rel_heights = c(.1,1))
  return(finalPlot)
}

figureS2cDf <- read.table(paste(plottingDataPath, 'figureS2_c.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_figure_s2_c(figureS2cDf)
saveFilePath = paste(plottingFilePath, 'figureS2_c.pdf')
ggsave(saveFilePath, plot=p,  width = 5, height = 4)

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
#Figure S2 (d) Fraction drivers in related genes

plot_figure_s2_d <- function(df){
  p <- ggplot(df[(df$nUnrelatedExpected > .01),],
            aes(x=nmut_IM341, colour=dominantSignature))+
  geom_smooth(aes(y=(nRelatedObserved)/(nRelatedObserved + nUnrelatedObserved),  linetype='unrelated'), method='loess', span=1, se=FALSE)+
  scale_x_log10()+
  theme_classic()+
  ylab('fraction drivers in related genes')+
  xlab('mutations in IM-341 genes')+
  ylim(0,1)+
  scale_colour_manual(values = c("#FF0000","#267574",
                                 'gray', "#ADFF2F", '#2A52BE', "#FFF600"))+
  ggtitle('S2(d.)')
  return(p)
}

figureS2dDf <- read.table(paste(plottingDataPath, 'figureS2_d.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_figure_s2_d(figureS2dDf)
saveFilePath = paste(plottingFilePath, 'figureS2_d.pdf')
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
#Figure S2 (e) DNDS significance values across 'related' and 'unrelated genes'

plot_figure_s2_e <- function(df){
  p <- ggplot(df[(df$CANCER_TYPE == 'Glioma') | (df$CANCER_TYPE == 'Colorectal Cancer') 
                 | (df$CANCER_TYPE == 'Endometrial Cancer'), ],
              aes(x=cancerTypeGene, y=dndsIsSignificantHyper, fill=CANCER_TYPE))+
    stat_summary(geom='bar')+
    ylim(0, 1)+
    theme(axis.text.x = element_text(angle=90))+
    ylab('Fraction genes signficant by DNDS')+
    emptyTheme+
    ggtitle('S2(e)')
  return(p)
}

figureS2eDf <- read.table(paste(plottingDataPath, 'figureS2_e.tsv', sep=''), sep='\t', header=TRUE)
p <- plot_figure_s2_e(figureS2eDf)
saveFilePath = paste(plottingFilePath, 'figureS2_e.pdf')
ggsave(saveFilePath, plot=p,  width = 4, height = 5, units = c("in"))

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
#Figure S2 (f) Gene mutation enrichment in hypermutated tumors

plot_figure_s2_f <- function(df){
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
  finalPlot <- plot_grid(ggplot()+ggtitle('S2(f)'), alignedPlot, nrow=2, rel_heights = c(.1,1))
  return(finalPlot)
}

figureS2fDf <- read.table(paste(plottingDataPath, 'figureS2_f.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_f(figureS2fDf) #FYI trying to display the variable p will cause an R viewport error (it saves fine though)
saveFilePath = paste(plottingFilePath, 'figureS2_f.pdf')
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
plot_figure_s2_g <- function(df){
  p <- ggplot(df, aes(x= deviation, fill=mutationType, color=mutationType))+
    geom_histogram(color='black')+
    geom_text_repel(data=df[df$deviation > max(df$deviation) - 5,], 
                    aes(label=allele, y=0), nudge_y=400, force=20)+
    theme_classic()+
    scale_colour_viridis_d(drop=FALSE)+
    scale_fill_viridis_d(drop=FALSE)+
    ylab('n distinct observed mutations')+
    xlab('standard deviations above mean n mutations\nat the same pentanucleotide context')+
    ggtitle('S2(g)')
  return(p)
}


#Figure S2 (g) mutation enrichment and pole pentanucleotides
figureS2gDf <- read.table(paste(plottingDataPath, 'figureS2_g.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_g(figureS2gDf) 
saveFilePath = paste(plottingFilePath, 'figureS2_g.pdf')
ggsave(saveFilePath, plot=p,  width = 6, height = 4, units = c("in"))


#
######
################
##########################
#####################################
##########################
################
#######
#
#Figure S2 (h) cumulative hotspots explained by related and unrelated mutations
plot_figure_s2_h <- function(df){
  plot_prevalence_curve <- function(df, title){
    p <- ggplot(df, aes(color=class, x=percentile, y=val))+
      theme_classic()+
      geom_path()+
      xlab('nth most recurrent hotspot')+
      theme(axis.ticks.x = element_blank(), axis.text.x=element_blank())+
      ylab('percent of all hotspot mutations')+
      ggtitle(title)
    return(p)
  }
  endoPlot <- plot_prevalence_curve(df[df$cancerType == 'Endometrial_Cancer',], 'Endometrial_Cancer')
  coloPlot <- plot_prevalence_curve(df[df$cancerType == 'Colorectal_Cancer',], 'Colorectal_Cancer')
  gliomaPlot <- plot_prevalence_curve(df[df$cancerType == 'Glioma',], 'Glioma')
  leg <- get_legend(endoPlot)
  #Remove legends and add a significance asterisk 
  endoPlot <- endoPlot + theme(legend.position = 'none') + geom_text(aes(x=.5, y=.75), label='*', color='black')
  coloPlot <- coloPlot + theme(legend.position = 'none') + geom_text(aes(x=.5, y=.75), label='*', color='black')
  gliomaPlot <- gliomaPlot + theme(legend.position = 'none') + geom_text(aes(x=.5, y=.75), label='ns', color='black')
  alignedPlot <- plot_grid(endoPlot, coloPlot, gliomaPlot, leg, ncol=4, rel_widths = c(1,1,1,.5))
  return(alignedPlot)
}

figureS2hDf <- read.table(paste(plottingDataPath, 'figureS2_h.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_h(figureS2hDf)
plotWithTitle <- plot_grid(ggplot()+ggtitle('S2(h.)'), p, nrow=2, rel_heights = c(.1, 1))
saveFilePath = paste(plottingFilePath, 'figureS2_h.pdf')
ggsave(saveFilePath, plot=plotWithTitle,  width = 10, height = 4, units = c("in"))

#
######
################
##########################
#####################################
##########################
################
#######
#
#Figure S2 (i) Comparing genes mutated in cancers of identical signature exposure

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

plot_figure_s2_i <- function(df){
  l <- list()
  for(comp in unique(df$cohort)){
    dfHere <- df[df$cohort == comp,]
    l[[comp]] <- plot_comp(dfHere, comp)
  }
  p <- plot_grid(plotlist=l)
  finalP <- plot_grid(ggplot()+ggtitle('S2(i)'), p, nrow=2, rel_heights = c(.1,1))
  return(finalP)
}

figureS2iDf <- read.table(paste(plottingDataPath, 'figureS2_i.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_i(figureS2iDf) #FYI trying to display the variable p will cause an R viewport error (it saves fine though)
saveFilePath = paste(plottingFilePath, 'figureS2_i.pdf')
ggsave(saveFilePath, plot=p,  width = 15, height = 15, units = c("in"))

#
######
################
##########################
#####################################
##########################
################
#######
#
#Figure S2 (j) fraction biallelic inactivation

plot_figure_s2_j <- function(df){
  p <- ggplot(df, aes(x=cancerType))+
    geom_bar(aes(y=1, fill='LOH', colour=type), position='dodge2', stat='identity')+
    geom_bar(aes(y=nComposite/total, fill='composite', colour=type), position='dodge2', stat='identity')+
    scale_colour_manual(values=c('black', 'white'))+
    theme(axis.text.x = element_text(angle=90))+
    emptyTheme+
    ylab('fraction of all biallelic inactivation')+
    ggtitle('S2(j.)')
  return(p)
}

figureS2jDf <- read.table(paste(plottingDataPath, 'figureS2_j.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_j(figureS2jDf) #FYI trying to display the variable p will cause an R viewport error (it saves fine though)
saveFilePath = paste(plottingFilePath, 'figureS2_j.pdf')
ggsave(saveFilePath, plot=p,  width = 5, height = 5, units = c("in"))


#
######
################
##########################
#####################################
##########################
################
#######

#Figure S2 (k) Permutation and composite mutation

plot_figure_s2_k <- function(df){
  pBars <- ggplot(df[df$pVal < .05,], aes(x=cancerType, fill=geneType))+
    geom_bar()+
    theme_classic()+
    theme(axis.text.x = element_text(angle=90))+
    ylab('n genes')+
    ggtitle('genes enriched for composite mutation\nvia permutation test')
  pCounts <- ggplot(df[(df$pVal < .05) & (df$geneType == 'oncogene'),], aes(x=1, y=nObs, fill=geneCancerType))+
                geom_bar(stat='identity')+
                coord_polar("y", start=0)+
                emptyTheme+
                theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
                ylab('N cases with composite')+
                xlab('')+
                ggtitle('Significant Oncogenes')
                
  alignedPlot <- plot_grid(pBars, pCounts, ncol=2, rel_widths=c(.75,1)) 
  finalPlot <- plot_grid(ggplot()+ggtitle('S2(k.)'), alignedPlot, nrow=2, rel_heights = c(.1,1))
  return(finalPlot)
}

figureS2kDf <- read.table(paste(plottingDataPath, 'figureS2_k.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_k(figureS2kDf) #FYI trying to display the variable p will cause an R viewport error (it saves fine though)
saveFilePath = paste(plottingFilePath, 'figureS2_k.pdf')
ggsave(saveFilePath, plot=p,  width = 9, height = 4, units = c("in"))


#
####
###########
######################
#############
#########
####
#Figure S2 (l) Double mutation violating the infinite sites hypothesis

plot_figure_s2_l <- function(df){

  fun.1 <- function(x) x^2 

  p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))+
    geom_point(data=df, aes(x=nAllele/nGene, y=nDoubleHit/nBiallelicLoss, size=nDoubleHit, colour=geneType))+
    geom_text_repel(data=df[df$nDoubleHit >= 4,], aes(x=nAllele/nGene, y=nDoubleHit/nBiallelicLoss, label=allele), force=5)+
    xlab('Fraction of all drivers in gene\ncaused by double hit allele')+
    ylab('Fraction of all biallelic inactivation of gene\nattributable to double hit')+
    ylim(0,1)+
    xlim(0,1)+
    stat_function(fun = fun.1, linetype='dashed')+
    theme_classic()+
    ggtitle('s2(l.)')
  return(p)
}


figureS2lDf <- read.table(paste(plottingDataPath, 'figureS2_l.tsv', sep=''), sep = '\t', header=TRUE)
p <- plot_figure_s2_l(figureS2lDf) #FYI trying to display the variable p will cause an R viewport error (it saves fine though)
saveFilePath = paste(plottingFilePath, 'figureS2_l.pdf')
ggsave(saveFilePath, plot=p,  width = 5, height = 5, units = c("in"))

#
####
###########
######################
#############
#########
####

