#Written by Noah Friedman

plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/'
library(ggplot2)
library(grid)
require(cowplot)
library(egg)
library(dplyr)
library(data.table); setDTthreads(6)
library(stringr)
library(ggrepel)
library(ggnewscale)

emptyTheme <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())
#
######
###############
########################
##################################
########################
###############
######
#

#multiple scales: https://eliocamp.github.io/codigo-r/2018/09/multiple-color-and-fill-scales-with-ggplot2/

#PLOT FIGURE 3A

make_multiple_samples_plot <- function(df){
  p <- ggplot(df, aes(y=reorder(branchId, -1*orderingVal)))+
    
    #scale 1: cancer_type
    geom_tile(aes(x=-.5, fill=cancerType))+
    scale_fill_viridis_d()+
    new_scale_fill()+
    
    #scale 2: which id you are
    geom_tile(aes(x=.5, fill=reorder(trunkId, -1*orderingVal)))+
    geom_tile(aes(x=1.5, fill=reorder(trunkId, -1*orderingVal)))+
    scale_fill_manual(
      values = c(rep_len(c("gray", "#d3d3d3"), length(unique(df$trunkId))))
    )+
    new_scale_fill()+
    
    #scale 3: second hit TSG
    geom_tile(aes(x=2.25, fill=secondHitTSG), width=.5)+
    scale_fill_manual(values = c('white', '#75816b'))+
    new_scale_fill()+
    
    #scale 4: denovo biallelic TSG
    geom_tile(aes(x=2.75, fill=denovoBiallelicInactivation), width=.5)+
    scale_fill_manual(values = c('white', '#A87D7D'))+
    new_scale_fill()+
    
    #scale 5: parallel evolution
    geom_tile(aes(x=3.25, fill=convergentEvolution), width=.5)+
    scale_fill_manual(values = c('white', '#4682b4'))+
    
    #RELATED vs UNRELATED POINTS   #alpha if statement hides stuff if it is zero https://stackoverflow.com/questions/55262551/ggplot2-point-size-by-numeric-do-not-display-point-when-value-0
    geom_point(aes(x=.25, size=nRelatedTrunk, colour='related', alpha = ((nRelatedTrunk == 0) | (branchNumber != 1)) ))+
    geom_point(aes(x=.75, size=nUnrelatedTrunk, colour='unrelated', alpha = ((nUnrelatedTrunk == 0) | (branchNumber != 1)) ))+
    
    geom_point(aes(x=1.25, size=nRelatedBranch, colour='related', alpha = nRelatedBranch == 0))+
    geom_point(aes(x=1.75, size=nUnrelatedBranch, colour='unrelated', alpha = nUnrelatedBranch == 0))+
    scale_alpha_manual(values = c(1,0))+
    
    #LABELS/prettify
    geom_segment(aes(x=1, xend=1, y=0.5, yend=length(unique(df$branchId)) + .5), colour='black')+
    emptyTheme+
    theme(legend.position = 'none')+
    theme(axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text.x = element_text(angle=90), axis.text.y = element_blank())+
    
    scale_x_continuous(breaks= c(-0.5,.5,1.5, 2.25, 2.75, 3.25), labels=c('cancer type', 'truncal drivers', 'branch drivers', 'second hit TSG', 'de-novo biallelic TSG mutation', 'convergent evolution'), position = "top")
    
  return(p)
}

df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/multipleSampleSummary.tsv', sep='\t', header=TRUE)
p <- make_multiple_samples_plot(df)
#Get legends for alignment
legendCancerType <- get_legend(ggplot(df, aes(x=1, y=1, fill=cancerType))+ geom_tile()+ scale_fill_viridis_d())
legendSecondHitTSG <- get_legend(ggplot(df, aes(x=1, y=1, fill=secondHitTSG))+ geom_tile()+ scale_fill_manual(values = c('white', '#75816b')))
legendDenovoBiallelicTSG <- get_legend(ggplot(df, aes(x=1, y=1, fill=denovoBiallelicInactivation))+ geom_tile()+ scale_fill_manual(values = c('white', '#A87D7D')))
legendConvergentEvolution <- get_legend(ggplot(df, aes(x=1, y=1, fill=convergentEvolution))+ geom_tile()+ scale_fill_manual(values = c('white', '#4682b4')))
legendRelatedUnrelated <- get_legend(ggplot(df, aes(x=1,y=1, size=nRelatedBranch, colour='related')) +geom_point())

combinedLegends <- plot_grid(legendCancerType, legendRelatedUnrelated, ncol=2)
combinedPlot <- plot_grid(p, ggplot(), combinedLegends, nrow=3, rel_heights = c(1, .05, .25))
ggsave('~/Desktop/plot.pdf', plot=combinedPlot,  width = 5, height = 10, units = c("in"))

df$trunkId
#
######
###############
########################
##################################
########################
###############
######
#
#PLOT FIGURE 3B



#
######
###############
########################
##################################
########################
###############
######
#
#PLOT FIGURE 3C

#Essential genes
plot_essential_gene_mutagenesis <- function(df){
  p <- ggplot(df)+
      #geom_smooth(aes(x=score, y=nDoubleTruncating/(geneSize*nCancerType), colour='truncating'))+
      #geom_smooth(aes(x=score, y=nDoubleMissense/(geneSize*nCancerType), colour='missense'))+
      #geom_smooth(aes(x=score, y=nDoubleSilent/(geneSize*nCancerType), colour='silent'))
    geom_histogram(x=score, y=nTruncating/(geneSize*nCancerType))
    #geom_smooth(aes(x=score, y=nTruncating/(geneSize*nCancerType), colour='truncating'))+
    #geom_smooth(aes(x=score, y=nMissense/(geneSize*nCancerType), colour='missense'))+
    #geom_smooth(aes(x=score, y=nSilent/(geneSize*nCancerType), colour='silent'))
  return(p)
}

df <- read.table('~/Desktop/WORK/dataForLocalPlotting/essentialGeneInfo.tsv', sep='\t', header=TRUE)
p <- plot_essential_gene_mutagenesis(df)
ggsave('~/Desktop/plot.pdf', plot=p,  width = 5, height = 5, units = c("in"))

#
######
###############
########################
##################################
########################
###############
######
#
#PLOT FIGURE 3D

#
######
###############
########################
##################################
########################
###############
######
#
#PLOT FIGURE 3E

