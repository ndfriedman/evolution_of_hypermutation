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

#adjust this as needed
plottingDataPath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure3/FIGURE3_PLOTTING_FILES/plotDataFiles/'
plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure3/FIGURE3_PLOTTING_FILES/figurePdfs/'

#
######
###############
########################
##################################
########################
###############
######
#

#TRUNCATING MUTATIONS IN ESSENTIAL GENES

plot_figure_3z1 <- function(df){
  p <- ggplot(df, aes(x=nNeutralTruncating))+
    geom_smooth(aes(y=nEssentialTruncatingGenes, color='Truncating mutations in essential genes'))+
    geom_smooth(aes(y=nEssentialDoubleTruncatingGenes, color='Double Truncating mutations in essential genes'))+
    geom_smooth(aes(y=nNeutralDoubleTruncatingGenes, color='Double Truncating mutations in neutral genes'))+
    theme_classic()+
    ylab('n genes with truncating mutations')+
    xlab('n truncating mutations')
  return(p)
}

#Figure 3z1
figure3z1DataFrame <- read.table(paste(plottingDataPath, 'figure_3z1.tsv', sep=''), sep = '\t', header=TRUE)
#We remove the top decile of most hypermutated 
plt3z1 <- plot_figure_3z1(figure3z1DataFrame)
saveFilePath = paste(plottingFilePath, 'figure3z1.pdf')
ggsave(saveFilePath,
       plot=plt3z1,  width = 5.5, height = 4, units = c("in"), limitsize = FALSE)

#
######
###############
########################
##################################
########################
###############
######
#

#SUMMARIZING CONTINUED TUMOR EVOLUTION PLOT

plot_figure_3z2 <- function(df){
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
  
  legendCancerType <- get_legend(ggplot(df, aes(x=1, y=1, fill=cancerType))+ geom_tile()+ scale_fill_viridis_d())
  legendSecondHitTSG <- get_legend(ggplot(df, aes(x=1, y=1, fill=secondHitTSG))+ geom_tile()+ scale_fill_manual(values = c('white', '#75816b')))
  legendDenovoBiallelicTSG <- get_legend(ggplot(df, aes(x=1, y=1, fill=denovoBiallelicInactivation))+ geom_tile()+ scale_fill_manual(values = c('white', '#A87D7D')))
  legendConvergentEvolution <- get_legend(ggplot(df, aes(x=1, y=1, fill=convergentEvolution))+ geom_tile()+ scale_fill_manual(values = c('white', '#4682b4')))
  legendRelatedUnrelated <- get_legend(ggplot(df, aes(x=1,y=1, size=nRelatedBranch, colour='related')) +geom_point())
  
  combinedLegends <- plot_grid(legendCancerType, legendRelatedUnrelated, ncol=2)
  combinedPlot <- plot_grid(p, ggplot(), combinedLegends, nrow=3, rel_heights = c(1, .05, .25))
  return(combinedPlot)
}

figure3z2DataFrame <- read.table(paste(plottingDataPath, 'figure_3z2.tsv', sep=''), sep = '\t', header=TRUE)
#We remove the top decile of most hypermutated 
plt3z2 <- plot_figure_3z2(figure3z2DataFrame)
saveFilePath = paste(plottingFilePath, 'figure3z2.pdf')
ggsave(saveFilePath,
       plot=plt3z2,  width = 5, height = 10, units = c("in"), limitsize = FALSE)


#
######
###############
########################
##################################
########################
###############
######
#
#CLONALITY PLOT

plot_figure_3z3 <- function(df){
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
    emptyTheme
  return(p) 
}


figure3z3DataFrame <- read.table(paste(plottingDataPath, 'figure_3z3.tsv', sep=''), sep = '\t', header=TRUE)
plt3z3 <- plot_figure_3z3(figure3z3DataFrame)
saveFilePath = paste(plottingFilePath, 'figure3z3.pdf')
ggsave(saveFilePath,
       plot=plt3z3,  width = 5, height = 10, units = c("in"), limitsize = FALSE)




