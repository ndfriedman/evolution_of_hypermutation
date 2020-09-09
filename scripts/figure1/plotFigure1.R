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

#TO BEAUTIFY OUR PLOTS WE ADD "EMPTY THEME"
emptyTheme <- theme(axis.line = element_blank(),
                    #axis.text.x = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())
#
###
######
###########
#######
###
#

plottingFilePath = '/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/'

#
###
######
###########
#######
###
#

#PLOT FIGURE 1A

plot_n_cases_figure <- function(df){
  bottomScaleFactor <- 10
  p <- ggplot(df, aes(x=reorder(label, fracHypermutatedOrdering)))+
    geom_bar(aes(y=nTotal, fill='Total Cases'), stat='identity')+
    #optional include the n high mut burden as well
    geom_bar(aes(y=-1*bottomScaleFactor*nHypermutated - 1*bottomScaleFactor*nHighMutBurden, fill='Hypermutated'), stat='identity')+
    
    geom_bar(aes(y=-1*bottomScaleFactor*nHypermutated, fill='Hypermutated'), stat='identity')+
    
    theme(axis.text.x = element_text(angle=90))+
    theme(axis.ticks.x = element_blank())+
    scale_fill_manual(values=c('#858585', 'light gray'))+
    #THIS IS MANUALLY SET YOU NEED TO CHANGE THE SCALE FACTOR IF YOU CHANGE IT
    scale_y_continuous(breaks= c(-4000, -2000, 0, 2000, 4000, 6000), labels=c('400', '200', '0', '2000', '4000', '6000'))+
    emptyTheme+
    guides(fill=guide_legend(title="Tumor Classification"))+
    xlab('Cancer Type')+
    ylab('N hypermutated                    N total')+
    coord_flip()
  return(p)
}

plot_percent_cases_figure <- function(df){
  bottomScaleFactor <- 10
  p <- ggplot(df, aes(x=reorder(label, fracHypermutatedOrdering)))+
    geom_bar(aes(y=fracHypermutated), stat='identity')+
    
    theme(axis.text.x = element_text(angle=90))+
    theme(axis.ticks.x = element_blank())+
    scale_fill_manual(values=c('#858585', 'light gray'))+
    emptyTheme+
    xlab('Cancer Type')+
    ylab('Fraction of cases hypermutated')+
    coord_flip()
  return(p)
}

#figure1bDataFrame <- read.table(paste(plottingFilePath, 'figure1bCancerTypeSummary.tsv', sep=''), sep='\t', header=TRUE)
#figure1bDataFrame <- read.table(paste(plottingFilePath, 'figure1bCancerTypeSummary.tsv', sep=''), sep='\t', header=TRUE)
figure1aDataFrame <- read.table('/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/figure1aCancerTypeSummary.tsv', sep='\t', header=TRUE)

p <- plot_percent_cases_figure(figure1aDataFrame)
#p <- plot_n_cases_figure(figure1bDataFrame)
ggsave('~/Desktop/plot.pdf', plot=p,  width = 6, height = 4)

#
###
######
###########
#######
###
#

#PLOT FIGURE 1B

plot_signatures_figure <- function(df){
  p <- ggplot(df, aes(x=1, fill=signature, y=nHyperHigh/sum(df$nHyperHigh)))+
    geom_bar(stat='identity')+
    emptyTheme+
    scale_fill_manual(values=c("#FF0000","#267574",
                               'gray', "#ADFF2F","#00dd5d", '#ffb347', '#2A52BE', "#FFF600"))+
    xlab('')+
    theme(axis.text.x = element_blank())+
    ylab('fraction of all hypermutated cases')
  return(p)
}


figure1bDataFrame <- read.table('/Users/friedman/Desktop/hypermutationProjectFinal/scripts/figure1/FIGURE1_PLOTTING_FILES/figure1bSignatureSummary.tsv',
                                sep='\t', header=TRUE)

p <- plot_signatures_figure(figure1bDataFrame)
ggsave('~/Desktop/plot.pdf',
       plot=p,  width = 3, height = 5)


#
###
######
###########
#######
###
#

#PLOT FIGURE 1C

plot_data <- function(df){
  #plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nHotspots))+
  plt <- ggplot(df, aes(x=reorder(cohort, orderingVal), y=nOncMuts))+
    #geom_boxplot(fatten = NULL, outlier.shape=NA)+
    geom_violin(bw = 1, aes(fill=factor(cohort,
                                  levels = c('normal_Colorectal', 'hyper_Colorectal', 'normal_Endometrial', 'hyper_Endometrial',
                                  'normal_Glioma', 'hyper_Glioma', 'normal_Other', 'hyper_Other'))))+
    stat_summary(fun.y = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.75, size = 1, linetype = "solid")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
    #scale_colour_manual(values =  c('black', "#267574", 'gray', '#ADFF2F', "#9acd32", '#2A52BE'), name="Dominant\nSignature")+
    ylab('N Driver Mutations')+
    #ylab('N hotspot mutations')+
    xlab('Cancer Type')+
    emptyTheme+
    coord_cartesian(ylim=c(0,50))+
    #geom_jitter(aes(colour=factor(cohort,
    #                              levels = c('normal_Colorectal', 'hyper_Colorectal', 'normal_Endometrial', 'hyper_Endometrial',
    #                                         'normal_Glioma', 'hyper_Glioma', 'normal_Other', 'hyper_Other'))),
    #            shape=16, position=position_jitter(0.1), alpha=0.75)+
    
    scale_fill_manual(values =c('orange', '#b36200', 'lavender', '#301934', '#add8e6', 'blue', 'gray', '#333333'), name='Cohort')
  
  return(plt)
}

figure1cDataFrame <- read.table(paste(plottingFilePath, 'figure1c_nOncMutByCohort.tsv', sep=''), sep='\t', header=TRUE)
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


#Make plot
plt <- plot_data(figure1cDataFrame)
ggsave('~/Desktop/plot.pdf',
       plot=plt,  width = 6, height = 4, units = c("in"), limitsize = FALSE)

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
    emptyTheme
  return(p)
}

dfObsExp <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/figure1d.tsv', sep = '\t', header=TRUE)
p <- plot_figure_1d(dfObsExp)
ggsave('~/Desktop/plot.pdf',
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
    emptyTheme+
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
    emptyTheme+
    ggtitle('MSI only')+
    xlim(0,200)
  return(p)
}

#plot figure 1e
df <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/indelRateInfo.tsv', sep = '\t', header=TRUE)
p <- plot_figure_1e(df)
ggsave('~/Desktop/plot.pdf',
       plot=p,  width = 3, height = 4, units = c("in"), limitsize = FALSE)



#
############
###############################
#######DEPRECATED PLOTS


make_obs_exp_plot <- function(includeLegend = TRUE){
  plt <- ggplot()+
    
    stat_summary_bin(data=dfObsExp[dfObsExp$dominantSignature == 'mean_MMR',], aes(x=nmut, y=obsHotspot, colour='Observed_MMR'), bins=10)+
    stat_summary_bin(data=dfObsExp[dfObsExp$dominantSignature == 'mean_MMR',], aes(x=nmut, y=expectedHotspot, colour='Expected_MMR'), bins=10)+
    geom_smooth(data=dfObsExp[dfObsExp$dominantSignature == 'mean_MMR',], aes(x=nmut, y=obsHotspot, colour='Observed_MMR'), method='loess', se=FALSE, span = 25)+
    geom_smooth(data=dfObsExp[dfObsExp$dominantSignature == 'mean_MMR',], aes(x=nmut, y=expectedHotspot, colour='Expected_MMR'), method='loess', se=FALSE, span =25)+
    
    stat_summary_bin(data=dfObsExp[dfObsExp$dominantSignature == 'mean_10',], aes(x=nmut, y=obsHotspot, colour='Observed_POLE'), bins=5)+
    stat_summary_bin(data=dfObsExp[dfObsExp$dominantSignature == 'mean_10',], aes(x=nmut, y=expectedHotspot, colour='Expected_POLE'), bins=5)+
    geom_smooth(data=dfObsExp[dfObsExp$dominantSignature == 'mean_10',], aes(x=nmut, y=obsHotspot, colour='Observed_POLE'), method='loess', se=FALSE, span = 25)+
    geom_smooth(data=dfObsExp[dfObsExp$dominantSignature == 'mean_10',], aes(x=nmut, y=expectedHotspot, colour='Expected_POLE'), method='loess', se=FALSE, span =25)+
    
    stat_summary_bin(data=dfObsExp[(dfObsExp$dominantSignature == 'mean_11') & (dfObsExp$nmut < 1150),], aes(x=nmut, y=obsHotspot, colour='Observed_TMZ'), bins=5)+
    stat_summary_bin(data=dfObsExp[(dfObsExp$dominantSignature == 'mean_11') & (dfObsExp$nmut < 1150),], aes(x=nmut, y=expectedHotspot, colour='Expected_TMZ'), bins=5)+
    geom_smooth(data=dfObsExp[(dfObsExp$dominantSignature == 'mean_11') & (dfObsExp$nmut < 1150),], aes(x=nmut, y=obsHotspot, colour='Observed_TMZ'), method='loess', se=FALSE, span = 25)+
    geom_smooth(data=dfObsExp[(dfObsExp$dominantSignature == 'mean_11') & (dfObsExp$nmut < 1150),], aes(x=nmut, y=expectedHotspot, colour='Expected_TMZ'), method='loess', se=FALSE, span =25)+
    
    xlab('N Nonsynonymous Mutations\nin IMPACT 341 Genes')+
    scale_x_continuous(breaks=c(0,100,200,400))+
    ylab('N hotspots per case')+
    #scale_x_log10()+
    #ggtitle('Observed and expected hotspot\nburden in hypermutated tumors')+
    emptyTheme+
    #scale_color_manual(values=c('gray', 'black'))+
    scale_color_manual(values=c("#9CA89C",
                                "gray", '#6699cc',"#267574","#ADFF2F", '#2A52BE', "#FFF600"))+
    labs(colour = "Number of Hotspots:")
  if(includeLegend == FALSE){
    plt <- plt + theme(legend.position = 'none')
  }
  return(plt)
}

#dfObsExp <- read.table(paste(plottingFilePath, 'figure1e_observedVsExpected.tsv', sep=''), sep = '\t', header=TRUE)
dfObsExp <- read.table('/Users/friedman/Desktop/WORK/dataForLocalPlotting/figure1d.tsv', sep = '\t', header=TRUE)



my_comparisons <- list( c("nEssential", "nEssentialExp"), c("nOncogene", "nOncogeneExp"), c("nTsg", "nTsgExp") )
p <- ggplot(df, aes(x=variable, y=value))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle=90))+
  ylab('N Mutations')+
  xlab('Mutation type')+
  emptyTheme+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(28, 28, 28))+
  stat_compare_means(comparisons = my_comparisons)+
  ggtitle('Observed vs expected truncating mutations\nbased on exome mutation rates')+
  scale_y_log10()
  #coord_cartesian(ylim=c(0,30))
  
ggsave('~/Desktop/plt.pdf', plot=p,  width = 4, height = 5, units = c("in"))

s = 1
p <- ggplot()+
  #geom_smooth(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_mmr'), method = 'lm', se=FALSE)+
  #geom_smooth(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_mmr'), method = 'lm', se=FALSE)+
  #geom_smooth(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_mmr'), method = 'lm', se=FALSE)+
  
  geom_smooth(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_pole'), method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_pole'),  method = 'lm', se=FALSE)+
  geom_smooth(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_pole'),  method = 'lm', se=FALSE)+
  
  #stat_summary_bin(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_mmr'))+
  #stat_summary_bin(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_mmr'))+
  #stat_summary_bin(data = df[(df$signature == 'MMR'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_mmr'))+
  
  stat_summary_bin(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nOncogene - nOncogeneExp, colour='oncogene_pole'), bins=5)+
  stat_summary_bin(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nTsg - nTsgExp, colour='tsg_pole'), bins=5)+
  stat_summary_bin(data = df[(df$signature == 'POLE'),], aes(x = tmb, y=nEssential - nEssentialExp, colour='essential_pole'), bins=5)+
  
  
  xlab('TMB')+
  ylab('Difference between observed and\nexpected n truncating muts')+
  scale_color_manual(values=c('#013220', '#037D50', '#FF6347',
                             '#FF8C00', '#82CFFD', '#0D4F8B'))+
  labs(colour = 'Signature and\ngene type')+
  emptyTheme

ggsave(paste(plottingFilePath, 'figure1f.pdf', sep=''),
       plot=p,  width = 4, height = 5, units = c("in"))



#
###
######

plt <- make_obs_exp_plot(includeLegend = FALSE)
legend <- get_legend(make_obs_exp_plot(includeLegend = TRUE))

nCasesHistogram <- ggplot(dfObsExp, aes(x=nmut))+
  geom_histogram(bins=100)+
  scale_y_log10()+
  emptyTheme+
  xlab('')+
  ylab('N cases')+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

alignedPlot <- plot_grid(nCasesHistogram, plt, nrow=2, rel_heights = c(.3, 1))
alignedPlotWithLegend <- plot_grid(alignedPlot, legend, ncol=2, rel_widths = c(1, .5))

ggsave(paste(plottingFilePath, 'figure1e.pdf', sep=''),
       plot=alignedPlotWithLegend,  width = 6, height = 5, units = c("in"))

