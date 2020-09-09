#written by Noah Friedman 
import sys
import argparse
import os
import pandas as pd
import numpy as np
import math

from collections import Counter
sys.path.append('/Users/friedman/Desktop/hypermutationProjectFinal/scripts/utilityScripts')

import analysis_utils
import maf_analysis_utils
#import imp
#analysis_utils = imp.load_source('analysis_utils', '/Users/friedman/Desktop/mnt/juno/work/taylorlab/friedman/myUtils/analysis_utils.py')
#maf_analysis_utils = imp.load_source('maf_analysis_utils', '/Users/friedman/Desktop/mnt/juno/work/taylorlab/friedman/myUtils/analysis_utils.py')

def create_facets_dict_key(row):
	return row['Tumor_Sample_Barcode'] + '_' + row['idCol']

def create_facets_clonality_dict(facetsDf):
	#todo make it incorportate patient info too
  	facetsDf['idCol'] = facetsDf.apply(lambda row: str(row['Chromosome']) + '_' + str(row['Start_Position']), axis=1)
  	facetsDf = data_compacting_and_cleaning_util.create_expected_mut_copies_col(facetsDf)
  	facetsDf['naiveClonalStatus'] = facetsDf.apply(lambda row: row['VAF']/row['purity'], axis=1)
  	d = dict()
  	for index, row in facetsDf.iterrows():
  		d[create_facets_dict_key(row)] = (row['naiveClonalStatus'], row['ccf_Mcopies_upper'])
  	return d

def analyze_clonality_across_muts_in_pole_cases(mutsDf, facetsDict, writeDir = ''):
	cases = set(mutsDf['Tumor_Sample_Barcode'])
	for case in cases:
		caseMuts = mutsDf[mutsDf['Tumor_Sample_Barcode'] == case]
		data = []
		if caseMuts.shape[0] > 0:
			for index, row in caseMuts.iterrows():
				facetsDictKey = create_facets_dict_key(row)
				if facetsDictKey in facetsDict:
					ccf, ccf_Mcopies = facetsDict[facetsDictKey]
					if not np.isnan(float(ccf_Mcopies)):
						data.append(ccf_Mcopies)
			if len(data) > 0:
				figureTitle = case + '; NMuts: ' + str(caseMuts.shape[0])
				histogram_util.plot_simple_normed_histogram(data, case, writeDir, figureTitle)
			print '_________'


#a heuristic for testing if a mutation is a double hit of the same location
#assumes that a mutation that occurs in a balanced region with 
#TODO make sure that we dont miscall variants in a small founder clone as doubles

#TODO CHANGE THIS TO CALL DOUBLE MUTATIONS BASED ON SAMPLE PURITY
def is_mut_double_hit(row, 
	flatGenome, #BOOLEAN to tell us if facets said the case's genome is flat
 doubleFactor=2): #factor by which a mutations vaf needs to be bigger than the median to be considered double.  

	if row['Chromosome'] == 'X': return False #DONT CALL A MUTATION ON X AS DOUBLE BY MISTAKE

	if (math.isnan(row['ccf_Mcopies']) or math.isnan(row['ccf_1copy'])) and not flatGenome: #Null values for ccf only matter if the genome isnt flat
		return False
	else:
		if not flatGenome and row['tcn'] - row['lcn'] != row['lcn']: #if the region isnt balanced we wont call it a double hit
			return False
		else:
			if row['t_var_freq'] >= 1.75*row['medianClonalVaf']:
				return True

#validate putative double hit mutations by confirming another SNP in the same gene is balanced
def annotate_double_hit_mutations(maf):

	def validate_double_hit(row):
		if row['isDoubleHit'] != True: return None
		if row['otherMutsInGeneExist'] != True: return 'no_other_muts'
		otherMutsMaf = maf[(maf['geneCase'] == row['geneCase']) & (maf['varUuid'] != row['varUuid'])]
		if len(list(set(otherMutsMaf['isBalanced']))) == 1 and list(set(otherMutsMaf['isBalanced']))[0] == True: return 'validated_by_other_snps'
		else: return 'other_SNPs_unbalanced'

	maf['geneCase'] = maf.apply(lambda row: str(row['Hugo_Symbol']) + '_' + str(row['Tumor_Sample_Barcode']), axis=1)
	geneCaseDict = dict(maf['geneCase'].value_counts())
	maf['otherMutsInGeneExist'] = maf['geneCase'].apply(lambda x: True if x in geneCaseDict and geneCaseDict[x] > 1 else False)
	maf['doubleHitValidates'] = maf.apply(lambda row: validate_double_hit(row), axis=1)
	return maf

#MARK mutations that occur at a balanced region of the genome 
#use this for simplified analyses of balanced regions
def mark_mutation_is_balanced(row):
	if row['Chromosome'] == 'X': return False #we dont include X mutations cause its complicated
	
	#if flatGenome: return True #cases with purity = NA, flat genomes are considered to be flat everywhere
	#if (math.isnan(row['ccf_Mcopies']) or math.isnan(row['ccf_1copy'])) #dont analyze regions where Mcopies or 1copy are null (flat genome cases have already been called true)
	#	return False

	if row['tcn'] - row['lcn'] != row['lcn']: #this means the region isnt balanced
		return False
	else: #this means it is
		return True

def mark_maf_with_ccf_for_flat_genomes(maf, gene):
	
	maf = maf[maf['oncogenic'].notnull()]
	print set(maf['oncogenic'])

	for case in set(maf['Tumor_Sample_Barcode']):
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
		caseMafGene = caseMaf[caseMaf['Hugo_Symbol'] == gene]
		print caseMafGene['ccf_Mcopies']




#ALTERNATIVE METHOD for ading clonality information to mafs
#this method probably does not work on non hypermutated cases


from sklearn.cluster import MeanShift
def assign_variants_to_clonal_cluster(vafs, ids):
    
    #mark which clusters returned by the clustering are clonal by iterating over clusters by mean vaf 
    #and returning the clusters once we have at least minClonalMuts mutations accumulated
    def assign_clonal_subclonal_clusters(clusterDf, minClonalMuts = 10):
        l = []
        for cluster in set(df['cluster']):
            clusterDf = df[df['cluster'] == cluster]
            l.append((np.nanmean(clusterDf['vaf']), clusterDf.shape[0], cluster))
        runningMutSum = 0
        clonalClusters = []
        for meanVaf, nMut, cluster in sorted(l, reverse=True):
            clonalClusters.append(cluster)
            runningMutSum += nMut
            if runningMutSum >= minClonalMuts:
                return clonalClusters
    
    a = np.array(vafs).reshape(-1, 1)
    if a.shape[0] <= 1: return dict() #if there arent any variants please return an empty dict
    
    clustering = MeanShift().fit(a)
    prediction = clustering.predict(a)
    
    #We make a dataframe 
    listOfDicts = []
    la = list(a)
    lp = list(prediction)
    for i in range(0, len(list(a))):
        listOfDicts.append({
            'vaf': la[i], 'cluster': lp[i], 'varUuid': ids[i]
        })
    df = pd.DataFrame(listOfDicts)
    
    minCMut = max(.1*df.shape[0], 10) #at least 10% of mutation in every case are called clonal
    clonalClusters = assign_clonal_subclonal_clusters(df, minClonalMuts = minCMut)
    if clonalClusters != None:
    	df['clonal'] = df['cluster'].apply(lambda x: True if x in clonalClusters else False)
    else:
    	print df.shape
    	return dict()
    return dict(zip(df['varUuid'], df['clonal']))

def add_clonal_calls_to_maf_based_on_vaf_hypermutators(maf, cases=None):

	#DOES ONE OF TWO THINGS:
	#if the ccf_mcopies not null > great than the number of ccf_mcopies null use those calls

	clonalThresh = .8 #the threshold for ccf_mcopies_lower for call a mutation clonal

	if cases == None: cases = set(maf['Tumor_Sample_Barcode'])
	maf['isClonal'] = None
	maf['clonalityMethod'] = None
	cntr = 0
	for case in list(cases):
		cntr += 1
		if cntr%5 == 0: 
			print cntr, len(cases)
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
		if caseMaf.shape[0] > 0:

			if caseMaf[caseMaf['ccf_Mcopies'].notnull()].shape[0] > caseMaf[caseMaf['ccf_Mcopies'].isnull()].shape[0]:
				maf['isClonal'] = maf.apply(lambda row: True if row['ccf_Mcopies_lower'] > clonalThresh and row['Tumor_Sample_Barcode'] == case
				 else False if row['Tumor_Sample_Barcode'] == case and row['ccf_Mcopies_lower'] <= clonalThresh
				 else None if row['Tumor_Sample_Barcode'] == case
				 else row['isClonal'], axis=1)

				maf['clonalityMethod'] = maf.apply(lambda row: 'FACETS_CCF' if row['Tumor_Sample_Barcode'] == case
				else row['clonalityMethod'], axis=1)
			else:
				caseIds = list(caseMaf['varUuid'])
				caseVafs = list(caseMaf[caseMaf['t_var_freq'].notnull()]['t_var_freq'])
				
				isClonalDict = assign_variants_to_clonal_cluster(caseVafs, caseIds)
				maf['isClonal'] = maf.apply(lambda row: 
						isClonalDict[row['varUuid']] if row['Tumor_Sample_Barcode'] == case and row['varUuid'] in isClonalDict 
						else row['isClonal'], axis=1)

				maf['clonalityMethod'] = maf.apply(lambda row: 'VAF_CLUSTERING' if row['Tumor_Sample_Barcode'] == case
				else row['clonalityMethod'], axis=1)
	return maf


def calculate_delta_vaf_across_mutation_pairs(maf, #MAF to analyze
	doDoubleHitCorrection=True, #does a correction to count all cases where nmut copies = 2 as two mutations at the same spot
	genes = None #genes to look at
	):
	

	#FIRST AND FOREMOST, its only relevant to do this analysis at diploid parts of the genome
	maf = maf[maf['tcn'] == 2]

	dDeltaVaf = {}
	if genes == None: 
		genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	cntr = 0
	nCases = len(set(maf['Tumor_Sample_Barcode']))
	for case in set(maf['Tumor_Sample_Barcode']):
		#COUNTING
		print 'analyzing case number ', cntr, ' out of ', nCases
		cntr += 1

		#SET SOME VARIABLES WE WILL NEED FOR THIS ANALYSIS
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
		maxTVarFreq = max(caseMaf['t_var_freq'])
		medianVaf = np.nanmedian(caseMaf['t_var_freq'])
		meanVaf = np.nanmean(caseMaf['t_var_freq'])

		caseMaf['isDouble'] = None #make a dummy column, fill it out if the do double correction is ons
		if doDoubleHitCorrection:
			flatGenome = False
			if caseMaf[caseMaf['ccf_Mcopies'].notnull()].shape[0] == 0: #CATCH ALL THE PURITY=NA flat genomes cases
				flatGenome = True

			caseMaf['isDouble'] = caseMaf.apply(lambda row: is_mut_double_hit(row, medianVaf, flatGenome, 1.75), axis=1)

		for gene in genes:
			geneMaf = caseMaf[caseMaf['Hugo_Symbol'] == gene]

			if geneMaf.shape[0] > 0:

				#DO NOT CONSIDER GENES WHERE THERE IS LOH or the gene is on CHR X
				if geneMaf[geneMaf['lcn'] == 0].shape[0] > 0 or geneMaf[geneMaf['Chromosome'] == 'X'].shape[0] > 0: 
					pass
				else:

					if geneMaf[geneMaf['isDouble'] == True].shape[0] > 0:
						pass
						#ALERT DO WE INCLUDE DOUBLE MUTS AT THE SAME LOCUS OR NOT???
						#diff = max(list(geneMaf['t_var_freq'])))
						#if gene in dDeltaVaf:
						#	dDeltaVaf[gene][0] = dDeltaVaf[gene][0] + [0] #append a zero for double mutations (note that this does not take into account that the double mutations may not occur at the same time!!?>)
						#	dDeltaVaf[gene][1] = dDeltaVaf[gene][1] + [(0.5*max(list(geneMaf['t_var_freq'])))/meanVaf] #take max vaf, divide it by two and compare the ratio to median vaf
						#else:
						#	dDeltaVaf[gene] = [[0],[(0.5*max(list(geneMaf['t_var_freq'])))/meanVaf]] #make it a list of lists
					###TODO add code so we dont do shit like when PTEN IS LOF'D already
					else:
						if geneMaf.shape[0] > 1:
							maxVafToMedianVafRatio = 1.0*max(list(geneMaf['t_var_freq']))/meanVaf

							differences = analysis_utils.calculate_all_pairwise_differences(np.array(list(geneMaf['t_var_freq'])))
							differencesNormed = [1.0*i/maxTVarFreq for i in differences]

							if gene in dDeltaVaf:
								dDeltaVaf[gene][0] = dDeltaVaf[gene][0] + [min(differencesNormed)] # we only take the closest pair if there is more than one pair
								dDeltaVaf[gene][1] = dDeltaVaf[gene][1] + [maxVafToMedianVafRatio]
							else:
								dDeltaVaf[gene] = [[min(differencesNormed)], [maxVafToMedianVafRatio]] # we only take the closest pair if there is more than one pair
	return dDeltaVaf


#TODO add functionality to do percentile for only the biggest mutation as well
def mark_mutations_by_vaf_mut_percentile(maf):

	#only do this at balanced regions of the genome
	maf['isBalanced'] = maf.apply(lambda row: mark_mutation_is_balanced(row), axis=1)
	maf = maf[maf['isBalanced'] == True]

	print maf.shape

	caseDict = {}
	for case in set(maf['Tumor_Sample_Barcode']):
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
		percentiles = dict(zip(caseMaf['varUuid'], caseMaf['t_var_freq'].rank(pct=True)))
		caseDict[case] = percentiles
	maf['caseVafPercentileRank'] = maf.apply(lambda row: caseDict[row['Tumor_Sample_Barcode']][row['varUuid']], axis=1)		
	return maf

def prepare_maf_for_pyclone(maf):

	#IF THE COPY NUMBER OF NAN THAT MEANS ITS BALANCED
	maf['tcn'] = maf['tcn'].apply(lambda x: 2 if math.isnan(x) else x)
	maf['lcn'] = maf['lcn'].apply(lambda x: 1 if math.isnan(x) else x)

	maf['mutation_id'] = maf.apply(lambda row: str(row['Hugo_Symbol']) + '_' + str(row['Chromosome']) + '_' + str(row['Start_Position']), axis=1)
	maf['ref_counts'] = maf['t_ref_count'].apply(lambda x: int(x))
	maf['var_counts'] = maf['t_alt_count'].apply(lambda x: int(x))
	maf['normal_cn'] = 2 #ALERT FIX THIS TO HAVE AN EDGE CASE FOR MEN ON THE X AND Y
	maf['minor_cn'] = maf['lcn'].apply(lambda x: int(x))
	maf['major_cn'] = maf.apply(lambda row: int(row['tcn'] - row['lcn']), axis=1)

	return maf[['Tumor_Sample_Barcode', 'mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn']]

#omnibus function that returns a maf of hypermutated cases and their purities
def assign_clonality_information_for_hypermutated_cases(mafWithClonalityInfo, facetsWhitelist, facetsBlacklist):
	idsToInclude = set(mafWithClonalityInfo[mafWithClonalityInfo['tcn'].notnull()]['Tumor_Sample_Barcode']) | facetsWhitelist - facetsBlacklist
	
	mafWithClonalityInfo = mafWithClonalityInfo[mafWithClonalityInfo['Tumor_Sample_Barcode'].isin(idsToInclude)]

	#NOTE WE HAVE TO DEAL WITH REGIONS THAT ARE NAN FROM FACETs

	mafWithClonalityInfo['FlatGenome'] = mafWithClonalityInfo['Tumor_Sample_Barcode'].apply(lambda x: True if x in facetsWhitelist else False)
	#ADD clonality info to the facets whitelist cases
	mafWithClonalityInfo['tcn'] = mafWithClonalityInfo.apply(lambda row: 2 if row['Tumor_Sample_Barcode'] in facetsWhitelist
		else row['tcn'], axis=1)
	mafWithClonalityInfo['lcn'] = mafWithClonalityInfo.apply(lambda row: 1 if row['Tumor_Sample_Barcode'] in facetsWhitelist
		else row['lcn'], axis=1)
	
	mafWithClonalityInfo['varUuid'] = mafWithClonalityInfo.apply(lambda row: str(row['Chromosome']) + '_' + str(row['Start_Position']), axis=1)

	print 'adding balanced annotation'
	mafWithClonalityInfo['isBalanced'] = mafWithClonalityInfo.apply(lambda row: True if row['tcn'] - row['lcn'] == row['lcn'] else False, axis=1)
	
	print 'adding clonal calls'
	mafWithClonalityInfo = add_clonal_calls_to_maf_based_on_vaf_hypermutators(mafWithClonalityInfo, cases=None)

	return mafWithClonalityInfo
	"""print 'adding median vaf'
	mafWithClonalityInfo = maf_analysis_utils.mark_cases_with_median_clonal_vaf_of_case(mafWithClonalityInfo)
	print mafWithClonalityInfo['medianClonalVaf']
	print 'adding double annotation'
	mafWithClonalityInfo['isDouble'] = mafWithClonalityInfo.apply(lambda row: is_mut_double_hit(row, row['FlatGenome'], 1.75), axis=1)"""




#clonality_analysis_util.create_clonality_whitelist_command(mafWithClonalityInfo, excludeIds=facetsWhitelist | facetsBlacklist)   
def create_clonality_whitelist_command(mafWithClonalityInfo, purityCutoff = .15, excludeIds=set([])):
	print 'copy the following to terminal and run'
	cntr = 0
	for i in set(mafWithClonalityInfo['Tumor_Sample_Barcode'])-set(mafWithClonalityInfo[mafWithClonalityInfo['tcn'].notnull()]['Tumor_Sample_Barcode'])-excludeIds:
		caseMaf = mafWithClonalityInfo[mafWithClonalityInfo['Tumor_Sample_Barcode'] == i]
		purity = 2*np.nanmedian(caseMaf['t_var_freq'])
		if purity > purityCutoff:
			print 'facetsInbound ' + i 
			cntr += 1
	print cntr 

#returns cases I have manaully facets whitelisted
#these are completely CNA flat cases with purity=NA
def get_facets_whitelist():
	return set(
		['P-0000069-T01-IM3', 'P-0001248-T01-IM3', 'P-0001703-T01-IM3', 'P-0003767-T01-IM5', 'P-0004255-T01-IM5',
		'P-0004260-T02-IM5', 'P-0004688-T01-IM5', 'P-0005021-T01-IM5', 'P-0006612-T01-IM5', 'P-0006753-T01-IM5',
		'P-0009157-T01-IM5', 'P-0010308-T01-IM5', 'P-0010393-T02-IM6', 'P-0010499-T01-IM5', 'P-0010671-T01-IM5',
		'P-0011345-T01-IM5', 'P-0011385-T01-IM5', 'P-0012397-T01-IM5', 'P-0012402-T01-IM5', 'P-0012445-T01-IM5',
		'P-0012670-T01-IM5', 'P-0012726-T01-IM5', 'P-0012881-T01-IM5', 'P-0013227-T01-IM5', 'P-0013350-T01-IM5',
		'P-0013400-T01-IM5', 'P-0013676-T01-IM5', 'P-0014258-T01-IM6', 'P-0014787-T01-IM6', 'P-0015885-T01-IM6',
		'P-0016773-T01-IM6', 'P-0017675-T01-IM5', 'P-0017839-T01-IM6', 'P-0017925-T01-IM6', 'P-0018005-T01-IM6',
		'P-0018616-T01-IM6', 'P-0019264-T01-IM6', 'P-0019360-T01-IM6', 'P-0019545-T01-IM6', 'P-0019658-T01-IM6',
		'P-0019871-T01-IM6', 'P-0021090-T01-IM6', 'P-0025554-T01-IM6', 'P-0026278-T01-IM6', 'P-0026523-T01-IM6',
		'P-0026962-T01-IM6', 'P-0028144-T01-IM6', 'P-0029690-T01-IM6', 'P-0029778-T01-IM6', 'P-0032113-T01-IM6',
		'P-0032181-T01-IM6', 'P-0032660-T01-IM6', 'P-0033425-T01-IM6', 'P-0033605-T01-IM6', 'P-0034448-T01-IM6',
		'P-0035281-T01-IM6', 'P-0036123-T01-IM6', 'P-0036500-T01-IM6', 'P-0036503-T01-IM6', 'P-0036568-T01-IM6',
		'P-0036860-T01-IM6', 'P-0037288-T01-IM6', 'P-0037582-T01-IM6', 'P-0003524-T01-IM5', 'P-0004051-T01-IM5',
		'P-0004865-T01-IM5', 'P-0006207-T01-IM5', 'P-0007831-T01-IM5', 'P-0007997-T01-IM5', 'P-0008345-T01-IM5',
		'P-0010504-T01-IM5', 'P-0010828-T01-IM5', 'P-0011540-T01-IM5', 'P-0013537-T01-IM5', 'P-0014388-T01-IM6',
		'P-0014780-T01-IM6', 'P-0016023-T01-IM6', 'P-0016099-T01-IM6', 'P-0016801-T01-IM6', 'P-0017713-T01-IM6',
		'P-0019464-T01-IM6', 'P-0019649-T01-IM6', 'P-0020143-T01-IM6', 'P-0020242-T01-IM6', 'P-0020331-T01-IM6',
		'P-0021077-T01-IM6', 'P-0021572-T01-IM6', 'P-0023271-T01-IM6', 'P-0023555-T01-IM6', 'P-0024633-T01-IM6',
		'P-0025073-T01-IM6', 'P-0026456-T01-IM6', 'P-0027438-T01-IM6', 'P-0032496-T01-IM6', 'P-0034308-T01-IM6',
		'P-0035916-T01-IM6', 'P-0000157-T01-IM3', 'P-0007843-T01-IM5', 'P-0011357-T01-IM5', 'P-0021897-T01-IM6',
		'P-0020295-T01-IM6', 'P-0008646-T01-IM5', 'P-0010803-T01-IM5', 'P-0011226-T01-IM5', 'P-0013876-T01-IM5',
		'P-0015288-T01-IM6', 'P-0017681-T01-IM6', 'P-0018781-T01-IM6', 'P-0019437-T01-IM6', 'P-0020151-T01-IM6',
		'P-0020757-T01-IM6', 'P-0025648-T01-IM6', 'P-0027375-T01-IM6', 'P-0029228-T01-IM6', 'P-0030260-T01-IM6',
		'P-0032548-T01-IM6', 'P-0032602-T01-IM6', 'P-0033410-T01-IM6', 'P-0035146-T01-IM6', 'P-0035147-T01-IM6',
		'P-0006170-T01-IM5', 'P-0005197-T01-IM5', 'P-0004379-T01-IM5', 'P-0004379-T02-IM6', 'P-0011570-T01-IM5',
		'P-0012171-T02-IM6', 'P-0012333-T01-IM5', 'P-0013557-T01-IM5', 'P-0016825-T01-IM6', 'P-0016972-T01-IM6',
		'P-0017862-T01-IM6', 'P-0018437-T01-IM6', 'P-0020271-T01-IM6', 'P-0024488-T01-IM6', 'P-0025242-T01-IM6',
		'P-0025662-T01-IM6', 'P-0029972-T01-IM6', 'P-0031195-T01-IM6', 'P-0032548-T01-IM6', 'P-0032602-T01-IM6',
		'P-0035146-T01-IM6', 'P-0035147-T01-IM6', 'P-0035770-T01-IM6', 'P-0036297-T01-IM6', 'P-0009964-T02-IM5',
		'P-0013787-T01-IM5', 'P-0015626-T01-IM6', 'P-0024219-T01-IM6', 'P-0032589-T01-IM6', 'P-0037227-T01-IM6',
		'P-0024733-T01-IM6'
		])

#get cases I manually reviewed in facets and decided the fit was weird
def get_facets_blacklist():
	return set([
			'P-0001649-T01-IM3', 'P-0004362-T01-IM5', 'P-0006368-T01-IM5', 'P-0009858-T01-IM5', 'P-0010649-T01-IM5',
			'P-0021600-T01-IM6', 'P-0024184-T01-IM6', 'P-0027380-T01-IM6', 'P-0036284-T01-IM6', 'P-0013462-T01-IM5',
			'P-0017986-T01-IM6', 'P-0003167-T01-IM5'
		])

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















