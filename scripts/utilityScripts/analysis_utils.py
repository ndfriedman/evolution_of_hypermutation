#written by Noah Friedman 
#a script containing functions used for analysis task
import sys
import argparse
import os
import pandas as pd
import numpy as np
import re
import math

import scipy.stats

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

import imp
#maf_analysis_utils = imp.load_source('maf_analysis_utils', pathPrefix + '/ifs/work/taylorlab/friedman/myUtils/maf_analysis_utils.py')
#import maf_analysis_utils


def path_fix(basePath):
	pathPrefix = ''
	if os.getcwd() == '/Users/friedman/Desktop/mnt': pathPrefix = '/Users/friedman/Desktop/mnt'
	elif os.getcwd() == '/Users/friedman/Desktop/WORK/hypermutationProjectJupyterScripts': pathPrefix = '/Users/friedman/Desktop/WORK/hypermutationProjectJupyterScripts'
	return pathPrefix + basePath

#loads in a big file with a progress bar that gives us a realistic estimate of how long its gonna take
def load_in_df_with_progress(filePath, nLinesFile, cSize=10000):
	chunks = []
	cntr = 0.0
	for chunk in pd.read_table(filePath, chunksize=cSize):
		cntr += 1
		print round(cntr * cSize/nLinesFile, 4)*100 ,'percent done'
		chunks.append(chunk)
	print 'performing big concat then returning'
	return pd.concat(chunks)

########USEFUL code area

#SET UP IMPACT SIGNATURES INFORMATION
"""impactSigs = pd.read_table(pathPrefix + '/ifs/res/taylorlab/impact_sigs/mixedpact_data_mutations_unfiltered.sigs.tab.txt')
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/dmp/mskimpact/data_clinical_sample.txt')
impactSigs['cancer_type'] = impactSigs['pid'].apply(lambda x: cDict[x] if x in cDict else None)""" 


#hypermuation project utilities

#function for making comparissons between gene incidences in different cancer types
#Used for figure 2c
def make_comparissons(maf, mode = 'gene', mutationType='msi',
        cancerType1 = 'Endometrial Cancer', cancerType2 = 'Colorectal Cancer'):
    tsgs = get_gene_and_cohort_list_utils.get_tsgs()
    oncogenes = get_gene_and_cohort_list_utils.get_oncogenes()
    indelClassifications = ['Frame_Shift_Del', 'Frame_Shift_Ins']  #TODO actually only include MSI indels not just random indels
    truncatingClassifications = ['']
    
    if mutationType == 'msi':
        maf = maf[maf['Variant_Classification'].isin(indelClassifications)]
        maf = maf[maf['correctedAllele'].notnull()]
    elif mutationType == 'pole':
        maf = maf[maf['Variant_Classification'].isin(['Nonsense_Mutation'])] 
    
    if mode == 'gene':
        maf['allele'] = maf['Hugo_Symbol']
    else:
        maf['allele'] = maf['Hugo_Symbol'] + '_' + maf['HGVSp_Short']
    
    c1Maf = maf[maf['cancerType'] == cancerType1]
    c2Maf = maf[maf['cancerType'] == cancerType2]
    
    nC1 = 1.0*len(set(c1Maf['Tumor_Sample_Barcode']))
    nC2 = 1.0*len(set(c2Maf['Tumor_Sample_Barcode']))
    
    listOfDicts = []
    for allele in set(maf['allele']):
        allele = str(allele)
        aMafC1 = c1Maf[c1Maf['allele'] == allele]
        aMafC2 = c2Maf[c2Maf['allele'] == allele]
        
        gene = ''
        if mode == 'gene':
            gene = allele
        else:
            gene = allele.split('_')[0]
        
        
            
        geneType = 'tsg' if gene in tsgs else 'oncogene' if gene in oncogenes else None
        c1Count = len(set(aMafC1['Tumor_Sample_Barcode']))
        c2Count = len(set(aMafC2['Tumor_Sample_Barcode']))
        listOfDicts.append({'Allele': allele, 'GeneType': geneType, 'Gene': gene,
                            'N_C1': c1Count, 'N_C2': c2Count, 'total_C1': nC1, 'total_C2': nC2,
                            'c1_cancerType': cancerType1, 'c2_cancerType': cancerType2,
                           'perCase_c1': c1Count/nC1, 'perCase_c2': c2Count/nC2})
    
    df = pd.DataFrame(listOfDicts)    
    df['n_NotPresent_c1'] = df['N_C1'].apply(lambda x: nC1 - x)
    df['n_NotPresent_c2'] = df['N_C2'].apply(lambda x: nC2 - x)
    
    #get fisher's test results
    df['p_proportions_z_score'] = df.apply(lambda row: proportions_ztest(np.array([row['N_C1'], row['N_C2']]),
                                                                         np.array([nC1, nC2]))[1], axis=1)
    
    df = df[df['p_proportions_z_score'].notnull()] #remove null z scores, otherwise the qVal wont calculate
    fdrDict = dict(zip(df['Allele'], fdrcorrection(df['p_proportions_z_score'])[1]))
    df['qVal'] = df['Allele'].apply(lambda x: fdrDict[x])
    
    wntGenes = get_gene_and_cohort_list_utils.get_pathway_genes('WNT')
    pi3kGenes = get_gene_and_cohort_list_utils.get_pathway_genes('PI3K') | set(['INPPL1', 'JAK1']) #manually add these guys
    df['pathway'] = df['Gene'].apply(lambda x: 'WNT' if x in wntGenes else 'PI3K' if x in pi3kGenes else 'OTHER')
    return df
    
    
#USED FOR WORKING WITH MSI INDELS
#returns the name of the gene and left most position spanning the msi indel
#uses re stuff
def get_left_aligned_allele_name(hgvsNames):
    
    positions = []
    geneName = hgvsNames[0].split('_p.')[0]
    for entry in hgvsNames:
        if len(entry.split('_p.')) == 2: #ignore weirdly formatted hgvs names
            variantNotation = entry.split('_p.')[1]
            number = variantNotation[1:]
            refAA = variantNotation[0]
            position = re.match('\d*', number).group(0)
            positions.append((position, refAA))
    
    if len(positions) == 0: return None #if all the hgvs names were ill formatted return None
    
    minEntry = sorted(positions)[0] #this is a sorted list of tuples the first thing is the position second is the reference aa
    return geneName + '_p.' + str(minEntry[1]) + str(minEntry[0])

#collapses all indels within 1 bp of each other for a start to be at the same location/name for matching
def standardize_allele_names(msiLengthInfo, observedMuts):
    
    neverObservedSites = set([]) #all the names of sites from criags msi file we cant match with the real maf
    msiSitesToNameMapping = {} #a dictionary mapping each msi site allele from craigs file to its corrected name
    mafMsiSiteToNameMapping = {} #a dictionary mapping each msi site allele from the maf to its corrected name
    
    cntr = 0.0
    for hgvs in set(msiLengthInfo['allele']):
        
        cntr += 1
        #if cntr%500 == 0: print 100*(cntr/len(set(msiLengthInfo['allele']))), 'percent done'
        
        startPos = msiLengthInfo[msiLengthInfo['allele'] == hgvs]['Start_Position']
        
        
        if startPos.shape[0] == 1:
            #we want all names given to indels near (within 1position) of the start position of the MSI site in Craig's file
            putativeVariantNames = list(set(observedMuts[(abs(observedMuts['Start_Position'] - int(startPos)) < 2)
                                                    & (observedMuts['Variant_Type'].isin(set(['INS', 'DEL'])))]['allele']))
            
            if len(putativeVariantNames) == 0:
                neverObservedSites.add(hgvs) #if it cant be matched in the MAF we add it to never observed sites
                #note some of there are likely to be actually matched but missed by my method
            else:
                trueVariantName = get_left_aligned_allele_name(putativeVariantNames)
                
                #NOW WE PROPERLY create the mappings
                for putativeVariantName in putativeVariantNames:
                    mafMsiSiteToNameMapping[putativeVariantName] = trueVariantName
                msiSitesToNameMapping[hgvs] = trueVariantName
        else:
            pass #ignore variants with multiple start position in the msi info file
        
    return neverObservedSites, msiSitesToNameMapping, mafMsiSiteToNameMapping


#
########
###################
###############################
###############################################
###############################
###################
#########
#

def get_ids_by_hypermutant_status(hypermutantIdDir='/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType='', hypermutantStatus = 'Hypermutated'):
	cancerTypeAdj = re.sub(' ', '_', cancerType)
	path = os.path.join(hypermutantIdDir, cancerTypeAdj + '.tsv')
	df = pd.read_table(path)
	if hypermutantStatus == 'all':
		return set(df['Tumor_Sample_Barcode'])
	else:
		return set(df[df['hypermutantClassification'] == hypermutantStatus]['Tumor_Sample_Barcode'])


def enumerate_related_unrelated_genes_for_hypermutation_analysis(allImpactMuts, cTypes=['Endometrial Cancer', 'Colorectal Cancer', 'Glioma'], pathPrefix='~/Desktop/mnt/'):
	genesImplicatedInCancerTypes= maf_analysis_utils.create_dictionary_mapping_genes_to_cancer_types_with_implication(allImpactMuts, pathPrefix=pathPrefix, cancerTypes=cTypes, t=0.04)
	hypermutationInitiatingGenes = set(['MSH6', 'MLH1', 'MSH2', 'PMS2', 'POLE'])
	relatedGenesDict = {}
	for cancerType, genes in genesImplicatedInCancerTypes.items():
	    if cancerType == 'Colorectal Cancer' or 'Endometrial Cancer':
	        genes = genes | hypermutationInitiatingGenes
	        relatedGenesDict[cancerType] = genes
	    else:
	        relatedGenesDict[cancerType] = genes
	return relatedGenesDict

#enumerates the ids of all hypermutated cases in impact
def enumerate_all_hypermutated_cases(pathPrefix = '/Users/friedman/Desktop/mnt'):
	hypermutatedCases = set([])
	d = pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds'
	for f in os.listdir(d):
		fPath = os.path.join(d, f)
		df = pd.read_table(fPath)
		hypermutants = set(df[df['hypermutantClassification'] == 'Hypermutated']['Tumor_Sample_Barcode'])
		hypermutatedCases = hypermutants | hypermutatedCases
	return hypermutatedCases


############

def map_cases_to_msi_sensor_class(df, msiSensorInfo='/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_clinical_sample.txt'):
	dfMSI = pd.read_table(msiSensorInfo, skiprows=[1,2,3])
	msiClassDict = dict(zip(dfMSI['#Sample Identifier'], dfMSI['MSI Type']))
	df['caseMsiClass'] = df['Tumor_Sample_Barcode'].apply(lambda x: msiClassDict[x] if x in msiClassDict else None)
	return df

#utility function to quickly tell us the chromosomes of all the IMPACT genes
def map_impact_genes_to_chromosome(impactRegionsFile = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/genelist.with_aa.interval_list'):
	df = pd.read_table(impactRegionsFile, header=None)
	df = df.rename(columns = {0: 'Chromosome', 1:'pos1', 2:'pos2', 3:'strand', 4:'geneTranscript'}) #rename the columns of the dataframe because it is ill formatted
	df['Gene'] = df['geneTranscript'].apply(lambda x: x.split(':')[0])
	dfReduced = df.drop_duplicates(subset=['Gene']) #drop all the duplicate rows
	return dict(zip(dfReduced['Gene'], dfReduced['Chromosome']))

def map_impact_genes_to_chromosome_arm(armFilePath = '/juno/work/taylorlab/friedman/myAdjustedDataFiles/IMPACTv6_gene_cytoband.txt'):
	armDf = pd.read_table(armFilePath)
	armDf['arm'] = armDf['Cytoband'].apply(lambda x: x.split('p')[0] + 'p' if 'p' in x
											else x.split('q')[0] + 'q' if 'q' in x
											else None)
	return dict(zip(armDf['Gene_Symbol'], armDf['arm']))

#KAPLAN MEIER ANALYSIS

def test_significance(df1, df2):
	results = logrank_test(df1['os_years'], df2['os_years'], df1['CENSOR'], df2['CENSOR'], alpha=.99)
	return results.p_value
	#results.print_summary()

def make_kaplan_meier_plots(dfs):
	def correct_censorship(dfs): #creates a sensorship column which invertys the pt vital status colum (0 is alive, 1 is dead)
		for i in range(len(dfs)):
			df = dfs[i]
			df['CENSOR'] = df['PT_VITAL_STATUS'].apply(lambda x: 1 if x == 0 else 0)
			dfs[i] = df
		return dfs

	dfs = correct_censorship(dfs)
	kmf = KaplanMeierFitter()
	kmf.fit(dfs[0]['os_years'], dfs[0]['CENSOR'], label='cond')
	ax = kmf.plot()
	kmf.fit(dfs[1]['os_years'], dfs[1]['CENSOR'], label='nonCond')
	ax = kmf.plot(ax=ax)
	fig = ax.get_figure()
	#fig.savefig('testFig.pdf')

	print test_significance(dfs[0], dfs[1])


#INFORMATION FUNCTIONS ##################################################

#function taken from https://stackoverflow.com/questions/46736258/deleting-diagonal-elements-of-a-numpy-array
#removes all the diagonals from a matrix
def skip_diag_strided(A):
    m = A.shape[0]
    strided = np.lib.stride_tricks.as_strided
    s0,s1 = A.strides
    return strided(A.ravel()[1:], shape=(m-1,m), strides=(s0+s1,s1)).reshape(m,-1)

#a function that takes a numpy array and returns all pairwise differences (excludes the diagonal)
def calculate_all_pairwise_differences(arr):
	if len(arr) == 0: return [] #return an empty list for something where there are issues
	diffMatrix = arr[:,np.newaxis] - arr #get the differences in magnitudes between elements of a matrix
	diffMatrixNoDiag = skip_diag_strided(diffMatrix) #get rid of the diagnoal of the matrix (the differences in ccf between the exact same variant which will be 0)
	diffs = [i for i in diffMatrixNoDiag.flatten() if i >= 0] #ignore negative values so we dont double count entries
	return diffs  

def mean_confidence_interval(data, confidence=0.95):
	data = [i for i in data if not math.isnan(i)] #remove nas
	a = 1.0 * np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
	return m, m-h, m+h

#TODO ---put this in analysis utils
def normalize_counter(cntrObj, mode='Round', nDigitsRound=2):
	cntrObj = cntrObj.copy() #operate on a copy of the counter
	total = sum(cntrObj.values(), 0.0)
	for key in cntrObj:
		cntrObj[key] /= total
		if mode == 'Round': cntrObj[key] = round(cntrObj[key], nDigitsRound)
	return cntrObj

#returns the mean value of a column specified by colname
def get_mean_of_df_col(df, colname, idColumn = 'Tumor_Sample_Barcode'):
	df = df.drop_duplicates(subset=[idColumn])
	return np.nanmean(np.asarray(list(df[colname])))

#returns the median of a df column
def get_median_of_df_col(df, colname, idColumn = 'Tumor_Sample_Barcode'):
	df = df.drop_duplicates(subset=[idColumn])
	return np.nanmedian(np.asarray(list(df[colname])))


def get_n_cases_with_mutation_present(maf, gene, idCol='Tumor_Sample_Barcode'):
	n = len(set(maf[maf['Hugo_Symbol'] == 'TP53'][idCol]))
	print n, 1.0*n/len(set(maf['Tumor_Sample_Barcode']))


def get_age_information(ageInformationPath = '/ifs/work/taylorlab/friedman/msk-impact/msk-impact/darwin/darwin_age.txt'):
	ageDf = pd.read_table(ageInformationPath)
	d = dict(zip(ageDf['PATIENT_ID'], ageDf['AGE']))
	return d

#OLD PATH: '/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_clinical_sample.txt'
def get_cancer_type_information(cancerTypeDfPath = '/juno/work/taylorlab/friedman/myAdjustedDataFiles/cancerTypeInfo_asOfNov192019.txt', mode='pid'):
	cancerTypeDf = pd.read_table(cancerTypeDfPath)
	if mode == 'pid':
		cancerTypeDf['#Sample Identifier'] = cancerTypeDf['#Sample Identifier'].apply(lambda x: x[:9])
	d = dict(zip(cancerTypeDf['#Sample Identifier'], cancerTypeDf['Cancer Type']))
	return d

def get_cancer_type_detailed_information(cancerTypeDfPath = '/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_clinical_sample.txt', mode='pid'):
	cancerTypeDf = pd.read_table(cancerTypeDfPath)
	if mode == 'pid':
		cancerTypeDf['#Patient Identifier'] = cancerTypeDf['#Patient Identifier'].apply(lambda x: x[:9])
	d = dict(zip(cancerTypeDf['#Patient Identifier'], cancerTypeDf['Cancer Type Detailed']))
	return d

#TODO MAKE THIS MORE SUSTAINABLE
def get_gene_length_info(bedFilePath = '/ifs/res/pwg/data/gencode/gencode.v19.all_gene_bounds.bed'):
	impactGenes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	bedDf = pd.read_table(bedFilePath)	
	bedDf = bedDf[bedDf['OR4F5'].isin(impactGenes)]
	bedDf['geneLength'] = bedDf.apply(lambda row: row['70008'] - row['69090'], axis=1)
	return dict(zip(bedDf['OR4F5'], bedDf['geneLength']))

#returns a dictionary mapping each gene in impact to the size of the cds targeted by the panel
def get_cds_size_targeted_by_impact(infoFilePath = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/impact_gene_reference_signatures.tsv'):
	df = pd.read_table(infoFilePath)
	return dict(zip(df['Hugo_Symbol'], df['cds_length']))

#gives a set of tumor surpressor genes
def get_tumor_supressor_genes():
	return set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])

#SIGNATURE SPECIFIC ANALYSIS UTILS #############################################################

#util to give the top N most epxressed signatures:
def get_n_top_signatures(row, n=2):
	#signatureCols = list(row.columns.values)
	row = row[['mean_' + str(i) for i in range(1,31)]]
	l = list(row)
	return list(reversed([str(i + 1) + ':' + str(l[i]) for i in np.argsort(l)[-n:]])) #I plus one to take into account the signatures ordering

#enumerates a set of all 96 possible trinucelotides/changes
def get_all_possible_quadNucs():
	allSpectra = []
	for firstLetter in ['A','C','G','T']:
		for change in ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']:
			for lastLetter in ['A','C','G','T']:
				allSpectra.append(firstLetter+change+lastLetter)
	return(set(allSpectra))

#STATISTICAL TESTS ####################################################

def test_significance(df1, df2):
	results = logrank_test(df1['os_years'], df2['os_years'], df1['CENSOR'], df2['CENSOR'], alpha=.99)
	return results.p_value

#does the mann whitney u test
def do_mann_whitney_test(dist1, dist2):
	return scipy.stats.mannwhitneyu(dist1, dist2).pvalue

#FISHER TEST PLEASE!!
def do_p_val_tests_with_comps(dist1, dist2, comparissonMessage=None):
	print comparissonMessage
	print 'Ns for each distribution: ', len(list(dist1)), len(list(dist2))
	print 'MEANS of distributions: ', np.nanmean(np.asarray(list(dist1))), np.nanmean(np.asarray(list(dist2)))
	print 'MEDIANS of distributions: ', np.nanmedian(np.asarray(list(dist1))), np.nanmedian(np.asarray(list(dist2)))
	print 'P val: ', do_mann_whitney_test(dist1, dist2)

#calcualtes a correlation and pearson p val for two columns of the same df
def do_pearson_correlation_one_df(df, col1, col2):
	df = df[np.isfinite(df[col1])]
	df = df[np.isfinite(df[col2])]
	c1 = np.asarray(df[col1])
	c2 = np.asarray(df[col2])
	#print c1, c2
	return scipy.stats.pearsonr(c1,c2)

#util function to pretty print tumor sample barcodes for cbioportal
def print_for_cbio_portal(s):
	for v in s:
		print v


###############################################
#DO STUFF WITH 

#RUN MAF TO MAF 



def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()


if __name__ == '__main__':
    main()















