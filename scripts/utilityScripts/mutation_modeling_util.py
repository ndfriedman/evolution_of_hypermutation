#written by Noah Friedman
#a collection of utilities designed to model mutation acquisition in silico 
import sys
import argparse
import os
import pandas as pd
import numpy as np
import math

from collections import Counter
#sys.path.append('/ifs/work/taylorlab/friedman/')
#import variant_consequence_analysis_util
#import mutationSigUtils


def add_zero_cols_to_counts_df(df):
    allBases = ['A', 'C', 'G', 'T']
    changes = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG'] #format: 'CA' means a change from C>A
    allQuadNucs = [firstBase + change + lastBase for firstBase in allBases for change in changes for lastBase in allBases] #enumerate all 96 quadnucs for signatures
    allAnnotationTypes = ['truncating', 'oncogenic', 'hotspot', 'nonSilent']
    allColHeaders = [quadNuc + '_' + annotation for quadNuc in allQuadNucs for annotation in allAnnotationTypes]
    missingColumns = set(allColHeaders) - set(df.columns.values)
    for c in missingColumns:
        df[c] = 0
    return df

#CURRENTLY DEPRECTED BUT IS A GOOD IDE A FOR WHAT SHIT WOULD LOOK LIKE IF WE RANDOMLY PICK SIGS
def get_trinucleotide_probabilities(trinucEntries):
		rowAsDict = trinucEntries.to_dict()
		trinucChoices = rowAsDict.keys()
		trinucCounts = rowAsDict.values()
		trinucProbabilities = [1.0*i/sum(trinucCounts) for i in trinucCounts]
		return trinucChoices, trinucProbabilities

def normalize_counter(cntrObj, mode='Round', nDigitsRound=2):
	total = sum(cntrObj.values(), 0.0)
	for key in cntrObj:
		cntrObj[key] /= total
		if mode == 'Round': cntrObj[key] = round(cntrObj[key], nDigitsRound)
	return cntrObj

####################################

###################################SAMPLING AREA

#note this could be done a lot faster as a matrix multiplication; todo implement it that way?
def get_quadnuc_fracs_given_decomposition(caseDecomposition, spectraD):
	quadNucs = spectraD['Signature.1'].keys()
	quadNucD = {} #initialize the fractions for each quadnuc
	for quadNuc in quadNucs:
		quadNucD[quadNuc] = 0

	for key, value in spectraD.items():
		decompositionFraction = caseDecomposition[key]
		for quadNuc, percent in value.items():
			quadNucD[quadNuc] = quadNucD[quadNuc] + decompositionFraction*percent

	#RETURNS a dictionary for each quadnuc, summing up to one
	return quadNucD
			
#NOTE THIS IS FOR IMPACT-468
def summarize_mutation_chances(allPossibleMutsSummary, mutType):
	allBases = ['A', 'C', 'G', 'T']
	changes = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG'] #format: 'CA' means a change from C>A
	quadnucs = [firstBase + change + lastBase for firstBase in allBases for change in changes for lastBase in allBases] #enumerate all 96 quadnucs for signatures
	mutTypes = ['_hotspot', '_nonSilent', '_oncogenic', '_truncating']
	d = {}
	for quadNuc in quadnucs:
		nAllMuts = 1.0*np.nansum(allPossibleMutsSummary[quadNuc + '_nonSilent'])
		nSelectedMuts = np.nansum(allPossibleMutsSummary[quadNuc + '_' + mutType])
		d[quadNuc] = nSelectedMuts/nAllMuts
	return d
		

###################################OVERALL MODELING/COMPARISSON

def get_average_oncogenic_mut_burdens_for_tcga(caseDecompositionDict, 
	spectraD, 
	quadNucDict,
	quadNucDictAdjusted,
	geneOncogenicOrHotspotProbDict,
	geneMutProbInfo,
	n=1000,
	mutMode='oncogenic'):

	oncogenicMutFracs = []
	for model in ['model1', 'model2', 'model3']:

		nNonImpactMuts, spectraChosen, quadNucsChosen, genesChosen, mutationsChosen = pick_n_tcga_mutations_given_spectra(
		caseDecompositionDict, 
		spectraD, 
		quadNucDict,
		quadNucDictAdjusted,
		geneOncogenicOrHotspotProbDict,
		geneMutProbInfo,
		n=n,
		mutMode=mutMode,
		modelName=model)

		nNonOncogenicMuts = len([i for i in mutationsChosen if i == 'nonOncogenic'])
		nNotOncogenic = nNonImpactMuts + nNonOncogenicMuts
		nOncogenic = n - nNotOncogenic
		fracOncogenic = nOncogenic*1.0/n
		oncogenicMutFracs.append(fracOncogenic)
	averageFrac = np.nanmean(oncogenicMutFracs)
	return averageFrac

def get_simulated_oncogenic_mut_frac_per_case(
	signaturesDf,
	spectraD, 
	quadNucDict,
	quadNucDictAdjusted,
	geneOncogenicOrHotspotProbDict,
	geneMutProbInfo,
	n=1000,
	mutMode='oncogenic',
	idColName='Tumor_Sample_Barcode'):
	
	mutFracDict = {}
	cases = set(signaturesDf[idColName])
	signatureColNames = ['Signature.' + str(i) for i in range(1,31)]
	cntr = 0
	for case in cases:
		if cntr%50 == 0: print cntr
		cntr += 1

		localDf = signaturesDf[signaturesDf[idColName] == case]
		caseSignatures = localDf[signatureColNames].iloc[0].to_dict()
		avgVal = get_average_oncogenic_mut_burdens_for_tcga(
			caseSignatures, 
			spectraD, 
			quadNucDict,
			quadNucDictAdjusted,
			geneOncogenicOrHotspotProbDict,
			geneMutProbInfo,
			n=n,
			mutMode=mutMode)

		mutFracDict[case] = avgVal

	return mutFracDict


def compare_observed_to_expected_mut_burden_data(actualMutInfoDf, simulatedMutInfo, idCol='Tumor_Sample_Barcode'):
	listOfDicts = []
	for caseId, value in simulatedMutInfo.items():
		localD = {}
		acutalMutCaseInfo = actualMutInfoDf[actualMutInfoDf[idCol] == caseId]
		if acutalMutCaseInfo.shape[0] > 0:
			actualNOncogenicSnps = acutalMutCaseInfo.iloc[0]['nOncogenicSnps']
			actualNmut = acutalMutCaseInfo.iloc[0]['Nmut']
			simulatedNOncogenicSNPs = actualNmut*value
			localD['Tumor_Sample_Barcode'] = caseId
			localD['Nmut'] = actualNmut
			localD['nOncogenicSnps'] = actualNOncogenicSnps
			localD['simulatedNOncogenicSNPs'] = simulatedNOncogenicSNPs
			if actualNOncogenicSnps > 0:
				localD['ratioObservedToExpected'] = math.log(1.0*actualNOncogenicSnps/simulatedNOncogenicSNPs, 2)
			else:
				localD['ratioObservedToExpected'] = None
			listOfDicts.append(localD)
	return pd.DataFrame(listOfDicts)

###################################################FUNCTIONS

#TODO change this TO include two different ways of treating silent mutation data


def calculate_quadnuc_based_oncogenic_susceptibility_dict(df):

	def calculate_oncogenic_susceptibility_of_quadnuc_given_panel(panelDf, quadNuc):

		quadNucTotalMuts = sum(panelDf[quadNuc])
		quadNucSilentMuts = sum(panelDf[quadNuc + '_silent'])
		quadNucOncogenicMuts = sum(panelDf[quadNuc + '_oncogenic'])
		nNonSilentMutsAtQuadNuc = quadNucTotalMuts - quadNucSilentMuts
		return 1.0*quadNucOncogenicMuts/nNonSilentMutsAtQuadNuc

	d = {}
	quadNucs = set([start + change + end for start in ['A', 'C', 'T', 'G'] 
		for change in ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
			for end in ['A', 'C', 'T', 'G']])

	print 'gendarameee'
	print len(quadNucs)

	impact368genes = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])
	impact410genes = set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	impact468genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 

	dfImpact368 = df[df['Hugo_Symbol'].isin(impact368genes)]
	dfImpact410 = df[df['Hugo_Symbol'].isin(impact410genes)]
	dfImpact468 = df[df['Hugo_Symbol'].isin(impact468genes)]

	for quadnuc in quadNucs:
		localD = {
		'IMPACT_368':calculate_oncogenic_susceptibility_of_quadnuc_given_panel(dfImpact368, quadnuc),
		'IMPACT_410':calculate_oncogenic_susceptibility_of_quadnuc_given_panel(dfImpact410, quadnuc),
		'IMPACT_468':calculate_oncogenic_susceptibility_of_quadnuc_given_panel(dfImpact468, quadnuc)
		}
		d[quadnuc] = localD
	return d

#calculates the probability that an oncogenic mutation occurs at a gene given a quad nuc
#NOTE MAYBE write differently so its not just a repeat of the above fucntion
def calculate_quadnuc_and_gene_based_oncogenic_susceptibility_dict(df):

	def calculate_oncogenic_susceptibility_of_quadnuc_at_gene_given_panel(panelDf, quadNuc, gene):
		quadNucTotalMuts = sum(panelDf[quadNuc])
		quadNucSilentMuts = sum(panelDf[quadNuc + '_silent'])
		quadNucGeneOncogenicMuts = sum(panelDf[panelDf['Hugo_Symbol'] == gene][quadNuc + '_oncogenic'])
		nNonSilentMutsAtQuadNuc = quadNucTotalMuts - quadNucSilentMuts
		return 1.0*quadNucGeneOncogenicMuts/nNonSilentMutsAtQuadNuc

	d = {}

	quadNucs = set([start + change + end for start in ['A', 'C', 'T', 'G'] 
		for change in ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
			for end in ['A', 'C', 'T', 'G']])

	impact368genes = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])
	impact410genes = set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	impact468genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']) 
	
	dfImpact368 = df[df['Hugo_Symbol'].isin(impact368genes)]
	dfImpact410 = df[df['Hugo_Symbol'].isin(impact410genes)]
	dfImpact468 = df[df['Hugo_Symbol'].isin(impact468genes)]

	cntr = 0
	for gene in impact468genes:
		cntr += 1
		if cntr%25==0: print 'considering the ', cntr ,'th gene'
		quadNucsD = {}
		for quadnuc in quadNucs:
			localD = {
			'IMPACT_368':calculate_oncogenic_susceptibility_of_quadnuc_at_gene_given_panel(dfImpact368, quadnuc, gene),
			'IMPACT_410':calculate_oncogenic_susceptibility_of_quadnuc_at_gene_given_panel(dfImpact410, quadnuc, gene),
			'IMPACT_468':calculate_oncogenic_susceptibility_of_quadnuc_at_gene_given_panel(dfImpact468, quadnuc, gene)
			}
			quadNucsD[quadnuc] = localD
		d[gene] = quadNucsD
	return d


#tries to calculates the fraction of MSI mutations the impact panel that ought to be oncogenic theoretically
#it calculates the fraction of repeat regions in the impact panel that occur in tumor suppresors and multiplies by 2/3 (on the theory that 1/3 of all indels will be inframe indels)
def calculate_pan_impact_likelihood_of_oncogenic_mmr_indel(repeatRegionInfo = '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/IMPACT_repeat_stats.txt'):
	tumorSuppressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
	repeatInfoDf = pd.read_table(repeatRegionInfo) 
	repeatInfoDf['nRepeatBases'] = repeatInfoDf.apply(lambda row: row['bp']*row['rep.pct'], axis=1)
	totalNRepeatBasesInTumorSuppressors = np.nansum(repeatInfoDf[repeatInfoDf['Hugo_Symbol'].isin(tumorSuppressors)]['nRepeatBases'])
	totalNBases = np.nansum(repeatInfoDf['nRepeatBases'])
	repeatTumorSuppressorFrac = totalNRepeatBasesInTumorSuppressors/totalNBases
	return repeatTumorSuppressorFrac*(.95)


####################################

#converts  the counts information to fraction chances
#creates an intermediate chances data frame witht the chances of each specific mutation happening at a gene by motif
def convert_counts_information_to_fraction(countsDf):

	allBases = ['A', 'C', 'G', 'T']
	changes = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG'] #format: 'CA' means a change from C>A
	allQuadNucs = [firstBase + change + lastBase for firstBase in allBases for change in changes for lastBase in allBases] #enumerate all 96 quadnucs for signatures

	countsDf = countsDf.fillna(0)

	allMutQuadNucCounts = {}
	for quadnuc in allQuadNucs:
		allMutQuadNucCounts[quadnuc] = np.nansum(countsDf[quadnuc + '_nonSilent'])

	listOfDicts = []
	for quadnuc in allQuadNucs:
		nQuadnucMut = allMutQuadNucCounts[quadnuc]
		geneHotspotDict = dict(zip(countsDf['Hugo_Symbol'], countsDf[quadnuc + '_hotspot']))
		geneOncogenicDict = dict(zip(countsDf['Hugo_Symbol'], countsDf[quadnuc + '_oncogenic']))
		geneTruncatingDict = dict(zip(countsDf['Hugo_Symbol'], countsDf[quadnuc + '_truncating']))
		for gene in set(countsDf['Hugo_Symbol']):

			listOfDicts.append({
				'gene': gene,  'quadNuc': quadnuc,
				'nHotspot': geneHotspotDict[gene], 'hotspotChance': (1.0*geneHotspotDict[gene])/nQuadnucMut,
				'nOncogenic': geneOncogenicDict[gene], 'oncogenicChance': (1.0*geneOncogenicDict[gene])/nQuadnucMut,
				'nTruncating': geneTruncatingDict[gene], 'truncatingChance': (1.0*geneTruncatingDict[gene])/nQuadnucMut
				})

	return pd.DataFrame(listOfDicts)

#OTHER WAY TO SET UP DATA
"""
quadNucsNmutDict = {}
for qn in set(dfMuts['quadNuc']):
    nQnMuts = sum(dfMuts[dfMuts['quadNuc'] == qn]['nNonSilent'])
    quadNucsNmutDict[qn] = nQnMuts

dfMuts['hotspotChance'] = dfMuts.apply(lambda row: 1.0*row['nHotspot']/quadNucsNmutDict[row['quadNuc']], axis=1)
dfMuts['oncogenicChance'] = dfMuts.apply(lambda row: 1.0*row['nOncogenic']/quadNucsNmutDict[row['quadNuc']], axis=1)
dfMuts['truncatingChance'] = dfMuts.apply(lambda row: 1.0*row['nTruncating']/quadNucsNmutDict[row['quadNuc']], axis=1)"""


#functions to get quad nuc based oncogenic fractions

#THIS function counts the number of oncogenic, hotspot etc mutations at each quad nuc which is a critixal tool for downstream quantifying the expected oncogenic mutation burdne
#USAGE: df = mutation_modeling_util.quantify_quadnuc_mut_type_susceptibility_per_mutation(allMutsPossible, mode='perGene')
#where allMutsPossible = analysis_utils.load_in_df_with_progress(filePath = pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/allMutsContextSummary.tsv', nLinesFile=4243342) 

def quantify_quadnuc_mut_type_susceptibility_per_mutation(simMafData, #a maf of all possible mutations in impact
	#Default is to only do this calculation with the IMPACT 341 panel
	genes= set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 
		'YES1']),
	mode='perQuadNuc'
):
	allBases = ['A', 'C', 'G', 'T']
	changes = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG'] #format: 'CA' means a change from C>A
	allQuadNucs = [firstBase + change + lastBase for firstBase in allBases for change in changes for lastBase in allBases] #enumerate all 96 quadnucs for signatures

	d = {}

	#mafSelect = simMafData[simMafData['Hugo_Symbol'].isin(genes)]
	mafSelect = simMafData

	genes = set(mafSelect['Hugo_Symbol'])

	#SILENT is synonymous mutations and intron variants
	#NOTE I HABVENT included intron splice region variants, UTRS, etc but realistically they represent <1% of all mutations so whatever decision I make isnt biasing things much
	silentConsequences = ['synonymous_variant', 'splice_region_variant,synonymous_variant','non_coding_transcript_exon_variant', 'intron_variant,non_coding_transcript_variant']

	mafNonSilent = mafSelect[~mafSelect['Consequence'].isin(silentConsequences)]
	mafOncogeneic = mafSelect[mafSelect['oncogenic'].notnull()]
	mafHotspot = mafSelect[mafSelect['is-a-hotspot'] == 'Y']
	mafTruncating = mafSelect[mafSelect['Consequence'] == 'stop_gained']

	listOfDicts = []
	for quadNuc in allQuadNucs:
		print quadNuc
		qMafNonSilent = mafNonSilent[mafNonSilent['quadNuc'] == quadNuc]
		qMafOncogenic = mafOncogeneic[mafOncogeneic['quadNuc'] == quadNuc]
		qMafHotspot = mafHotspot[mafHotspot['quadNuc'] == quadNuc]
		qMafTruncating = mafTruncating[mafTruncating['quadNuc'] == quadNuc]
		if mode =='perQuadNuc':
			listOfDicts.append({
				'quadNuc': quadNuc, 'nNonSilent': qMafNonSilent.shape[0],
				'nOncogenic': qMafOncogenic.shape[0],
				'nHotspot ': qMafHotspot.shape[0],
				'nTruncating': qMafTruncating.shape[0]
			})
		elif mode =='perGene':
			qSilentGeneCounts = dict(qMafNonSilent['Hugo_Symbol'].value_counts())
			qOncogenicGeneCounts = dict(qMafOncogenic['Hugo_Symbol'].value_counts())
			qHotspotGeneCounts = dict(qMafHotspot['Hugo_Symbol'].value_counts())
			qTruncatingCounts = dict(qMafTruncating ['Hugo_Symbol'].value_counts())
			for gene in genes:
				silentCount = 0
				if gene in qSilentGeneCounts: silentCount = qSilentGeneCounts[gene]
				oncogenicCount = 0
				if gene in qOncogenicGeneCounts: oncogenicCount = qOncogenicGeneCounts[gene]
				hotspotCount = 0
				if gene in qHotspotGeneCounts: hotspotCount = qHotspotGeneCounts[gene]
				truncatingCount = 0
				if gene in qTruncatingCounts: truncatingCount = qTruncatingCounts[gene]
				listOfDicts.append({
					'quadNuc': quadNuc, 'nNonSilent': silentCount,
					'nOncogenic': oncogenicCount,
					'nHotspot ': hotspotCount,
					'nTruncating': truncatingCount,
					'gene': gene
				})
	return pd.DataFrame(listOfDicts)

######################################

#LATEST VERSION JUNE 2 20202
def get_expected_mut_chance_given_quadnuc_fractions(quadNucFractions, susceptibilityDict):
	p = 0
	for key, value in quadNucFractions.items():
		p += value * susceptibilityDict[key]
	return p


def get_expected_oncogenic_val_given_quadnuc_fractions(quadNucFractions, susceptibilityDict, impactVersion):
	p = 0
	for key, value in quadNucFractions.items():
		p += value * susceptibilityDict[key][impactVersion]
	return p

#a version for a different structure
def get_expected_oncogenic_val_given_quadnuc_fractions_v2(quadNucFractions, susceptibilityDict, suffix=''):
	p = 0
	for key, value in quadNucFractions.items():
		p += value * susceptibilityDict[key + '_' + suffix]
	return p


#SAME FUNCTION AS ABOVE 
def get_expected_oncogenic_val_given_quadnuc_fractions_and_gene(quadNucFractions, susceptibilityDict, impactVersion, gene):
	p = 0
	for key, value in quadNucFractions.items():
		p += value * susceptibilityDict[gene][key][impactVersion]
	return p

def compare_observed_and_expected_gene_mutation_rates(susceptibilityDf, observedDf):
	for gene in set(susceptibilityDf['Hugo_Symbol']):
		print gene

###&&&&&&&&*************New code 2-18-2019!!!
def get_chance_of_oncogenic_mutation_at_quadnuc(df, quadnuc):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	allMutsAtQuadNuc = df[df['quadNuc'] == quadnuc]
	nMutsAtQuadNuc = sum(allMutsAtQuadNuc['totalNmut'])
	nOncogenicMutAtQuadNuc = sum(allMutsAtQuadNuc['nOncogenicMut'])
	return 1.0*nOncogenicMutAtQuadNuc/nMutsAtQuadNuc

def get_chance_of_oncogenic_mutation_in_gene_at_quadnuc(df, quadnuc, gene):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	allMutsAtQuadNuc = df[df['quadNuc'] == quadnuc]
	nMutsAtQuadNuc = sum(allMutsAtQuadNuc['totalNmut'])
	nOncogenicMutAtQuadNuc = sum(allMutsAtQuadNuc[allMutsAtQuadNuc['Hugo_Symbol'] == gene]['nOncogenicMut'])
	return 1.0*nOncogenicMutAtQuadNuc/nMutsAtQuadNuc

#given a decomposition of a Sample expressed as a dict (ie 50%sig1, 25% sig2, 25%sig3) randomly picks a signature to have caused a mutation
def pick_spectra_given_decomposition(decompositionDict):
	return np.random.choice(
			decompositionDict.keys(), p=decompositionDict.values()
		)

#One liner to pick a quadnuc given a signature that theoretically genrted it
def pick_quadnuc_given_spectra(spectraDict):
	return np.random.choice(
			spectraDict.keys(), p=spectraDict.values()
		)

#one liner to return a gene that is mutated given a quadNuc
def pick_gene_given_quadnuc(d, quadNuc): 
	return np.random.choice(
			d[quadNuc][0], p=d[quadNuc][1]
		)

#one liner to pick a mutation at random given a 
def pick_mutation_given_gene_and_quadNuc(gene, quadNuc, d, mode='oncogenic'):
	return np.random.choice(
			d[gene + '_' + quadNuc][mode][0],
			p= d[gene + '_' + quadNuc][mode][1]
		)


def asses_simulated_mutation_rates_by_signatures(sigMagnitudesDict, k, spectraD, n, quadNucDict, geneOncogenicOrHotspotProbDict, mMode='oncogenic'):
    d = dict()
    signatures = ['Signature.' + str(i) for i in range(1,31)]
    for sig in signatures:
        sigMagnitudesDict = {key:value for (key, value) in [('Signature.' + str(i), 0) for i in range(1,31)]}
        sigMagnitudesDict[sig] = 1
        d[sig] = mutation_modeling_util.do_k_simulations_of_mutations(sigMagnitudesDict, k, spectraD, n, quadNucDict, geneOncogenicOrHotspotProbDict, mutMode=mMode)
    return d


#given a maf of all possible mutations in impact summarizes for each gene the number of oncogenic/not oncogenic mutations that occur at each quadnuc and the 
def initiate_gene_mut_mapping(allMutsMaf): #TO improve runtime we precalculate a dict mapping 
															#gene: maf of all mutations and maf of oncogenic mutations
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	listOfDicts = []
	allQuadNucs = set(allMutsMaf['quadNuc'])
	for gene in set(allMutsMaf['Hugo_Symbol']):
		print 'analyzing gene ', gene
		geneMuts = allMutsMaf[allMutsMaf['Hugo_Symbol'] == gene]
		chromosome = geneMuts['Chromosome'].iloc[0]
		impactPanel = None

		#Set the proper impact information
		impact368Panel = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])
		impact410Panel = set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
		impact468Panel = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
		if gene in impact368Panel:
			impactPanel = 'IMPACT_368'
		elif gene in impact410Panel:
			impactPanel = 'IMPACT_410'
		elif gene in impact468Panel:
			impactPanel = 'IMPACT_468'
		else:
			impactPanel = None

		for quadNuc in allQuadNucs:
			geneMutsAtQuadNuc = geneMuts[geneMuts['quadNuc'] == quadNuc]
			nMutsInGeneAtQuadNuc = geneMutsAtQuadNuc.shape[0]
			nOncogenicMutsInGeneAtQuadNuc = geneMutsAtQuadNuc[geneMutsAtQuadNuc['oncogenic'].isin(oncogenicMutColNames)].shape[0]
			localD = {'totalNmut': nMutsInGeneAtQuadNuc, 'nOncogenicMut': nOncogenicMutsInGeneAtQuadNuc,
			'Hugo_Symbol': gene,'quadNuc': quadNuc, 'impactPanel': impactPanel, 'Chromosome': chromosome}
			listOfDicts.append(localD)
			
	return pd.DataFrame(listOfDicts)


#an option intialization function that loads in the simulated mutation mafs and summarizes them one at a time
def load_in_and_summarize_simulated_mut_information(mafDirPath = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/simulatedMafs'):
	chromosomes = ['X'] + [str(i) for i in range(1,23)]
	chromosomes.reverse() #reverse the order to do small stuff first for debugging purposes
	listOfDataFrames = []

	for chromosome in chromosomes:
		fileName = 'chr' + chromosome + '_allmut_mafAnno_trinuc.maf'
		path = os.path.join(mafDirPath, fileName)
		print 'loading in file ', path
		maf = pd.read_table(path)
		print 'adding quadNuc to chr', chromosome
		maf['quadNuc'] = maf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
		print 'sumarizing chr', chromosome
		mutMappingDf = initiate_gene_mut_mapping(maf)
		listOfDataFrames.append(mutMappingDf)
		
	concatDf = pd.concat(listOfDataFrames)
	return concatDf 


#TODO: write a function for reading probabilitiy data from disk and converting it into a dict of the appropriate format

#TODO MAKE A SPECIf

def get_impact_version_given_case_name(caseName):
	if 'IM3' in caseName: return 'IMPACT_368'
	elif 'IM5' in caseName: return 'IMPACT_410'
	elif 'IM6' in caseName: return 'IMPACT_468'
	else: return None

		

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--fastaPath', help='Path for full FASTA file', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/b37.fasta')
	parser.add_argument('--bedFilePath', help='Path of full Bed file', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/genelist.with_aa.interval_list')   
	parser.add_argument('--scratchBedPath', help='pathForTemporaryBedFile', default='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/tempScriptFiles/tempBed')
	
	args = parser.parse_args()
	
	enumerate_all_possible_snps(sequenceDfPath = args.bedFilePath,
		genes = ['TP53'],
		scratchBedFilePath = args.scratchBedPath)

if __name__ == '__main__':
    main()



#DEPRECATED MUTATION MODELING CODE


#CONSTANTS FOR MUTATION MODELING
"""
ONCOGENIC_MUTATION_CHANCE_IN_IMPACT = 0.0403674462019 #fraction of all possible non silent mutations in impact that are oncogenic
IMPACT_INDEL_FRAC = 0.156 #the fraction of mutations that are indels in our cohort
NON_MSI_INDEL_ONCOGENICITY = 0.3497
MSI_INDEL_ONCOGENICITY = 0.4856

#FUNCTION that uses a maf and an oncogenicity dictionary and computes observed and expected get_observed_expected_values_for_cohort
#TAKES about 5 minutes to run on the entire impact cohort
def get_observed_expected_values_for_cohort(cases, #the tumor sample barcodes we will consider
	signatures, #the decomposed signatures for impact
	oncogenicSusceptibilityDict, # a dictionary mapping quadnucs to their susceptibility to oncogenic mutations
	nmutDict, #a dictionary mapping tumor sample barcode to nmut
	nmutOncogenicDict, #a dictionary mapping tumor sample barcodes to the number of oncogenic mutations in the case
	perCaseUniqueOncogenicInfo, #a dictionary mapping cases to their n of unique genes oncogenically mutated
	caseIndelFracDict, #a dictionary mapping each cases to the fraction of mutations in the case that are indels
	minNmut = 10 #the number of mutations below which we do not believe we can adequately provide a estimate of the expected mutation burden
	):

	listOfDicts = []

	#allow for interactive run from spyder
	pathPrefix = ''
	if os.getcwd() == '/Users/friedman/Desktop/mnt':
		pathPrefix = os.getcwd()

	cases = set(signatures[signatures['Tumor_Sample_Barcode'].isin(cases)]['Tumor_Sample_Barcode'])
	#Save time by ignoring cases that dont have enough mutations
	#cases = set(signatures[(signatures['Nmut'] >= minNmut) & (signatures['Tumor_Sample_Barcode'].isin(cases))]['Tumor_Sample_Barcode'])

	signatureColumnNames = ['Signature.' + str(i) for i in range(1,31)]
	renameDict = dict([('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)])
	signatures = signatures.rename(columns=renameDict)
	signatures['dominantSignature'] = signatures.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=None, prefix='Signature.', notEnoughMuts= True), axis=1)

	#TODO MAKE A SPECIfic nmut dict based on the appropriate impact version

	cntr = 0
	nCases = 0
	for case in cases:
		cntr += 1
		if cntr%100 == 0: print 'analyzing ', cntr, 'th case out of ', len(cases), ' cases' 
		
		impactVersion = get_impact_version_given_case_name(case)
		if impactVersion == None:
			pass
		elif case not in nmutDict:
			pass
		else:
			nCases +=1

			caseSigs = signatures[signatures['Tumor_Sample_Barcode'] == case]
			dominantSignature = caseSigs['dominantSignature'].iloc[0]
			msiClass = caseSigs['caseMsiClass'].iloc[0]
			nmutInCase = nmutDict[case]
			signatureColsOnly = caseSigs[signatureColumnNames]

			decompositionDict = signatureColsOnly.to_dict(orient='records')[0]
			quadNucFractions = get_quadnuc_fracs_given_decomposition(decompositionDict, spectraPath = pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
			perSnpOncogenicRate = get_expected_oncogenic_val_given_quadnuc_fractions(quadNucFractions, oncogenicSusceptibilityDict, impactVersion)
			
			#set up indel oncogenicity based on msi stable/instable
			indelOncogenicity = NON_MSI_INDEL_ONCOGENICITY
			if msiClass == 'Instable':
				indelOncogenicity = MSI_INDEL_ONCOGENICITY

			#treat cases with insufficient impact mutations differently--assume a constant indel rate 
			#also assume a constant change of an oncogenic snp
			if dominantSignature == 'insufficientMutBurden':
				oncogenicRatePerMut = IMPACT_INDEL_FRAC*indelOncogenicity + (1- IMPACT_INDEL_FRAC)*ONCOGENIC_MUTATION_CHANCE_IN_IMPACT
			else:
				caseIndelFrac = 0
				if case in caseIndelFracDict: #take the observed fraction of indels in a case
					caseIndelFrac = caseIndelFracDict[case]
				oncogenicRatePerMut = caseIndelFrac*indelOncogenicity + (1- caseIndelFrac)*perSnpOncogenicRate
				
			expectedNmutOnc = nmutInCase*oncogenicRatePerMut
			nOncogenicMuts = 0
			nUniqueOncogenicMuts = 0
			if case in nmutOncogenicDict:
				nOncogenicMuts = nmutOncogenicDict[case]
			if case in perCaseUniqueOncogenicInfo:
				nUniqueOncogenicMuts = perCaseUniqueOncogenicInfo[case]
			listOfDicts.append({'Nmut': nmutInCase, 'Expected Oncogenic': expectedNmutOnc, 'Observed Oncogenic': nOncogenicMuts, 'Unique Oncogenic': nUniqueOncogenicMuts,
				'ObservedMinusExpected': nOncogenicMuts - expectedNmutOnc, 'ObservedUniqueMinusExpected': nUniqueOncogenicMuts - expectedNmutOnc,
				'Tumor_Sample_Barcode': case, 'IMPACT_VERSION': impactVersion, 'MSI_STATUS': msiClass})
	df = pd.DataFrame(listOfDicts)
	df['nCases'] = nCases
	return df

def get_observed_expected_values_for_genes(cases, #the tumor sample barcodes of cases we ought to consider
	nmutDict,
	signatures, #the decomposed signatures of the cases we ought to consider
	oncogenicSusceptibilityDictByGene, #a dictionary mapping quadnuc muts in each gene to oncogenic probabilities
	nmutOncogenicInGeneDict #a dictionary that tells us how many oncogenic mutations we observe in a gene in a cohort
	):
	
	listOfDicts = []
	#TELL US WHICH GENES ARE IN WHICH PANEL AND WHAT THEY DO
	impactPanelGenesDict = {
	'IMPACT_368': set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1']),
	'IMPACT_410': set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2']),
	'IMPACT_468': set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'AMER1', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'GTF2I', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MLL4', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SETD8', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	}
	tumorSuppressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
	
	signatureColumnNames = ['Signature.' + str(i) for i in range(1,31)]
	renameDict = dict([('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)])
	signatures = signatures.rename(columns=renameDict)
	signatures['dominantSignature'] = signatures.apply(lambda row: mutationSigUtils.get_dominant_signature(row.to_dict(), cols=None, prefix='Signature.', notEnoughMuts= True), axis=1)

	cntr = 0
	for gene in impactPanelGenesDict['IMPACT_468']:
		cntr += 1
		print cntr

		indelOncogenicity = 0
		tumorSuppressor = False
		if gene in tumorSuppressors:
			indelOncogenicity = .95 #this is the proportion of tumor-suppresor indels that are oncogenic
			tumorSuppressor = True
		#TODO consider only cases that are with the right impavt psnel
		nExpectedMuts = 0.0
		for case in cases:

			if case in set(signatures['Tumor_Sample_Barcode']): #we can only look at cases that are in our signatures data

				impactVersion = get_impact_version_given_case_name(case)
				if impactVersion == None:
					pass
				elif case not in nmutDict:
					pass
				else: #get information about the case from 
					if gene in impactPanelGenesDict[impactVersion]:
						caseSigs = signatures[signatures['Tumor_Sample_Barcode'] == case]
						signatureColsOnly = caseSigs[signatureColumnNames]
						dominantSignature = caseSigs['dominantSignature'].iloc[0]
						nmutInCase = nmutDict[case]

						decompositionDict = signatureColsOnly.to_dict(orient='records')[0]
						quadNucFractions = get_quadnuc_fracs_given_decomposition(decompositionDict, spectraPath = '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
				
						nExpectedMuts = nExpectedMuts + nmutInCase*get_expected_oncogenic_val_given_quadnuc_fractions_and_gene(quadNucFractions, oncogenicSusceptibilityDictByGene, impactVersion, gene)
					else:
						pass #there is no point calculating observed expected for genes that arent in the panel
		print nExpectedMuts

		nmutOncogenic = 0
		if gene in nmutOncogenicInGeneDict:
			nmutOncogenic = nmutOncogenicInGeneDict[gene]
		listOfDicts.append({
				'gene': gene,
				'nExpectedMuts': nExpectedMuts,
				'isTumorSupressor': tumorSuppressor,
				'nObservedMuts': nmutOncogenic
		})

	df = pd.DataFrame(listOfDicts)
	return df"""


