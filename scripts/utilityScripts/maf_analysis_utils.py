#written by Noah Friedman
#funcions for performing analyses and functions on mafs 
import sys
import argparse
import os
import pandas as pd
import numpy as np
import re

from collections import Counter
sys.path.append('/ifs/work/taylorlab/friedman/')

pathPrefix = ''
if os.getcwd() == '/Users/friedman/Desktop/mnt':
	pathPrefix = '/Users/friedman/Desktop/mnt'

import imp
#analysis_utils = imp.load_source('analysis_utils', pathPrefix + '/ifs/work/taylorlab/friedman/myUtils/analysis_utils.py')
#ignore the path prefix babadookery
#analysis_utils = imp.load_source('analysis_utils', '/Users/friedman/Desktop/mnt/ifs/work/taylorlab/friedman/myUtils/analysis_utils.py')


#OMNIBUS maf prepping function (adds tid, pid etc)
def prep_maf(df,
			sampleNameCol='Tumor_Sample_Barcode',
			pid=True, tid=True, quadNuc=True, cancerType=True,
			secondSig=None):

	if tid: df['tid'] = df[sampleNameCol].apply(lambda x: x[:13])
	if pid: df['pid'] = df[sampleNameCol].apply(lambda x: x[:9])
	if quadNuc:
		df['quadNuc'] = df.apply(lambda row: 
   			mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'])
    		if row['Variant_Type'] == 'SNP' else None, axis=1)
	if cancerType:
		cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = pathPrefix +'/ifs/work/taylorlab/friedman/msk-impact/msk-impact/data_clinical_sample.txt')
		df['cancer_type'] = df['pid'].apply(lambda x: cDict[x] if x in cDict else None)
	
	if secondSig != None: #mark the second most common signature
		pass
	return df

def get_nmut_mb_from_impact_id(case, nmut):
	captureAreaDict = {'IM3': 896665, 'IM5':1016478, 'IM6':1139322}
	nmut_Mb = None
	imVersion = None
	if 'IM3' in case: imVersion = 'IM3'
	elif 'IM5' in case: imVersion = 'IM5'
	elif 'IM6' in case: imVersion = 'IM6'
	if imVersion != None:
		nmut_Mb = (1000000.0*nmut)/captureAreaDict[imVersion]
	return nmut_Mb

####little utilities
def fix_mll_genes(maf):
    maf['Hugo_Symbol'] = maf['Hugo_Symbol'].apply(lambda x:
        'KMT2A' if x == 'MLL'
        else 'KMT2B' if x == 'MLL2'
        else 'KMT2C' if x == 'MLL3'
        else x)   
    return maf

##quick command line based tools for counting information in a maf
def count_n_oncogenic_muts_in_maf(mafPath):
	cmd = 'grep -n Oncogenic ' + mafPath + ' | wc -l'
	process = os.popen(cmd)
	return int(process.read())

#quick way to get nmuts for maf
def get_per_case_mut_info(nmutDfPath = '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/nmutInfo_impact_filtered.tsv'):
    df = pd.read_table(nmutDfPath)
    return dict(zip(df['Tumor_Sample_Barcode'], df['Nmut']))

####CNA MAF ANALYSIS UTILITIES
def mark_cases_with_flat_genomes(maf):
	flatGenomeSet = set()
	for case in set(maf['Tumor_Sample_Barcode']):
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]

		if caseMaf[caseMaf['ccf_Mcopies'].notnull()].shape[0] == 0: #CATCH ALL THE PURITY=NA flat genomes cases
			flatGenomeSet.add(case)
	maf['FlatGenome'] = maf['Tumor_Sample_Barcode'].apply(lambda x: True if x in flatGenomeSet else False)
	return maf

def mark_cases_with_median_vaf_of_case(maf):
	vafMapping = dict()
	for case in set(maf['Tumor_Sample_Barcode']):
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
		medianVaf = np.nanmedian(caseMaf['t_var_freq'])
		vafMapping[case] = medianVaf
	maf['medianVaf'] = maf['Tumor_Sample_Barcode'].apply(lambda x: vafMapping[x])
	return maf

#marks each mutation with the median clonal VAF of the case (a pseudo purity value)
def mark_cases_with_median_clonal_vaf_of_case(maf):
	vafMapping = dict()
	cntr = 0
	for case in set(maf['Tumor_Sample_Barcode']):
		if cntr %100 == 0: print cntr
		cntr +=1
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
		clonalThresh = .8 #the threshold to use for ccf mcopies for clonal vs subclonal
		caseMaf['clonal'] = caseMaf['ccf_Mcopies'].apply(lambda x: 1 if x > clonalThresh else 0)
		caseMafClonal = caseMaf[caseMaf['clonal'] == 1]
		medianVaf = np.nanmedian(caseMafClonal['t_var_freq'])
		vafMapping[case] = medianVaf
	maf['medianClonalVaf'] = maf['Tumor_Sample_Barcode'].apply(lambda x: vafMapping[x])
	return maf

#BEST most efficient way to quickly count up the number of times a value occurs
def get_n_oncogenic_muts_per_case_dict(maf):
	oncoKbOncogenicAnnotations = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	oncoMaf = maf[maf['oncogenic'].isin(oncoKbOncogenicAnnotations)]
	return dict(oncoMaf['Tumor_Sample_Barcode'].value_counts())

#returns the fraction of mutations in a cohort that are indels
def get_indel_frac_for_cohort(maf):
	indelClassificationNames = set(['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins'])
	return 1.0*maf[maf['Variant_Classification'].isin(indelClassificationNames)].shape[0]/maf.shape[0]

#a utility function designed to add information about hotspots (which hotspot, what 4nuc etc) to a df with one entry per case
def add_hotspot_maf_info_to_df(df, mafDf,
 allHotspots, hotspotsToFocusOn, #a set of all hospots and a set of hotspots to focus on (those not to call 'other')
 idCol = 'pid', removeDups=True):
	hotspotsWeCareAbout = mafDf[mafDf['Amino_Acid_Change'].isin(allHotspots)]
	if removeDups:
		hotspotsWeCareAbout = hotspotsWeCareAbout.drop_duplicates(subset=[idCol], keep=False)
	hotspotLabelDict = dict(zip(hotspotsWeCareAbout[idCol], hotspotsWeCareAbout['Amino_Acid_Change']))

	hotspotFourNucDict = dict(zip(hotspotsWeCareAbout[idCol], hotspotsWeCareAbout['quadNuc']))
	df['Amino_Acid_Change'] = df[idCol].apply(lambda x: hotspotLabelDict[x] if x in hotspotLabelDict else None)
	df['Amino_Acid_Change_Adj'] = df[idCol].apply(lambda x:
		hotspotLabelDict[x] if (x in hotspotLabelDict and hotspotLabelDict[x] in hotspotsToFocusOn) else 'otherHotspot' if x in hotspotLabelDict else None)
	df['quadNuc'] = df[idCol].apply(lambda x: hotspotFourNucDict[x] if x in hotspotFourNucDict else None)
	return df

###############ADDING COLUMNS

#a utility function that we use to add a column to a dataframe that summarizes True/False for the column X which is true if its hotspot/oncogenic/truncating
def add_mut_effect_summary_col(mafDf):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	mafDf['mutKnowOncogenicOrTruncating'] = mafDf.apply(lambda row: 
		True if row['is-a-hotspot'] == 'Y'
		else True if row['oncogenic'] in oncogenicMutColNames
		else True if row['Consequence'] == 'stop_gained'
		else True if row['Consequence'] == 'frameshift_variant' #TODO is this something that should be included
		else False, 
		axis=1)
	return mafDf

#UTILITIES FOR ASSIGING mutations to the signature that most likely caused them
#creates the reference four nucleotide context for signatures

#df['quadNuc'] = df.apply(lambda row: maf_analysis_utils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
def create_reference_four_nuc(refTri, refAllele, altAllele, variantType):
	#properly invert when needed
	def invert_allele_and_ref_tri(altAllele):
		nucleotideDict = {'A': 'T', 'G': 'C', 'C':'G', 'T':'A'}
		return nucleotideDict[altAllele]

	if variantType != 'SNP': return None  #there is no reference trinuc for non snps
	if not isinstance(refTri, basestring): return None #if the ref tri is not a string we better return none
	if len(refTri) < 3: return None # if the ref tri is less than length 3 (at the end of an exon), we cant do anything
	refTri = str(refTri) #just make sure the ref tri is a string here to avoid funny business
	refAlleleFromReftri = refTri[1]
	alt = altAllele
	if refAlleleFromReftri != refAllele:
		alt = invert_allele_and_ref_tri(altAllele)
	quadNuc = refTri[:2] + alt + refTri[2]
	return quadNuc



#Marks mutations in cases with multiple mutations as shared
def mark_private_vs_shared_mutations(maf):

	def combineDicts(bigDict, littleDict): 
		for key, value in littleDict.items():
			bigDict[key] = value
		return bigDict

	maf['varUuid'] = maf.apply(lambda row: 
		str(row['Start_Position']) + '_' + row['Hugo_Symbol'] + '_' + row['Tumor_Seq_Allele2'] + '_' + row['pid'], axis=1)
	d = dict()
	for pid in set(maf['pid']):

		patientMaf = maf[maf['pid'] == pid]
		
		nMutsDict = dict(patientMaf['varUuid'].value_counts())

		d = combineDicts(d, nMutsDict)
	maf['isSharedMut'] = maf['varUuid'].apply(lambda x: True if d[x] == 2 else False)
	return maf

#marks private mutations by how many alterations already happened in the gene
#NOT a fully fleshed out function

def mark_mutations_by_gene_mut_type(maf, implicationDictStrong, implicationDictWeak):

	maf['geneMutType'] = maf.apply(lambda row:
                                            'strongly_recurrent' if row['Hugo_Symbol'] in implicationDictStrong[row['cancer_type']]
                                             else 'weakly_recurrent' if row['Hugo_Symbol'] in implicationDictWeak[row['cancer_type']]
    										 else 'not_recurrent', axis=1
      )
	return maf

def mark_mutations_by_nth_alteration_in_gene(maf, oncogenicOnly = True):
	for case in set(maf['Tumor_Sample_Barcode']):
		caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]

		print case

		beforeMaf = caseMaf[caseMaf['isSharedMut'] == True]
		beforeMafOncogenic = beforeMaf[beforeMaf['oncogenic'].notnull()]
		afterMaf = caseMaf[caseMaf['isSharedMut'] == False]
		afterMafOncogenic = afterMaf[afterMaf['oncogenic'].notnull()]

		for gene in set(afterMafOncogenic['Hugo_Symbol']):
			oncAfterGeneMaf = afterMafOncogenic[afterMafOncogenic['Hugo_Symbol'] == gene]
			oncBeforeGeneMaf = beforeMafOncogenic[beforeMafOncogenic['Hugo_Symbol'] == gene]
			print gene, oncBeforeGeneMaf.shape[0], oncAfterGeneMaf.shape[0] 
		print '________________'
	return 0


#a function that enumerates all activating mutations in a gene across a cohort maf
def enumerate_activating_muts_across_cohort(gene, mafDf = None):

	#TODO appropriately allow the user to load the maf df if there is no default
	#if mafDf.sh:
	#	mafDf = pd.read_table('/ifs/work/taylorlab/friedman/myAdjustedDataFiles/annotatedOncoPlusHotspotMafAllImpact_trinuc')

	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	geneOncogenicOrHotspotMuts = mafDf[(mafDf['Hugo_Symbol'] == gene) &((mafDf['is-a-hotspot'] == 'Y') |(mafDf['oncogenic'].isin(oncogenicMutColNames)))]

	#colsToKeep = ['Tumor_Sample_Barcode', 'is-a-hotspot', 'oncogenic', 'Ref_Tri', 'Tumor_Seq_Allele2', 'Reference_Allele']
	#return geneOncogenicOrHotspotMuts[[colsToKeep]]
	return geneOncogenicOrHotspotMuts

#a function that given a maf of mutaitons and quadnuc spectra enriched for a signature returns the number of mutations that occur at the enriched spectra
def summarize_signature_attribution_for_case(mafDf, enrichedSigMotifs):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])

	listOfDicts = []
	cases = set(mafDf['Tumor_Sample_Barcode'])
	for case in cases:
		localD = {}
		caseDf = mafDf[mafDf['Tumor_Sample_Barcode'] == case]

		#get a maf for the current cases hotspot and oncogenic muts 
		hotspotMutDf = caseDf[caseDf['is-a-hotspot'] == 'Y']
		oncogenicMutDf = caseDf[caseDf['oncogenic'].isin(oncogenicMutColNames)]
		nHotspotMutations = hotspotMutDf.shape[0]
		nOncogenicMutations = oncogenicMutDf.shape[0]

		#get data about how many of each occurs at a motif
		nOncogenicMutationsAtEnrichedMotif = oncogenicMutDf[oncogenicMutDf['quadNuc'].isin(enrichedSigMotifs)].shape[0]
		nHotpsotMutationsAtEnrichedMotif = hotspotMutDf[hotspotMutDf['quadNuc'].isin(enrichedSigMotifs)].shape[0]
		nMut = caseDf.shape[0]

		#append info to the local dictionary (row of the future df)
		localD['Tumor_Sample_Barcode'] = case
		localD['nHotspots'] = nHotspotMutations
		localD['Nmut'] = nMut
		localD['nOncogenicMutations'] = nOncogenicMutations
		localD['nOncogenicMutationsAtEnrichedMotif'] = nOncogenicMutationsAtEnrichedMotif
		localD['nHotpsotMutationsAtEnrichedMotif'] = nHotpsotMutationsAtEnrichedMotif  

		listOfDicts.append(localD)

	df = pd.DataFrame(listOfDicts)
	return df      

#returns a dataframe with the cases that have mutations in a gene and the type of mutation in that gene
#note the way we return it is designed for the R long format
def asses_per_case_mut_info_for_gene(mafDf, gene, quadnucSet):

	def classify_mut_residue_data(quadNucs, qNucSet):
		mutType = None
		if len(quadNucs) == 1: #if there is only one quadnuc mutation figure out which type it is
			v = quadNucs.pop()
			if v in qNucSet: return 'favoredMutation'
			else: return 'notFavoredMutation'
		else: #if there are more than one mutations in the gene we need to classify whether the mutations are mixed, from the favored process only or from both
			nFavoredMutations = 0
			nNotFavoredMutations = 0
			for v in quadNucs:
				if v in qNucSet:
					nFavoredMutations += 1
				else:
					nNotFavoredMutations += 1
				if nFavoredMutations >= 1 and nNotFavoredMutations >= 1:
					return 'mixed'
				elif nFavoredMutations >= 1:
					return 'mutlipleFavoredMutation'
				else:
					return 'multipleNotFavoredMutation'

	cases = set(mafDf['Tumor_Sample_Barcode'])
	geneMuts = mafDf[mafDf['Hugo_Symbol'] == gene]
	#we only put information in the dataframe we return if there are snp muts
	listOfDicts = []
	for case in cases:
		caseMuts = geneMuts[geneMuts['Tumor_Sample_Barcode'] == case]
		caseQuadNucs = set(caseMuts['quadNuc'])
		if len(caseQuadNucs) != 0:
			if None not in caseQuadNucs:
				classification = classify_mut_residue_data(caseQuadNucs, quadnucSet)
				localD = dict()
				localD['gene'] = gene
				localD['Tumor_Sample_Barcode'] = case
				localD['mutClassification'] = classification
				localD['nGeneMut'] = len(caseQuadNucs)
				listOfDicts.append(localD)
	return pd.DataFrame(listOfDicts)

#a utility that asseses the SNP burden for each case
def asses_snp_burden_across_cohort(maf):
	cases = set(maf['Tumor_Sample_Barcode'])
	cntr = 0
	listOfDicts = []
	for case in cases:
		if cntr%500 == 0: print cntr, len(cases)
		localDict = dict()
		cntr += 1

		caseMuts = maf[maf['Tumor_Sample_Barcode'] == case]
		caseSnps = caseMuts[caseMuts['Variant_Type'] == 'SNP']
		caseIndels = caseMuts[(caseMuts['Variant_Type'] == 'INS') | (caseMuts['Variant_Type'] == 'DEL')]
		localDict['Tumor_Sample_Barcode'] = case
		localDict['nSnps'] = caseSnps.shape[0]
		localDict['nIndels'] = caseIndels.shape[0]
		localDict['nMuts'] = caseMuts.shape[0]
		oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
		localDict['nOncogenicMutations'] = caseMuts[caseMuts['oncogenic'].isin(oncogenicMutColNames)].shape[0]
		listOfDicts.append(localDict)

	return pd.DataFrame(listOfDicts)

def summarize_nmut_info_across_cohort(maf):
	listOfDicts = []
	nmutDict = dict(maf['Tumor_Sample_Barcode'].value_counts())
	cntr = 0
	for case, nmut in nmutDict.items():
		cntr += 1
		if cntr%1000 == 0:
			print cntr

		nmut_Mb = get_nmut_mb_from_impact_id(case, nmut)
		listOfDicts.append({
				'Nmut_Mb': nmut_Mb,
				'Nmut': nmut,
				'Tumor_Sample_Barcode': case
			})
	return pd.DataFrame(listOfDicts)

#script to create a NMUT/MB estimate based on a filtered impact maf
def calculate_nmut_mb_info_from_filtered_maf(maf, write=True):
	exonicVarClass = set(["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Translation_Start_Site"])
	maf = maf[maf['Variant_Classification'].isin(exonicVarClass)]
	
	listOfDicts = []
	cntr = 0
	captureAreaDict = {'IM3': 896665, 'IM5':1016478, 'IM6':1139322}
	for case in set(maf['Tumor_Sample_Barcode']):
		cntr += 1
		if cntr%1000 == 0: print cntr
		nmut_Mb = get_nmut_mb_from_impact_id(case)
		listOfDicts.append({'Tumor_Sample_Barcode': case, 'Nmut_Mb': nmut_Mb})
	df = pd.DataFrame(listOfDicts)

	#ADD CANCER TYPE INFO AND PROPERLY RENAME IT
	renameMapping = {'Pleural Mesothelioma, Epithelioid Type': 'Mesothelioma','Breast Invasive Ductal Carcinoma': 'Breast Cancer','Bladder Urothelial Carcinoma': 'Bladder Cancer','Upper Tract Urothelial Carcinoma': 'Bladder Cancer',
	'Colon Adenocarcinoma': 'Colorectal Cancer', 'Glioblastoma Multiforme': 'Glioma','Adenocarcinoma of the Gastroesophageal Junction': 'Esophagogastric Cancer','Pancreatic Neuroendocrine Tumor': 'Pancreatic Cancer',
	'Endometrial Carcinoma': 'Endometrial Cancer','Stomach Adenocarcinoma': 'Esophagogastric Cancer','Rectal Adenocarcinoma': 'Colorectal Cancer','High-Grade Serous Ovarian Cancer': 'Ovarian Cancer','Breast Invasive Lobular Carcinoma': 'Breast Cancer',
	'Oligodendroglioma': 'Glioma','Serous Ovarian Cancer': 'Ovarian Cancer','Prostate Adenocarcinoma': 'Prostate Cancer','Breast Invasive Carcinoma, NOS': 'Breast Cancer','Esophageal Adenocarcinoma': 'Esophagogastric Cancer',
	'Invasive Breast Carcinoma': 'Breast Cancer','Pancreatic Adenocarcinoma': 'Pancreatic Cancer','Uterine Endometrioid Carcinoma': 'Endometrial Cancer','Colorectal Adenocarcinoma': 'Colorectal Cancer','Mucinous Adenocarcinoma of the Colon and Rectum': 'Colorectal Cancer'
	}

	df['pid'] = df['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
	cDict = analysis_utils.get_cancer_type_information(cancerTypeDfPath = analysis_utils.path_fix('/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_clinical_sample.txt'))
	df['cancer_type'] = df['pid'].apply(lambda x: cDict[x] if x in cDict else None)
	df['cancer_type'] = df['cancer_type'].apply(lambda x: renameMapping[x] if x in renameMapping else x)

	if write:
		df.to_csv(analysis_utils.path_fix('/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/tmbInfo.tsv'), index=False, sep='\t')

	return df


#a utility that tells us the number of mutations per case in the filtered maf
def generate_filtered_mut_per_case_dict(filteredMafDir = '/ifs/work/taylorlab/friedman/mskImpactAsOfMarch2019/dmp/mskimpact/data_mutations_extended.txt'):
	df = pd.read_table(filteredMafDir)
	print df
	return 0

#a utility function that enumerates the top N genes with the most oncogenic mutations in distinct samples
#also enumerates fraction of cases with muts in each case
#returns three things:
#i. A counter of N cases with oncogenic muts
#ii. A dict with gene name: cohort ranking
#iii. A dict mapping gene to fraction of cohort with an oncogenic mutation in that case
def enumerate_top_n_oncogenic_mutated_genes_across_cohort(cohortMaf, n=None):

	#Add the pid column if its not already there
	if 'pid' not in cohortMaf.columns.values:
		cohortMaf['pid'] = cohortMaf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])

	oncoKbOncogenicAnnotations = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	oncogenicMutations = cohortMaf[cohortMaf['oncogenic'].isin(oncoKbOncogenicAnnotations)]
	oncogenicMutations['patientGeneMutated'] = oncogenicMutations.apply(lambda row: row['pid'] + '_' + row['Hugo_Symbol'], axis=1)
	oncogenicMutationsSansPatientDuplicates = oncogenicMutations.drop_duplicates(subset=['patientGeneMutated']) #NOTE this analysis isnt perfect as it treats independent primaries as genetically related

	nPatients = len(set(cohortMaf['pid']))	
	occurenceCounter = None
	if n!= None:
		occurenceCounter = Counter(oncogenicMutationsSansPatientDuplicates['Hugo_Symbol']).most_common(n)
	else:
		occurenceCounter = Counter(oncogenicMutationsSansPatientDuplicates['Hugo_Symbol']).most_common()

	fractionalDict = {key: value for (key, value) in occurenceCounter}
	for key, value in fractionalDict.items():
		fractionalDict[key] = 1.0*value/nPatients

	rankingDict = dict()
	cntr = 1 #counter is 1 indexed is that a problem
	for gene in occurenceCounter:
		rankingDict[gene[0]] = cntr
		cntr += 1

	return occurenceCounter, rankingDict, fractionalDict

#MAKEs a summary of oncogenic mutation info for a cohort for each case, namely the cases mutation burden, n related and unrelated mutations, n second hit mutations, and n second hit mutations in related genes

#We use two thresholds for ???
THRESH_FOR_RECURRENT_MUTATION_BIG = 0.1
THRESH_FOR_RECURRENT_MUTATION_SMALL = 0.01
def summarize_oncogenic_mutation_info(maf, pathPrefix):
	#ENUMERATE specific cancer types we will look at
	maf = maf[maf['oncogenic'].notnull()] #currently only look at oncogenic mutations

	cTypes = set([re.sub('_', ' ', x.strip('.tsv')) for x in os.listdir(pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds')])
	cTypes.remove('Cancer of Unknown Primary')
	print 'making dict'

	#temp remove
	#cTypes = set(['Endometrial Cancer'])
	maf = maf[maf['cancer_type'].isin(cTypes)]

	cancerTypeImplicationDictBig = create_dictionary_mapping_genes_to_cancer_types_with_implication(maf, pathPrefix, cancerTypes=cTypes, t=THRESH_FOR_RECURRENT_MUTATION_BIG)
	#make a dictionary of 'semi-related genes'
	cancerTypeImplicationDictSmall = create_dictionary_mapping_genes_to_cancer_types_with_implication(maf, pathPrefix, cancerTypes=cTypes, t=THRESH_FOR_RECURRENT_MUTATION_SMALL) 
	for key, value in cancerTypeImplicationDictBig.items():
		curCTypeGenes = cancerTypeImplicationDictSmall[key]
		cancerTypeImplicationDictSmall[key] = curCTypeGenes - value

	#TEMP NOAH ALERT ALERT
	#nMutDict = get_per_case_mut_info(nmutDfPath = pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/nmutInfo_impact_filtered.tsv')

	#print cancerTypeImplicationDict
	maf['stronglyRelated'] = maf.apply(lambda row: None if row['cancer_type'] not in cTypes 
		else True if row['Hugo_Symbol'] in cancerTypeImplicationDictBig[row['cancer_type']] 
		else False, axis=1)
	maf['weaklyRelated'] = maf.apply(lambda row: None if row['cancer_type'] not in cTypes 
		else True if row['Hugo_Symbol'] in cancerTypeImplicationDictSmall[row['cancer_type']] 
		else False, axis=1)

	maf['geneAndCase'] = maf.apply(lambda row: row['Hugo_Symbol'] + '_' + row['Tumor_Sample_Barcode'], axis=1) #used to let us know how many times something is mutated
	mafNoMultiples = maf.drop_duplicates(subset=['geneAndCase'])

	nMutStronglyRelatedDict = dict(maf[maf['stronglyRelated'] == True]['Tumor_Sample_Barcode'].value_counts())
	nMutWeaklyRelatedDict = dict(maf[maf['weaklyRelated'] == True]['Tumor_Sample_Barcode'].value_counts())
	nMutUnrelatedDict = dict(maf[(maf['stronglyRelated'] == False) & (maf['weaklyRelated'] == False )]['Tumor_Sample_Barcode'].value_counts())
	nStronglyRelatedGenesAlteredDict = dict(mafNoMultiples[mafNoMultiples['stronglyRelated'] == True]['Tumor_Sample_Barcode'].value_counts())
	nWeaklyRelatedGenesAlteredDict = dict(mafNoMultiples[mafNoMultiples['weaklyRelated'] == True]['Tumor_Sample_Barcode'].value_counts())

	#LONG OVER THE TOP WAY TO DO THINGS TO SET UP THE DF TODO MAKE THIS CODE BETTER
	listOfDicts = []
	for case in set(maf['Tumor_Sample_Barcode']):
		
		nmut = 0
		if case in nMutDict: nmut = nMutDict[case]
		
		nMutStronglyRelated = 0
		if case in nMutStronglyRelatedDict: nMutStronglyRelated = nMutStronglyRelatedDict[case]

		nMutWeaklyRelated = 0
		if case in nMutWeaklyRelatedDict: nMutWeaklyRelated = nMutWeaklyRelatedDict[case]

		nMutUnrelated = 0
		if case in nMutUnrelatedDict: nMutUnrelated = nMutUnrelatedDict[case]

		nStronglyRelatedGenesAltered = 0
		if case in nStronglyRelatedGenesAlteredDict: nStronglyRelatedGenesAltered = nStronglyRelatedGenesAlteredDict[case]

		nWeaklyRelatedGenesAltered = 0
		if case in nWeaklyRelatedGenesAlteredDict: nWeaklyRelatedGenesAltered = nWeaklyRelatedGenesAlteredDict[case]

		#TODO add cancer type
		listOfDicts.append({
			'Tumor_Sample_Barcode': case, 'Nmut': nmut, 
			'nMutStronglyRelated': nMutStronglyRelated, 'nMutWeaklyRelated': nMutWeaklyRelated, 'nMutOncogenic': nMutStronglyRelated + nMutWeaklyRelated + nMutUnrelated,
			'nMutUnrelated': nMutUnrelated, 
			'nStronglyRelatedGenesAltered': nStronglyRelatedGenesAltered, 'nWeaklyRelatedGenesAltered': nWeaklyRelatedGenesAltered,
			'nStronglyRelatedMultipleMutations': nMutStronglyRelated - nStronglyRelatedGenesAltered, 'nWeaklyRelatedMultipleMutations': nMutWeaklyRelated - nWeaklyRelatedGenesAltered
		})

	df = pd.DataFrame(listOfDicts)
	df['nStronglyRelatedSingleAlterations'] = df.apply(lambda row: row['nMutStronglyRelated'] - row['nStronglyRelatedMultipleMutations'], axis=1)
	df['nWeaklyRelatedSingleAlterations'] = df.apply(lambda row: row['nMutWeaklyRelated'] - row['nWeaklyRelatedMultipleMutations'], axis=1)
	return df

######################LOH/CNA analysis tools

def do_gene_loh_summary(maf, genes=None):
	listOfDicts = []

	if genes == None: 
		genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])

	for case in set(maf['Tumor_Sample_Barcode']):
		localD = {}
		localD = {'Tumor_Sample_Barcode': case}
		for gene in genes:

			mutAtGene = False
			oneMutOncogenic = False
			lohAtGene = False
			caseMaf = maf[maf['Tumor_Sample_Barcode'] == case]
			geneMaf = caseMaf[caseMaf['Hugo_Symbol'] == gene]
			if geneMaf.shape[0] > 0:
				mutAtGene = True
			if geneMaf[geneMaf['oncogenic'].isin(oncogenicMutColNames)].shape[0] > 0:
				oneMutOncogenic = True
			if geneMaf[geneMaf['lcn'] == 0].shape[0] > 0:
				lohAtGene = True
			lohPlusOncogenic = False
			if oneMutOncogenic & lohAtGene:
				lohPlusOncogenic = True

			
			localD[gene + '_mut'] = mutAtGene
			localD[gene + '_loh'] = lohAtGene
			localD[gene + '_oncogenicMut'] = oneMutOncogenic
			localD[gene + '_lohPlusOncogenic'] = lohPlusOncogenic
		listOfDicts.append(localD)

	return pd.DataFrame(listOfDicts)

#a function that takes a cohort maf and comes up with a mutation recurrence rnaking for the entirety of the impact panel
def enumerate_gene_mut_ranking_for_cohort(cohortMaf, mode='oncogneic'):

	genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'KMT2A', 'KMT2B', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])

	mafToAnalyze = cohortMaf
	if mode == 'oncogneic':
		mafToAnalyze = cohortMaf[cohortMaf['oncogenic'].notnull()]
	elif mode == 'hotspot':
		mafToAnalyze = cohortMaf[cohortMaf['is-a-hotspot'] == 'Y']
	elif mode == 'vus':
		mafToAnalyze = cohortMaf[~cohortMaf['oncogenic'].notnull()]
	mutCounts = dict(mafToAnalyze['Hugo_Symbol'].value_counts())

	for gene in genes:
		if gene not in mutCounts:
			mutCounts[gene] = 0

	mutCounts = Counter(mutCounts).most_common(len(genes))
	runningSum = 0
	d = {}
	for key, val in mutCounts:
		runningSum += 1.0/len(mutCounts)
		d[key] = runningSum

	return d
	


#a function that given a cohort maf enumerates which tumor suppressors and oncogenes are reccurently mutated (mutated in greater than 'thresh' fraction of the cohort)
def enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(cohortMaf, thresh=.1):
    tumorSupressors = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
    occurenceCounter, rankingDict, fractionalDict = enumerate_top_n_oncogenic_mutated_genes_across_cohort(cohortMaf, n=50)
    recurrentTumorSupressors = []
    recurrentOncogenes = []
    for key, value in fractionalDict.items():
        if value > thresh:
            if key in tumorSupressors:
                recurrentTumorSupressors.append(key)
            else:
                recurrentOncogenes.append(key) 
    return set(recurrentTumorSupressors), set(recurrentOncogenes)

#creates a dictionary that maps cancer types to the implication of mutations in that cancer type
def create_dictionary_mapping_genes_to_cancer_types_with_implication(maf, pathPrefix='', cancerTypes=None, t=0.05):
	
	d = {}
	for cType in cancerTypes:
		print cType
		cTypeAllIds = analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir=pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType=cType, hypermutantStatus = 'all')
		cTypeAllMuts = maf[maf['Tumor_Sample_Barcode'].isin(cTypeAllIds)]
		cTypeNormalIds = analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir=pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType=cType, hypermutantStatus = 'Normal')
		cTypeNormalMuts = maf[maf['Tumor_Sample_Barcode'].isin(cTypeNormalIds)]
		tumorS, oncoG = enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(cTypeNormalMuts, thresh=t)
		recurrentGenes = tumorS | oncoG
		d[cType] = recurrentGenes
	return d

#An extended version of this dictionary having more than one result
#ALERT
def create_dictionary_mapping_genes_to_cancer_types_with_implication_multiple_thresh(maf, pathPrefix='', cancerTypes=None, t1=0.01, t2 = 0.1):
	
	d = {}
	for cType in cancerTypes:
		print cType
		cTypeAllIds = analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir=pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType=cType, hypermutantStatus = 'all')
		cTypeAllMuts = maf[maf['Tumor_Sample_Barcode'].isin(cTypeAllIds)]
		cTypeNormalIds = analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir=pathPrefix + '/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType=cType, hypermutantStatus = 'Normal')
		cTypeNormalMuts = maf[maf['Tumor_Sample_Barcode'].isin(cTypeNormalIds)]
		tumorS, oncoG = enumerate_recurrently_mutated_tumor_supressors_and_oncogenes(cTypeNormalMuts, thresh=t)
		recurrentGenes = tumorS | oncoG
		d[cType] = recurrentGenes
	return d



def annotate_whether_indel_is_at_msi_site(msiMafPath = '/ifs/work/taylorlab/bielskic/MSI/b37_dmp_microsatellites.vep.IMPACT468.maf'):
	df = pd.read_table(msiMafPath)
	return 0


def mark_genes_in_maf_that_have_multiplets(maf, oncogenicOnly = True):
	oncogenicMutColNames = set(['Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic'])
	if oncogenicOnly: #if we only do the 
		maf = maf[maf['oncogenic'].isin(oncogenicMutColNames)]

    #only mark IMPACT genes
	genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'KMT2A', 'KMT2B', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
    #TODO COME UP WITH A REAL AND EFFECTIVE WAY TO MARK DOUBLETD

def main():

	parser = argparse.ArgumentParser(description='Arg parser for this script')
	parser.add_argument('--argument', help='stub for parser argument', default='')

	args = parser.parse_args()

if __name__ == '__main__':
    main()















