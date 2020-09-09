#written by noah friedman friedman@mskcc.org

import sys
import os
import subprocess
import pandas as pd
import argparse
import re


#run code to get mutational signatures
def run_mutational_signatures_code(scriptDir, outputFilePath, sMaf, signaturesFilePath, triuncOnly=False, tMaf=None):
	def run_triunc_command(scriptDir, sMaf, tMaf):
		#temp
		scriptDir = '/home/friedman/friedman/mutation-signatures'
		triuncScriptPath = os.path.join(scriptDir, 'make_trinuc_maf.py')
		cmd = 'python {sPath} {sourceMafPath} {TargetMafPath}'.format(sPath = triuncScriptPath, sourceMafPath = sMaf, TargetMafPath = tMaf)
		print 'executing Triunc command: ', cmd
		subprocess.Popen(cmd, shell=True).wait()
		print 'Triunc command completed'
		return tMaf

	def run_signatures_code(scriptDir, sMaf, outputFile, signaturesFilePath):
		signaturesScriptPath = os.path.join(scriptDir, 'main.py')
		cmd = 'python {sPath} {sigFilePath} {sourceMafPath} {targetFilePath} --spectrum_output spectrumNoah.txt'.format(sPath = signaturesScriptPath, sigFilePath = signaturesFilePath, sourceMafPath = sMaf, targetFilePath = outputFile)
		print 'executing Signatures command: ', cmd
		subprocess.Popen(cmd, shell=True).wait()
		print 'signatures command completed', cmd

	triuncMaf = ''
	print sMaf
	print 'Sirt'
	if 'triunc' in sMaf: triuncMaf = sMaf
	else:
		triuncMaf = run_triunc_command(scriptDir, sMaf, tMaf)
		if triuncOnly:
			print 'triunc only mode, returning'
			return
	print triuncMaf
	print 'RUMi'
	targetSignaturesFile = outputFilePath
	run_signatures_code(scriptDir, triuncMaf, targetSignaturesFile, signaturesFilePath)

#sanity check that the code works by summing up the counts
def sanity_check_mutation_counts(signaturesFilePath):
	#util function for getting sample names from the key name
	#takes the sample name, and a key, 0 or 1 referring to which value we want
	def extract_sample_name(sampleNameComparison, index):
		if '!' in sampleNameComparison:
			sampleName = sampleNameComparison.split('!')[index]
		elif '&' in sampleNameComparison:
			sampleName = sampleNameComparison.split('&')[index]
		else:
			print 'Major error! sample name improperly formatted'
			sys.exit()
		return sampleName


	def create_sanity_check_dictionary(signaturesDf): #create a dictionary that we will use to do the sanity check
		sanityCheckDict = dict()
		for index, row in signaturesDf.iterrows(): 
			sampleNameComparison = row['Sample Name'] 
			sampleName = extract_sample_name(sampleNameComparison, 0)
			if sampleName in sanityCheckDict:
				curListOfRows = sanityCheckDict[sampleName]
				curListOfRows.append(row)
				sanityCheckDict[sampleName] = curListOfRows
			else:
				sanityCheckDict[sampleName] = [row]

		return sanityCheckDict

	#for each entry in the dictionary we create another dictionary and sanity check it
	def perform_sanity_check_on_dict(sanityCheckDict):
		for key, value in sanityCheckDict.items():
			subDict = dict()
			for v in value:
				sampleNameComparison = v['Tumor_Sample_Barcode']
				sampleName = extract_sample_name(sampleNameComparison, 1)
				if sampleName in subDict:
					l = subDict[sampleName]
					l.append(v['Number of Mutations'])
					subDict[sampleName] = l
				else:
					subDict[sampleName] = [v['Number of Mutations']]
			prevS = -1
			for key1, v in subDict.items():
				v1, v2 = v
				s = v1 + v2
				if prevS != -1 and s != prevS:
					print 'error on key:', key
					print v1, ' + ', v2, '!=', prevS
					sys.exit()
				prevS = s
			print 'sums validated for ', key
	
	#makes logic easier by duplicating and row data--ie if we have 1&2 we also create 2&1 which is identical but makes logic and code easier to manage
	def add_extra_and_rows(df):
		for index, row in df.iterrows():
			if '&' in row['Sample Name']:
				v1, v2 = row['Sample Name'].split('&')
				reversedSampleName = v2 + '&' + v1
				row['Sample Name'] = reversedSampleName
				df = df.append(row)
		return df

	print 'performing sanity check of mutation counts'
	signaturesDf = pd.read_table(signaturesFilePath)
	signaturesDf = signaturesDf[['Sample Name', 'Number of Mutations']]
	signaturesDf = add_extra_and_rows(signaturesDf)
	sanityCheckDict = create_sanity_check_dictionary(signaturesDf)
	perform_sanity_check_on_dict(sanityCheckDict)


def subset_triunc_signature_fractions(signature, nucleotides, spectrumFile='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'):
	mutationSigTable = pd.read_table(spectrumFile)
	mutationSigTable['signature'] = mutationSigTable.index
	curSignature = mutationSigTable[mutationSigTable['signature'] == 'Signature.' + signature]
	s = 0
	for nuc in nucleotides:
		s += float(curSignature[nuc])
	return s

#UTILITIES FOR ASSIGING mutations to the signature that most likely caused them
#creates the reference four nucleotide context for signatures

#df['quadNuc'] = df.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
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

#Creates a string of the form XX(X>Y)ZZ that is strand specific
def create_strand_specific_pentanucleotide_change(refPenta, refAllele, altAllele, variantType):

	complementDict = {'C': 'G', 'G': 'C', 'T': 'A', 'A': 'T'}

	if variantType != 'SNP': return None
	if not isinstance(refPenta, basestring): return None
	if len(refPenta) < 5: return None

	if refAllele != refPenta[2]:
		altAllele = complementDict[altAllele]
	#if refAllele != refPenta[2]: #this means the ref penta is corrected
		#refPenta = refPenta[::-1] #reverse the reference pentanucleotide so that it matches
		#refPenta = ''.join([complementDict[x] for x in refPenta])#and complement it
	return refPenta[:2] + '(' + refPenta[2] + '>' + altAllele + ')' + refPenta[3:]


#assigns the signature that was most likely to create a mutation
def assign_most_likely_mutation(spectrumDict, row, 
	n = 1,  #the number of most likely mutations we would like to return
	signatures=None, #if the argument signatures is specified we use that as the possible signatures
	signaturesPrefix = 'Signature.', #prefix for spec
	returnSpectrumProb=False #if true we return the most likely spectrum prob, not just the spectrum
	): 
	fourNuc = create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'])
	if signatures == None: 
		signatures = [signaturesPrefix + str(i) for i in range(1,31)]
	row = row[signatures]
	rowAsDict = row.to_dict()
	l = []
	for key, value in rowAsDict.items():
		key = re.sub('mean_', 'Signature.', key) #make sure we have the proper key to index into the spectrum dict
		l.append((key, value*spectrumDict[key][fourNuc]))
	signatures = [i[0] for i in l]
	values = [i[1] for i in l]
	sortedL = [x for _,x in sorted(zip(values, signatures), reverse=True)]
	if returnSpectrumProb:
		returnL = [(sig, spectrumDict[sig][fourNuc]) for sig in sortedL[:n]]
		return returnL
	else:
		return sortedL[:n] #returns the n most ocmmon signatures

#Expands a df to include columns for the probabilities of each spectrum
#note the signature columns must not yet be combined together at this stage
def expand_df_to_include_spectrum_probs(df, spectrumDicts, artifSpectra):

	def multiply_function(row, signatureCols, spectrumDs, curTriNuc, artificialSpectra=None):
		curTriNucSpecDict = {key:val for key,val in [(signatureCols[i], spectrumDs['Signature.' + str(i + 1)][curTriNuc]) for i in range(0,30)]}
		rowAsDict = row[signatureCols].to_dict() 
		if row['Nmut'] < 10 and artificialSpectra != None: 
			rowAsDict = artificialSpectra

		l = [rowAsDict[key]*curTriNucSpecDict[key] for key in rowAsDict.keys()]
		return sum(l)

	signatureColumns = ['mean_' + str(i) for i in range(1,31)]
	trinucKeys = spectrumDicts['Signature.1'].keys()
	
	for key in trinucKeys:
		df[key + '_prob'] = df.apply(lambda row: multiply_function(row, signatureColumns, spectrumDicts, key, artifSpectra), axis=1)

	return df

def get_spectrum_probability(spectrumDict, row, dominantSignature):
	fourNuc = create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'])
	return spectrumDict[dominantSignature][fourNuc]

# a little utility function to convert the mutations spectrum file to a dictionary of dictionaries for quicker lookup and access
def convert_spectrum_file_to_dict_of_dicts(spectrumFile='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'):
	d = {}
	df = pd.read_table(spectrumFile)
	for index, row in df.iterrows():
		localD = row.to_dict()
		d[str(index)] = localD
	return d

#a function used to enumerate which spectra are enriched for each signature. This is useful for understanding which spectra we can associate with a signature
def get_enriched_spectra_for_signatures(spectrumFile='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt', spectraSignificanceThresh=.05, pathPrefix='',
	signaturesToIgnore= #ignore signatures we dont care about 
	set(['Signature.5','Signature.8','Signature.9','Signature.12','Signature.16','Signature.19','Signature.22','Signature.23','Signature.24','Signature.25','Signature.27','Signature.28','Signature.29','Signature.30'])):
	
	def filter_entries_by_thresh(spectraDict, thresh): #filters keys for the dictionary sp
		newDictofSets = {}
		for key, subDict in spectraDict.items():
			s = set()
			for subDictKey, value in subDict.items():
				if value > thresh:
					s.add(subDictKey)
			newDictofSets[key] = s
		return newDictofSets


	d = convert_spectrum_file_to_dict_of_dicts(pathPrefix + spectrumFile)
	for keyToDrop in signaturesToIgnore:
		d.pop(keyToDrop, None)

	filteredD = filter_entries_by_thresh(d, spectraSignificanceThresh)
	filteredD['Signature.MMR'] = filteredD['Signature.6'] | filteredD['Signature.15'] | filteredD['Signature.20'] | filteredD['Signature.21'] | filteredD['Signature.26']
	filteredD['Signature.APOBEC'] = filteredD['Signature.2'] | filteredD['Signature.13']
	#pop all the signatures we just combined
	filteredD.pop('Signature.6'), filteredD.pop('Signature.15'), filteredD.pop('Signature.20'), filteredD.pop('Signature.21'), filteredD.pop('Signature.26'), filteredD.pop('Signature.2'), filteredD.pop('Signature.13') 
	return filteredD

def annotate_mutations_with_signatures_in_case(outputFilename, outputDir, mutationsFileToAnnotate='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/maf2mafAnnotatedMay16filteredMafWithIsHotspot.maf',
	signaturesFile='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/mutationSigFiles/may16unfiltered30sigs.txt'
	):

	def convert_signatures_file_to_dict(signaturesFile):
		dictToDict = {}
		df = pd.read_table(signaturesFile)
		dfDictList = df.to_dict(orient='records')
		for row in dfDictList:
			barcode = row['Sample Name']
			del row['Number of Mutations']
			del row['Sample Name']
			dictToDict[barcode] = row
		return dictToDict

	barcodeToSignaturesDict = convert_signatures_file_to_dict(signaturesFile)
	mutationsDf = pd.read_table(mutationsFileToAnnotate)
	dfMutationsList = mutationsDf.to_dict(orient='records')
	listOfDicts = []
	for row in dfMutationsList:
		barcode = row['Tumor_Sample_Barcode']
		if barcode in barcodeToSignaturesDict:
			signaturesInfo = barcodeToSignaturesDict[barcode]
			rowCopy = row.copy()
			rowCopy.update(signaturesInfo)
			listOfDicts.append(rowCopy)
	df = pd.DataFrame(listOfDicts)
	writePath = os.path.join(outputDir, outputFilename)
	print 'writing file to ', writePath
	df.to_csv(writePath, sep='\t', index=False)


def create_limited_spectrum_file(signaturesToInclude, oldSpectrumFile='/ifs/work/taylorlab/friedman/myUtils/newSignatures.txt', outputDir='/ifs/work/taylorlab/friedman/myAdjustedDataFiles/spectrumFiles'):
	spectrumDf = pd.read_table(oldSpectrumFile)
	spectrumDf = spectrumDf.ix[signaturesToInclude]
	writePath = os.path.join(outputDir, 'bladderSignatures.txt')
	print 'writing file to ', writePath
	spectrumDf.to_csv(writePath, index=True, sep='\t')

#function used for defining hotspot mutation percentage for a cohort
#multiplies average pan cohort mutation fraction by each of the trinuc contexts
#then sums all of them
def get_spectrum_mutation_frac_for_cohort(signatureFracs, spectrumFilePath='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'):
	spectrumFile = pd.read_table(spectrumFilePath)
	for i in range(len(signatureFracs)):
		spectrumFile.iloc[i] = spectrumFile.iloc[i].apply(lambda x: 1.0*x*signatureFracs[i])
	return spectrumFile.sum(axis=0)

#####################UTILITIES FOR MERGING MUTATIONAL SIGNATURE COLUMNS

def merge_signature_columns(df, mode='Stratton', drop=True,
	smokingMerge=False, mmrAgingMerge=False, tmzMerge= False, confidence=True, mean=True, prefix='mean_'):
	if mode == 'Stratton':
		if confidence: df['confidence_APOBEC'] = df.apply(lambda row: max(row['confidence_2'], row['confidence_13']), axis=1)
		if mean: df[prefix + 'APOBEC'] = df.apply(lambda row: row[prefix + '2'] + row[prefix + '13'], axis=1)
		if confidence: df['confidence_MMR'] = df.apply(lambda row: max(row['confidence_6'], row['confidence_15'], row['confidence_20'], row['confidence_21'], row['confidence_26']), axis=1)
		if mean: df[prefix + 'MMR'] = df.apply(lambda row: row[prefix + '6'] + row[prefix + '15'] + row[prefix + '20'] + row[prefix + '21'] + row[prefix + '26'], axis=1)	
		if mean:
			if smokingMerge: #smoking merge, if specified merges the smoking signature, mutyh, aflatoxin and chewing tobacco (which are all usually smoking)
				df[prefix + 'SMOKING'] = df.apply(lambda row: row[prefix + '4'] + row[prefix + '18'] + row[prefix + '24'] + row[prefix + '29'], axis=1)
			if mmrAgingMerge: #mmr aging merge: combines signature 1 and signature mmr
				df[prefix + 'AGING/MMR'] = df.apply(lambda row: row[prefix + 'MMR'] + row[prefix + '1'], axis=1)
			if tmzMerge:
				df[prefix + 'TMZ/23'] = df.apply(lambda row: row[prefix + '11'] + row[prefix + '23'], axis=1)
		dropCols = []
		if mean: dropCols += [prefix + '2', prefix + '13', prefix + '6', prefix + '15', prefix + '20', prefix + '21', prefix + '26']
		if smokingMerge: dropCols += [prefix + '4', prefix + '18', prefix + '24', prefix + '29']
		if confidence: dropCols += ['confidence_2', 'confidence_13', 'confidence_6', 'confidence_15', 'confidence_20', 'confidence_21', 'confidence_26']
		if drop: #drop cols if asked
			df = df.drop(dropCols, axis=1)
		return df 

	#TODO implement mean merges for SBS mode
	elif mode == 'SBS':
		df['confidence_APOBEC'] = df.apply(lambda row: max(row['confidence_SBS2'], row['confidence_SBS13']), axis=1)
		df['confidence_MMR'] = df.apply(lambda row: max(row['confidence_SBS6'], row['confidence_SBS15'], row['confidence_SBS20'], row['confidence_SBS21'], row['confidence_SBS26'], row['confidence_SBS44']), axis=1)
		df['confidence_BRCA'] = df.apply(lambda row: max(row['confidence_SBS3'], row['confidence_SBS39']), axis=1)
		df['confidence_UV'] = df.apply(lambda row: max(row['confidence_SBS7a'], row['confidence_SBS7b'], row['confidence_SBS7c'], row['confidence_SBS7d']), axis=1)
		df['confidence_POLE'] = df.apply(lambda row: max(row['confidence_SBS10a'], row['confidence_SBS10b']), axis=1)
		df['confidence_Sig17'] = df.apply(lambda row: max(row['confidence_SBS17a'], row['confidence_SBS17b']), axis=1)
		df = df.drop(['confidence_SBS2', 'confidence_SBS13', 
			'confidence_SBS6', 'confidence_SBS15', 'confidence_SBS20', 'confidence_SBS21', 'confidence_SBS26', 'confidence_SBS44', 
			'confidence_SBS3', 'confidence_SBS39',
			'confidence_SBS7a', 'confidence_SBS7b', 'confidence_SBS7c', 'confidence_SBS7d',
			'confidence_SBS10a', 'confidence_SBS10b',
			'confidence_SBS17a', 'confidence_SBS17b'], axis=1)
		return df
	else:
		print 'error improper mode specified'
		sys.exit()

#utility function to do the ninja work required to give me the adjusted signature names
def get_adjusted_signature_column_names(mode = 'Stratton'):
	if mode == 'Stratton':
		cols = ['mean_' + str(i) for i in range(30)]
		removeCols = ['mean_2', 'mean_13', 'mean_6', 'mean_15', 'mean_20', 'mean_21', 'mean_26']
		for c in removeCols: cols.remove(c)
		cols.append('mean_APOBEC')
		cols.append('mean_MMR')
		return cols
	elif mode == 'SBS':
		print 'uh oh noah was too lazy to implement me'
	else:
		print 'error improper mode specified'
		sys.exit()

#returns the dominant signautre for a row of a df expressed as a dict with the specified signatures under consideration
def get_dominant_signature(rowAsDict, cols=None, prefix='mean', notEnoughMuts= True):
	if notEnoughMuts == True:
		if rowAsDict['Nmut'] < 10: return 'insufficientMutBurden'
	if cols == None:
		cols = get_adjusted_signature_column_names()
	tupList = []
	for key, value in rowAsDict.items():
		if prefix in key: tupList.append((key, value))
	sortedSigs = sorted(tupList, key = lambda tup: tup[1], reverse=True)
	return sortedSigs[0][0]

	#TODO return magnitude etc to help with classification

#rowwise lambda function to find hte domniant signaturefor each case
def find_nth_most_predominant_signature(row, n=1, mode='Name', signaturesToConsider=None):
	if signaturesToConsider == None:
		signaturesToConsider = ['mean_1','mean_10', 'mean_11', 'mean_12', 'mean_14', 'mean_16', 'mean_17', 'mean_18', 'mean_19',
		'mean_22', 'mean_23', 'mean_24', 'mean_25', 'mean_27', 'mean_28','mean_29', 'mean_3', 'mean_30',
		'mean_4', 'mean_5', 'mean_7', 'mean_8', 'mean_9', 'mean_APOBEC', 'mean_MMR']
 	sigs = row[signaturesToConsider]
 	sortedTupleList = sorted(list(sigs.to_dict().items()), key=lambda x: x[1], reverse=True)
 	if mode == 'Name':
 		return sortedTupleList[n-1][0]
 	elif mode == 'Magnitude':
 		return sortedTupleList[n-1][1]

#a rowwise function to find the second most common signature in a case
def find_second_most_common_signature(row, primarySig, returnMode, 
                            sigNamesToSpecify = set(['mean_1', 'mean_3', 'mean_4', 'mean_7', 'mean_10', 'mean_11','mean_14', 'mean_17', 'mean_MMR', 'mean_APOBEC']), #a set of signatures we actually mark on the chart
                            signatureColumns=None
                            ):
    colNames = row.to_dict().keys()
    if signatureColumns == None:
    	signatureColumns = [i for i in list(row.keys()) if 'mean' in i]
    	if len(signatureColumns) == 0: signatureColumns = [i for i in list(row.keys()) if 'Signature' in i]

    rowSigsOnly = row[signatureColumns]
    rowAsDict = rowSigsOnly.to_dict()
    items = rowAsDict.items()
    sortedItems = sorted(items, key=lambda x: x[1], reverse=True)

    l = [list(t) for t in zip(*sortedItems)][1]
    if sum(l) > 1.1: 
    	print sum(l), sortedItems
    	print '______________________'

    if sortedItems[0][0] == primarySig:
        if returnMode == 'name':
            sigName = sortedItems[1][0]
            if sigName in sigNamesToSpecify:
                return sigName
            else:
                return 'other'
        else:
            return sortedItems[1][1]
    else:
        if returnMode == 'name':
            sigName = sortedItems[0][0]
            if sigName in sigNamesToSpecify:
                return sigName
            else:
            	#return sortedItems[0][0]
                return 'other' 
        else:
            return sortedItems[0][1]

def count_quad_nuc_sprectra_from_maf(maf):

	def Merge(dict1, dict2):  #UTILITY function to merge dicitonaries
		z = dict2.copy()
		z.update(dict1)
		return z

	listOfDataframes = []

	if ['Ref_Tri'] not in maf.columns.values:
		print 'error, the maf we use needs to have a column called ref tri'
		return
	allBases = ['A', 'C', 'G', 'T']
	changes = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG'] #format: 'CA' means a change from C>A
	allQuadNucs = [firstBase + change + lastBase for firstBase in allBases for change in changes for lastBase in allBases] #enumerate all 96 quadnucs for signatures
	maf['quadNuc'] = maf.apply(lambda row: create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)
	allCases = set(maf['Tumor_Sample_Barcode'])
	for quadNuc in allQuadNucs:
		qMaf = maf[maf['quadNuc'] == quadNuc]

		#Convluted (but hopefully somewhat more efficient) code to count the number of times a trinucleotide occurs in each case
		valueCountsDict = dict(qMaf['Tumor_Sample_Barcode'].value_counts()) #cases where the trincleotide occurs
		zeroBarcodes = allCases - set(valueCountsDict.keys())
		casesWhereQuadNucDoesNotOccur = {key:value for (key, value) in [(barcode, 0) for barcode in zeroBarcodes]}#create another dictionary to mark other cases as 0
		allCasesDictionary = Merge(valueCountsDict, casesWhereQuadNucDoesNotOccur)
		df = pd.DataFrame(allCasesDictionary.items())
		df = df.rename(columns={0:'Tumor_Sample_Barcode', 1:quadNuc}) #fix the column names of this df we have constructed
		listOfDataframes.append(df) #store all these columns so we can merge them later

	df = listOfDataframes[0]
	for df_ in listOfDataframes[1:]:
		df = df.merge(df_, on='Tumor_Sample_Barcode')
	return df

def main():

	parser = argparse.ArgumentParser(description='Noahs script!')
	parser.add_argument('--inputMaf', help='maf to run signatures/triunc on', default='/ifs/work/taylorlab/friedman/clinicalData/msk-impact/msk-impact/data_mutations_extended.txt')
	parser.add_argument('--outputDir', help='output directory', default='/juno/work/taylorlab/friedman/myAdjustedDataFiles')
	parser.add_argument('--outputFilename', help='output filename', default=None)
	#CHANGE TO MUT SIGNATURES 2019 path!!!!
	parser.add_argument('--mutationalSignaturesScriptPath', help='path to the mutational signatures script', default='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures')
	parser.add_argument('--trinucOnly', help='mode for whether we just generate the triunc only file', default=False)
	parser.add_argument('--mode', help='mode for whether we just generate the triunc only file', default='trinucOnly')
	parser.add_argument('--spectrumFilePath', help='path to the spectrum file', default='/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')  #signaturesFilePath = '/ifs/work/taylorlab/friedman/myUtils/newSignatures.txt'


	args = parser.parse_args()

	print args.mode

	trinucOnly = False
	if args.mode == 'trinucOnly':
		args.mode = 'runMutSigs' #change mode to move with logic (alert confusing logic and dsylexic variable naming)
		trinucOnly = True

	if args.mode == 'annotateMutations':
		annotate_mutations_with_signatures_in_case(args.outputFilename, args.outputDir, mutationsFileToAnnotate=args.inputMaf)

	elif args.mode == 'createSpectrumFile':
		signaturesToInclude = ['SBS1', 'SBS2', 'SBS5', 'SBS10a', 'SBS10b', 'SBS13', 'SBS31', 'SBS35']
		create_limited_spectrum_file(signaturesToInclude)

	elif args.mode == 'runMutSigs':
		args.mutationalSignaturesOutputPath = os.path.join(args.outputDir, ''.join([args.outputFilename, 'mutationalSignatuesOutput.txt']))
		outputPath = os.path.join(args.outputDir, args.outputFilename)
		run_mutational_signatures_code(args.mutationalSignaturesScriptPath, args.mutationalSignaturesOutputPath, args.inputMaf, args.spectrumFilePath, triuncOnly=trinucOnly, tMaf=outputPath)

	else: print 'invalid mode specified', args.mode


if __name__ == '__main__':
    main()


#
#EXAMPLE RUN: python friedman/myUtils/mutationSigUtils.py --inputMaf /juno/work/ccs/gongy/megatron_Jan5th/Result/cohort_level/mut_somatic.maf --outputDir /juno/work/taylorlab/friedman/myAdjustedDataFiles --outputFilename exomeRecaptureSignatures.tsv --mode runMutSigs
