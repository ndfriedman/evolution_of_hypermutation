#written by noah friedman
#a utility used to attribute mutations to the process that caused them


import sys
import pandas as pd
from collections import Counter
import numpy as np

#meant to be run as a lambda mutation on a maf to generate attributions
#NOTE your maf must have a column called 'quadNuc' describining the nucleotide change and motif
def attribute_mutations_to_signatures(maf, mutAttributionDf):
	dictOfDicts = {}
	for case in set(mutAttributionDf['Sample Name']):
		dictOfDicts[case] = dict(mutAttributionDf[mutAttributionDf['Sample Name'] == case].iloc[0])
	possibleQuadNucs = [firstLetter + change + lastLetter for firstLetter in ['A', 'C', 'T', 'G']
	for change in ['CA', 'CG', 'CT', 'TA', 'TG', 'TC'] for lastLetter in ['A', 'C', 'T', 'G']]

	print 'attributing mutations'
	maf['mutAttribution'] = maf.apply(lambda row: dictOfDicts[row['Tumor_Sample_Barcode']][row['quadNuc']]
	 if row['quadNuc'] in possibleQuadNucs and row['Tumor_Sample_Barcode'] in dictOfDicts else None, axis=1)
	return maf

def do_smoking_sig_correction(row, colname):
	if row['signature'] == 'Signature.4':
		return 1
	else:
		return row[colname]

#combines the columns in the signature decomposition
def combine_signature_decomposition_columns(signatureDecomposition):
	d = {}
	dropCols = []
	#change the combo rules format to a dict of lists
	for key, value in signatureComboRules.items():
		if key in signatureDecomposition.columns.values: #if the column to be combined isnt there it means they have been combined upstream
			dropCols.append(key)
			if value in d:
				d[value] = d[value] + [key]
			else:
				d[value] = [key]

	#now sum the decomposition columns by iterating over our new dictionary
	for key, value in d.items():
		signatureDecomposition[key] = sum([signatureDecomposition[v] for v in value])

	#drop the now redudant signature columns
	signatureDecomposition = signatureDecomposition.drop(dropCols, axis=1)
	return signatureDecomposition 

def summarize_quad_nuc_mut_attribution(signatureDecomposition, spectrumDf):
	#figure out signature attribution for cases
	listOfDicts = []
	cntr = 0
	for index, row in signatureDecomposition.iterrows():
		if cntr%100 == 0: print 'analyzing case number ', cntr, ' out of ', len(set(signatureDecomposition['Sample Name'])), ' cases'
		d = {}
		sampleId = row['Sample Name']
		rowAsDict = dict(row)
		del row['Sample Name']
		signaturesPresentInCase = set([v[0] for v in rowAsDict.items() if v[1] == 1]) #get the signatures present in the case
		
		#TEMP COMMENTED THIS OUT
		#if agingIsAlwaysPresent == 'agingIsAlwaysPresent': #add the 
		#	signaturesPresentInCase = signaturesPresentInCase | set(['Signature.1'])
		
		#summarize and add to the list of dictionaries which I turn into a dataframe
		spectrumDfCase = spectrumDf[spectrumDf['signature'].isin(signaturesPresentInCase)]
		d['Sample Name'] = sampleId
		d['signaturesPresentInCase'] = '|'.join(signaturesPresentInCase)
		for quadnuc in quadNucColumns:
			sigs = '|'.join(set(spectrumDfCase[quadnuc]) - set([''])) #join all signatures present with a pipe and get rid of all empty strings
			d[quadnuc] = sigs
		listOfDicts.append(d)
		cntr += 1

	return pd.DataFrame(listOfDicts)


#maps which signatures get mashed together
#please change this to your hearts content
signatureComboRules = {
	'Signature.2': 'Signature.APOBEC',
	'Signature.13': 'Signature.APOBEC',

	'Signature.4': 'Signature.SMOKING',
	'Signature.18': 'Signature.SMOKING',
	'Signature.24': 'Signature.SMOKING',
	'Signature.29': 'Signature.SMOKING',

	'Signature.6': 'Signature.MMR',
	'Signature.15': 'Signature.MMR',
	'Signature.20': 'Signature.MMR',
	'Signature.21': 'Signature.MMR',
	'Signature.26': 'Signature.MMR'
}


def main():
	print 'usage:\n python attribute_mutations_to_signatures.py path/to/mut/spectrum/file thresholdForMutationMotif thresholdNMutForSignature path/to/file/listing/signatures/ output/file/path agingIsAlwaysPresent? doSmokingCorrection[optional]'

	#READ in input
	spectrumFilePath = sys.argv[1]
	inclusionThreshold = float(sys.argv[2])
	thresholdNMutForSignature = float(sys.argv[3])
	signaturesFilePath = sys.argv[4]
	outputFilePath = sys.argv[5]
	agingIsAlwaysPresent = sys.argv[6]
	smokingMode = None
	if len(sys.argv) >= 8:
		smokingMode = sys.argv[7]


	#format spectrum file
	spectrumDf = pd.read_table(spectrumFilePath)
	spectrumDf.reset_index(level=0, inplace=True)
	spectrumDf = spectrumDf.rename(columns={'index': 'signature'})
	quadNucColumns = set(spectrumDf.columns.values) - set(['signature']) #all columns of the spectra file except the signature column

	#BINARIZE THE SPECTRUMS
	for col in quadNucColumns:
		spectrumDf[col] = spectrumDf[col].apply(lambda x: 1 if x >= inclusionThreshold else 0)

	#CHANGE ALL C>A mutations to 1 for signature 4 if the smoking correction is specified
	nucleotides = ['C', 'T', 'G', 'A']
	if smokingMode == 'doSmokingCorrection':
		smokingCols = [firstLetter + 'CA' + lastLetter for firstLetter in nucleotides for lastLetter in nucleotides]
		for col in smokingCols:
			spectrumDf[col] = spectrumDf.apply(lambda row: do_smoking_sig_correction(row, col), axis=1)

	#DO A GROUPBY TO MASH SPECTRA TOGETHER
	spectrumDf['signature'] = spectrumDf['signature'].apply(lambda x: signatureComboRules[x] if x in signatureComboRules else x)
	spectrumDf = spectrumDf.groupby(['signature']).sum()
	#WE have to reset the index again after the group by
	spectrumDf.reset_index(level=0, inplace=True)
	spectrumDf = spectrumDf.rename(columns={'index': 'signature'})
	#re-binarize everything
	for col in quadNucColumns:
		spectrumDf[col] = spectrumDf[col].apply(lambda x: 1 if x > 1 else x)


	#read in signatures input, combine it then binarize it
	signatureDecomposition = pd.read_table(signaturesFilePath)
	#DO a correction if the column names are mean_1, etc
	print signatureDecomposition.columns.values
	if 'mean_1' in signatureDecomposition.columns.values:
		renameDict = dict([(i, 'Signature.' + str(i.strip('mean_'))) for i in signatureDecomposition.columns.values if 'mean_' in i])
		renameDict['Nmut'] = 'Number of Mutations'
		renameDict['Tumor_Sample_Barcode'] = 'Sample Name'
		signatureDecomposition = signatureDecomposition.rename(columns=renameDict)
		dropCols = [i for i in signatureDecomposition.columns.values if 'confidence' in i] + ['impact_version', 'Nmut_Mb']
		if 'dominantSignature' in signatureDecomposition.columns.values: dropCols += ['dominantSignature']
		signatureDecomposition = signatureDecomposition.drop(columns=dropCols)

	signatureDecomposition = combine_signature_decomposition_columns(signatureDecomposition)
	signatureColNames = set(signatureDecomposition.columns.values) - set(['Sample Name', 'Number of Mutations'])

	print 'attributing n mutations to signature'

	for signatureCol in signatureColNames:
		#binarize the signatures to get the signatures that are present by getting nmut attributable then binarizing it based on our threshold
		signatureDecomposition[signatureCol] = signatureDecomposition.apply(lambda row:
			1 if row['Number of Mutations']*row[signatureCol] > thresholdNMutForSignature else 0, axis=1)


	print 'figuring out attibution stuff for cases'

	#we change the spectrum df here to be either an empty string or signatureName for later processing
	for col in quadNucColumns:
		spectrumDf[col] = spectrumDf.apply(lambda row: row['signature'] if row[col] == 1 else '', axis=1)

	outputDf = summarize_quad_nuc_mut_attribution(signatureDecomposition, spectrumDf)

	print 'writing output to ', outputFilePath
	outputDf.to_csv(outputFilePath, index=False, sep='\t')


if __name__ == '__main__':
    main()









