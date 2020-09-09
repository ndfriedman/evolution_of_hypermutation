#written by Noah Friedman
#THIS IS THE MAIN SCRIPT TO TAKE THE GENE MUT MAFs AND CONVERT THEM INTO 
import os
import pandas as pd
from collections import Counter
import sys

sys.path.append('/juno/work/taylorlab/friedman/myUtils')
import mutationSigUtils

cntr = 0
mainDir = '/juno/work/taylorlab/friedman/myAdjustedDataFiles/allPossibleIMPACTSnps'
listOfDicts = []

for file in os.listdir(mainDir):
	if 'trinuc' in file:
		localDict = {}
		gene = file.split('_')[0]

		try:
			geneDf = pd.read_csv(os.path.join(mainDir, file),
	        	header=0,
	        	usecols=["Ref_Tri", "Hugo_Symbol", "Variant_Classification", "is-a-hotspot", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type", "oncogenic", "quadNuc"], sep='\t')

			if geneDf.shape[0] > 1:
				#geneDf['quadNuc'] = geneDf.apply(lambda row: mutationSigUtils.create_reference_four_nuc(row['Ref_Tri'], row['Reference_Allele'], row['Tumor_Seq_Allele2'], row['Variant_Type']), axis=1)

				localDict = {}

				nonSilentVariantTypes = ["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Translation_Start_Site"]

				nonSilentGeneDf = geneDf[geneDf['Variant_Classification'].isin(nonSilentVariantTypes)]
				truncatingDf = nonSilentGeneDf[nonSilentGeneDf['Variant_Classification'] == 'Nonsense_Mutation']
				hotspotDf = nonSilentGeneDf[nonSilentGeneDf['is-a-hotspot'].notnull()]
				oncogenicDf = nonSilentGeneDf[nonSilentGeneDf['oncogenic'].notnull()]

				nonSilentQuadNucDict = dict(nonSilentGeneDf ['quadNuc'].value_counts())
				truncatingQuadNucDict = dict(truncatingDf['quadNuc'].value_counts())
				oncogenicQuadNucDict = dict(oncogenicDf['quadNuc'].value_counts())
				hotspotQuadNucDict = dict(hotspotDf['quadNuc'].value_counts())

				localDict['Hugo_Symbol'] = gene

				for key,value in nonSilentQuadNucDict.items():
					if key[1] == key[2]:
						print key, file
					localDict[key+'_nonSilent'] = value

				for key,value in truncatingQuadNucDict.items():
					if key[1] == key[2]:
						print key, file
					localDict[key+'_truncating'] = value

				for key,value in oncogenicQuadNucDict.items():
					if key[1] == key[2]:
						print key, file
					localDict[key+'_oncogenic'] = value

				for key,value in hotspotQuadNucDict.items():
					if key[1] == key[2]:
						print key, file
					localDict[key+'_hotspot'] = value

				listOfDicts.append(localDict)
				#TODO TOMORROW: make this work for quad nucs etc
		except:
			#AVOID CORRUPTED FILES
			print gene, 'not found'
			pass


		print cntr, gene
		cntr += 1

df = pd.DataFrame(listOfDicts)
#TEMP
print 'writing data to csv', '/juno/work/taylorlab/friedman/myAdjustedDataFiles/allPossibleIMPACTMutationsSummary.tsv'
df.to_csv('/juno/work/taylorlab/friedman/myAdjustedDataFiles/allPossibleIMPACTMutationsSummary.tsv', index=False, sep='\t')















