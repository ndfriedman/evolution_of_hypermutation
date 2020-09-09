#written by Noah Friedman
#this works on the whole hypermutated cohort and clearly articulates expected data for a variety of contexts


import sys
import argparse
import os
import pandas as pd
import numpy as np
import re

sys.path.append('/Users/friedman/Desktop/hypermutationProjectFinal/scripts/utilityScripts')
import analysis_utils 
import mutationSigUtils
import mutation_modeling_util
import configuration_util

filePathDict = configuration_util.get_all_files_path_dict()

def calculate_oncogenic_mut_susceptibility_of_genes_by_signature(oncogenicSDict, suffix='_hotspot_rate'):
    listOfDicts = []
    sigNames = ['Signature.' + str(i) for i in range(1,31)]
    for i in range(1,31):
        curSig = 'Signature.' + str(i)
        d = {}
        for s in sigNames:
            d[s] = 0
        d[curSig] = 1
        #PRETEND we got a case with 100% signature i on the decomposition
        quadNucFractions = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(d, spectraPath = pathPrefix + '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt')
        #v = mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions(quadNucFractions, oncogenicSDict, 'IMPACT_468')
        #ALERT NOAH I CHANGED THIS HERE
        v = mutation_modeling_util.get_expected_oncogenic_val_given_quadnuc_fractions_v2(quadNucFractions, oncogenicSDict, suffix)

        listOfDicts.append({'Signature_Name': curSig, 'ExpectedFracOfMutsOncogenic': v})
    return pd.DataFrame(listOfDicts)

def quantify_quadnuc_oncogenic_susceptibility_per_mutation(simMafData, relatedGenes=None):
    allBases = ['A', 'C', 'G', 'T']
    changes = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG'] #format: 'CA' means a change from C>A
    allQuadNucs = [firstBase + change + lastBase for firstBase in allBases for change in changes for lastBase in allBases] #enumerate all 96 quadnucs for signatures
   
    d = {}
    
    for quadNuc in allQuadNucs:
        nPossibleMuts = sum(simMafData[quadNuc]) - sum(simMafData[quadNuc + '_silent'])
        nPossibleHotspotMuts = sum(simMafData[quadNuc + '_oncogenic'])
        d[quadNuc + '_oncogenic_rate'] = (1.0*nPossibleHotspotMuts)/(1.0*nPossibleMuts)
        
    return d

#fills in impossible columns with zeroes 
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

im3GenesOnly = True
im3Genes = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])

print 'loading signatures info'
impactSigsPath = '/juno/work/taylorlab/friedman/myAdjustedDataFiles/impactSignatureCalls_Nov20_2019_not_merged.tsv'
print 'loading spectra d'
spectraPath = '/ifs/work/taylorlab/friedman/noahFirstProject/signature_sig_copy/mutation-signatures/Stratton_signatures30.txt'
print 'loading expected mut rate df'
countsDf = pd.read_table('/juno/work/taylorlab/friedman/myAdjustedDataFiles/allPossibleIMPACTMutationsSummary.tsv')

countsDf = add_zero_cols_to_counts_df(countsDf)

allHypermutatedCases = analysis_utils.enumerate_all_hypermutated_cases(pathPrefix='')
#TODO add a mode that supports working with more than just hypermutated cases

impactSigs = pd.read_table(filePathDict['IMPACT_SIGNATURE_DECOMPOSITIONS'])
impactSigs['pid'] = impactSigs['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
impactSigs = impactSigs.rename(columns= dict([('mean_' + str(i), 'Signature.' + str(i)) for i in range(1,31)]))

spectraD = mutationSigUtils.convert_spectrum_file_to_dict_of_dicts(spectrumFile=filePathDict['SIGNATURE_SPECTRUM'])

if im3GenesOnly:
    countsDf = countsDf[countsDf['Hugo_Symbol'].isin(im3Genes)]
genes = set(countsDf['Hugo_Symbol'])

expectedMutRateDf = mutation_modeling_util.convert_counts_information_to_fraction(countsDf)

print 'calculating chances of mutation'

#save runtime by presubsetting the dataframe by genes
gdfDict = {}
for gene in genes:
	gdf = expectedMutRateDf[expectedMutRateDf['gene'] == gene]
	gdfDict[gene] = gdf

listOfDfs = []
cntr = 0.0


#FOR NOW JUST DO THIS FOR ENDOMETRIAL COLON AND GLIOMA CAUSE I DONT HAVE TIME
"""hyperEndometrial = analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir='/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType='Endometrial Cancer', hypermutantStatus = 'Hypermutated')
hyperColorectal = analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir='/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType='Colorectal Cancer', hypermutantStatus = 'Hypermutated')
hyperGlioma = analysis_utils.get_ids_by_hypermutant_status(hypermutantIdDir='/ifs/work/taylorlab/friedman/hypermutationAnalysisProj/projectDataAndConfigFiles/hypermutationStatusIds', cancerType='Glioma', hypermutantStatus = 'Hypermutated')

allHypermutatedCases = hyperEndometrial | hyperColorectal | hyperGlioma
"""

allCases = set(impactSigs['Tumor_Sample_Barcode'])
print 'beginning iteration over cases, note each iteration takes between 5-10 and seconds per case'
for case in list(allCases):
    if cntr %10 == 0:
        print cntr/len(allHypermutatedCases)*100, ' percent complete'

    cntr += 1
    caseSigs = impactSigs[impactSigs['Tumor_Sample_Barcode'] == case]
    quadNucFracs = mutation_modeling_util.get_quadnuc_fracs_given_decomposition(caseSigs.iloc[0], spectraD)

    expectedCopy = expectedMutRateDf.copy()
    expectedCopy['hotspotChance'] =  expectedCopy.apply(lambda row: quadNucFracs[row['quadNuc']]* row['hotspotChance'], axis=1)
    expectedCopy['oncogenicChance'] =  expectedCopy.apply(lambda row: quadNucFracs[row['quadNuc']]* row['oncogenicChance'], axis=1)
    expectedCopy['truncatingChance'] =  expectedCopy.apply(lambda row: quadNucFracs[row['quadNuc']]* row['truncatingChance'], axis=1)

    #THE Groupby operation gives us the information we want
    hotspotInfo = expectedCopy.groupby("gene").hotspotChance.sum().reset_index()
    oncogenicInfo = expectedCopy.groupby("gene").oncogenicChance.sum().reset_index()
    truncatingInfo = expectedCopy.groupby("gene").truncatingChance.sum().reset_index()

    #merge them all
    dfMerged = reduce(lambda left,right: pd.merge(left,right,on='gene'), [hotspotInfo, oncogenicInfo, truncatingInfo])
    dfMerged['case'] = case
    listOfDfs.append(dfMerged)

df = pd.concat(listOfDfs)
print 'writing data'


df.to_csv('/juno/work/taylorlab/friedman/hypermutationAnalysisProj/mutSimulation/expectedMutationTables/allCasesExpectedGeneMutInfo.tsv', index=False, sep='\t')

#oncogenicSusceptibilityDict = mutation_modeling_util.calculate_quadnuc_based_oncogenic_susceptibility_dict(simDfSummary)
#allMutsPossible = analysis_utils.load_in_df_with_progress(filePath = pathPrefix + '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/allMutsContextSummary.tsv', nLinesFile=4243342) 









