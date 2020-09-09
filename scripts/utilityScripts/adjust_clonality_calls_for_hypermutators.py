#written by noah friedman
#a script that adds clonality calls to maf

import sys
import argparse
import os
import pandas as pd
import numpy as np
import sys
from collections import Counter

sys.path.append('/Users/friedman/Desktop/hypermutationProjectFinal/scripts/utilityScripts')
import analysis_utils 
import mutationSigUtils 
import maf_analysis_utils
import clonality_analysis_util
import get_gene_and_cohort_list_utils
import configuration_util
filePathDict = configuration_util.get_all_files_path_dict()



mafToAnnotatePath = sys.argv[1]
writePath = sys.argv[2]
mode = sys.argv[3]

mafToAnnotate = pd.read_table(mafToAnnotatePath)
hypermutationStatusDir = filePathDict['HYPERMUTATION_STATUS_IDS']

hypermutantIds = []
if mode == 'IMPACT':
	hypermutantIds = get_gene_and_cohort_list_utils.get_all_hypermutant_ids(hypermutantIdDir=hypermutationStatusDir)

facetsWhitelist = clonality_analysis_util.get_facets_whitelist()
facetsBlacklist = clonality_analysis_util.get_facets_blacklist()

hypermutatedMaf = mafToAnnotate[mafToAnnotate['Tumor_Sample_Barcode'].isin(hypermutantIds)]
notHypermutatedMaf = mafToAnnotate[~mafToAnnotate['Tumor_Sample_Barcode'].isin(hypermutantIds)]
hypermutatedMaf = clonality_analysis_util.assign_clonality_information_for_hypermutated_cases(
	hypermutatedMaf, facetsWhitelist, facetsBlacklist)

print 'combiningMafs'
combinedMaf = pd.concat([hypermutatedMaf, notHypermutatedMaf])

print 'writing clonal annotated file to ', writePath
combinedMaf.to_csv(writePath)

