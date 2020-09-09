#written by noah friedman
#a script that is intended to add annotations to mafs as specified
import sys
import os
import pandas as pd

#usage python mafToBeAnnotatedPath sourceMafForAnnotationPath all column names we want to include from source maf space separated
mafToBeAnnotatedPath = sys.argv[1]
sourceMafForAnnotationPath = sys.argv[2]
writePath = sys.argv[3]
colsToMap = sys.argv[4:]

mafToBeAnnotated = pd.read_table(mafToBeAnnotatedPath)
sourceMafForAnnotation = pd.read_table(sourceMafForAnnotationPath)

print 'adding varUuid annotation to maf to be annotated'
mafToBeAnnotated['varUuid'] = mafToBeAnnotated.apply(lambda row:
	str(row['Tumor_Sample_Barcode']) + '_' + str(row['Start_Position']) + '_' + str(row['Tumor_Seq_Allele2']), axis=1)
print 'adding varUuid annotation to sourceMafForAnnotation'
sourceMafForAnnotation['varUuid'] = sourceMafForAnnotation.apply(lambda row:
	str(row['Tumor_Sample_Barcode']) + '_' + str(row['Start_Position']) + '_' + str(row['Tumor_Seq_Allele2']), axis=1)

colDicts = {}
for col in colsToMap:
	colDicts[col] = dict(zip(sourceMafForAnnotation['varUuid'], col))

print 'adding nnotations to maf'
for key, d in colDicts.items():
	mafToBeAnnotated[key] = mafToBeAnnotated['varUuid'].apply(lambda x:
		d[x] if x in d else None)

print 'writing maf to', writePath
mafToBeAnnotated.to_csv(writePath)
