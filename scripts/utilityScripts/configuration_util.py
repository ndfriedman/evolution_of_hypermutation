#written by noah friedman
#a utility that returns paths for different types of files we use
#note this script and all others must be run from the base path of the project
import pandas as pd
import os

fileDirPath = '/Users/friedman/Desktop/hypermutationProjectFinal/fileDirectory.txt'

def parse_file_path_dir():
	f = open(fileDirPath, 'r')
	lines = f.readlines()
	headers = lines[0].strip('\n').split('\t')

	listOfDicts = []
	for line in lines:
		if line[0] != '#' and line[0] != '\n': #ignore lines that start with a # and that are empty
			d = {}
			lineSplit = line.strip('\n').split('\t')
			nSplitVals = len(lineSplit)
			for i in range(len(headers)):
				if i + 1 > nSplitVals:
					d[headers[i]] = None
				else:
					d[headers[i]] = lineSplit[i]
			listOfDicts.append(d)
	df = pd.DataFrame(listOfDicts)
	#SET UP ALL PATHS APPROPRIATELY
	basePath = get_project_base_path()
	df['path'] = df['path'].apply(lambda x: os.path.join(basePath, str(x)))
	return df

#gets the project's basePath by reading a file 
def get_project_base_path():
	return '/Users/friedman/Desktop/hypermutationProjectFinal/'

def get_all_files_path_dict():
	configurationDf = parse_file_path_dir()
	return dict(zip(configurationDf['fileDescriptor'], configurationDf['path']))

if __name__ == '__main__':
	# run!
	print get_all_files_path_dict()
