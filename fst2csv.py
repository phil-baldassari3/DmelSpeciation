#this script is for quickly converting the vcftools fst output files to csv files so that they can be more easily used in R
#put this in the directory that it is needed in and set the directory path

import os,sys

directory = #'path/to/directory/'

for file in os.listdir(directory):
    if file.endswith('.fst'):
    	#opening and reading file
    	fst_file = open(file)
		read_fst = fst_file.read()
		
		#replacing spaces with commas
		fst_csv = read_fst.replace(' ', ',')
		
		#create new file and write in the replaced data
        with open('{path}'.format(path = directory) + '{file}'.format(file = file) + '.csv', 'w') as csv:
            csv.write(fst_csv)
    else:
        continue       
		