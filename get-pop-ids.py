'''
script to loop through seq folder and create a pop.txt file that contains the population ID for each fly.
The nature of the for loop is not alphbetical so you will need to sort to a new file in terminal
'''
#Put this script in directory with the seq files

import os
import re


#set directories
seq_directory = #path/to/seq/directory/

vcf_directory = #path/to/new/directory/ <--- make sure you end with "/"


with open('{path}'.format(path = vcf_directory) + 'unsorted_pop.txt', 'w') as pop:
        for file in os.listdir(seq_directory):
                if file.endswith('.seq'):
                        pattern = '-|0|1|2|3|4|5|6|7|8|9'
                        string = file
                        result = re.split(pattern, string)[0]
                        pop.write(result + '\n')
                else:
                        continue       
