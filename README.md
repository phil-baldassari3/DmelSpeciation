# Kulathinal Lab DrosSpec Scripts

Here is Phil's repository for scripts for the DrosSpec Project

This repo contains:

check_SNPs.py --->
This python script is used to loop through seq files in a directory and contains a function that will print out the a given locus from each file based on the indices passed into the function.  Note: since python starts counting characters from zero, you will need to think about subtracting one to get the exact locus you want to print.

check-looper.py --->
This python script is the same script that is located in each of the directories on my machine that contains the seq files that I am working with.  It loops through the directory and has multiple functions.  This script contains the same function in the check_SNPs.py script.  It also contains a function to print the length of each file and to tell which file if any has a different length.

seq2fasta.py --->
This python script converts the seq files to individual fasta files by looping through the seq directory, creating new files in the fasta directory (you set the directories), appending ">" plus the title of each file to the beginning line, and adding a line break to the end of each file.  Then it changes the file extensions of the new files to ".fas"
