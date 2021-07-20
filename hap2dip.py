#using modified search+replace.py to convert haploid vcf into homozygous diploid by replacing 0 or 1 with 0/0 or 1/1
#yes, I know this is probably an easier way to do this, but I need to get this done
#this script is only going to convert my specific files for the Dros_Speciation project

###Converting zim-cos_ChrX_0.05.vcf
###########################
# Read in the file
with open('zim-cos_ChrX_0.05.vcf', 'r') as Xfile:
  Xfiledata = Xfile.read()
# Replace the target string
Xfiledata = Xfiledata.replace('	0	0', '	0/0	0/0')
Xfiledata = Xfiledata.replace('	1	1', '	1/1	1/1')
Xfiledata = Xfiledata.replace('	1	0', '	1/1	0/0')
Xfiledata = Xfiledata.replace('	0	1', '	0/0	1/1')
# Write the file out again
with open('dipzim-cos_ChrX_0.05.vcf', 'w') as Xfile1:
  Xfile1.write(Xfiledata)

print("Ploidy of ChrX changed.  Please wait...")
############################

# Read in the file
with open('dipzim-cos_ChrX_0.05.vcf', 'r') as Xfile1:
  Xfiledata1 = Xfile1.read()
# Replace the target string
Xfiledata1 = Xfiledata1.replace('/0/', '/')
Xfiledata1 = Xfiledata1.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_ChrX_0.05.vcf', 'w') as Xfile2:
  Xfile2.write(Xfiledata1)

print("Corrected erroneous triploids in ChrX.  Please wait...")
###########################

# Read in the file
with open('dipzim-cos_ChrX_0.05.vcf', 'r') as Xfile3:
  Xfiledata2 = Xfile3.read()
# Replace the target string
Xfiledata2 = Xfiledata2.replace('0' + '\n' + 'X', '0/0' + '\n' + 'X')
Xfiledata2 = Xfiledata2.replace('1' + '\n' + 'X', '1/1' + '\n' + 'X')
# Write the file out again
with open('dipzim-cos_ChrX_0.05.vcf', 'w') as Xfile4:
  Xfile4.write(Xfiledata2)

print("Corrected remaining haploids in ChrX.  Please wait...")
###########################
# Read in the file
with open('dipzim-cos_ChrX_0.05.vcf', 'r') as Xfile5:
  Xfiledata3 = Xfile5.read()
# Replace the target string
Xfiledata3 = Xfiledata3.replace('/0/', '/')
Xfiledata3 = Xfiledata3.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_ChrX_0.05.vcf', 'w') as Xfile6:
  Xfile6.write(Xfiledata3)

print("Corrected erroneous triploids in ChrX. Done!")
print("Please wait...")

############################################################################################################
###Converting zim-cos_Chr2L_0.05.vcf
###########################
# Read in the file
with open('zim-cos_Chr2L_0.05.vcf', 'r') as 2Lfile:
  2Lfiledata = 2Lfile.read()
# Replace the target string
2Lfiledata = 2Lfiledata.replace('	0	0', '	0/0	0/0')
2Lfiledata = 2Lfiledata.replace('	1	1', '	1/1	1/1')
2Lfiledata = 2Lfiledata.replace('	1	0', '	1/1	0/0')
2Lfiledata = 2Lfiledata.replace('	0	1', '	0/0	1/1')
# Write the file out again
with open('dipzim-cos_Chr2L_0.05.vcf', 'w') as 2Lfile1:
  2Lfile1.write(2Lfiledata)

print("Ploidy of Chr2L changed.  Please wait...")
############################

# Read in the file
with open('dipzim-cos_Chr2L_0.05.vcf', 'r') as 2Lfile1:
  2Lfiledata1 = 2Lfile1.read()
# Replace the target string
2Lfiledata1 = 2Lfiledata1.replace('/0/', '/')
2Lfiledata1 = 2Lfiledata1.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr2L_0.05.vcf', 'w') as 2Lfile2:
  2Lfile2.write(2Lfiledata1)

print("Corrected erroneous triploids in Chr2L.  Please wait...")
###########################

# Read in the file
with open('dipzim-cos_Chr2L_0.05.vcf', 'r') as 2Lfile3:
  2Lfiledata2 = 2Lfile3.read()
# Replace the target string
2Lfiledata2 = 2Lfiledata2.replace('0' + '\n' + '2L', '0/0' + '\n' + '2L')
2Lfiledata2 = 2Lfiledata2.replace('1' + '\n' + '2L', '1/1' + '\n' + '2L')
# Write the file out again
with open('dipzim-cos_Chr2L_0.05.vcf', 'w') as 2Lfile4:
  2Lfile4.write(2Lfiledata2)

print("Corrected remaining haploids in Chr2L.  Please wait...")
###########################
# Read in the file
with open('dipzim-cos_Chr2L_0.05.vcf', 'r') as 2Lfile5:
  2Lfiledata3 = 2Lfile5.read()
# Replace the target string
2Lfiledata3 = 2Lfiledata3.replace('/0/', '/')
2Lfiledata3 = 2Lfiledata3.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr2L_0.05.vcf', 'w') as 2Lfile6:
  2Lfile6.write(2Lfiledata3)

print("Corrected erroneous triploids in Chr2L. Done!")
print("Please wait...")

############################################################################################################
###Converting zim-cos_Chr2R_0.05.vcf
###########################
# Read in the file
with open('zim-cos_Chr2R_0.05.vcf', 'r') as 2Rfile:
  2Rfiledata = 2Rfile.read()
# Replace the target string
2Rfiledata = 2Rfiledata.replace('	0	0', '	0/0	0/0')
2Rfiledata = 2Rfiledata.replace('	1	1', '	1/1	1/1')
2Rfiledata = 2Rfiledata.replace('	1	0', '	1/1	0/0')
2Rfiledata = 2Rfiledata.replace('	0	1', '	0/0	1/1')
# Write the file out again
with open('dipzim-cos_Chr2R_0.05.vcf', 'w') as 2Rfile1:
  2Rfile1.write(2Rfiledata)

print("Ploidy of Chr2R changed.  Please wait...")
############################

# Read in the file
with open('dipzim-cos_Chr2R_0.05.vcf', 'r') as 2Rfile1:
  2Rfiledata1 = 2Rfile1.read()
# Replace the target string
2Rfiledata1 = 2Rfiledata1.replace('/0/', '/')
2Rfiledata1 = 2Rfiledata1.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr2R_0.05.vcf', 'w') as 2Rfile2:
  2Rfile2.write(2Rfiledata1)

print("Corrected erroneous triploids in Chr2R.  Please wait...")
###########################

# Read in the file
with open('dipzim-cos_Chr2R_0.05.vcf', 'r') as 2Rfile3:
  2Rfiledata2 = 2Rfile3.read()
# Replace the target string
2Rfiledata2 = 2Rfiledata2.replace('0' + '\n' + '2R', '0/0' + '\n' + '2R')
2Rfiledata2 = 2Rfiledata2.replace('1' + '\n' + '2R', '1/1' + '\n' + '2R')
# Write the file out again
with open('dipzim-cos_Chr2R_0.05.vcf', 'w') as 2Rfile4:
  2Rfile4.write(2Rfiledata2)

print("Corrected remaining haploids in Chr2R.  Please wait...")
###########################
# Read in the file
with open('dipzim-cos_Chr2R_0.05.vcf', 'r') as 2Rfile5:
  2Rfiledata3 = 2Rfile5.read()
# Replace the target string
2Rfiledata3 = 2Rfiledata3.replace('/0/', '/')
2Rfiledata3 = 2Rfiledata3.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr2R_0.05.vcf', 'w') as 2Rfile6:
  2Rfile6.write(2Rfiledata3)

print("Corrected erroneous triploids in Chr2R. Done!")
print("Please wait...")

############################################################################################################
###Converting zim-cos_Chr3L_0.05.vcf
###########################
# Read in the file
with open('zim-cos_Chr3L_0.05.vcf', 'r') as 3Lfile:
  3Lfiledata = 3Lfile.read()
# Replace the target string
3Lfiledata = 3Lfiledata.replace('	0	0', '	0/0	0/0')
3Lfiledata = 3Lfiledata.replace('	1	1', '	1/1	1/1')
3Lfiledata = 3Lfiledata.replace('	1	0', '	1/1	0/0')
3Lfiledata = 3Lfiledata.replace('	0	1', '	0/0	1/1')
# Write the file out again
with open('dipzim-cos_Chr3L_0.05.vcf', 'w') as 3Lfile1:
  3Lfile1.write(3Lfiledata)

print("Ploidy of Chr3L changed.  Please wait...")
############################

# Read in the file
with open('dipzim-cos_Chr3L_0.05.vcf', 'r') as 3Lfile1:
  3Lfiledata1 = 3Lfile1.read()
# Replace the target string
3Lfiledata1 = 3Lfiledata1.replace('/0/', '/')
3Lfiledata1 = 3Lfiledata1.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr3L_0.05.vcf', 'w') as 3Lfile2:
  3Lfile2.write(3Lfiledata1)

print("Corrected erroneous triploids in Chr3L.  Please wait...")
###########################

# Read in the file
with open('dipzim-cos_Chr3L_0.05.vcf', 'r') as 3Lfile3:
  3Lfiledata2 = 3Lfile3.read()
# Replace the target string
3Lfiledata2 = 3Lfiledata2.replace('0' + '\n' + '3L', '0/0' + '\n' + '3L')
3Lfiledata2 = 3Lfiledata2.replace('1' + '\n' + '3L', '1/1' + '\n' + '3L')
# Write the file out again
with open('dipzim-cos_Chr3L_0.05.vcf', 'w') as 3Lfile4:
  3Lfile4.write(3Lfiledata2)

print("Corrected remaining haploids in Chr3L.  Please wait...")
###########################
# Read in the file
with open('dipzim-cos_Chr3L_0.05.vcf', 'r') as 3Lfile5:
  3Lfiledata3 = 3Lfile5.read()
# Replace the target string
3Lfiledata3 = 3Lfiledata3.replace('/0/', '/')
3Lfiledata3 = 3Lfiledata3.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr3L_0.05.vcf', 'w') as 3Lfile6:
  3Lfile6.write(3Lfiledata3)

print("Corrected erroneous triploids in Chr3L. Done!")
print("Please wait...")

############################################################################################################
###Converting zim-cos_Chr3R_0.05.vcf
###########################
# Read in the file
with open('zim-cos_Chr3R_0.05.vcf', 'r') as 3Rfile:
  3Rfiledata = 3Rfile.read()
# Replace the target string
3Rfiledata = 3Rfiledata.replace('	0	0', '	0/0	0/0')
3Rfiledata = 3Rfiledata.replace('	1	1', '	1/1	1/1')
3Rfiledata = 3Rfiledata.replace('	1	0', '	1/1	0/0')
3Rfiledata = 3Rfiledata.replace('	0	1', '	0/0	1/1')
# Write the file out again
with open('dipzim-cos_Chr3R_0.05.vcf', 'w') as 3Rfile1:
  3Rfile1.write(3Rfiledata)

print("Ploidy of Chr3R changed.  Please wait...")
############################

# Read in the file
with open('dipzim-cos_Chr3R_0.05.vcf', 'r') as 3Rfile1:
  3Rfiledata1 = 3Rfile1.read()
# Replace the target string
3Rfiledata1 = 3Rfiledata1.replace('/0/', '/')
3Rfiledata1 = 3Rfiledata1.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr3R_0.05.vcf', 'w') as 3Rfile2:
  3Rfile2.write(3Rfiledata1)

print("Corrected erroneous triploids in Chr3R.  Please wait...")
###########################

# Read in the file
with open('dipzim-cos_Chr3R_0.05.vcf', 'r') as 3Rfile3:
  3Rfiledata2 = 3Rfile3.read()
# Replace the target string
3Rfiledata2 = 3Rfiledata2.replace('0' + '\n' + '3R', '0/0' + '\n' + '3R')
3Rfiledata2 = 3Rfiledata2.replace('1' + '\n' + '3R', '1/1' + '\n' + '3R')
# Write the file out again
with open('dipzim-cos_Chr3R_0.05.vcf', 'w') as 3Rfile4:
  3Rfile4.write(3Rfiledata2)

print("Corrected remaining haploids in Chr3R.  Please wait...")
###########################
# Read in the file
with open('dipzim-cos_Chr3R_0.05.vcf', 'r') as 3Rfile5:
  3Rfiledata3 = 3Rfile5.read()
# Replace the target string
3Rfiledata3 = 3Rfiledata3.replace('/0/', '/')
3Rfiledata3 = 3Rfiledata3.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr3R_0.05.vcf', 'w') as 3Rfile6:
  3Rfile6.write(3Rfiledata3)

print("Corrected erroneous triploids in Chr3R.")
print(Done!")
