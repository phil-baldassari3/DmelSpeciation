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
###Converting zim-cos_Chr_2L_0.05.vcf
###########################
# Read in the file
with open('zim-cos_Chr_2L_0.05.vcf', 'r') as _2Lfile:
  _2Lfiledata = _2Lfile.read()
# Replace the target string
_2Lfiledata = _2Lfiledata.replace('	0	0', '	0/0	0/0')
_2Lfiledata = _2Lfiledata.replace('	1	1', '	1/1	1/1')
_2Lfiledata = _2Lfiledata.replace('	1	0', '	1/1	0/0')
_2Lfiledata = _2Lfiledata.replace('	0	1', '	0/0	1/1')
# Write the file out again
with open('dipzim-cos_Chr_2L_0.05.vcf', 'w') as _2Lfile1:
  _2Lfile1.write(_2Lfiledata)

print("Ploidy of Chr2L changed.  Please wait...")
############################

# Read in the file
with open('dipzim-cos_Chr_2L_0.05.vcf', 'r') as _2Lfile1:
  _2Lfiledata1 = _2Lfile1.read()
# Replace the target string
_2Lfiledata1 = _2Lfiledata1.replace('/0/', '/')
_2Lfiledata1 = _2Lfiledata1.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr_2L_0.05.vcf', 'w') as _2Lfile2:
  _2Lfile2.write(_2Lfiledata1)

print("Corrected erroneous triploids in Chr2L.  Please wait...")
###########################

# Read in the file
with open('dipzim-cos_Chr_2L_0.05.vcf', 'r') as _2Lfile3:
  _2Lfiledata2 = _2Lfile3.read()
# Replace the target string
_2Lfiledata2 = _2Lfiledata2.replace('0' + '\n' + '2L', '0/0' + '\n' + '2L')
_2Lfiledata2 = _2Lfiledata2.replace('1' + '\n' + '2L', '1/1' + '\n' + '2L')
# Write the file out again
with open('dipzim-cos_Chr_2L_0.05.vcf', 'w') as _2Lfile4:
  _2Lfile4.write(_2Lfiledata2)

print("Corrected remaining haploids in Chr2L.  Please wait...")
###########################
# Read in the file
with open('dipzim-cos_Chr_2L_0.05.vcf', 'r') as _2Lfile5:
  _2Lfiledata3 = _2Lfile5.read()
# Replace the target string
_2Lfiledata3 = _2Lfiledata3.replace('/0/', '/')
_2Lfiledata3 = _2Lfiledata3.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr_2L_0.05.vcf', 'w') as _2Lfile6:
  _2Lfile6.write(_2Lfiledata3)

print("Corrected erroneous triploids in Chr2L. Done!")
print("Please wait...")

############################################################################################################
###Converting zim-cos_Chr_2R_0.05.vcf
###########################
# Read in the file
with open('zim-cos_Chr_2R_0.05.vcf', 'r') as _2Rfile:
  _2Rfiledata = _2Rfile.read()
# Replace the target string
_2Rfiledata = _2Rfiledata.replace('	0	0', '	0/0	0/0')
_2Rfiledata = _2Rfiledata.replace('	1	1', '	1/1	1/1')
_2Rfiledata = _2Rfiledata.replace('	1	0', '	1/1	0/0')
_2Rfiledata = _2Rfiledata.replace('	0	1', '	0/0	1/1')
# Write the file out again
with open('dipzim-cos_Chr_2R_0.05.vcf', 'w') as _2Rfile1:
  _2Rfile1.write(_2Rfiledata)

print("Ploidy of Chr2R changed.  Please wait...")
############################

# Read in the file
with open('dipzim-cos_Chr_2R_0.05.vcf', 'r') as _2Rfile1:
  _2Rfiledata1 = _2Rfile1.read()
# Replace the target string
_2Rfiledata1 = _2Rfiledata1.replace('/0/', '/')
_2Rfiledata1 = _2Rfiledata1.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr_2R_0.05.vcf', 'w') as _2Rfile2:
  _2Rfile2.write(_2Rfiledata1)

print("Corrected erroneous triploids in Chr2R.  Please wait...")
###########################

# Read in the file
with open('dipzim-cos_Chr_2R_0.05.vcf', 'r') as _2Rfile3:
  _2Rfiledata2 = _2Rfile3.read()
# Replace the target string
_2Rfiledata2 = _2Rfiledata2.replace('0' + '\n' + '2R', '0/0' + '\n' + '2R')
_2Rfiledata2 = _2Rfiledata2.replace('1' + '\n' + '2R', '1/1' + '\n' + '2R')
# Write the file out again
with open('dipzim-cos_Chr_2R_0.05.vcf', 'w') as _2Rfile4:
  _2Rfile4.write(_2Rfiledata2)

print("Corrected remaining haploids in Chr2R.  Please wait...")
###########################
# Read in the file
with open('dipzim-cos_Chr_2R_0.05.vcf', 'r') as _2Rfile5:
  _2Rfiledata3 = _2Rfile5.read()
# Replace the target string
_2Rfiledata3 = _2Rfiledata3.replace('/0/', '/')
_2Rfiledata3 = _2Rfiledata3.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr_2R_0.05.vcf', 'w') as _2Rfile6:
  _2Rfile6.write(_2Rfiledata3)

print("Corrected erroneous triploids in Chr2R. Done!")
print("Please wait...")

############################################################################################################
###Converting zim-cos_Chr_3L_0.05.vcf
###########################
# Read in the file
with open('zim-cos_Chr_3L_0.05.vcf', 'r') as _3Lfile:
  _3Lfiledata = _3Lfile.read()
# Replace the target string
_3Lfiledata = _3Lfiledata.replace('	0	0', '	0/0	0/0')
_3Lfiledata = _3Lfiledata.replace('	1	1', '	1/1	1/1')
_3Lfiledata = _3Lfiledata.replace('	1	0', '	1/1	0/0')
_3Lfiledata = _3Lfiledata.replace('	0	1', '	0/0	1/1')
# Write the file out again
with open('dipzim-cos_Chr_3L_0.05.vcf', 'w') as _3Lfile1:
  _3Lfile1.write(_3Lfiledata)

print("Ploidy of Chr3L changed.  Please wait...")
############################

# Read in the file
with open('dipzim-cos_Chr_3L_0.05.vcf', 'r') as _3Lfile1:
  _3Lfiledata1 = _3Lfile1.read()
# Replace the target string
_3Lfiledata1 = _3Lfiledata1.replace('/0/', '/')
_3Lfiledata1 = _3Lfiledata1.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr_3L_0.05.vcf', 'w') as _3Lfile2:
  _3Lfile2.write(_3Lfiledata1)

print("Corrected erroneous triploids in Chr3L.  Please wait...")
###########################

# Read in the file
with open('dipzim-cos_Chr_3L_0.05.vcf', 'r') as _3Lfile3:
  _3Lfiledata2 = _3Lfile3.read()
# Replace the target string
_3Lfiledata2 = _3Lfiledata2.replace('0' + '\n' + '3L', '0/0' + '\n' + '3L')
_3Lfiledata2 = _3Lfiledata2.replace('1' + '\n' + '3L', '1/1' + '\n' + '3L')
# Write the file out again
with open('dipzim-cos_Chr_3L_0.05.vcf', 'w') as _3Lfile4:
  _3Lfile4.write(_3Lfiledata2)

print("Corrected remaining haploids in Chr3L.  Please wait...")
###########################
# Read in the file
with open('dipzim-cos_Chr_3L_0.05.vcf', 'r') as _3Lfile5:
  _3Lfiledata3 = _3Lfile5.read()
# Replace the target string
_3Lfiledata3 = _3Lfiledata3.replace('/0/', '/')
_3Lfiledata3 = _3Lfiledata3.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr_3L_0.05.vcf', 'w') as _3Lfile6:
  _3Lfile6.write(_3Lfiledata3)

print("Corrected erroneous triploids in Chr3L. Done!")
print("Please wait...")

############################################################################################################
###Converting zim-cos_Chr_3R_0.05.vcf
###########################
# Read in the file
with open('zim-cos_Chr_3R_0.05.vcf', 'r') as _3Rfile:
  _3Rfiledata = _3Rfile.read()
# Replace the target string
_3Rfiledata = _3Rfiledata.replace('	0	0', '	0/0	0/0')
_3Rfiledata = _3Rfiledata.replace('	1	1', '	1/1	1/1')
_3Rfiledata = _3Rfiledata.replace('	1	0', '	1/1	0/0')
_3Rfiledata = _3Rfiledata.replace('	0	1', '	0/0	1/1')
# Write the file out again
with open('dipzim-cos_Chr_3R_0.05.vcf', 'w') as _3Rfile1:
  _3Rfile1.write(_3Rfiledata)

print("Ploidy of Chr3R changed.  Please wait...")
############################

# Read in the file
with open('dipzim-cos_Chr_3R_0.05.vcf', 'r') as _3Rfile1:
  _3Rfiledata1 = _3Rfile1.read()
# Replace the target string
_3Rfiledata1 = _3Rfiledata1.replace('/0/', '/')
_3Rfiledata1 = _3Rfiledata1.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr_3R_0.05.vcf', 'w') as _3Rfile2:
  _3Rfile2.write(_3Rfiledata1)

print("Corrected erroneous triploids in Chr3R.  Please wait...")
###########################

# Read in the file
with open('dipzim-cos_Chr_3R_0.05.vcf', 'r') as _3Rfile3:
  _3Rfiledata2 = _3Rfile3.read()
# Replace the target string
_3Rfiledata2 = _3Rfiledata2.replace('0' + '\n' + '3R', '0/0' + '\n' + '3R')
_3Rfiledata2 = _3Rfiledata2.replace('1' + '\n' + '3R', '1/1' + '\n' + '3R')
# Write the file out again
with open('dipzim-cos_Chr_3R_0.05.vcf', 'w') as _3Rfile4:
  _3Rfile4.write(_3Rfiledata2)

print("Corrected remaining haploids in Chr3R.  Please wait...")
###########################
# Read in the file
with open('dipzim-cos_Chr_3R_0.05.vcf', 'r') as _3Rfile5:
  _3Rfiledata3 = _3Rfile5.read()
# Replace the target string
_3Rfiledata3 = _3Rfiledata3.replace('/0/', '/')
_3Rfiledata3 = _3Rfiledata3.replace('/1/', '/')
# Write the file out again
with open('dipzim-cos_Chr_3R_0.05.vcf', 'w') as _3Rfile6:
  _3Rfile6.write(_3Rfiledata3)

print("Corrected erroneous triploids in Chr3R.")
print("Done!")
