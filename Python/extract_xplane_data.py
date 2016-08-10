import os
import glob
import fnmatch
import os

txtfiles = []
for root, dirnames, filenames in os.walk('.\\'):
    for filename in fnmatch.filter(filenames, '*_data*.txt'):
        txtfiles.append(os.path.join(root, filename))
        
        
print txtfiles

for txtfile in txtfiles:
    datafile = open(txtfile,'r')
    data = datafile.read()

    datafile.close()