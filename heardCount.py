import os
import glob
import fnmatch
import os

txtfiles = []
for root, dirnames, filenames in os.walk('G:\Sensorimotor_RT_8-3-2015\MEG_Subject23'):
    for filename in fnmatch.filter(filenames, '*replay*.txt'):
        txtfiles.append(os.path.join(root, filename))
##print txtfiles

heardCount = []
missedCount = []
for txtfile in txtfiles:
    datafile = open(txtfile,'r')
    data = datafile.read()
    numHeard = data.count('heard')
    numMissed = data.count('missed')
    
    heardCount.append(numHeard)
    missedCount.append(numMissed)
    datafile.close()
sumHeard = sum(heardCount)
sumCount =  sum(missedCount)

print float(sumHeard)/(sumCount+sumHeard)