#!/usr/bin/python
import argparse
import numpy as np
import sys

"""
Joint two Xhi/EPS output file and order all values
Author:  Claudio Attaccalite
"""

def exit_error(stringa):
    print(stringa)
    sys.exit(1)
#
# parse command line
#
parser = argparse.ArgumentParser(prog='join_xhi',description='Join two Xhi/EPS output and order all values',\
         epilog="Copyright Claudio Attaccalite 2015")
parser.add_argument('--files' , nargs='*', default=None)
parser.add_argument('-sortcol', type=int , default=0)
args = parser.parse_args()

print("\n Joint Xhi/EPS files \n\n")


if (args.files == None): exit_error('type "join_xhi.py --help" for help ',)

print(" Files to join: %s  \n" % str(args.files)) 

# get the number of column from the first file

data=np.genfromtxt(args.files[0],comments="#")
ncol=data.shape[1]
alldata=np.zeros([0,ncol],float)

print(" Number of columns: %d \n" % (ncol))

for filename in args.files:
    try:
        data=np.genfromtxt(filename,comments="#")
    except:
        exit_error(" Error opening file : " + filename)

    if (ncol != data.shape[1]): exit_error("Files have different number of columns ")
    alldata=np.concatenate([alldata,data])

alldata = sorted(alldata, key=lambda a_entry: a_entry[args.sortcol]) 

np.savetxt('xhi_joint.out',alldata,fmt='%2.15e')

