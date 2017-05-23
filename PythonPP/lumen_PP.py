#!/usr/bin/python3
import argparse
import numpy as np
import re
import sys

"""
Analise ypp output to extract Kerr, two-photon absorption by means of Richardson extrapolation
Author:  C. Attaccalite and M. Grüning
"""
#
# parse command line
#
parser = argparse.ArgumentParser(prog='lumen_PP',description='Analise ypp output to extract Kerr and two-photon absorption',epilog="Copyright C. Attaccalite and M. Grüning 2017")
parser.add_argument('-nx', help="number of harmonics", type=int , default=4, dest="nX")
parser.add_argument('-nr', help="number of intensities",type=int , default=2, dest="nR")
parser.add_argument('-J', help="job identifies",type=str , default=None, dest="jobname")
args = parser.parse_args()

print("\n * * * Analize ypp output to extract Kerr and two-photon absorption * * * \n\n")

args = parser.parse_args()

if args.jobname == None:
    file_begin="o.YPP-X_probe"
else:
    file_begin="o-"+args.jobnem+".YPP-X_probe"

#
# Read the number of frequencies
#
xhi0=open(file_begin+"_int_1_order_0","r")
lines=xhi0.read()
pattern=r'Number of freqs  :\s*(\d*)'
try:
    match = re.search(pattern, lines, re.MULTILINE)
    nfreqs= int(match.group(1))
except:
    exit_error("Error reading nfreqs")
xhi0.close()

print("Number of frequency step: %d \n " % nfreqs)

XHI=np.zeros([args.nR,args.nX,nfreqs,7],dtype=float)

for iR in range(0,args.nR):
    for iX in range(0,args.nX):
        file_name=file_begin+"_int_"+str(iR+1)+"_order_"+str(iX)
        print("Reading %s " % file_name)
        XHI[iR,iX,:,:]=np.genfromtxt(file_name,comments="#")
#
# Apply Richardson to correct XHI2(2w) and XHI3(3w)
# 


#
# Extract Kerr and Two-photon absorption 
# 

