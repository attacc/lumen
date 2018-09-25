#!/usr/bin/python3
import argparse
import numpy as np
import re
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt

"""
Calculate current using finite differences from the polarization
Author:  C. Attaccalite
"""
#
# parse command line
#
parser = argparse.ArgumentParser(prog='lumen_PP',description='Current from the polarization',epilog="Copyright C. Attaccalite")
parser.add_argument('-f', help="polarization file",type=str , default=None, dest="polname")
args = parser.parse_args()

print("\n * * * Calculate Current from the Polarization * * * \n\n")

args = parser.parse_args()

if args.polname == None:
    print('type "current.py --help" for help ',)
    exit(0)


fs2aut=41.341373336561361  # convertion femptosecond to atomic units of time

data=np.genfromtxt(args.polname,comments="#")

data_der=np.empty_like(data)

delta  =(data[1,0]-data[0,0])*fs2aut
npoints=len(data[:,0])
ndata  =len(data[0,:])-1
print("Number of points: %d " % npoints)
print("Number of column: %d " % ndata)
print("DeltaT          : %f " % delta)


for ip in range(1,npoints-1):
    for idt in range(1,ndata):
        data_der[ip,idt]=(data[ip+1,idt]-data[ip-1,idt])/(2.0*delta)

for idt in range(1,ndata):
    data_der[0,idt]        =(data[1,idt]-data[0,idt])/(delta)
    data_der[npoints-1,idt]=(data[npoints-1,idt]-data[npoints-2,idt])/(delta)

data_der[:,0]=data[:,0]

#for idt in range(1,ndata):
# interp= InterpolatedUnivariateSpline(data[:,0], data[:,idt], k=3)
# derf = interp.derivative()
# data_der[idt]=derf(data[0])

np.savetxt('der_polarization.dat',data_der,fmt='%2.15e')
