#!/usr/bin/python3
import argparse
import numpy as np
import re
import sys
from scipy.interpolate import InterpolatedUnivariateSpline

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

data=np.genfromtxt(args.polname,comments="#")

f1 = InterpolatedUnivariateSpline(data[0], data[1], k=3)
f2 = InterpolatedUnivariateSpline(data[0], data[2], k=3)
f3 = InterpolatedUnivariateSpline(data[0], data[3], k=3)

df1 = f1.derivative()
df2 = f2.derivative()
df3 = f3.derivative()

data_der=np.empty_like(data)

data_der[0]=data[0]
data_der[1]=df1(data[0])
data_der[2]=df2(data[0])
data_der[3]=df3(data[0])

np.savetxt('der_polarization.dat',data,fmt='%2.15e')
