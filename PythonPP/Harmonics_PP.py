#!/usr/bin/python3
import argparse
import numpy as np
import re
import sys

"""
Refine the harmonic generation from ypp by Richardson extrapolation
Author:  C. Attaccalite and M. Grüning
"""
#
# Parameters
#
N_Col = 7
Def_Harm = 4
Def_Int = 2
Def_header="""

 ZZth Harmonic Generation
 refined with Richardson extrapolation

E[eV] Im[XhiZZ_x]  Re[XhiZZ_x]  Im[XhiZZ_y]  Re[XhiZZ_y]  Im[XhiZZ_z]  Re[XhiZZ_z]
"""
#
# Useful functions
#
def Eliminate_lower_harmonics(inp,power):
    intensities = np.shape(inp)[0]
    levels = min(intensities,power)
    out = np.zeros((levels,levels),dtype=float)
    out[0,:]=inp[:levels]
    for m in range(1,levels):
        for n in range(levels-1,0,-1):
            factor=2**(power-m)
            out[m,n-1]=(factor*out[m-1,n-1]-out[m-1,n])/(factor-1)
    return out[levels-1,0]
#
#
# parse command line
#
parser = argparse.ArgumentParser(prog='Harmonics_PP',description='Refine the harmonic generation from ypp by Richardson extrapolation',epilog="Copyright C. Attaccalite and M. Grüning 2017")
parser.add_argument('-nx', help="number of harmonics", type=int , default=Def_Harm, dest="nX")
parser.add_argument('-nr', help="number of intensities",type=int , default=Def_Int, dest="nR")
parser.add_argument('-J', help="job identifier",type=str , default=None, dest="jobname")
args = parser.parse_args()

print("\n * * * Refine the harmonic generation from ypp by Richardson extrapolation * * * \n\n")
#
# Assign input values
#
args = parser.parse_args()

N_harm = args.nX + 1
N_int = args.nR 

if args.jobname == None:
    file_begin="o.YPP-X_probe"
else:
    file_begin="o-"+args.jobname+".YPP-X_probe"
#
# Read Number of frequencies from file...
#
xhi0=open(file_begin+"_int_1_order_0","r")
lines=xhi0.read()
xhi0.close()
#
pattern=r'Number of freqs  :\s*(\d*)'
try:
    match = re.search(pattern, lines, re.MULTILINE)
    nfreqs= int(match.group(1))
except:
    print("Error reading nfreqs !!")
    sys.exit(1)

print("Number of frequency step: %d \n " % nfreqs)
#
#  Read Harmonic generation susceptibilities from files...
#  
print("\nReading ypp response functions... \n")

XHI=np.zeros([N_int,N_harm,nfreqs,N_Col],dtype=float)

for iR in range(0,N_int):
    for iX in range(2,N_harm):
        file_name=file_begin+"_int_"+str(iR+1)+"_order_"+str(iX)
        print("Reading %s " % file_name)
        try:
            XHI[iR,iX,:,:]=np.genfromtxt(file_name,comments="#")
        except:
            print("Error reading file "+file_name+" !! ")
            sys.exit(1)
#
# Refine the Harmonics by eliminating lower spurious harmonics  
#

print("\nWriting refined response functions... \n")

XHI_out=np.zeros([nfreqs,N_Col,N_harm],dtype=float)
for iX in range(2,N_harm):
    XHI_out[:,0,iX]=XHI[0,1,:,0]  # copy energies
    for dim in range(1,N_Col):
        for freq in range(nfreqs):
            XHI_out[freq,dim,iX]=Eliminate_lower_harmonics(XHI[:,iX,freq,dim],iX)
#
    xhi_filename_out=file_begin+"_refined_order_"+str(iX)
    print("Writing %s " % xhi_filename_out)
    xhi_header=Def_header.replace('ZZ',str(iX))
    np.savetxt(xhi_filename_out,XHI_out[:,:,iX],header=xhi_header)
