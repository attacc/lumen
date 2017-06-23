#!/usr/bin/python3
import argparse
import numpy as np
import re
import sys

"""
Extract two-photon absorption and intensity-dependent refractive index from ypp 
Author:  C. Attaccalite and M. Grüning
"""
#
# Units
#
SVCMm12VMm1=29.98*10**3.0
AU2VMm1    =5.14220632*10**11.0
#
# Parameters
#
N_Col = 7
Def_Harm = 4
Def_Int = 2
Def_header="""

 Intensity dependent refractive index and two-photon absorption

E[eV] Im[Xhi3_x]  Re[Xhi3_x]  Im[Xhi3_y]  Re[Xhi3_y]  Im[Xhi3_z]  Re[Xhi3_z]
"""
#
# Useful functions
#
def Extract_third_order(inp):
    intensities = np.shape(inp)[0]
    out = np.zeros((intensities-1),dtype=complex)
    for II in range(intensities-1):
        out[II] = 4.0*(inp[II] - inp[II+1])/3.0
    return out

def Eliminate_fifth_order(inp):
    out = (4.0*inp[1] - inp[0])/3.0
    return out
#
# parse command line
#
parser = argparse.ArgumentParser(prog='NLAbsorption_PP',description='Extract two-photon absorption and intensity-dependent refractive index from ypp',epilog="Copyright C. Attaccalite and M. Grüning 2017")
parser.add_argument('-nr', help="number of intensities",type=int , default=Def_Int, dest="nR")
parser.add_argument('-J', help="job identifier",type=str , default=None, dest="jobname")
args = parser.parse_args()

print("\n * * * Extract two-photon absorption and intensity-dependent refractive index * * * \n\n")
#
# Assign input values
#
args = parser.parse_args()

N_int = args.nR

if (N_int<2):
    print("Error: At least two intensities needed !!")
    sys.exit(1)

if args.jobname == None:
    file_begin="o.YPP-X_probe"
else:
    file_begin="o-"+args.jobname+".YPP-X_probe"
#
# Read Number of frequencies from file...
#
xhi3=open(file_begin+"_int_1_order_1","r")
lines=xhi3.read()
xhi3.close()
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
# ... and Electric field
#
real_sn=r'[+\-]?\d*\.\d*[E][+\-]\d\d?'
pattern=r'Efield Denominator =\s*('+real_sn+')\s*('+real_sn+')'
try:
    match = re.search(pattern, lines, re.MULTILINE)
    Divide_Efield=float(match.group(1))+1j*float(match.group(2))
except:
    print("Error reading efield denominator !!")
    sys.exit(1)
print("Efield denominator (in au): %s " % (str(Divide_Efield)))
Divide_Efield=Divide_Efield*SVCMm12VMm1/AU2VMm1
#
#  Read first susceptibilities from files...
#  
print("\nReading ypp response function... \n")

XHI=np.zeros([N_int,nfreqs,N_Col],dtype=float)

for iR in range(0,N_int):
        file_name=file_begin+"_int_"+str(iR+1)+"_order_1"
        print("Reading %s " % file_name)
        try:
            XHI[iR,:,:]=np.genfromtxt(file_name,comments="#")
        except:
            print("Error reading file "+file_name+" !! ")
            sys.exit(1)
#
# Extract TPA and intensity-dependent refractive index
#
print("\nWriting X3(\omega)... \n")
XHI_out=np.zeros([nfreqs,N_Col],dtype=float)
X_in=np.zeros([N_int],dtype=complex)
X_out=np.zeros([N_int-1],dtype=complex)
XHI_out[:,0]=XHI[0,:,0]  # copy energies
for dim in range(3):
    for freq in range(nfreqs):
        X_in[:] = (1j*XHI[:,freq,2*dim+1]+XHI[:,freq,2*(dim+1)])
        X_out=Extract_third_order(X_in)*Divide_Efield**2
        if (N_int>2):
            X_out[1] = 4*X_out[1]
            X_tmp=Eliminate_fifth_order(X_out)
        else:
            X_tmp = X_out[0]
        XHI_out[freq,2*dim+1]= X_tmp.imag
        XHI_out[freq,2*(dim+1)]= X_tmp.real
        #
xhi_filename_out=file_begin+"_intensity_dependent_NL"
print("Writing %s " % xhi_filename_out)
xhi_header=Def_header
np.savetxt(xhi_filename_out,XHI_out[:,:],header=xhi_header)
#
