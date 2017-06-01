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

xhi0=open(file_begin+"_int_1_order_0","r")
lines=xhi0.read()

#
# Read the number of frequencies
#

pattern=r'Number of freqs  :\s*(\d*)'
try:
    match = re.search(pattern, lines, re.MULTILINE)
    nfreqs= int(match.group(1))
except:
    print("Error reading nfreqs !!")
    sys.exit(1)
#
# Read efield intensity
#
#pattern=r'Number of freqs  :\s*(\d*)'
#try:
#    match = re.search(pattern, lines, re.MULTILINE)
#    nfreqs= int(match.group(1))
#except:
#    print("Error reading nfreqs !!")
#    sys.exit(1)
#
xhi0.close()


print("Number of frequency step: %d \n " % nfreqs)

XHI=np.zeros([args.nR,args.nX,nfreqs,7],dtype=float)

for iR in range(0,args.nR):
    for iX in range(0,args.nX):
        file_name=file_begin+"_int_"+str(iR+1)+"_order_"+str(iX)
        print("Reading %s " % file_name)
        try:
            XHI[iR,iX,:,:]=np.genfromtxt(file_name,comments="#")
        except:
            print("Error reading file "+file_name+" !! ")
            sys.exit(1)
#
# Apply Richardson to correct XHI2(2w) 
# XHI2(2w: w, w )
#
# Remove any possible linear dependence from the field
# intensity
#
# XHI2 = 2/E^2 [ P(E) - 2 * P(E/2) ] = 2 XHI2(E) - XHI2(E/2)

XHI2=np.zeros([nfreqs,7],dtype=float)

XHI2=2.0*XHI[0,2,:,:]-XHI[1,2,:,:]

#
# Apply Richardson to correct XHI3(3w) 
# XHI3(3w: w, w, w )
#
XHI3=np.zeros([nfreqs,7],dtype=float)
#
# Kerr and Two-photon absorption 
# XHI3(w: w, -w , w)
#
KERR=np.zeros([nfreqs,7],dtype=float)
#
E_square=1.0
#
if args.nR == 2:
    #
    # Remove any linear dependence 
    # from the field intensity in XH3
    #
    # XHI3 = 1/3 [ 4 P(E)/E^3 - 8/E^3 P(E/2) ] = 1/3 [4 XHI3(E) - XHI3(E/2) ]
    #
    XHI3=1.0/3.0*(4.0*XHI[0,3,:,:]-XHI[1,3,:,:])
    #
    # And for the Kerr Kerr
    #
    # KERR = 4/3 [ XHI(E) - XHI(E/2) ] /E^2
    #
    KERR[:,1:]=4.0/3.0*(XHI[0,1,:,1:]-XHI[1,1,:,1:])/E_square
    #
elif args.nR == 3:
    #
    # Remove any linear and quadratic dependence 
    # from the field intensity in XH3
    #
    # P1 = 1/3 [ 4 P(E)   - 8 P(E/2) ] 
    # P2 = 1/3 [ 4 P(E/2) - 8 P(E/4) ] 
    # XHI3 = 2 * [ P1 -4 * P2]/E^3 = 1/3 [ 8 XHI3(E) - 6 * XHI3(E/2) + XHI3(E/4) ]
    #
    XHI3    = 1.0/3.0*(8.0* XHI[0,3,:,:] - 6.0 * XHI[1,3,:,:] + XHI[2,3,:,:])
    #
    # KERR = 2/3 * [ 4 XHI(E) - 12 XHI(E/2) + 8 XHI(E/4) ]/E^2
    #
    KERR    = 2.0/3.0*(4.*XHI[0,1,:,:] - 12.0*XHI[1,1,:,:] + 8.0*XHI[2,1,:,:])/E_square
    #

XHI2[:,0]=XHI[0,1,:,0]
XHI3[:,0]=XHI[0,1,:,0]
KERR[:,0]=XHI[0,1,:,0]


xhi2_header="""

 Second Harmonic Generatio\n")
 corrected with Richardson extrapolation\n")

E[eV] Im[Xhi2_x]  Re[Xhi2_x]  Im[Xhi2_y]  Re[Xhi2_y]  Im[Xhi2_z]  Re[Xhi2_z]
"""
np.savetxt("xhi2.dat",XHI2,header=xhi2_header)


xhi3_header="""

 Third Harmonic Generation
 corrected with Richardson extrapolation

E[eV] Im[Xhi3_x]  Re[Xhi3_x]  Im[Xhi3_y]  Re[Xhi3_y]  Im[Xhi3_z]  Re[Xhi3_z]
"""
np.savetxt("xhi3.dat",XHI3,header=xhi3_header)


kerr_header="""

 Kerr and Two-Photon absorption (TPA)
 obtained by Richardson extrapolation

E[eV] Im[Xhi3_x]  Re[Xhi3_x]  Im[Xhi3_y]  Re[Xhi3_y]  Im[Xhi3_z]  Re[Xhi3_z]
"""
np.savetxt("kerr.dat",KERR,header=kerr_header)

