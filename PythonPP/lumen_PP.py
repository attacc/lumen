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
# Units
#
SVCMm12VMm1=29.98*10**3.0
AU2VMm1    =5.14220632*10**11.0

#

# parse command line
#
parser = argparse.ArgumentParser(prog='lumen_PP',description='Analise ypp output to extract Kerr and two-photon absorption',epilog="Copyright C. Attaccalite and M. Grüning 2017")
parser.add_argument('-nx', help="number of harmonics", type=int , default=5, dest="nX")
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

print("Number of frequency step: %d \n " % nfreqs)
xhi0.close()

#
# Read Efield denominators
#
real_sn=r'[+\-]?\d*\.\d*[E][+\-]\d\d?'
#
Divid_Efield=np.zeros([args.nX+1],dtype=complex)
Scale_factor=np.zeros([args.nX+1],dtype=float)

for i_order in range(0,args.nX+1):
    xhi_file=open(file_begin+"_int_1_order_"+str(i_order),"r")
    lines=xhi_file.read()
    pattern=r'Efield Denominator =\s*('+real_sn+')\s*('+real_sn+')'
    try:
        match = re.search(pattern, lines, re.MULTILINE)
        Divid_Efield[i_order]=float(match.group(1))+1j*float(match.group(2))
    except:
        print("Error reading efield denominator !!")
        sys.exit(1)
    print("Efield denominator %d: %s " % (i_order,str(Divid_Efield[i_order])))
    xhi_file.close()

#
# Put units
#
Divid_Efield[0]=Divid_Efield[0]*SVCMm12VMm1/AU2VMm1
for iX in range(2,args.nX+1):
    Divid_Efield[iX]=Divid_Efield[iX]*(SVCMm12VMm1/AU2VMm1)**(iX-1.0)

print("\nReading ypp response functions... \n")

XHI=np.zeros([args.nR,args.nX+1,nfreqs,7],dtype=float)

for iR in range(0,args.nR):
    for iX in range(0,args.nX+1):
        file_name=file_begin+"_int_"+str(iR+1)+"_order_"+str(iX)
        print("Reading %s " % file_name)
        try:
            XHI[iR,iX,:,:]=np.genfromtxt(file_name,comments="#")
        except:
            print("Error reading file "+file_name+" !! ")
            sys.exit(1)
#
# Define Polarization
#
# Polarizations
# P(E), P(E/2),  P(E/4)
#
P  =np.zeros([nfreqs,3,args.nX+1],dtype=complex)
P_2=np.zeros([nfreqs,3,args.nX+1],dtype=complex)
if args.nR == 3:
    P_4=np.zeros([nfreqs,3,args.nX+1],dtype=complex)
#
# P(w=0, E ) = P[:,:,0]
# P(  w, E ) = P[:,:,1]
# P(2*w, E ) = P[:,:,2]
# P(3*w, E ) = P[:,:,3]
# P(4*w, E ) = P[:,:,4]
# P(5*w, E ) = P[:,:,5]
#
for iX in range(0,args.nX+1):
    #
    # P(w; E)
    # 
    P[:,0,iX]=1j*XHI[0,iX,:,1]+XHI[0,iX,:,2] # x
    P[:,1,iX]=1j*XHI[0,iX,:,3]+XHI[0,iX,:,4] # y
    P[:,2,iX]=1j*XHI[0,iX,:,5]+XHI[0,iX,:,6] # z
    P[:,:,iX]=P[:,:,iX]/Divid_Efield[iX]

Scale_factor[0]=1.0/2.0
for iX in range(1,args.nX+1):
    Scale_factor[iX]=(1.0/2.0)**(iX-1.0)

for iX in range(0,args.nX+1):
    #
    # P(w; E/2)
    # 
    P_2[:,0,iX]=1j*XHI[1,iX,:,1]+XHI[1,iX,:,2] # x
    P_2[:,1,iX]=1j*XHI[1,iX,:,3]+XHI[1,iX,:,4] # y
    P_2[:,2,iX]=1j*XHI[1,iX,:,5]+XHI[1,iX,:,6] # z
    P_2[:,:,iX]=P_2[:,:,iX]/(Divid_Efield[iX]/Scale_factor[iX])
    #

if args.nR == 3:
    Scale_factor[0]=1.0/4.0
    for iX in range(1,args.nX+1):
        Scale_factor[iX]=(1.0/4.0)**(iX-1.0)

    for iX in range(0,args.nX+1):
        P_4[:,0,iX]=1j*XHI[2,iX,:,1]+XHI[2,iX,:,2] # x
        P_4[:,1,iX]=1j*XHI[2,iX,:,3]+XHI[2,iX,:,4] # y
        P_4[:,2,iX]=1j*XHI[2,iX,:,5]+XHI[2,iX,:,6] # z
        P_4[:,:,iX]=P_4[:,:,iX]/(Divid_Efield[iX]/Scale_factor[iX])
#
# Apply Richardson to correct XHI2(2w) 
# XHI2(2w: w, w )
#
# Remove any possible linear dependence from the field
# intensity
#
# XHI2 = 2/E^2 [ P(E) - 2 * P(E/2) ]
#
XHI2    =np.zeros([nfreqs,3],dtype=complex)
XHI2_out=np.zeros([nfreqs,7],dtype=float)
#
XHI2 = 2.0*(P[:,:,2]-2.0*P_2[:,:,2])*Divid_Efield[2]
#
# Apply Richardson to correct XHI3(3w) 
# XHI3(3w: w, w, w )
#
XHI3    =np.zeros([nfreqs,3],dtype=complex)
XHI3_out=np.zeros([nfreqs,7],dtype=float)
#
# Kerr and Two-photon absorption 
# XHI3(w: w, -w , w)
#
KERR    =np.zeros([nfreqs,3],dtype=complex)
KERR_out=np.zeros([nfreqs,7],dtype=float)  # to write on file
#


if args.nR == 2:
    #
    # Remove any linear dependence 
    # from the field intensity in XH3 and Kerr
    #
    # XHI3 = 4/3 [ P(E) - 2 P(E/2) ] / E^3 
    #
    XHI3=4.0/3.0*(P[:,:,3]-2.0*P_2[:,:,3])*Divid_Efield[3]
    #
    KERR=4.0/3.0*(P[:,:,1]-2.0*P_2[:,:,1])*Divid_Efield[3]
    #
    #
elif args.nR == 3:
    #
    # Remove any linear and quadratic dependence 
    # from the field intensity in XH3
    #
    # XHI3 = 8/3 * [ P(E) - 6 P(E/2) + 8.0 * P (E/4) ] /E^3
    #
    XHI3    = 8.0/3.0*(P[:,:,3] - 6.0*P_2[:,:,3] + 8.0*P_4[:,:,3])*Divid_Efield[3]
    #
    KERR    = 8.0/3.0*(P[:,:,1] - 6.0*P_2[:,:,1] + 8.0*P_4[:,:,1])*Divid_Efield[3]
    #
#
# Copy energy coloumn
#
KERR_out[:,1]=KERR[:,0].imag
KERR_out[:,2]=KERR[:,0].real
KERR_out[:,3]=KERR[:,1].imag
KERR_out[:,4]=KERR[:,1].real
KERR_out[:,5]=KERR[:,2].imag
KERR_out[:,6]=KERR[:,2].real
KERR_out[:,0]=XHI[0,1,:,0]  # copy energies

XHI2_out[:,1]=XHI2[:,0].imag
XHI2_out[:,2]=XHI2[:,0].real
XHI2_out[:,3]=XHI2[:,1].imag
XHI2_out[:,4]=XHI2[:,1].real
XHI2_out[:,5]=XHI2[:,2].imag
XHI2_out[:,6]=XHI2[:,2].real
XHI2_out[:,0]=XHI[0,1,:,0]  # copy energies

XHI3_out[:,1]=XHI3[:,0].imag
XHI3_out[:,2]=XHI3[:,0].real
XHI3_out[:,3]=XHI3[:,1].imag
XHI3_out[:,4]=XHI3[:,1].real
XHI3_out[:,5]=XHI3[:,2].imag
XHI3_out[:,6]=XHI3[:,2].real
XHI3_out[:,0]=XHI[0,1,:,0]  # copy energies

xhi2_header="""

 Second Harmonic Generatio\n")
 corrected with Richardson extrapolation\n")

E[eV] Im[Xhi2_x]  Re[Xhi2_x]  Im[Xhi2_y]  Re[Xhi2_y]  Im[Xhi2_z]  Re[Xhi2_z]
"""
np.savetxt("xhi2.dat",XHI2_out,header=xhi2_header)


xhi3_header="""

 Third Harmonic Generation
 corrected with Richardson extrapolation

E[eV] Im[Xhi3_x]  Re[Xhi3_x]  Im[Xhi3_y]  Re[Xhi3_y]  Im[Xhi3_z]  Re[Xhi3_z]
"""
np.savetxt("xhi3.dat",XHI3_out,header=xhi3_header)


kerr_header="""

 Kerr and Two-Photon absorption (TPA)
 obtained by Richardson extrapolation

E[eV] Im[Xhi3_x]  Re[Xhi3_x]  Im[Xhi3_y]  Re[Xhi3_y]  Im[Xhi3_z]  Re[Xhi3_z]
"""
np.savetxt("kerr.dat",KERR_out,header=kerr_header)

