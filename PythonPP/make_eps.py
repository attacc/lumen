#!/usr/bin/python3
import argparse
import numpy as np
import re
import sys
import math


def lorentzian(freq, pole, smear):
    return 1.0/math.pi*(smear/((freq-pole)**2+smear**2))


def bose(E,T):
    #
    # Bose function
    #
    if T < 1e-10:
        return 0.0
    else:
        return math.exp(-E/T)/(1.0-math.exp(-E/T))

ha2ev   = 27.211396132
au2kelvin = 3.1577513e5


"""
Build dielectric constant from excitons
Author:  C. Attaccalite
"""
#
# parse command line
#
parser = argparse.ArgumentParser(prog='make_eps',description='Build epsilon from exciton decomposition',epilog="Copyright C. Attaccalite 2018")
parser.add_argument('-f', help="exciton file name",type=str , default="o.exc_E_sorted", dest="exc_fname")
parser.add_argument('-sm', help="smearing (default 0.1 eV)",type=float , default=0.1, dest="smearing")
parser.add_argument('-ns', help="energy steps (default 1000)",type=int , default=1000, dest="ensteps")
parser.add_argument('-er', help="energy range (default [0,10] ev)",type=float , default=[0,10], dest="enrange",nargs=2)
parser.add_argument('-t', help="temperature(default 0.0 [Kelvin])",type=float , default=0.0, dest="temp")
parser.add_argument('-mu', help="chemical potential for the Bose [eV]",type=float , default=-1.0, dest="mu")
args = parser.parse_args()

print("\n * * * Build dielecrict constant from exciton analisys * * * \n\n")

args = parser.parse_args()

data=np.genfromtxt(args.exc_fname,comments="#")

delta_e=(args.enrange[1]-args.enrange[0])/float(args.ensteps)

temp = args.temp/au2kelvin*ha2ev

if args.mu != -1.0:
    print("Temperature [eV] :"+str(temp))
    print("Chemical potential [eV] :"+str(args.mutemp))

eps =np.zeros([args.ensteps,2],dtype=float)
jdos=np.zeros([args.ensteps,2],dtype=float)
lum =np.zeros([args.ensteps,2],dtype=float)

for iw in range(args.ensteps):
    w=args.enrange[0]+delta_e*float(iw)
    eps [iw,0]=w
    jdos[iw,0]=w
    lum [iw,0]=w
    for ixc in range(data.shape[0]):
        eps[iw,1]=eps[iw,1]+lorentzian(w, data[ixc][0], args.smearing)*data[ixc][1]
        jdos[iw,1]=jdos[iw,1]+lorentzian(w, data[ixc][0], args.smearing)
        if args.mu != -1.0:
            exc_en=data[ixc][0]-data[0][0]+args.mu
            lum[iw,1]=lum[iw,1]+lorentzian(w, data[ixc][0], args.smearing)*data[ixc][1]*bose(exc_en,temp)


np.savetxt('eps_excitons.out',eps,fmt='%2.15e')
np.savetxt('jdos_excitons.out',jdos,fmt='%2.15e')
if args.mu != -1.0:
    np.savetxt('lum_excitons.out',lum,fmt='%2.15e')

