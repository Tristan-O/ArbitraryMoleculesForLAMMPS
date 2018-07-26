#******************General Imports**********************************
import os
import sys
import subprocess
import numpy as np
import scipy as sp
import shutil
import datetime
import math
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.interpolate import interp1d
from scipy.integrate import quad

# The definition of the integrand for the pair interaction
# A is the prefactor of the gaussian and B is the length scale
def Ugauss(x,A,B):
    return A*np.exp(-B*x**2)

def integrand(x,A,B): 
    #print y, dw, a
    return 4*math.pi*x**2*(1-np.exp(-B*Ugauss(x,A,B)))
    
def ExcludedVolume():  
    
    # GAUSSIAN POTENTIAL
    ''' Gaussian Interaction Length Scale'''
    B = 4.0
    ''' Interaction Strength '''
    y = [0.01,0.1,1.0,10,100,1000,10000,2.2448,22.4483]
    ''' Actual Fluid Volume Fraction '''
    # -This is called the packing fraction in the literature!
    # SEE: https://aip.scitation.org/doi/pdf/10.1063/1.4870524
    # Equation of state and the freezing point in the hard-sphere model, 2014
    # -The Hard-sphere has a freezing point around a packing fraction of 0.492
    # -The Gaussian fluid exhibits unusual behavior at T* ~ 4E-3 (A~1000) and rho* ~ 0.4
    # where rho* is dimensionless number density and T* is kbT/A.
    VfracReal = [0.3,0.4,0.45,0.5,0.55,0.6]
    LengthSimulation = np.linspace(71.1,71.2,11)
    
    excludedvol = []
    for i in y:
        #print i
        temp = quad(integrand,0,np.inf,args=(i,B))
        excludedvol.append(temp[0])
    
    DensityGaussianFluid = []
    temp1 = []
    for j in excludedvol:
        for i in VfracReal:
            temp2 = i/j
            temp1.append(temp2)
        DensityGaussianFluid.append(temp1)
        temp1 = []
    
    ''' The loop below calculates the solvent number in a simulation cell for a given interaction strength (A)'''
    numbersolvent = []
    VolumeSimulation = []
    # choose an excluded volume to use below 
    A = 22.44839
    print 'A = ' + str(A)
    print 'y[8] = ' + str(y[8])
    excvol = excludedvol[8] # make sure corresponds to A
    FluidVolFrac = 0.45
    for i in LengthSimulation:
        densitydesired = FluidVolFrac/excvol
        temp = densitydesired*i**3
        vol = i**3
        numbersolvent.append(temp)
        VolumeSimulation.append(vol)
            
    
    print "Interaction Strength"
    print y
    print "ExcludedVolume"
    print excludedvol
    print "DensityGaussianFluid"
    print DensityGaussianFluid
    
    plt.loglog(y,excludedvol,'-')
    plt.title('Gaussian Fluid Excluded Volume')
    plt.ylabel('exc. vol.')
    plt.xlabel("interaction strength")
    plt.savefig('ExcVolvsStrength.png')
    plt.savefig('ExcVolvsStrength')
    plt.close()
    
    plt.loglog(excludedvol,DensityGaussianFluid,'-')
    plt.legend((VfracReal))
    plt.text(2,1,'varying solvent packing fractions')
    plt.title('Gaussian Fluid Density')
    plt.ylabel("density")
    plt.xlabel("exc. vol.")
    plt.savefig('DensityGaussianFluid.png')
    plt.savefig('DensityGaussianFluid')
    plt.close()
    
    plt.plot(LengthSimulation,numbersolvent,'-')
    plt.title('Number of Gaussian Beads $A=$'+str(A)+" " +"$B=$"+str(B)+" " +"$\phi=$"+str(FluidVolFrac))
    plt.ylabel("number beads")
    plt.xlabel("side length simulation cell")
    plt.savefig('NumberGaussianBeads.png')
    plt.savefig('NumberGaussianBeads')
    plt.close()
    
    AlignedNumber = np.column_stack((LengthSimulation,numbersolvent))
    np.savetxt("NumberGaussVSLengthCell.txt", AlignedNumber, header="A="+str(A)+" " +"B="+str(B)+" " +"phi="+str(FluidVolFrac))
    
    
if __name__ == "__main__":					# Code that runs the program if called from the command line
	ExcludedVolume()		