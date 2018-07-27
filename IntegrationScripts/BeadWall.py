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

# The definition of the integrand for fluctuations
def integrand(x,y,dw,a,x0): 
    #print y, dw, a
    return np.exp(-1*(x**2)/2/a**2)*(1-np.tanh((x+y-x0)/dw))*x**2

def findForce(r, U):
    #F= -dU/dr
    if len(r) != len(U):
        print Error
        sys.exit()

    F = [ -(U[1]-U[0]) / (r[1]-r[0]) ]
    for i in range(2,len(U)-1):
        F.append( -(U[i]-U[i+1]) / (r[i]-r[i+1]) )
        # if abs(F[-1]) <= 0.01:
        #     print 'Zero'
        #     print i
        #     print U[i]
        #     print r[i]
    F.append( -(U[-1]-U[-2]) / (r[-1]-r[-2]) )
    F.append( F[-1] ) #I dont think this is the best idea but...
    print('F[-1] = ', F[-1])

    return F
    
def BeadWall():  
    '''
    Remember to account for the C parameter that multiplies everything in the hamiltonian!
    '''
    
    # USER INPUTS
    ''' Interaction Strength '''
    Excee = 1000
    ''' Characteristic width of segment '''
    a = 0.25
    ''' Interface "thickness" '''
    # large dw, thicker interface
    dw = 0.5
    ''' Position of Interface '''
    Rc = 10.
    ''' Distance Maximum '''
    xmax = 2.**(1./6.)*Rc+3 # ***********************************************************************************************************************************
    NumberXpts = 10000
    
    ''' LJ Parameters '''
    # LJ 12-6
    E = 1E-4
    sigma = Rc / (2.**(1./6.))# 16 
    rcut126 = (2.**(1./6.))*sigma
   
    y = np.linspace(0.0001,xmax,NumberXpts)
    y1 = np.linspace(0.0001,Rc,NumberXpts)
    Integrand = []
    LJ126 = []

    # prefactor1 = Excee/2.0
    # prefactor2 has a factor of 2 when integration was changed
    # from (-Inf,Inf) to [0,Inf)
    # prefactor2 = 2/math.sqrt(2.0*math.pi*a**2)
    #print "prefactor1*prefactor2"
    #print prefactor1*prefactor2
    prefactor = 2.*Excee/( (2.*a**2.)**1.5 * math.pi**0.5)
    # Performing the integration at varying distance from the interface
    for i in y:
        temp = quad(integrand,0,np.inf,args=(i,dw,a, Rc))
        #print temp[0]
        temp1 = prefactor*temp[0]
        #print temp1
        Integrand.append(temp1)

    for i in y1:
        ''' LJ-12/6 Interaction '''
        r126min = 2.0**(0.16667)
        #print r126min
        if i > rcut126:
            break
        else:
            temp3 = 4*E*(((i)/sigma)**(-12) - ((i)/sigma)**(-6)) - 4*E*((((rcut126)/sigma)**(-12) - (((rcut126)/sigma)**(-6))))
        LJ126.append(temp3)
        

        
    # Calculate second excluded volume
    int1 = []
    int3 = []
    
    for i in Integrand:
        temp1 = (1-np.exp(-1*i))
        #print temp1
        int1.append(temp1)
    for i in LJ126:
        temp3 = (1-np.exp(-1*i))
        int3.append(temp3)
        
        
    IntBeadWall = sp.integrate.simps(int1,y)
    print "Exc. Vol. Bead-Wall"
    print IntBeadWall 
    IntLJ126 = sp.integrate.simps(int3,y1)
    print "Exc. Vol. LJ-12/6"
    print IntLJ126

    index = np.arange(1,np.shape(y1)[0]+1, dtype=int)
    with open('BeadWallInteraction.txt', 'w') as f:
        f.write("#Excee="+str(Excee)+" " +"a="+str(a)+" " +"dw="+str(dw)+" " +"Rc="+str(Rc)+"\n")
        f.write("TANH\n")
        f.write("N "+str(NumberXpts)+ " R " + str(y[0]) + " " + str(y[-1]) + "\n")
        f.write('\n')      
        force = findForce(y,Integrand)
        for i in range(len(index)):
            f.write('%i %f %E %E\n' % (index[i], y[i], Integrand[i], force[i]))

    with open('BeadWallLJ126.txt', 'w') as f:
        f.write("#E="+str(E)+" " +"sigma="+str(sigma)+"\n")
        f.write("LJ126\n")
        f.write("N "+str(len(y1))+" R " + str(y1[0]) + " " + str(y1[-1]) + "\n")
        f.write('\n')

        force = findForce(y1,LJ126)
        for i in range(len(index)):
            f.write('%i %f %E %E\n' % (index[i], y1[i], LJ126[i], force[i]))

    
    plt.plot(y,Integrand,'b-', y1,LJ126, 'r-')
    plt.ylim(0,Integrand[0]*1.5)
    plt.xlim((0,xmax))
    plt.savefig('TANH_ColloidInt_And_LJ126_Colloid-AtomInt.png')
    plt.close()
    
    
if __name__ == "__main__":                  # Code that runs the program if called from the command line
    BeadWall()