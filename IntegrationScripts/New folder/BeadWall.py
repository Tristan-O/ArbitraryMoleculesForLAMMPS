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
    xmax = 2.**(1./6.)*Rc # ***********************************************************************************************************************************
    NumberXpts = 10000
    
    ''' LJ Parameters '''
    # LJ 9-3
    E1 = 1E-2
    sigma1 = 17.
    rcut93 = (0.4)**(1/6)*sigma1
    # LJ 12-6
    E2 = 1E-4
    sigma2 = Rc / (2.**(1./6.))# 16 ********* I DON'T KNOW IF I WAS SUPPOSED TO ACTUALLY SHIFT THIS POTENTIAL OVER OR JUST CHANGE SIGMA SO IT GOES TO ZERO AT THE RIGHT SPOT *************
    rcut126 = (2.**(1./6.))*sigma2
#shift = 
    #rcut126 = sigma2
    #rcut126 = 12
    # Colloid 
    # E3 = 4*pi^2*rho_wall*rho_colloid*epsilon*sigma^6 (where epsilon and simga are LJ parameters)
    E3 = 1E9
    sigma3 = 8.
    rcutcolloid = xmax
    # Harmonic
    E4 = 120
    rcut = Rc + 1.25
    sigma = 1. 
    
    ''' Yukawa '''
    E5 = 1E1
    kappa = 0.08
    rcutyukawa = Rc + 5
    
    
    
    
    y = np.linspace(0.0001,xmax,NumberXpts)
    Integrand = []
    LJ93 = []
    LJ126 = []
    Colloid = []
    Harmonic = []
    Yukawa = []
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
        
        ''' Harmonic Interaction '''
        temp5 = E4*(i-rcut)**2
        if i > rcut:
            Harmonic.append(0)
        else:
            Harmonic.append(temp5)
        
   
    y1 = np.linspace(0.0001,xmax,NumberXpts)
    for i in y1:
        ''' LJ-9/3 Interaction '''
        r93min = (2/5.)**(0.16667)
        r1 = (sigma1/i)**3
        r2 = (sigma1/i)**9
        r1min = (sigma1/rcut93)**3
        r2min = (sigma1/rcut93)**9
        #print r93min
        if i > rcut93:
            temp2 = 0.
        else: 
            temp2 = E1*(2*r2/15 - r1) - E1*(2*r2min/15 - r1min)
        LJ93.append(temp2)
        
        ''' LJ-12/6 Interaction '''
        r126min = 2.0**(0.16667)
        #print r126min
        if i > rcut126:
            temp3= 0.
        else:
            temp3 = 4*E2*(((i)/sigma2)**(-12) - ((i)/sigma2)**(-6)) - 4*E2*((((rcut126)/sigma2)**(-12) - (((rcut126)/sigma2)**(-6))))
        LJ126.append(temp3)
        
        ''' Colloid Interaction '''
        if i > rcutcolloid:
            temp4 = 0.
        else:
            i2 = i + sigma3/2
            temp4 = E3*( sigma**6./7560*((6*sigma3/2-i2)/i2**7 + (i2+8.*sigma3/2)/(i2+2.*sigma3/2)**7) - 1/6*((2.*sigma3/2.*(i2+sigma3/2)+i2*(i2+2.*sigma3/2)*(np.log(i2)-np.log(i2+2.*sigma3/2)))/(i2*(i2+2.*sigma3/2))) )
        Colloid.append(temp4)
        
        ''' Yukawa Interaction '''
        if i > rcutyukawa:
            temp5 = 0.
        else:
            temp5 = E5*np.exp(-1*kappa*i)/i - E5*np.exp(-1*kappa*rcutyukawa)/rcutyukawa
        Yukawa.append(temp5)
        
    # Calculate second excluded volume
    int1 = []
    int2 = []
    int3 = []
    int4 = []
    for i in Integrand:
        temp1 = (1-np.exp(-1*i))
        #print temp1
        int1.append(temp1)
    for i in Harmonic:
        temp2 = (1-np.exp(-1*i))
        #print temp2
        int2.append(temp2)
    for i in LJ126:
        temp3 = (1-np.exp(-1*i))
        int3.append(temp3)
    for i in Yukawa:
        temp5 = (1-np.exp(-1*i))
        int4.append(temp5)
        
        
    IntBeadWall = sp.integrate.simps(int1,y)
    print "Exc. Vol. Bead-Wall"
    print IntBeadWall
    IntHarmonic = sp.integrate.simps(int2,y)
    print "Exc. Vol. Harmonic"
    print IntHarmonic
    IntLJ126 = sp.integrate.simps(int3,y1)
    print "Exc. Vol. LJ-12/6"
    print IntLJ126
    IntYukawa = sp.integrate.simps(int4,y1)
    print "Exc. Vol. Yukawa"
    print IntYukawa
        
        
    y1 = y1
    y2 = y1
    y3 = y1 + Rc

    index = np.arange(1,np.shape(y1)[0]+1, dtype=int)


    # AlignedWallBead = np.column_stack((index,y,Integrand,findForce(y, Integrand)))

    # AlignedLJ93 = np.column_stack((index,y1,LJ93,findForce(y1, LJ93)))
    # AlignedLJ126 = np.column_stack((index,y2,LJ126,findForce(y2, LJ126)))
    # AlignedColloid = np.column_stack((index,y3,Colloid,findForce(y3, Colloid)))
    # AlignedHarmonic = np.column_stack((index,y,Harmonic,findForce(y, Harmonic)))
    # AlignedYukawa = np.column_stack((index,y1,Yukawa, findForce(y1,Yukawa)))

    #np.savetxt("WallBeadInteraction.txt", AlignedWallBead, header="#Excee="+str(Excee)+" " +"a="+str(a)+" " +"dw="+str(dw)+" " +"x0="+str(x0)+"\nTANH_PAIR_POTENTIAL\nN "+str(NumberXpts)+"\n")
    with open('BeadWallInteraction.txt', 'w') as f:
        f.write("#Excee="+str(Excee)+" " +"a="+str(a)+" " +"dw="+str(dw)+" " +"Rc="+str(Rc)+"\n")
        f.write("TANH\n")
        f.write("N "+str(NumberXpts)+ " R " + str(y[0]) + " " + str(y[-1]) + "\n")
        f.write('\n')      
        force = findForce(y,Integrand)
        for i in range(len(index)):
            f.write('%i %f %E %E\n' % (index[i], y[i], Integrand[i], force[i]))

    #np.savetxt("WallBeadLJ93.txt", AlignedLJ93, header="#E1="+str(E1)+" " +"sigma="+str(sigma))
    #np.savetxt("WallBeadLJ126.txt", AlignedLJ126, header="E2="+str(E2)+" " +"sigma="+str(sigma)+"\nLJ12-6_PAIR_POTENTIAL\nN "+str(NumberXpts)+"\n")
    with open('BeadWallLJ126.txt', 'w') as f:
        f.write("#E2="+str(E2)+" " +"sigma="+str(sigma)+"\n")
        f.write("LJ126\n")
        f.write("N "+str(NumberXpts)+" R " + str(y2[0]) + " " + str(y2[-1]) + "\n")
        f.write('\n')

        force = findForce(y2,LJ126)
        for i in range(len(index)):
            f.write('%i %f %E %E\n' % (index[i], y2[i], LJ126[i], force[i]))
    # np.savetxt("WallBeadColloid.txt", AlignedColloid, header="E3="+str(E3)+" " +"sigma="+str(sigma))
    # np.savetxt("WallBeadHarmonic.txt", AlignedHarmonic, header="E4="+str(E4)+" " +"rcut="+str(rcut))
    # np.savetxt("WallBeadYukawa.txt", AlignedHarmonic, header="E4="+str(E4)+" " +"rcut="+str(rcut))
    
    plt.plot(y,Integrand,'b-', y1,LJ126, 'r-')
    plt.ylim(0,Integrand[0]*1.5)
    plt.xlim((0,xmax))
    plt.savefig('TANH_ColloidInt_And_LJ126_Colloid-AtomInt.png')
    plt.close()
    # plt.title('Interface-Bead Interaction Integrand')
    # plt.ylabel('$W/k_bT$')
    # plt.xlabel("$h [b_s/6^{0.5}]$")
    # plt.savefig('Interface-Bead Interaction.png')
    # plt.savefig('Interface-Bead Interaction')
    # plt.close()
    
    # plt.plot(y,Integrand,'-',y1,LJ93,'b-',y2,LJ126,'r-',y3,Colloid,'g-',y,Harmonic,'y-',y1,Yukawa,'k-')
    # # plt.ylim(-0.1*Excee,(Excee+0.1*Excee))
    # # plt.xlim((0,xmax))
    # plt.title('Interface-Bead Interaction')
    # plt.ylabel('$W/k_bT$')
    # plt.xlabel("$h [b_s/6^{0.5}]$")
    # #plt.setp(linewidth=3.0)
    # plt.legend(('Compressibility'+" "+"$\epsilon= $"+str(Excee)+"\n"+" "+"a="+str(a)+" "+"$d_w=$"+str(dw),'LJ-9/3'+" "+"$\epsilon= $"+str(E1),'LJ-12/6'+" "+"$\epsilon= $"+str(E2),'colloid'+" "+"$\epsilon= $"+str(E3),'harmonic'+" "+"$\epsilon= $"+str(E4)+"\n"+"rcut="+str(rcut),'yukawa'+" "+"$\epsilon= $"+str(E5)+"\n"+"rcut="+str(rcutyukawa)+"kappa="+str(kappa)))
    # plt.savefig('Interface-Bead Interaction & LJ.png')
    # plt.savefig('Interface-Bead Interaction & LJ')
    # plt.close()
    
    # plt.plot(y,Integrand,'-',y,findForce(y,Integrand))
    # # plt.ylim(-0.1*Excee,(1*Excee+0.01*Excee))
    # # plt.xlim((0,xmax))
    # plt.title('Interface-Bead Interaction')
    # plt.ylabel('$W/k_bT$')
    # plt.xlabel("$h [b_s/6^{0.5}]$")
    # #plt.setp(linewidth=3.0)
    # plt.legend(('Compressibility'+" "+"$\epsilon= $"+str(Excee)+"\n"+" "+"a="+str(a)+" "+"$d_w=$"+str(dw),'harmonic'+" "+"$\epsilon= $"+str(E4)+"\n"+"rcut="+str(rcut)))
    # plt.savefig('Interface-Bead Interaction & Harmonic Approximation.png')
    # plt.savefig('Interface-Bead Interaction & Harmonic Approximation')
    # plt.close()
    
    
if __name__ == "__main__":                  # Code that runs the program if called from the command line
    BeadWall()
#     # saving output
#     AlignedWallBead = np.column_stack((y,Integrand))
#     AlignedLJ93 = np.column_stack((y1,LJ93))
#     AlignedLJ126 = np.column_stack((y2,LJ126))
#     AlignedColloid = np.column_stack((y3,Colloid))
#     AlignedHarmonic = np.column_stack((y,Harmonic))
#     AlignedYukawa = np.column_stack((y1,Yukawa))
#     np.savetxt("WallBeadInteraction.txt", AlignedWallBead, header="Excee="+str(Excee)+" " +"a="+str(a)+" " +"dw="+str(dw)+" " +"x0="+str(x0))
#     np.savetxt("WallBeadLJ93.txt", AlignedLJ93, header="E1="+str(E1)+" " +"sigma="+str(sigma))
#     np.savetxt("WallBeadLJ126.txt", AlignedLJ126, header="E2="+str(E2)+" " +"sigma="+str(sigma))
#     np.savetxt("WallBeadColloid.txt", AlignedColloid, header="E3="+str(E3)+" " +"sigma="+str(sigma))
#     np.savetxt("WallBeadHarmonic.txt", AlignedHarmonic, header="E4="+str(E4)+" " +"rcut="+str(rcut))
#     np.savetxt("WallBeadYukawa.txt", AlignedHarmonic, header="E4="+str(E4)+" " +"rcut="+str(rcut))
    
#     plt.plot(y,Integrand,'-')
#     plt.title('Interface-Bead Interaction')
#     plt.ylabel('$W/k_bT$')
#     plt.xlabel("$h [b_s/6^{0.5}]$")
#     plt.savefig('Interface-Bead Interaction.png')
#     plt.savefig('Interface-Bead Interaction')
#     plt.close()
    
#     plt.plot(y,Integrand,'-',y1,LJ93,'b-',y2,LJ126,'r-',y3,Colloid,'g-',y,Harmonic,'y-',y1,Yukawa,'k-')
#     plt.ylim(-0.1*Excee,(Excee+0.1*Excee))
#     plt.xlim((0,xmax))
#     plt.title('Interface-Bead Interaction')
#     plt.ylabel('$W/k_bT$')
#     plt.xlabel("$h [b_s/6^{0.5}]$")
#     #plt.setp(linewidth=3.0)
#     plt.legend(('Compressibility'+" "+"$\epsilon= $"+str(Excee)+"\n"+" "+"a="+str(a)+" "+"$d_w=$"+str(dw),'LJ-9/3'+" "+"$\epsilon= $"+str(E1),'LJ-12/6'+" "+"$\epsilon= $"+str(E2),'colloid'+" "+"$\epsilon= $"+str(E3),'harmonic'+" "+"$\epsilon= $"+str(E4)+"\n"+"rcut="+str(rcut),'yukawa'+" "+"$\epsilon= $"+str(E5)+"\n"+"rcut="+str(rcutyukawa)+"kappa="+str(kappa)))
#     plt.savefig('Interface-Bead Interaction & LJ.png')
#     plt.savefig('Interface-Bead Interaction & LJ')
#     plt.close()
    
#     plt.plot(y,Integrand,'-',y,Harmonic,'y-')
#     plt.ylim(-0.1*Excee,(1*Excee+0.01*Excee))
#     plt.xlim((0,xmax))
#     plt.title('Interface-Bead Interaction')
#     plt.ylabel('$W/k_bT$')
#     plt.xlabel("$h [b_s/6^{0.5}]$")
#     #plt.setp(linewidth=3.0)
#     plt.legend(('Compressibility'+" "+"$\epsilon= $"+str(Excee)+"\n"+" "+"a="+str(a)+" "+"$d_w=$"+str(dw),'harmonic'+" "+"$\epsilon= $"+str(E4)+"\n"+"rcut="+str(rcut)))
#     plt.savefig('Interface-Bead Interaction & Harmonic Approximation.png')
#     plt.savefig('Interface-Bead Interaction & Harmonic Approximation')
#     plt.close()
    
    
# if __name__ == "__main__":					# Code that runs the program if called from the command line
# 	BeadWall()		