# MD Simulation Module for Colloids and Polymers 
''' Updated 2018.02.20 to initialize CG Bead-Spring Polymers with an Explicit Solvent'''

# NOTE: 
#   - Currently only handles polymers with a harmonic bond. No dihedral or improper angles recorded.
#   - Currently changes the polymer ends to a different atom. Will help with post-processing.
 
# Edited by Tristan to allow for the inclusion of two colloids with the polymers and solvents

# import scipy and numpy 
import numpy as np
import scipy as sp
#import atomwrite        #  M.S.Shell's Implementation
import os
import math
import sys

    



def AtomArange(Lx, Ly, Lz, P, DOP, ns, Rc):

    """ Arranges the polymer atoms, N, onto a cubic lattice of side length, L.
            Input:
                N - Number of atoms
                L - Desired side length
            Output:
                Pos(N,3) - position vector
"""
    # N is the total number atoms in the system (DOP*P + ns)
    N = DOP*P + ns + 2*bool(Rc)
    Rc = 1.5*Rc #generate atoms 1.5x the radius away from the colloids, and check that atoms arent being pushed out past other molecules
    numpolyatoms = DOP*P
    # the polymer volume fraction
    Pvolfrac = (N/Lx/Ly/Lz)
    print "Segment Volume Fraction [LJ units, (bseg/sqrt(6))^3]"
    print Pvolfrac
    # make the position array.
    Pos = np.zeros((N,3), float)
    # compute integer grid # of locations for cubic lattice.
    NLat = int(N**(1./3.) + 1.)
    # LatScale scales the spacing between point!
    LatScale = 0.9
    LatSpac = LatScale / NLat
    # LatSpac = LatScale*Ly / NLat
    # LatSpac = LatScale*Lz / NLat
    print "The Lattice Spacing"
    print LatSpac
    # makes array of evenly spaced values.
    rx = np.linspace(-0.5*LatScale*Lx,0.5*LatScale*Lx, math.ceil(numpolyatoms**(1./3.)), dtype=float)
    ry = np.linspace(-0.5*LatScale*Ly,0.5*LatScale*Ly, math.ceil(numpolyatoms**(1./3.)), dtype=float)
    rz = np.linspace(-0.5*LatScale*Lz,0.5*LatScale*Lz, math.ceil(numpolyatoms**(1./3.)), dtype=float)

    #https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/26127012#26127012
    center_colloid1 = [-0.5*0.5*Lx, 0.0, 0.0]
    center_colloid2 = [+0.5*0.5*Lx, 0.0, 0.0]
    num_pts = int(4.*math.pi*Rc**2.) # Number of area units on sphere
    indices = np.arange(0, num_pts, dtype=float) + 0.5

    phi = np.arccos(1 - 2*indices/num_pts)
    theta = math.pi * (1 + 5**0.5) * indices

    points_colloid1 = np.transpose(np.array([(Rc+0.1) * np.cos(theta) * np.sin(phi), (Rc+0.01) * np.sin(theta) * np.sin(phi), (Rc+0.01) * np.cos(phi)]))
    points_colloid2 = np.transpose(np.array([(Rc+0.1) * np.cos(theta) * np.sin(phi), (Rc+0.01) * np.sin(theta) * np.sin(phi), (Rc+0.01) * np.cos(phi)]))
    print np.shape(points_colloid1)
    for i in range(len(points_colloid1)):
        points_colloid1[i] += center_colloid1
        points_colloid2[i] += center_colloid2



    # if Colloids are to be included, dont want atoms inside them, so I will exclude their volume by saying the box 
    # loop through x, y, z positions in lattice until done
    #   for every atom in the system except colloids
    # Tristan Edit: Loop through the polymer beads first

    i = 0
    if bool(Rc):
        N += 2 # add two colloids to number atoms
        Pos[-1] = np.array(center_colloid1, float)
        Pos[-2] = np.array(center_colloid2, float)
        if math.sqrt((Pos[-1][0]-Pos[-2][0])**2.+(Pos[-1][1]-Pos[-2][1])**2.+(Pos[-1][2]-Pos[-2][2])**2.) < Rc:
            print "Box too small, colloids are inside each other"
            print Pos[-1]
            print Pos[-2]
            sys.exit()


    for x in rx:
        for y in ry:
            for z in rz:
                Pos[i] = np.array([x,y,z], float)
                # add a random offset to help initial minimization.
                Offset = 1 * ( np.random.rand(3) - 0.5 )
                Pos[i] = Pos[i] + Offset

                if bool(Rc):
                    r_to_1 = math.sqrt((Pos[i][0]-Pos[-1][0])**2.+(Pos[i][1]-Pos[-1][1])**2.+(Pos[i][2]-Pos[-1][2])**2.)
                    r_to_2 = math.sqrt((Pos[i][0]-Pos[-2][0])**2.+(Pos[i][1]-Pos[-2][1])**2.+(Pos[i][2]-Pos[-2][2])**2.)
                    if r_to_1 < Rc:
                        #if atom is inside of colloid 1, move radially from colloid center to surface
                        #can do this easily just by scaling the offsets by the intial radial distance in xyz
                        minDist = Rc
                        minPoint = []
                        minpnt_indx = 0
                        pnt_indx = 0
                        for point in points_colloid1:
                            temp = math.sqrt((Pos[i][0]-point[0])**2.+(Pos[i][1]-point[1])**2.+(Pos[i][2]-point[2])**2.)
                            if temp < minDist:
                                minDist = temp
                                minPoint = point[:]
                                minpnt_indx = pnt_indx
                            pnt_indx += 1
                        # print '--------'
                        # print minDist
                        # print Pos[i]
                        # print minPoint
                        Pos[i] = minPoint[:]
                        points_colloid1 = np.delete(points_colloid1, minpnt_indx, 0)
                        # print minPoint

                    if r_to_2 < Rc:
                        minDist = Rc
                        minPoint = []
                        minpnt_indx = 0
                        pnt_indx = 0
                        for point in points_colloid2:
                            temp = math.sqrt((Pos[i][0]-point[0])**2.+(Pos[i][1]-point[1])**2.+(Pos[i][2]-point[2])**2.)
                            if temp < minDist:
                                minDist = temp
                                minPoint = point[:]
                                minpnt_indx = pnt_indx
                            pnt_indx += 1
                        Pos[i] = minPoint[:]
                        points_colloid2 = np.delete(points_colloid2, minpnt_indx, 0)
                        # print minPoint
                i += 1
                # if done placing atoms, return.
                if i >= numpolyatoms:
                    break
            if i >= numpolyatoms:
                break
            # So that polymer beads are not connected across periodic boundaries, see Tristan's notes 7-16-18
            rz = np.flip(rz, 0)
        if i >= numpolyatoms:
            break
        # So that polymer beads are not connected across periodic boundaries, see Tristan's notes 7-16-18
        ry = np.flip(ry, 0)

    #Tristan Edit: Loop through the solvent molecules now with different x,y,z values than the polymer ones
    # Reduce chance of overlap further with changing latscale:
    LatScale *= 2./3.
    rx = np.linspace(-0.5*LatScale*Lx,0.5*LatScale*Lx, math.ceil(ns**(1./3.)), dtype=float)
    ry = np.linspace(-0.5*LatScale*Ly,0.5*LatScale*Ly, math.ceil(ns**(1./3.)), dtype=float)
    rz = np.linspace(-0.5*LatScale*Lz,0.5*LatScale*Lz, math.ceil(ns**(1./3.)), dtype=float)

    for x in rx:
        for y in ry:
            for z in rz:
                Pos[i] = np.array([x,y,z], float)
                # add a random offset to help initial minimization.
                Offset = 0.1 * (np.random.rand(3) - 0.5)
                Pos[i] = Pos[i] + Offset

                if bool(Rc):
                    r_to_1 = math.sqrt((Pos[i][0]-Pos[-1][0])**2.+(Pos[i][1]-Pos[-1][1])**2.+(Pos[i][2]-Pos[-1][2])**2.)
                    r_to_2 = math.sqrt((Pos[i][0]-Pos[-2][0])**2.+(Pos[i][1]-Pos[-2][1])**2.+(Pos[i][2]-Pos[-2][2])**2.)
                    if r_to_1 < Rc:
                        #if atom is inside of colloid 1, move radially from colloid center to surface
                        #can do this easily just by scaling the offsets by the intial radial distance in xyz
                        minDist = Rc
                        minPoint = []
                        minpnt_indx = 0
                        pnt_indx = 0
                        for point in points_colloid1:
                            temp = math.sqrt((Pos[i][0]-point[0])**2.+(Pos[i][1]-point[1])**2.+(Pos[i][2]-point[2])**2.)
                            if temp < minDist:
                                minDist = temp
                                minPoint = point[:]
                                minpnt_indx = pnt_indx
                            pnt_indx += 1
                        # try:
                        Pos[i] = minPoint[:]
                        points_colloid1 = np.delete(points_colloid1, minpnt_indx, 0)
                        # except ValueError:
                        #     print Pos[i]
                        #     print Pos[-1]
                        #     print minDist
                        #     print minPoint
                        #     print temp
                        #     print minPoint
                        #     print r_to_1
                        #     print r_to_2
                        #     return
                        # print minPoint
                        # print temp

                    if r_to_2 < Rc:
                        minDist = Rc
                        minPoint = []
                        minpnt_indx = 0
                        pnt_indx = 0
                        for point in points_colloid2:
                            temp = math.sqrt((Pos[i][0]-point[0])**2.+(Pos[i][1]-point[1])**2.+(Pos[i][2]-point[2])**2.)
                            if temp < minDist:
                                minDist = temp
                                minPoint = point[:]
                                minpnt_indx = pnt_indx
                            pnt_indx += 1

                        # try:
                        Pos[i] = minPoint[:]
                        points_colloid2 = np.delete(points_colloid2, minpnt_indx, 0)
                        # except ValueError:
                        #     print Pos[i]
                        #     print Pos[-2]
                        #     print minDist
                        #     print minPoint
                        #     print temp
                        #     print minPoint
                        #     print r_to_1
                        #     print r_to_2
                        #     return
                        # print minPoint
                i += 1
                # if done placing atoms, return.
                if i >= numpolyatoms+ns:
                    break
            if i >= numpolyatoms+ns:
                break
        if i >= numpolyatoms+ns:
            break


    #Checking just to make sure my program doesn't have an error:
    if bool(Rc):
        for pos in Pos[:-2]:
            r_to_1 = math.sqrt((pos[0]-Pos[-1][0])**2.+(pos[1]-Pos[-1][1])**2.+(pos[2]-Pos[-1][2])**2.)
            r_to_2 = math.sqrt((pos[0]-Pos[-2][0])**2.+(pos[1]-Pos[-2][1])**2.+(pos[2]-Pos[-2][2])**2.)
            if r_to_1 < Rc or r_to_2 < Rc:
                print '*************************************************************************************************'
                print 'WARNING: ATOM AT ' + str(pos) + ' IS TOO CLOSE TO A COLLOID!'
                print '\tDistance to colloids is:'
                print r_to_1
                print r_to_2
                print 'When both should be greater than ' + str(Rc)
                print '*************************************************************************************************'

    # print 'Checking Pos'
    # for i in range(len(Pos)):
    #     pos1 = Pos[i]
    #     for j in range(i, len(Pos)):
    #         pos2 = Pos[j]
    #         if math.sqrt((pos1[0]-pos2[0])**2.+(pos1[1]-pos2[1])**2.+(pos1[2]-pos2[2])**2.) == 0.0 and i!=j:
    #             print i
    #             print j
    #             print math.sqrt((pos1[0]-pos2[0])**2.+(pos1[1]-pos2[1])**2.+(pos1[2]-pos2[2])**2.)


        
    
    data(Pos,Lx,Ly,Lz,P, DOP, ns, Rc)
    
    
def InitVis(Pos,Lx,Ly,Lz,P, DOP, ns, Rc):    
    """ Creates a visualization pdb file for use with Chimera
    
    Inputs: 
        Pos - atom positions
        L   - length of the box
    Outputs:
        Pos - the atom positions
""" 

    fn = "anim"
    fnNum = "0"
    fnExt = ".pdb"
    i = 1
    # checks to see if the file already exists
    filename = "%s%s%s" % (fn,fnNum,fnExt)
    #Names = np.arange(len(Pos)) + 1
    #Names.astype(str)
    
    while os.path.isfile(filename) == True:
        fnNum = "%s" % i
        filename = str(fn + fnNum + fnExt)
        i = i + 1
    print "PDB file visualization name"
    print filename
        
    #atomwrite.pdbfile(str(filename), Lx, Compressed = False).write(Pos)
    data(Pos,Lx,Ly,Lz,P,DOP, ns, Rc)
    return Pos

def data(Pos, Lx,Ly,Lz, P, DOP, ns, Rc):
    """ Creates the input data file for LAMMPS simulation.
    Inputs:
        Pos - atom positions
        L   - Box size
        P   - Number of polymers desired
    Outputs:
        Data File 
"""
    SizePos  = np.shape(Pos)
    print "Size of the Position"
    print SizePos
    NumAtoms = SizePos[0]
    print "The number of segments"
    print NumAtoms
    NumAtomsP = DOP
    print "The number of segments per polymer"
    print NumAtomsP
    NumBonds  = (DOP - 1)*P
    fn = "polymer"
    fnNum = "0"
    fnExt = ".data"
    
    i = 1
    filename = "%s%s%s" % (fn,fnNum,fnExt)
    
    # checks to see if file exists
    while os.path.isfile(filename) == True:
        fnNum = "%s" % i
        filename = str(fn + fnNum + fnExt)
        i = i + 1
    print filename
    
    # creates the new file with name "filename"
    f = file(str(filename), "w")
    
    #LAMMPS DESCRIPTION:
    f.write("This is the LAMMPS data file for %s\n" % filename)
    
    #header - box bounds **************************************************************
    f.write("\n%i atoms\n" % NumAtoms)
    f.write("%i bonds\n" % NumBonds)
    f.write("0 angles\n")
    f.write("0 dihedrals\n")
    f.write("0 impropers\n")
    
    ''' Currently set up for just 4 atom types: Polymer-ends = 1 ; polymer-segments = 2 ; solvent molecules = 3 '''
    f.write("\n4 atom types\n")
    '''Harmonic, Gaussian, Tabulated Tanh (colloid-bead/solvent), Tabulated LJ12-6 (colloid-colloid)'''
    f.write("2 bond types\n")
    f.write("0 angle types\n")
    f.write("0 dihedral types\n")
    f.write("0 improper types\n")
    
    #BOX SIZE: *************************************************************************
    f.write("\n%2.3f %2.3f xlo xhi" % (-1.*Lx/2. , 1.*Lx/2.))
    f.write("\n%2.3f %2.3f ylo yhi" % (-1.*Ly/2. , 1.*Ly/2.))
    f.write("\n%2.3f %2.3f zlo zhi\n" % (-1.*Lz/2. , 1.*Lz/2.))
    
    # MASSES: **************************************************************************
    # may need to edit the masses for hydrophobic versus hydrophilic components.
    f.write("\nMasses\n")
    f.write("\n")
    f.write("1 1.0\n") # all same mass
    f.write("2 1.0\n") # all same mass
    f.write("3 1.0\n") # all same mass
    f.write("4 1.0\n")
    
    # ATOMS: ***************************************************************************
    f.write("\nAtoms\n")
    f.write("\n")
    i = 1 # number of molecules counter
    j = 1 # molecule counter
    k = 0 # atom position
    numpolyatoms = P*DOP
    
    ''' This Section labels the polymer end groups. This was originally setup for TR-NEM polymer systems in which the polymer 
            end groups were hydrophobic. '''
    # # At the end of each polymer I want to insert a specific number of solvent molecules
    # num_ns_per_poly = ns / numpolyatoms
    # for ns_indx in range(math.ceil(ns % numpolyatoms)):
    #     # if the number of solvents don't exactly work out to be the number of


    while i <= NumAtoms: 
        
        k  = k + 1      
        f.write("%i " % (i)) # ID or also called atom index
        # if k == 1:
            # #insert a number of solvents before the next polymer
            # for ns_indx in range(num_ns_per_poly):
            #     ##f.write("%i " % (i)) # ID or also called atom index
            #     f.write("3 ") # the molecule is a solvent molecule
                
            #     f.write("%2.6f %2.6f %2.6f " % (Pos[i-1,0], Pos[i-1,1], Pos[i-1,2])) # positions
            #     f.write("%i " % j) # molecular tag
            #     f.write("1.0 ") # diameter
            #     f.write("1.0 ") # density
            
            #     f.write("0 0 0\n")
            #     j = j + 1 # increase the molecular tag index
            #     i = i + 1

        if i <= numpolyatoms:
            if k == NumAtomsP or k == 1:
                f.write("1 ") # the atom is an end-group
                if k == NumAtomsP: 
                    k = 0
            else:
                f.write("2 ") # the atom is in the middle of the polymer  
            
        
            f.write("%2.6f %2.6f %2.6f " % (Pos[i-1,0], Pos[i-1,1], Pos[i-1,2])) # positions
            f.write("%i " % j) # molecular tag
            f.write("1.0 ") # diameter
            f.write("1.0 ") # density
        
            f.write("0 0 0\n")
        
            #cheak if on new chain 
            chk = 1.*i/NumAtomsP
            if float.is_integer(chk):
                j = j + 1
            i = i + 1
        elif i>numpolyatoms and i <= NumAtoms - 2*bool(Rc): # num Colloids is 2, bool(Rc) is 0 if Rc is None or 0
            
            f.write("3 ") # the molecule is a solvent molecule
            
            f.write("%2.6f %2.6f %2.6f " % (Pos[i-1,0], Pos[i-1,1], Pos[i-1,2])) # positions
            f.write("%i " % j) # molecular tag
            f.write("1.0 ") # diameter
            f.write("1.0 ") # density
        
            f.write("0 0 0\n")
            j = j + 1 # increase the molecular tag index
            i = i + 1

        else:
            f.write("4 ") # the molecule is a colloid
            
            f.write("%2.6f %2.6f %2.6f " % (Pos[i-1,0], Pos[i-1,1], Pos[i-1,2])) # positions
            f.write("%i " % j) # molecular tag
            f.write("1.0 ") # diameter
            f.write("1.0 ") # density
        
            f.write("0 0 0\n")
            j = j + 1 # increase the molecular tag index
            i = i + 1



    # BONDS: *******************************************************************************
    f.write("\nBonds\n")
    f.write("\n")
    i = 1
    j = 1
    k = 1
    
    while k < numpolyatoms:
        if j == NumAtomsP:
            j = 1
            k = k + 1
        f.write("%i " % i)
        if j == NumAtomsP - 1 or j == 1:
            f.write("1 ") # the bond is hydrophobic-to-philic
        else:
            f.write("2 ") # the bond is hydrophobic-to-philic
        f.write("%i %i\n" % (k, k+1)) 
        j = j + 1
        i = i + 1
        k = k + 1
    
    f.close()
    
    return Pos
    
if __name__ == "__main__":
    import argparse as ap
    parser = ap.ArgumentParser(description='generates input LAMMPS congfiguration file and PDB visualization file')
    parser.add_argument('-N','--DOP',type=int,default=100, help='specify the polymer number degree of polymerication, number of statistical segments')
    parser.add_argument('-np','--numberpolymer',type=int, default=1, help='specify the number of polymers in the system')
    parser.add_argument('-Lx','--LboxSidex', type=float, default=10, help='specify the x-axis box side length for the simulation, LJ units, e.g. in units of statistical segment/sqrt(6)')
    parser.add_argument('-Ly','--LboxSidey', type=float, default=10, help='specify the y-axis box side length for the simulation, LJ units, e.g. in units of statistical segment/sqrt(6)')
    parser.add_argument('-Lz','--LboxSidez', type=float, default=10, help='specify the z-axis box side length for the simulation, LJ units, e.g. in units of statistical segment/sqrt(6)')
    parser.add_argument('-ns','--numbersolvent', type=int, default=0, help='specify the number of CG solvent molecules in the system')
    parser.add_argument('-Rc','--colloidRadius', type=float, default=0.0, help='Include two colloids in simulation')
    args=parser.parse_args()
    
    AtomArange(args.LboxSidex,args.LboxSidey,args.LboxSidez, args.numberpolymer, ars.DOP, args.numbersolvent, args.colloidRadius)