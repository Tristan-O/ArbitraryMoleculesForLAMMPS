# MD Simulation Module for Colloids and Polymers 
''' Updated 2018.02.20 to initialize CG Bead-Spring Polymers with an Explicit Solvent'''

# NOTE: 
#   - Currently only handles polymers with a harmonic bond. No dihedral or improper angles recorded.
#   - Currently changes the polymer ends to a different atom. Will help with post-processing.
 
# Edited by Tristan to allow for the inclusion of two colloids with the polymers and solvents

# import scipy and numpy 
import numpy as np
import scipy as sp
import atomwrite        #  M.S.Shell's Implementation
import os
import math

    



def AtomArange(P,DOP,Lx, Ly, Lz,ns,wc):

    """ Arranges the polymer atoms, N, onto a cubic lattice of side length, L.
            Input:
                N - Number of atoms
                L - Desired side length
                wc - with colloid (bool)
            Output:
                Pos(N,3) - position vector
"""
    # N is the total number atoms in the system (DOP*P + ns) without the colloids
    N = DOP*P + ns
    # the polymer volume fraction
    Pvolfrac = (N/Lx/Ly/Lz)
    print "Segment Volume Fraction [LJ units, (bseg/sqrt(6))^3]"
    print Pvolfrac
    # make the position array.
    if wc: 
        Pos = np.zeros((N+2,3), float)
        # position colloids (last two) at -L/2 and L/2
        Pos[-1] = np.array([-Lx/2.,0.,0.], float)
        Pos[-2] = np.array([Lx/2., 0.,0.], float)
        # Make edges of box now the colloid
        L/=2
    else:
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
    rx = np.linspace(-0.5*LatScale*Lx,0.5*LatScale*Lx, math.ceil(N**(1./3.)), dtype=float)
    ry = np.linspace(-0.5*LatScale*Ly,0.5*LatScale*Ly, math.ceil(N**(1./3.)), dtype=float)
    rz = np.linspace(-0.5*LatScale*Lz,0.5*LatScale*Lz, math.ceil(N**(1./3.)), dtype=float)

    # if Colloids are to be included, dont want atoms inside them, so I will exclude their volume by saying the box 
    # loop through x, y, z positions in lattice until done
    #   for every atom in the system except colloids
    i = 0

    for x in rx:
        for y in ry:
            for z in rz:
                Pos[i] = np.array([x,y,z], float)
                # add a random offset to help initial minimization.
                Offset = 0.1 * (np.random.rand(3) - 0.5)
                Pos[i] = Pos[i] + Offset
                i += 1
                # if done placing atoms, return.
                if i >= N:
                    break
            if i >= N:
                break
        if i >= N:
            break
    print Pos
    
    InitVis(Pos,Lx,Ly,Lz,P, DOP, ns, wc)
    
    
def InitVis(Pos,Lx,Ly,Lz,P, DOP, ns, wc):    
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
    data(Pos,Lx,Ly,Lz,P,DOP, ns, wc)
    return Pos

def data(Pos, Lx,Ly,Lz, P, DOP, ns, wc):
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
    
    ''' Currently set up for just 4 atom types: Polymer-ends = 1 ; polymer-segments = 2 ; solvent molecules = 3, colloid = 4 (wc=0 or 1)'''
    f.write("\n%i atom types\n" % (3+1*int(wc)))
    '''Harmonic, Gaussian, Tabulated Tanh (colloid-bead/solvent), Tabulated LJ12-6 (colloid-colloid)'''
    f.write("%i bond types\n" % (2+2*int(wc)))
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
    if wc:
        f.write("4 10.0\n") # Colloid bigger mass
    
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


    while i <= (NumAtoms-2*int(wc)): 
        
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

        if i<= numpolyatoms:
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
        elif i>numpolyatoms:
            
            f.write("3 ") # the molecule is a solvent molecule
            
            f.write("%2.6f %2.6f %2.6f " % (Pos[i-1,0], Pos[i-1,1], Pos[i-1,2])) # positions
            f.write("%i " % j) # molecular tag
            f.write("1.0 ") # diameter
            f.write("1.0 ") # density
        
            f.write("0 0 0\n")
            j = j + 1 # increase the molecular tag index
            i = i + 1

    if wc:
        f.write("%i 4 " % (NumAtoms+1))#this is a colloid
        f.write("%2.6f %2.6f %2.6f " % (Pos[-1,0], Pos[-1,1], Pos[-1,2])) # positions
        f.write("%i " % j) # molecular tag
        f.write("7.0 ") # diameter
        f.write("1.0 ") # density
        f.write("0 0 0\n")

        j += 1
        f.write("%i 4 " % (NumAtoms+2))#this is a colloid
        f.write("%2.6f %2.6f %2.6f " % (Pos[-2,0], Pos[-2,1], Pos[-2,2])) # positions
        f.write("%i " % j) # molecular tag
        f.write("7.0 ") # diameter
        f.write("1.0 ") # density
        f.write("0 0 0\n")



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
    parser.add_argument('-wc','--withcolloid', type=int, default=0, help='specify if two large Colloids are in the system')
    args=parser.parse_args()
    
    AtomArange(args.numberpolymer, args.DOP, args.LboxSidex,args.LboxSidey,args.LboxSidez, args.numbersolvent, int(bool(args.withcolloid)))