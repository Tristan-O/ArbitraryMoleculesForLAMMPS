# DerjaguinTest1

These are tools I am developing for using LAMMPS. Polymer.py gives a user the ability to generate arbitrarily complex polymer structures and write those out to LAMMPS data files. This currently only has support for generating atoms and bonds, but I am planning to give it the ability to add angles, dihedrals and more soon. Simulation_Tools.py provides some tools for analyzing data output by LAMMPS and also contains some pre-defined polymer structures that may be useful for a user (currently beadspring, multiblock, spine and ring). Look in exampleUsage/ for an example of how to use simulationTools.py and Polymer.py to generate and analyze LAMMPS data, or continue reading.

### How to use Polymer.py and Simulation_Tools.py

#### Simulation_Tools
To get start you will want to import Simulation_Tools.py:
```Python
import Simulation_Tools as sim
```
This will give you access to all of the tools in Simulation_Tools. For instance, it will allow you to generate some of the prebuilt polymer structures:
```Python
Lx = 100
Ly = 60
Lz = 60
n_p = 6 # number of polymers
n_s = None # number of solvent atoms
N=10 # number atoms per segment
Rc = 0
NSegments=5 # number of segments per ring

sim.make_dat_file(Lx, Ly, Lz, n_p, n_s=n_s, N=N, NSegments=NSegments, style='ring', outputFile='./polymer0.data')
```
Running this script will generate a LAMMPS data file containing 6 ring polymers with 5 blocks of 10 atoms each, and no solvent atoms. The data file will be called 'polymer0.data'. These polymers will be contained in a box with dimensions 100x60x60 (units to be specified in the accompanying input file for LAMMPS that the user must define).
make_dat_file has a number of style arguments that allow the user to specify what style of polymer they wish to generate. They are:
1. 'beadspring' - requires N to be specified as int - the number of monomers in your polymer
2. 'multiblock' - requires N to specified as int or list - with N as int this is equivalent to creating a beadspring of length N+2. If N is a list it should be a list of ints corresponding to the lengths of each block, not including the ends.
3. 'spine' - requires either N to be specified as a list of 3 ints: the end length, the arm length, and the length between the arms OR an int. Specifying N as int is equivalent to specifying it as a list with the same in three times. to be specified as int 
4. 'star'
5. 'ring'
