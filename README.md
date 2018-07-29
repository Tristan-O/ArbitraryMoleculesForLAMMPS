# DerjaguinTest1

These are tools I am developing for use in conjunction with LAMMPS. Polymer.py gives a user the ability to generate arbitrarily complex polymer structures and write those out to LAMMPS data files. This currently only has support for generating atoms and bonds, but I am planning to give it the ability to add angles, dihedrals and more soon. Simulation_Tools.py provides some tools for analyzing data output by LAMMPS and also contains some pre-defined polymer structures that may be useful for a user (currently beadspring, multiblock, spine and ring). Look in exampleUsage/ for an example of how to use simulationTools.py and Polymer.py to generate and analyze LAMMPS data, or continue reading.

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
N=10 # number monomers per segment
Rc = 0
NSegments=5 # number of segments per ring

sim.make_dat_file(Lx, Ly, Lz, n_p, n_s=n_s, N=N, NSegments=NSegments, style='ring', outputFile='./polymer0.data')
```
Running this script will generate a LAMMPS data file containing 6 ring polymers with 5 blocks of 10 monomers each, and no solvent atoms. The data file will be called 'polymer0.data'. These polymers will be contained in a box with dimensions 100x60x60 (units to be specified in the accompanying input file for LAMMPS that the user must define).
##### make_dat_file(Lx,Ly,Lz, n_p, n_s=None, N=None, NArms=None, NSegments=None, style=None, outputFile='./polymer0.data')
Has a number of style arguments that allow the user to specify what style of polymer they wish to generate. They are:
1. 'beadspring' - Requires N to be specified as int - the number of monomers in your polymer
2. 'multiblock' - Requires N to specified as int or list - with N as int this is equivalent to creating a beadspring of length N+2. If N is a list it should be a list of ints corresponding to the lengths of each block, not including the ends.
3. 'spine' - Requires either N to be specified as a list of 3 ints: the end length, the arm length, and the length between the arms OR an int. Specifying N as int is equivalent to specifying it as a list with the same in three times. 
4. 'star' - Requires EITHER N,NArms to be specified as int (this is the case where every arm has the same length) OR N as list of ints representing the length of each arm.
5. 'ring' - Requires N,NSegments to be specified as ints. Each segment will be a different type and will be N monomers long.

##### clean_log(infileloc, outfileloc, header='Steps Temp KinE PotE TotE Press Vol')
This function parses an output file of LAMMPS and returns a numpy array of the thermodynamic parameters enclosed, provided you specify them. NOTE header must contain ALL of the parameters that the LAMMPS file contains or this function will not work.

##### do_stats(params_list, infileloc, outfileloc=None)
Don't use this function, it is extremely slow, I am working on improving its speed.

##### def plot_params(outfile, data, stats_data='Steps Temp KinE PotE TotE Press Vol', numcols=None, numrows=None, t_lo=None, t_hi=None, txt='')
Produces an array of plots of all of the thermodynamic parameters produced by clean_log() described above and saves it to a file. NOTE this assumes that the timestep is always increasing in your output file, so if you used reset_timestep in your LAMMPS input file it may not work as intended, but you can alter the data argument to account for that.
* outfile is the filename that the plots will be written to.
* 'data' is assumed to be of the same format that was output by clean_log()
* stats_data will be the titles of the plots that it produces. It should be of the same format as the header argument in clean_log().
* numcols, numrows are the number of columns and rows that you wish to organize your plots in, or if the default is producing errors just ensure numcols * numrows is at least the number of plots you expect.
* t_lo,t_hi are optional bounds on the plots, useful if you wish to omit some warmup data
* txt is some text that will be put at the bottom of your plot.

##### analyze_dump(infile, style='beadspring', POLY_END_TYPE = 1, POLY_MID_TYPES = [2], COLLOID_TYPE = 4, TYPE_IGNORES = [3]
WIP - Analyzes an atom-style dump file and produces a histogram of the distances between monomers.



#### Polymer

Polymer organizes your atoms/monomers/whatever in an object-oriented approach. The hierarchy is as follows:
* Box
* Polymer
* Node, Chain, Loop, LoopSegment

The idea is that you can define a Box, and you can define a Polymer's structure by defining the relationships between Nodes, Chains, and LoopSegments. Then once those relationships have been defined (in a relatively simple and straightforward way), you can add that Polymer definition to your Box any number of times. Doing this is what gives atoms a spatial position in the box. Most cases will only need Nodes and Chains, but there is functionality to include polymers with closed-loops as well through objects called "Loop"s and "LoopSegment"s.

Every script using Polymer directly will need the following lines:
```Python
import Polymer as poly #Import the definitions

poly0 = poly.Polymer(TYPE) #First Node in polymer is an atom of type TYPE
node0 = poly0.get_node_by_id(0) #You need to access the first Node to begin to build on it
... # Further define your polymer here
box = poly.Box([Lx,Ly,Lz]) #Initialize a Box with those dimensions
box.add_polymer(poly0) #Box adds this polymer to the box
box.write_box('MYOUTPUTFILE') #Everything that this box knows is written to a LAMMPS-style data file
```
This code by itself will add a "polymer" of just one monomer to your box. Further defining your polymer is done like so:

```Python
import Polymer as poly #Import the definitions
currentType = 1 #notice me!
poly0 = poly.Polymer(currentType) #First Node in polymer is an atom of type 1
node0 = poly0.get_node_by_id(0) #You need to access the first Node to begin to build on it

# Further define your polymer here
# For instance, if I wanted to generate a star polymer (looks like an *):
NArms = 5 #the number of arms I want my star to have
N = 10 #the length of each arm
for i in range(NArms): #loop through NArms times
			currrentType += 1 #give every arm a different type
			chain = poly0.add_chain(N, currentType) #add a Chain object to the polymer
			chain.attach_parent_node(node0) #make that Chain stem from the first Node in the polymer, node0
                                      #this is referred to as making that Node its parent

box = poly.Box([Lx,Ly,Lz]) #Initialize a Box with those dimensions
box.add_polymer(poly0) #Box adds this polymer to the box
box.write_box('MYOUTPUTFILE') #Everything that this box knows is written to a LAMMPS-style data file
```
Note that I added specific examples instead of variable names in the code that was the same as above.

You can easily follow the logic behind this, with maybe a little help from the comments in it.

You can also add more Nodes to the ends of those Chains. Any Node can have as many child Chains as your python compiler will allow,
but a Chain can only have one parent Node and one child Node. Also important to note is that a Node cannot have another Node as its child.

Here is some slightly more complicated code to show how to make a spine polymer (looks something like below):  
  o   o   o
  |   |   |
o-o-o-o-o-o-o ...  
    |   |   |
    o   o   o
    
```Python
Nend = 3 #number monomers on the ends
Nbtwn = 5 #number monomers between the arms
Narm = 4 #number of monomers on each arm
numArms = 10 #number of arms to add
spineType = 1
armType = 2

endchain = poly0.add_chain(Nend-1, spineType) #-1 here because the first Node takes away the need for one
endchain.attach_parent_node(node0)
newNode = poly0.add_node(currAtomType)
endchain.attach_child_node(newNode)

for n in range(NArms-1):
	spinechain = poly0.add_chain(Nbtwn, spineType)
  spinechain.attach_parent_node(newNode)
  armchain = poly0.add_chain(Narm, armType)
	armchain.attach_parent_node(newNode)
	newNode = poly0.add_node(spineType)
	spinechain.attach_child_node(newNode)
			
armchain = poly0.add_chain(Narm, armType)
armchain.attach_parent_node(newNode)
endchain = poly0.add_chain(Nend,currAtomType)
endchain.attach_parent_node(newNode)
```
The final example is to generate polymers with loops. This is also very easy, using the Loop and LoopSegment objects defined in Polymer.py. An example of generating a ring polymer is below.
```Python
loop0 = Loop(node0) #first you define a Loop object to own every segment of the loop that you add
                   #passing in a Node serves to tell the Loop where it should begin generating on the polymer
                   #segments of these Loops are connected via normal nodes
                   
for i in range(5): #complete this loop with 5 LoopSegments
  loopseg = poly0.add_loop_segment(10,i+1, loop) #give this LoopSegment 10 atoms and a type i+1
                                        #and tell it that it is part of the loop called "loop0"
  loopseg.attach_parent_node(node0) #attach this LoopSegment to node0
  node0 = poly0.add_node(1) #redefine node0 to be a new Node
  loopseg0.attach_child_node(node0) #attach this new Node as a child of the current LoopSegment
loopseg0 = poly0.add_loop_segment(10,i+4, loop0)
loopseg0.attach_parent_node(node0)
chain = poly0.add_chain(5,2)
chain.attach_parent_node(node0)
loopseg0.attach_child_node(poly0.get_node_by_id(0)) #"come back full circle" as they say, make the original Node
                                                #of this Polymer the child of this LoopSegment to close the loop
```

Throughout these examples, I have been using the following functions:
##### Box(boxDims, pointSpacing=0.5)
Constructor of a Box object, this takes the dimensions of the box you will be working with as a list: [Lx, Ly, Lz]  will tell the box to span from -Lx/2:Lx/2, etc...
pointSpacing is the approximate spacing between the atoms that you desire. I have found that 0.5 works well but if you want your atoms to start closer together or farther apart you may change that parameter.

##### Box.add_polymer(polymer, startPos=[0,0,0], minDistance=None)
This function must be called on a box object. It takes a Polymer object to add to the box and a place to start that polymer at, and a minimum distance that the atoms must be spaced from each other, which will default to 0.5* the pointSpacing argument that was provided in the constructor of this Box.

##### Box.add_solvents(numSolvents, solventType, minDistance=0.5)
This function allows you to add discrete individual atoms not bonded to any others into your box. It takes the number of solvents you wish to add and the type of atom you want to assign to them, as well as a minimum distance as Box.add_polymer() has. It generates a 3D lattice of points to place the solvent atoms on, with some random displacement to help with minimization.

##### Box.write_box(outputFileLoc)
This function writes the box to a LAMMPS data file. This should only be done when all polymers/solvents have been added to it.

##### Polymer(initAtomType)
Constructor of a polymer, initializes a Polymer object with a single Node of the type provided by initAtomType.

##### Polymer.add_chain(chainLen, atomType)
Adds a chain to this Polymer of length chainLen and all of the atoms it contains will have type atomType. Returns a Chain.

##### Polymer.add_loop_segment(segLen, atomType, loop)
Adds a loop segment belonging to the Loop "loop". Similar to Polymer.add_chain(), returns a LoopSegment.

##### Polymer.add_node(atomType)
Adds a Node to a polymer. Returns a Node.

##### Polymer.get_chain_by_id(chainID), Polymer.get_node_by_id(nodeID), Polymer.get_loop_by_id(loopsegID)
returns the object with the specified ID, raises NameError if it is not in Polymer.

##### Chain.attach_parent_node(node)
Tells a Chain that it stems from node. There can only be one parent of a chain, and it cannot be reset.

##### Chain.attach_child_node(node)
Tells a Chain that there is another Node attached to the end of it. There can only be one child of a chain, and it cannot be reset.

##### Loop()
Initializes a Loop object, for use when you want to add a closed loop to your Polymer.

##### LoopSegment.attach_parent_node(node), LoopSegment.attach_child_node(node)
See Chain.attach_parent_node, Chain.attach_child_node

There are more functions but they should not be called by the user, they are used internally.
