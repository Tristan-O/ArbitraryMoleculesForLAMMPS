# Molecule.py
Object oriented programming approach to creating LAMMPS data files.

This script gives a user the ability to define molecules in an object oriented way and constructs a LAMMPS style input file based on that input.
Using Multidimensional Scaling techniques and some elementary graph theory, any arbitrary  molecular structure is intelligently constructed to fit into a box.

## Specifics of Use
- Uses SMACOF algorithm (courtesy of sklearn) to assign spatial positions to atoms such that they are relatively close to each other if the user dictated that they were bonded.
- Partitions the box and places molecules into partitions using a similar approach to Octree.
- Optionally displays 3D plots that show the partitioning process and the atom placing algorithm.

## Required Dependencies
- numpy
- sklearn
- (optional) matplotlib, mpl_toolkits
  - Without these, you cannot use the `debug=True` option to display the live plots described above, however it will still create the data file.

## How to use
There are two objects you must initialize first:
```python
import Molecule as m
box = m.Box([Lx,Ly,Lz]) #Lx,Ly,Lz are the box dimensions in each cartesian dimension
mol = m.Molecule()
```

From there, you can define the structure of the molecule m:
```python
atomType = 1
bondType = 1
atom1 = mol.add_atom( atomType ) #the argument here is atom type
atom2 = mol.add_atom( atomType )
atom3 = mol.add_atom( atomType )

mol.bond_atoms( bondType, atom1, atom2 ) #bond_atoms does not care which atoms you pass in first and second
mol.bond_atoms( bondType, atom2, atom3 )
mol.bond_atoms( bondType, atom3, atom1 )

angleType = 1
mol.angle_atoms( angleType, atom1, atom2, atom3 ) #this defines an angle. Here order is important, 
                                     # it is the same order that will be written to the input file.
```

Here we have defined a molecule to be a closed loop of three atoms, and it will have an angle as well of angle type 1.

To actually construct the data file now that we have our molecule well-defined:

```python
box.add_molecule( mol ) #this function clones mol, so if I wanted more than one
                        # copy of the same molecule in a box, I would simply call 
                        # this function again on box and pass in mol again.

box.define_atom_type( atomType, diameter=1, mass=1, density=1 ) #you MUST define 
                        # the atom types of every atom in your system this way

outfile = 'polymer0.data'
box.write( outfile )
```

This will create your input file. For more complicated molecules it does not get much harder. For example, lets define a 'ring molecule' of 20 atoms (a chain of atoms bonded to create a circular shape):

```python
mol = m.Molecule()
a_previous = mol.add_atom( 1 )  #first atom added to our molecule
a_1 = a_previous #keep track of our first atom in a_1
for i in range( 20 ):  #iterate 20 times
	a_current = mol.add_atom( 2 )  #add another atoms
	mol.bond_atoms( 1, a_previous,a_current )  #bond that atom to the previous atom
	a_previous = a_current  #define the most recent atom to be now previous
mol.bond_atoms( 1, a_previous, a_1 ) #finally, close the loop on our ring.

#now add an instance of this molecule to our box 8 times:
for _ in range(8):
  box.add_molecule( mol ) #remember Box.add_molecule( mol ) creates a clone of mol before adding it

box.define_atom_type( 1 ) #defined to have default values for diameter, mass and density
box.define_atom_type( 2 )
box.write( 'polymer0.data' )
```

Additionally, there is support for dihedrals and impropers:
```python
Molecule.dihedral_atoms( dihedralType, atom1, atom2, atom3, atom4 )
Molecule.improper_atoms( improperType, atom1, atom2, atom3, atom4 )
```
They work in exactly the same manner as angles (see ```mol.angle_atoms( ... )``` above) except that they call for four atoms instead of three. The order DOES matter for these atoms, unlike simple bonds.

Finally, you can define your own sections (say, Bond Coeffs)

```python
...

header = 'Bond Coeffs'
line1 = [1, 0, 4]
line2 = [2, 0, 5]
box.define_other_section( header,line1 )
box.define_other_section( header,line2 )

...
```

Which will append a section:
```
Bond Coeffs

1 0 4
2 0 5
```
To your file

#### debug
There is one thing I mentioned above that I did not go over in How to use. This script can produce nice 3D plots for debugging if you are uncertain your molecule is begin defined correctly and don't want to go through the generate data file to check, however this vastly slows down the process of building the molecules. If you have matplotlib and are having no issues with Tkinter, you can use ```debug=True``` in the call to initialize the Box object to enable these plots.

### Additional Reference
This was my reference for how a LAMMPS style input file is formatted:
https://lammps.sandia.gov/doc/2001/data_format.html

SMACOF and MDS are somewhat described here:
https://en.wikipedia.org/wiki/Stress_majorization
https://scikit-learn.org/stable/modules/generated/sklearn.manifold.MDS.html#sklearn.manifold.MDS

Octree is here (see Quadtree, the 2D equivalent):
http://gameprogrammingpatterns.com/spatial-partition.html

And I stole the code for live plotting from:
https://stackoverflow.com/questions/5179589/continuous-3d-plotting-i-e-figure-update-using-python-matplotlib
https://stackoverflow.com/questions/44881885/python-draw-3d-cube

### All User Functions
```python
class Box
	__init__(boxDims, debug=False)
    define_atom_type(atomType, mass=1., diameter=1., density=1.)
	define_other_section(header, line)
	add_molecule(m)
	write_box(outputFileLoc)

class Molecule
	__init__()
	add_atom(atomType)
	bond_atoms(bondType, atom1, atom2)
	angle_atoms(angleType, atom1, atom2, atom3)
	dihedral_atoms(dihedralType, atom1, atom2, atom3, atom4)
	improper_atoms(improperType, atom1, atom2, atom3, atom4)
```