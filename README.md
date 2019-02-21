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
box = m.Box([Lx,Ly,Lz])
mol = m.Molecule()
```

From there, you can define the structure of the molecule m:
```python
atomType = 1
atom1 = mol.add_atom( atomType ) #the argument here is atom type
atom2 = mol.add_atom( atomType )
atom3 = mol.add_atom( atomType )

mol.bond_atoms( atom1, atom2 ) #bond_atoms does not care whhich atoms you pass in first and second
mol.bond_atoms( atom2, atom3 ) #bond type is determined by what atom types are bonded
mol.bond_atoms( atom3, atom1 )

angleType = 1
mol.angle_atoms( angleType, atom1, atom2, atom3 ) #this defines an angle. Here order is important, 
                                     # it is the same order that will be written to the input file.
```

Here we have defined a molecule to be a closed loop of three atoms, and it will have an angle as well of angle type 1.

To actually construct the box now that we have our molecule well-defined:

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
for i in range( 20 ):  #iterate 20 times
	a_current = mol.add_atom( 2 )  #add another atoms
	mol.bond_atoms( a_previous,a_current )  #bond that atom to the previous atom
	a_previous = a_current  #define the most recent atom to be now previous
mol.bond_atoms( a_previous, a_1 ) #finally, close the loop on our ring.

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

This was my reference for how a LAMMPS style input file is formatted:
https://lammps.sandia.gov/doc/2001/data_format.html
