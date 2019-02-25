# Create a star molecule (i.e a molecule with one atom at the center and several arms that branch out)
print 'I will now make a star-shaped molecule'
import Molecule as m

mol = m.Molecule()
box = m.Box([1,1,1], debug=True)

center_atom_type = 1
center_atom = mol.add_atom( center_atom_type )

N_arms = 10 # 10 arms on our star
arm_length = 5 # 5 atoms on each arm, not including central node
arm_atom_type = 2
center_bond_type = 1
arm_bond_type = 2

for _ in range(N_arms): # loop over number of arms
	prev_atom = center_atom # keep track of previous atom
	for j in range(arm_length): # loop over the atoms in a single arm
		new_atom = mol.add_atom( arm_atom_type ) # add a new atom along this arm

		# bond the previously added atom and the new atom
		# bond type matters
		if j == 0:
			mol.bond_atoms( center_bond_type, prev_atom, new_atom )
		else:
			mol.bond_atoms( arm_bond_type, prev_atom, new_atom )
		prev_atom = new_atom # set the most recently added atom to be the previously added one
							 # in preparation for the next iteration of the loop

print 'We can check the number of atoms in the molecule is correct to ensure we did the setup correctly'
print 'Should be 1 + N_arms * arm_length = 1 + %i * %i = %i' %(N_arms, arm_length, 1+N_arms*arm_length)
print 'molecule has %i atoms' % len(mol.atomList)

#define atom types
box.define_atom_type( center_atom_type, diameter=5 ) # make the center atom large 
box.define_atom_type( arm_atom_type ) # give this atom type the default parameters

box.add_molecule( mol ) # add this molecule to our box. If we wanted two copies of this
						# molecule in the box we would just call this function again.

box.write_box( 'starMolecule.data' ) # generate the data file corresponding to this setup.