#Still not sure how to handle closed loops with polymer generation, I guess
#For future: consider making atomType class so that mass, radius etc... can be specified for those
import math
import numpy as np
from numpy import linalg as LA



class Box:
	#This can keep track of polymers, atomIDs and also be used to generate solvents
	#IDEA: if validating doesn't work out, I can generate a lattice of points that
	#	   the polymer can be generated on and remove available points as atoms go into them
	def __init__(self, boxDims, pointDensity=1):
		if (not isinstance(boxDims, list) and not isinstance(boxDims, tuple)) or any(boxDims[i] <= 0 for i in range(3)):
			raise TypeError('Object: Box constructor malformed')
		self.boxDims = boxDims[:]
		for i in range(3):
			self.boxDims[i] = self.boxDims[i] / 2.
		self.currentAtomID = 1
		self.atomList = []
		self.polymerAtomList = []#this is a temporary sort of thing, all of the atoms in this get added to the atomList after polymer has been added
		self.bondList = []
		self.pointDensity = float(pointDensity) #must be a float because I sometimes scale by this factor
		self.numPolymers = 0
		self.occupiedPoints = np.array([])

	def add_polymer(self, polymer, startPos=[0,0,0], minRadius=0.5):
		'''Adds a list of atoms for this polymer, which have Position, Type and ID information
		startPos is the starting point for this polymer to be built
		boxDims is the three side lengths of the box
		This assumes that the origin of the box is in the center'''
		self.numPolymers += 1
		if not isinstance(polymer, Polymer):
			raise TypeError('Object: Box add_polymer call malformed')
		print 'Adding polymer'
		if any(abs(startPos[i]) > self.boxDims[i] for i in range(3)):
			raise ValueError('Starting position of polymer is outside of box')

		startPos = np.array(startPos)
		polymer.nodeList[0].atom.pos = startPos
			
		cov = [[1,0,0],[0,1,0],[0,0,1]] #directions are independent (covariance is 0)
		forwardVec = np.random.multivariate_normal([0.,0.,0.], cov)
		forwardVec = forwardVec / LA.norm(forwardVec)
		self.generate(polymer.nodeList[0], startPos, startPos+forwardVec, minRadius)

		#clone atoms and bonds in this box so that changes to them in the polymer object will not affect the ones present in the box:
		atomList = []

		for atom in self.polymerAtomList:
			clone = atom.clone()
			clone.moleculeID = self.numPolymers
			atomList.append(clone) #cloned atoms do not have any bonds
			self.atomList.append(clone)

		bondList = []
		temp = []
		for atom in self.polymerAtomList: #clone bonds over
			for bond in atom.bondList:
				if bond not in bondList:
					for clone in atomList:
						if bond.atom1.cloneID == clone.cloneID:
							atom1 = clone
						elif bond.atom2.cloneID == clone.cloneID:
							atom2 = clone
					temp.append( Bond(self, atom1, atom2) )#adds this bond to both atom1 and atom2, which are clones, and these bonds are now children of this box
					bondList.append(bond) #since bond is known to two atoms by definition, we don't want it to get written twice
		self.bondList = temp

	def add_solvents(self, numSolvents, solventType):
		print 'Adding solvents'
		pass

	def write_box(self, outputFileLoc):
		'''Write everything in this box to LAMMPS data file'''
		print 'Writing box to file'
		with open(str(outputFileLoc), "wb") as f:
	    
		    #LAMMPS DESCRIPTION:
			f.write("This is the LAMMPS data file %s\n" % outputFileLoc)
		    
		    #header - box bounds **************************************************************
			f.write("\n%i atoms\n" % len(self.atomList))
			f.write("%i bonds\n" % len(self.bondList))
			f.write("0 angles\n")
			f.write("0 dihedrals\n")
			f.write("0 impropers\n")
		    
			atomTypes = []
			for i in self.atomList:
				if i.atomType not in atomTypes:
					atomTypes.append( i.atomType )
			atomTypes.sort()

			bondTypes = []
			for i in self.bondList:
				if i.bondType not in bondTypes:
					bondTypes.append( i.bondType )

			f.write("\n%i atom types\n" % len(atomTypes))
			f.write("%i bond types\n" % len(bondTypes))
			f.write("0 angle types\n")
			f.write("0 dihedral types\n")
			f.write("0 improper types\n")
		    
		    #BOX SIZE: *************************************************************************
			Lx = self.boxDims[0]
			Ly = self.boxDims[1]
			Lz = self.boxDims[2]
			f.write("\n%2.3f %2.3f xlo xhi" % (-1.*Lx , 1.*Lx))
			f.write("\n%2.3f %2.3f ylo yhi" % (-1.*Ly , 1.*Ly))
			f.write("\n%2.3f %2.3f zlo zhi\n" % (-1.*Lz , 1.*Lz))
		    
		    # MASSES: **************************************************************************
		    # may need to edit the masses for hydrophobic versus hydrophilic components.
			f.write("\nMasses\n")
			f.write("\n")
			print atomTypes
			for i in atomTypes:
				f.write("%i 1.0\n" % i) # all same mass
		    
		    # ATOMS: ***************************************************************************
			f.write("\nAtoms\n")
			f.write("\n")
			i=0
			for atom in self.atomList:
				i+=1 # Atom ID
				atom.atomID = i
				f.write("%i %i " % (atom.atomID, atom.atomType))
				f.write("%2.6f %2.6f %2.6f " % (atom.pos[0], atom.pos[1], atom.pos[2])) # positiona
				f.write("%i " % atom.moleculeID) # molecular tag ***************######################################################################################
				f.write("1.0 ") # diameter
				f.write("1.0 ") # density
				f.write("0 0 0\n")

		    # BONDS: *******************************************************************************
			f.write("\nBonds\n")
			f.write("\n")

			bondTypes.sort()
			bondDict = {}
			i=0
			for bondType_ in bondTypes:
				i+=1
				bondDict[bondType_] = i
			i=0
			for bond in self.bondList:
				i+=1
				f.write("%i " % i)
				f.write("%i " % bondDict[bond.bondType]) # the bond is hydrophobic-to-philic
				f.write("%i %i\n" % (bond.atom1.atomID, bond.atom2.atomID)) 

			print bondDict

	def addSemiRandPos(self, atom, prevPos1, prevPos2, minRadius):
		spread = 0.5
		foundGoodPos = False
		notGoodCounter = 0
		while not foundGoodPos:
			notGoodCounter += 1
			foundGoodPos = True

			if notGoodCounter >= 500:
				spread += 0.5
				print "Having trouble placing atom... Going to try increasing the allowed random spread to: {}".format(spread)
			# randDirection = np.random.multivariate_normal(forwardVec, cov)
			# randDirection = randDirection * self.pointDensity / LA.norm(randDirection)
			#find forward direction
			forwardVec = prevPos1[:] - prevPos2[:] + spread*( np.random.random(3)-0.5 )
			#scale forward direction to maintain desired density
			forwardVec *= self.pointDensity/LA.norm(forwardVec)
			atom.pos = prevPos1 + forwardVec
			for i in self.occupiedPoints:
				dist = LA.norm(i-atom.pos)
				if dist < minRadius:
					foundGoodPos = False
			for i in range(3):
				if abs(atom.pos[i]) > 0.9*self.boxDims[i]:#if atom is outside of box, change to move in a random direction
					foundGoodPos = False
					prevPos2 = np.random.random(3) + prevPos1
		prevPos2 = prevPos1[:]
		prevPos1 = atom.pos[:]
		try:
			self.occupiedPoints = np.vstack([self.occupiedPoints, atom.pos[:]])
		except ValueError: #Thrown if occupiedpoints is empty
			self.occupiedPoints = atom.pos[:]
		return (prevPos1, prevPos2)

	def generate(self, obj, prevPos1, prevPos2, minRadius):
		#recursively add polymer stuff
		#assign positions to atoms of obj
		if obj == None:
			return

		if isinstance(obj, Node):
			atomList = [ obj.atom ]
		else:
			atomList = obj.atomList

		for atom in atomList:
			prevPos1, prevPos2 = self.addSemiRandPos(atom, prevPos1, prevPos2, minRadius)
			self.polymerAtomList.append(atom)
			# for bond in atom.bondList:
			# 	if bond not in self.bondList:
			# 		self.bondList.append(bond)

		#test if there are children remaining in the tree (polymer is like a tree data structure)
		if isinstance(obj, Node):#if there are, call this function again on each of those children
			for chainTuple in obj.attachedChains:
				if chainTuple[1] == 'myChild': # if the node is the parent of this chain
					self.generate(chainTuple[0], prevPos1[:], prevPos2[:], minRadius)
					#important that these lists are not passed by reference if there's ever two chains attached to one node
		elif isinstance(obj, Chain):
			self.generate(obj.nodeChild, prevPos1[:], prevPos2[:], minRadius)


class Polymer:
	def __init__(self, initAtomType):
		if not isinstance(initAtomType, int):
			raise TypeError('Object: Polymer constructor malformed')
		self.currentChainID = 0 #ID of the most recently added chain
		self.currentNodeID = 0
		self.chainList = []
		self.nodeList = [ Node(self, initAtomType, self.currentNodeID) ]
		self.bondList = []
		self.currentNodeID += 1

	def add_chain(self, chainLen, atomType):
		self.chainList.append( Chain(self, chainLen, atomType, self.currentChainID) )
		self.currentChainID += 1
		return self.chainList[-1]

	def add_node(self, atomType):
		self.nodeList.append( Node(self, atomType, self.currentNodeID) )
		self.currentNodeID += 1
		return self.nodeList[-1]

	def get_chain_by_id(self, chainID):
		for i in self.chainList:
			if i.chainID == chainID:
				return i
		else: #if for loop doesn't get broken out of
			raise NameError('Could not find chain with ID %i' % chainID)

	def get_node_by_id(self, nodeID):
		for i in self.nodeList:
			if i.nodeID == nodeID:
				return i
		else: #if for loop doesn't get broken out of
			raise NameError('Could not find node with ID %i' % nodeID)


class Chain:
	def __init__(self, parentPolymer, chainLen, atomType, chainID):
		if not isinstance(parentPolymer, Polymer) or not isinstance(chainLen, int) or not isinstance(atomType, int) or not isinstance(chainID, int):
			raise TypeError('Object: Chain constructor malformed')
		self.parent = parentPolymer
		self.chainLen = int(chainLen)
		self.atomType = int(atomType)
		self.chainID = int(chainID)
		self.nodeParent = None
		self.nodeChild = None
		self.atomList = []
		for i in range(chainLen):
			#Add atoms
			self.atomList.append( Atom(self, atomType) )
			
		for i in range(chainLen-1):
			#Add chainLen-1 bonds between atoms
			Bond(self.parent, self.atomList[i], self.atomList[i+1])
	
	def attach_parent_node(self, node):
		#Attach this chain to its parent node
		#Bond here too
		if not isinstance(node, Node):
			raise TypeError('Object: Chain expected an instance of Node')
		node.attach( self, 'myChild' ) #this chain is a child of node
		self.nodeParent = node
		atom1 = node.atom
		atom2 = self.atomList[0]
		bond = Bond(self.parent, atom1, atom2)

	def attach_child_node(self, node):
		#Attach this chain to its parent node
		#Bond here too
		if not isinstance(node, Node):
			raise TypeError('Object: Chain expected an instance of Node')
		node.attach( self, 'myParent' ) #this chain is the parent of node
		self.nodeChild = node
		atom1 = node.atom
		atom2 = self.atomList[-1]
		bond = Bond(self.parent, atom1, atom2)


class Node:
	def __init__(self, parentPolymer, atomType, nodeID):
		if not isinstance(parentPolymer, Polymer) or not isinstance(nodeID, int) or not isinstance(atomType, int):
			raise TypeError('Object: Node constructor malformed')
		self.parent = parentPolymer
		self.attachedChains = []
		self.nodeID = int(nodeID)
		self.atom = Atom(self, int(atomType))

	def attach(self, chain, relationship):
		for chainTuple in self.attachedChains:
			if chain in chainTuple:
				raise IndexError('This chain is already in a relationship with this node') 
		self.attachedChains.append( (chain, relationship) )


class Bond:
	def __init__(self, parent, atom1, atom2):
		if (not isinstance(parent, Polymer) and not isinstance(parent, Box)) or not isinstance(atom1, Atom) or not isinstance(atom2, Atom):
			raise TypeError('Object: Bond constructor malformed')
		self.parent = parent
		self.atom1 = atom1
		self.atom2 = atom2
		#Define a string unique to each bond TYPE, ex. '1-1'
		#Want this string to be in ascending order of atom types, i.e. no '3-2', only '2-3'
		self.bondType = '%i-%i' % tuple(sorted( (atom1.atomType, atom2.atomType) )) #sorted returns a list
		atom1.bondList.append(self)
		atom2.bondList.append(self)
		self.parent.bondList.append(self)


class Atom:
	currentCloneID = 0
	def __init__(self, parent, atomType):
		if (not isinstance(parent, Chain) and not isinstance(parent, Node)) or not isinstance(atomType, int):
			raise TypeError('Object: Atom constructor malformed')
		self.parent = parent #parent is chain or node
		self.atomType = atomType
		self.pos = np.array([])
		self.bondList = []
		self.moleculeID = -1

	def clone(self):
		clone = Atom(self.parent, self.atomType)
		clone.pos = self.pos[:]
		clone.bondList = [] #empty on purpose
		# for bond in clone.bondList:
		# 	if bond.atom1 == self:
		# 		cloneBond = Bond(bond.parent, clone, bond.atom2)
		# 	elif bond.atom2 == self:
		# 		cloneBond = Bond(bond.parent, bond.atom1, clone)
		# 	clone.bondList.append(cloneBond)
		self.cloneID = Atom.currentCloneID
		clone.cloneID = Atom.currentCloneID
		Atom.currentCloneID += 1
		#print Atom.currentCloneID
		return clone




if __name__ == '__main__':
	print 'Example polymer generation'
	box = Box([100,100,100], pointDensity=01.) #Initialize a 10x10x10 box
	poly0 = Polymer(1) #First node in polymer is an atom of type 1

	node0 = poly0.get_node_by_id(0) #Get the first node object in the polymer 
	chain0 = poly0.add_chain(10, 2) #Generate a chain of len 10 and atom type 2
	chain0.attach_parent_node(node0) #Attach chain0 to node0

	chain1 = poly0.add_chain(5, 3) #Generate a chain of len 10 and atom type 3
	chain1.attach_parent_node(node0) #Attach chain0 to node0 at the the 0 end of the chain (end can be 0 or 1)

	node1 = poly0.add_node(4)
	chain1.attach_child_node(node1)

	chain3 = poly0.add_chain(10,5)
	chain3.attach_parent_node(node1)

	chain4 = poly0.add_chain(200,5)
	chain4.attach_parent_node(node1)

	box.add_polymer(poly0, startPos=[-5,-5,-5])
	box.add_polymer(poly0, startPos=[ 5, 5, 5]) #create that polymer one time in the box, starting at the position (5,5,5)
	box.add_solvents(20, 5) #generate solvents in the box
	#box.validate_atom_positions()
	box.write_box('./polymer0.data')
	import subprocess as prcs
	simulate_str = prcs.check_output("./lmp_serial -sf omp -pk omp 4 -in polymer.in")