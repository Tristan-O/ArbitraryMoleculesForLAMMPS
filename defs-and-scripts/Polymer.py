#Still not sure how to handle closed loops with polymer generation, I guess
import math
import numpy as np
from numpy import linalg as LA



class Box:
	#This can keep track of polymers, atomIDs and also be used to generate solvents
	#IDEA: if validating doesn't work out, I can generate a lattice of points that
	#	   the polymer can be generated on and remove available points as atoms go into them
	def __init__(self, boxDims, pointDensity=1):
		if (not isinstance(boxDims, list) and not isinstance(boxDims, tuple)) or any(boxDims[i] <= 0 for i in range(3)):
			print isinstance(boxDims, list)
			print isinstance(boxDims, tuple)
			print (boxDims[i] <= 0 for i in range(3))
			raise TypeError('Object: Box constructor malformed')
		self.boxDims = boxDims
		self.currentAtomID = 1
		self.atomList = []
		self.bondList = []
		self.pointDensity = float(pointDensity) #must be a float because I sometimes scale by this factor
		#create a lattice with some random offset
		#first pad the walls so points aren't created directly on the wall
		boxDims = 0.9*np.array(boxDims)
		x = np.arange(-boxDims[0], boxDims[0], pointDensity)
		y = np.arange(-boxDims[1], boxDims[1], pointDensity)
		z = np.arange(-boxDims[2], boxDims[2], pointDensity)
		self.lattice = np.transpose(np.array([x,y,z], float))
		for i in self.lattice:
			i += 0.1*pointDensity*( np.random.random()-0.5 )

		self.occupiedPoints = np.array([])

	def add_polymer(self, polymer, startPos=[0,0,0], minRadius=0.5):
		'''Adds a list of atoms for this polymer, which have Position, Type and ID information
		startPos is the starting point for this polymer to be built
		boxDims is the three side lengths of the box
		This assumes that the origin of the box is in the center'''
		if not isinstance(polymer, Polymer):
			raise TypeError('Object: Box add_polymer call malformed')
		print 'Adding polymer'
		if any(abs(startPos[i]) > self.boxDims[i] for i in range(3)):
			raise ValueError('Starting position of polymer is outside of box')

		startPos = np.array(startPos)
		polymer.nodeList[0].atom.pos = startPos

		def generate(obj, prevPos1, prevPos2):
			#recursively add polymer stuff
			#assign positions to atoms of obj
			print 'Generating...'
			if obj == None:
				return

			if isinstance(obj, Node):
				atomList = [ obj.atom ]
			else:
				atomList = obj.atomList

			for atom in atomList:
				#I dont think this is a good idea because of the way I have bonds stored to contain two literal atoms...
				#atom = append(Atom(atom, atom.atomType)) #create new atom object with same properties of old one
				#generate a vector is random direction, but with a normal distribution around a certain direction
				spread = 1.
				cov = [[spread,0,0],[0,spread,0],[0,0,spread]] #directions are independent (covariance is 0)
				forwardVec = prevPos1-prevPos2
				#scale forwardVec to keep the desired density 
				forwardVec = forwardVec*( self.pointDensity / ( LA.norm(forwardVec) if LA.norm(forwardVec)!=0.0 else 1.0 ) )

				foundGoodPos = False
				while not foundGoodPos:
					foundGoodPos = True
					randDirection = np.random.multivariate_normal(forwardVec, cov)
					for i in self.occupiedPoints:
						dist = LA.norm(self.occupiedPoints-randDirection)
						if dist < minRadius:
							foundGoodPos = False

				atom.pos = forwardVec + prevPos1
				try:
					self.occupiedPoints = np.vstack([self.occupiedPoints, atom.pos])
				except ValueError:
					self.occupiedPoints = atom.pos
				prevPos1 = prevPos2
				prevPos2 = forwardVec


			#test if there are children remaining in the tree (polymer is like a tree data structure)
			if isinstance(obj, Node):#if there are, call this function again on each of those children
				for chainTuple in obj.attachedChains:
					if chainTuple[1] == 'myChild': # if the node is the parent of this chain
						generate(chainTuple[0], prevPos1, prevPos2)
			elif isinstance(obj, Chain):
				generate(obj.nodeChild, prevPos1, prevPos2)
			
		cov = [[1,0,0],[0,1,0],[0,0,1]] #directions are independent (covariance is 0)
		forwardVec = np.random.multivariate_normal([0.,0.,0.], cov)
		forwardVec = forwardVec / LA.norm(forwardVec)
		generate(polymer.nodeList[0], startPos, forwardVec) #need the -1 here or your starting position will be negated
		

	def add_solvents(self, numSolvents, solventType):
		print 'Adding solvents'

	def validate_atom_positions(self, minRadius=1.):
		'''Make sure no atoms are too close to each other'''
		print 'Validating Atoms'
		invalid = True
		while invalid:
			invalid = False
			for i in range(currentAtomID):
				for j in range(i, currentAtomID):
					#box walls are padded when adding atoms so I don't need to worry about periodic boundaries
					pos1 = self.atomList[i].pos
					pos2 = self.atomList[j].pos
					R = LA.norm(pos1-pos2)
					if R < minRadius:
						invalid = True
						radialVec = (pos1-pos2)*(R-minRadius)/R
						self.atomList[i].pos = pos1 + radialVec
						break


	def write_box(self, outputFileLoc):
		'''Write everything to LAMMPS data file'''
		pass


class Polymer:
	def __init__(self, initAtomType):
		if not isinstance(initAtomType, int):
			raise TypeError('Object: Polymer constructor malformed')
		self.currentChainID = 0 #ID of the most recently added chain
		self.currentNodeID = 0
		self.chainList = []
		self.nodeList = [ Node(self, self.currentNodeID, initAtomType) ]
		self.bondList = []
		self.currentNodeID += 1

	def add_chain(self, chainLen, atomType):
		self.chainList.append( Chain(self, chainLen, atomType, self.currentChainID) )
		self.currentChainID += 1
		return self.chainList[-1]

	def add_node(self, atomType):
		self.nodeList.apppend( Node(self, self.currentNodeID) )
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
			self.parent.bondList.append( Bond(self, self.atomList[i], self.atomList[i+1]) )
	
	def attach_parent_node(self, node):
		#Attach this chain to its parent node
		#Bond here too
		if not isinstance(node, Node):
			raise TypeError('Object: Chain expected an instance of Node')
		node.attach( self, 'myChild' ) #this chain is a child of node
		self.nodeParent = node
		atom1 = node.atom
		atom2 = self.atomList[0]
		self.parent.bondList.append( Bond(self, atom1, atom2) )

	def attach_child_node(self, node):
		#Attach this chain to its parent node
		#Bond here too
		if not isinstance(node, Node):
			raise TypeError('Object: Chain expected an instance of Node')
		node.attach( self, 'myParent' ) #this chain is the parent of node
		self.nodeChild = node
		atom1 = node.atom
		atom2 = self.atomList[-1]
		self.parent.bondList.append( Bond(self, atom1, atom2) )

class Node:
	def __init__(self, parentPolymer, nodeID, atomType):
		if not isinstance(parentPolymer, Polymer) or not isinstance(nodeID, int) or not isinstance(atomType, int):
			raise TypeError('Object: Node constructor malformed')
		self.parent = parentPolymer
		self.attachedChains = []
		self.nodeID = int(nodeID)
		self.atom = Atom(self, int(atomType))

	def attach(self, chain, relationship):
		self.attachedChains.append( (chain, relationship) )


class Bond:
	def __init__(self, parent, atom1, atom2):
		if (not isinstance(parent, Chain) and not isinstance(parent, Node)) or not isinstance(atom1, Atom) or not isinstance(atom2, Atom):
			raise TypeError('Object: Bond constructor malformed')
		self.parent = parent
		self.atom1 = atom1
		self.atom2 = atom2
		#Define a string unique to each bond TYPE, ex. '1-1'
		#Want this string to be in ascending order of atom types, i.e. no '3-2', only '2-3'
		self.bondType = '%i-%i' % tuple(sorted( (atom1.atomType, atom2.atomType) )) #sorted returns a list
		


class Atom:
	def __init__(self, parent, atomType):
		if (not isinstance(parent, Chain) and not isinstance(parent, Node)) or not isinstance(atomType, int):
			raise TypeError('Object: Atom constructor malformed')
		self.parent = parent #parent is chain or node
		self.atomType = atomType
		self.pos = np.array([])




if __name__ == '__main__':
	print 'Example polymer generation'
	box = Box([10,10,10]) #Initialize a 10x10x10 box
	poly0 = Polymer(1) #First node in polymer is an atom of type 1
	node0 = poly0.get_node_by_id(0) #Get the first node object in the polymer 
	chain0 = poly0.add_chain(10, 2) #Generate a chain of len 10 and atom type 2
	chain0.attach_parent_node(node0) #Attach chain0 to node0 at the the 0 end of the chain (end can be 0 or 1)

	box.add_polymer(poly0, startPos=[5,5,5]) #create that polymer one time in the box, starting at the position (5,5,5)
	box.add_solvents(20, 5) #generate solvents in the box
	#box.validate_atom_positions()
	box.write_box('./polymer0.data')
	print box.occupiedPoints