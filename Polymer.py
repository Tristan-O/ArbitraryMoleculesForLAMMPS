#TODO: Add functionality for Angles, Dihedrals, Bond Coeffs
#Ways to improve:
#	- Position generation is done as you add polymers, but it could be done when you try to write it out, 
#	  and then the box would know how many polymers, how much "volume" they will take, could even possibly
#	  do something fancy with splitting the polymers into smaller boxes or something
#	- Look at scipy, they have a whole spatial module that could help with placing atoms. See KDTree

import math
import numpy as np
from numpy import linalg as LA

class Box:
	#This can keep track of polymers, atomIDs and also be used to generate solvents
	#IDEA: if validating doesn't work out, I can generate a lattice of points that
	#	   the polymer can be generated on and remove available points as atoms go into them
	def __init__(self, boxDims, pointSpacing=0.5):
		if (not isinstance(boxDims, list) and not isinstance(boxDims, tuple)) or any(boxDims[i] <= 0 for i in range(3)):
			raise TypeError('Object: Box constructor malformed')
		self.boxDims = boxDims[:]
		for i in range(3):
			self.boxDims[i] = self.boxDims[i] / 2.
		self.currentAtomID = 1
		self.atomList = []
		self.polymerAtomList = []#this is a temporary sort of thing, all of the atoms in this get added to the atomList after polymer has been added
		self.bondList = []
		self.pointSpacing = float(pointSpacing) #must be a float because I sometimes scale by this factor
		self.numPolymers = 0

		numPointsX = 2.*self.boxDims[0]/self.pointSpacing
		numPointsY = 2.*self.boxDims[1]/self.pointSpacing
		numPointsZ = 2.*self.boxDims[2]/self.pointSpacing
		self.latticePoints = np.zeros((numPointsX,numPointsY,numPointsZ,2)) # make 4D array
			# 3Dims to represent points
			# 1Dim to represent spatial point-tuple and atomType of atom occupying that point
			# This will keep track of what points atoms have occupied (including points occupied by an atom's volume)
		for i in range(numPointsx):
			for j in range(numPointsY):
				for k in range(numPointsZ):
					self.latticePoints[i,j,k,0] = (numPointsX*self.pointSpacing-self.boxDims[0],
													numPointsY*self.pointSpacing-self.boxDims[1],
													numPointsZ*self.pointSpacing-self.boxDims[2])
					self.latticePoints[i,j,k,1] = -1 # -1 is not a valid atom type in LAMMPS, this represents an open points (not occupied by any atom)

		self.atomTypes = {}

	def define_atom_type(self, atomType, mass=1., diameter=1., density=1.):
		if not isinstance(atomType, int): #in case the user gives a string or float or something
			atomType = int(atomType)
		self.atomTypes[atomType] = {'mass':mass, 'diameter':diameter, 'density':density}
		return

	def add_solvents(self, numSolvents, solventType, minDistance=0.5):
		MAX_NUM_FOR_CHECKING = 100000
		print 'Adding solvents'
		if numSolvents < 1:
			raise ValueError('Cannot generate less than one solvent')
		elif numSolvents > MAX_NUM_FOR_CHECKING:
			print "WARNING: Too many solvents to do rigorous checking of placement (computational expense)"

		atomList = []
		for solventNum in range(numSolvents):
			atomList.append( Atom(self, solventType) )

		#Assign positions to solvents
		LatScale = 0.9
		Lx = self.boxDims[0]
		Ly = self.boxDims[1]
		Lz = self.boxDims[2]
		rx = np.linspace(-LatScale*Lx, LatScale*Lx, math.ceil(numSolvents**(1./3.)), dtype=float)
		ry = np.linspace(-LatScale*Ly, LatScale*Ly, math.ceil(numSolvents**(1./3.)), dtype=float)
		rz = np.linspace(-LatScale*Lz, LatScale*Lz, math.ceil(numSolvents**(1./3.)), dtype=float)

		i=0
		for x in rx:
			for y in ry:
				for z in rz:
					Pos = np.array([x,y,z], float)
					# add a random offset to help initial minimization.
					spread = 0.5
					foundGoodPos = False
					notGoodCounter = 0
					while not foundGoodPos:
						notGoodCounter += 1
						foundGoodPos = True

						if notGoodCounter >= 5000:
							spread += 0.5
							print "Having trouble placing atom... Going to try increasing the allowed random spread to: {}".format(spread)
							if spread > 5*self.pointSpacing:#don't want atoms getting too far apart
								spread = 0.
			
						Pos += spread*( np.random.random(3)-0.5 )
						if numSolvents < MAX_NUM_FOR_CHECKING:
							for j in self.occupiedPoints:
								dist = LA.norm(j-Pos)
								if dist < minDistance:
									foundGoodPos = False
									break
						for j in range(3):
							if abs(Pos[j]) > self.boxDims[j]:#if atom is outside of box, change to move in a random direction
								foundGoodPos = False
							# elif spread > self.boxDims[j]:
							# 	print 'Something went wrong with placing a solvent... I will try again'
							# 	spread = 0
					atomList[i].pos = Pos[:]
					self.atomList.append(atomList[i])
					i += 1
					# if done placing atoms, return.
					if i >= numSolvents:
						break
				if i >= numSolvents:
					break
			if i >= numSolvents:
				break	

		for atom in atomList:
			try:
				self.occupiedPoints = np.vstack([self.occupiedPoints, atom.pos[:]])
			except ValueError: #Thrown if occupiedPoints is empty
				self.occupiedPoints = Pos[:]		

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
		    
			for i in self.atomList:
				if i.atomType not in self.atomTypes.keys():
					self.define_atom_type(i.atomType) #initialize default values for this atomType

			bondTypes = []
			for i in self.bondList:
				if i.bondType not in bondTypes:
					bondTypes.append( i.bondType )

			f.write("\n{} atom types\n".format(len(self.atomTypes)))
			f.write("{} bond types\n".format(len(bondTypes)))
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
			for i in sorted(self.atomTypes.keys()):
				f.write("{} {}\n".format(i, self.atomTypes[i]['mass']))
		    
		    # ATOMS: ***************************************************************************
			f.write("\nAtoms\n")
			f.write("\n")
			i=0
			for atom in self.atomList:
				i+=1 # Atom ID
				atom.atomID = i
				f.write("%i %i " % (atom.atomID, atom.atomType))
				f.write("%2.6f %2.6f %2.6f " % (atom.pos[0], atom.pos[1], atom.pos[2])) # positions, using old string format because I dont want to get rid of it
				f.write("{} ".format(atom.moleculeID)) # molecular tag ***************######################################################################################
				f.write("{} ".format(self.atomTypes[atom.atomType]['diameter'])) # diameter
				f.write("{} ".format(self.atomTypes[atom.atomType]['density'])) # density
				f.write("0 0 0\n")

		    # BONDS: *******************************************************************************
			if len(self.bondList)>0:
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

				print "Bond types and their corresponding atom type definitions:"
				print bondDict

	def add_semi_rand_pos(self, atom, prevPos1, prevPos2, minDistance):
		spread = 0.5
		foundGoodPos = False
		notGoodCounter = 0
		while not foundGoodPos:
			notGoodCounter += 1
			foundGoodPos = True

			if notGoodCounter >= spread*5000:
				spread += 0.5
				print "Having trouble placing atom... Going to try increasing the allowed random spread to: ",spread
			# randDirection = np.random.multivariate_normal(forwardVec, cov)
			# randDirection = randDirection * self.pointSpacing / LA.norm(randDirection)
			#find forward direction
			forwardVec = prevPos1[:] - prevPos2[:] + spread*( np.random.random(3)-0.5 )
			#scale forward direction to maintain desired density
			forwardVec *= self.pointSpacing/LA.norm(forwardVec)
			atom.pos = prevPos1 + forwardVec
			for i in self.occupiedPoints:
				dist = LA.norm(i-atom.pos)
				if dist < minDistance:
					foundGoodPos = False
			for i in range(3):
				if abs(atom.pos[i]) > 0.9*self.boxDims[i]:#if atom is outside of box, change to move in a random direction
					foundGoodPos = False
					prevPos2 = np.random.random(3) + prevPos1
				if spread > self.boxDims[i]:
					print 'Something went wrong with placing... I will try again'
					spread = 0

		prevPos2 = prevPos1[:]
		prevPos1 = atom.pos[:]
		try:
			self.occupiedPoints = np.vstack([self.occupiedPoints, atom.pos[:]])
		except ValueError: #Thrown if occupiedpoints is empty
			self.occupiedPoints = atom.pos[:]
		return (prevPos1, prevPos2)

	def generate(self, obj, prevPos1, prevPos2, minDistance, loop=None):
		#recursively add polymer stuff
		#assign positions to atoms of obj
		if obj == None:
			return

		if isinstance(obj, Node):
			atomList = [ obj.atom ]
		else:
			atomList = obj.atomList

		for atom in atomList:
			prevPos1, prevPos2 = self.add_semi_rand_pos(atom, prevPos1, prevPos2, minDistance)
			self.polymerAtomList.append(atom)

		#test if there are children remaining in the tree (polymer is like a tree data structure)
		if isinstance(obj, Node):#if there are, call this function again on each of those children
			for child in obj.children:
				if isinstance(child, Chain) or isinstance(child, Node):
					self.generate(child, prevPos1[:], prevPos2[:], minDistance)
					#important that these lists are not passed by reference if there's ever two chains attached to one node
				elif isinstance(child, Loop): #elif (objChild.loop.nodeParentSuper == obj and not objChild.loop.built) or objChild.loop == loop: #Only begin to generate this loop if you are starting at the loop's beginning or in the process of generating it
					self.generate_closed_loop(child, prevPos1[:], minDistance)
				else:
					raise ValueError('Object: Node has a child that is not of the correct type')
		elif isinstance(obj, Chain):
			self.generate(obj.child, prevPos1[:], prevPos2[:], minDistance)

	def generate_closed_loop(self, loop, prevPos, minDistance):
		#choose random direction, this will point to circle center
		#loop.center will be this
		Positions = []
		circleRadius = loop.numAtoms*self.pointSpacing/2./math.pi
		randDirection = np.random.random(3)-0.5
		randDirection /= LA.norm(randDirection)

		A = prevPos / LA.norm(prevPos)
		C = np.cross(A,randDirection)#random vector perpendicular to A
		
		B = np.cross(A,C)
		B /= LA.norm(B)

		for i in range(loop.numAtoms): 
			theta = 2.*math.pi*(i+1.)/(loop.numAtoms+1.)
			Pos = circleRadius*(math.cos(theta)-1)*A + circleRadius*math.sin(theta)*B + prevPos
			foundGoodPos = False
			while not foundGoodPos:
				foundGoodPos = True
				for j in self.occupiedPoints:
					dist = LA.norm(j-Pos)
					if dist < minDistance:
						foundGoodPos = False
						Pos += 2.*(np.random.random(3)-0.5)
				for j in range(3):
					if Pos[j] > 0.9*self.boxDims[j]:#if atom is outside of box, change to move in a random direction
						foundGoodPos = False
						Pos[j]= 0.9*self.boxDims[j] - np.random.random()
					elif Pos[j] < -0.9*self.boxDims[j]:#if atom is outside of box, change to move in a random direction
						foundGoodPos = False
						Pos[j]= -0.9*self.boxDims[j] + np.random.random()
			loop.atomList[i].pos = Pos
			self.occupiedPoints = np.vstack([self.occupiedPoints, Pos])
			self.polymerAtomList.append(loop.atomList[i])

		for obj in loop.loopElems:
			if isinstance(obj, Chain): #technically this check is unneccessary, since a Chain can only have one child and that child by continuation is in the Loop
				children = [obj.child]
				pos1 = obj.atomList[-1].pos[:] #atom.pos was set earlier in this function
			elif isinstance(obj, Node):
				children = obj.children
				pos1 = obj.atom.pos[:]
			for child in children:
				if child == loop:
					continue
				if child not in loop.loopElems and child != None:
					pos2 = prevPos - circleRadius*A #center of the circle
					self.generate(child, pos1[:], pos2[:], minDistance)#child is told to generate away from center of loop circle, this might be a problme though if the circle clips the edge of the box...

	def add_polymer(self, polymer, startPos=[0,0,0], minDistance=None):
		'''Adds a list of atoms for this polymer, which have Position, Type and ID information
		startPos is the starting point for this polymer to be built
		boxDims is the three side lengths of the box
		This assumes that the origin of the box is in the center'''
		#Make sure polymer is ready to be written
		if not isinstance(polymer, Polymer):
			raise TypeError('Object: Box add_polymer call malformed')
		if minDistance == None:
			minDistance=self.pointSpacing/2.#this used to be an argument to this function and I'm too lazy to replace it with self.pointSpacing everywhere so I'm just redefining it here
		for i in range(3):
			if abs(startPos[i]) >= 0.9*self.boxDims[i]:
				s='xyz'
				raise ValueError('Object: Box cannot add polymer at '+s[i]+'='+str(startPos[i])+' when the box dimensions are +/-'+str(self.boxDims[i]))
		foundGoodPos = False
		while foundGoodPos:
			foundGoodPos = True
			for j in self.occupiedPoints:
				dist = LA.norm(j-startPos)
				if dist < minDistance:
					foundGoodPos = False
					startPos += (np.random.random(3)-0.5)

		self.numPolymers += 1
		self.polymerAtomList = []
		print 'Adding polymer'

		startPos = np.array(startPos)
			
		forwardVec = np.random.random(3) - 0.5
		forwardVec = forwardVec / LA.norm(forwardVec)
		self.generate(polymer.nodeList[0], startPos, startPos+forwardVec, minDistance)

		#clone atoms and bonds in this box so that changes to them in the polymer object will not affect the ones present in the box:
		atomList = []
		temp = []

		for atom in self.polymerAtomList:
			clone = atom.clone()
			clone.moleculeID = self.numPolymers
			clone.parent = self
			atomList.append(clone)
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
					self.bondList.append( Bond(atom1, atom2) )#adds this bond to both atom1 and atom2, which are clones, and these bonds are now children of this box
					bondList.append(bond) #since bond is known to two atoms by definition, we don't want it to get written twice
		#self.bondList = temp


class Polymer:
	def __init__(self, node0):
		if not isinstance(node0, Node):
			raise TypeError('Object: Polymer constructor malformed')
		self.chainList = []
		self.nodeList = []
		self.add_node(node0)
		self.set_parent(node0)
		self.loopList = []
		self.bondList = []

	def add_chain(self, chain):
		self.chainList.append( chain )

	def add_node(self, node):
		self.nodeList.append( node )
 
	def add_loop(self, loop):
		self.loopList.append( loop )

	def set_parent(self, obj):
		if obj.parent == self:
			return
		elif isinstance(obj, Chain):
			children = [obj.child]
			self.chainList.append(obj)
		elif isinstance(obj, Node):
			children = obj.children
			self.nodeList.append(obj)
		elif isinstance(obj, Loop):
			children = obj.loopElems
			self.loopList.append(obj)
		else:
			return

		obj.parent = self
		for child in children:
			self.set_parent(child)

	def get_chain_by_id(self, chainID):
		return self.chainList[chainID]

	def get_node_by_id(self, nodeID):
		return self.nodeList[nodeID]


class Chain:
	def __init__(self, chainLen, atomType):
		if not isinstance(chainLen, int) or not isinstance(atomType, int):
			raise TypeError('Object: Chain constructor malformed')
		self.chainLen = int(chainLen)
		self.atomType = int(atomType)
		self.chainID = None
		self.child = None
		self.atomList = []
		self.parent = None
		for i in range(chainLen):
			#Add atoms
			self.atomList.append( Atom(self, atomType) )
			
		for i in range(chainLen-1):
			#Add chainLen-1 bonds between atoms
			Bond( self.atomList[i], self.atomList[i+1] )
	
	def add_child(self, child):
		#Bond here too
		if not isinstance(child, Node):
			raise TypeError('Object: Chain expected an instance of Node')
		if self.child != None:
			raise IndexError('Cannot reset child of Chain')
		self.child = child
		atom = self.atomList[-1]
		Bond(child.atom, atom)
		if self.parent != None:
			self.parent.set_parent(child)


class Loop:
	def __init__(self, loopElemList):
		self.numAtoms = 0
		self.loopElems = loopElemList #loop below requires two elements
		self.atomList = [ ]
		for i in range(len(loopElemList)):
			#self.loopElems.append(loopElemList[i])
			if len(loopElemList) > 1 and i != len(loopElemList) -1:
				self.loopElems[i].add_child(self.loopElems[i+1])

			if isinstance(loopElemList[i], Chain):
				self.numAtoms += loopElemList[i].chainLen
				for j in loopElemList[i].atomList:
					self.atomList.append(j)
			elif isinstance(loopElemList[i], Node):
				self.numAtoms += 1
				self.atomList.append(loopElemList[i].atom)
		self.Positions = []


class Node:
	def __init__(self, atomType):
		if not isinstance(atomType, int):
			raise TypeError('Object: Node constructor malformed')
		self.children = []
		self.atom = Atom(self, int(atomType))
		self.parent = None

	def add_child(self, child):
		if not (isinstance(child, Chain) or isinstance(child, Node) or isinstance(child, Loop)):
			raise ValueError('Object: Node can only take an instance of Chain, Node or Loop as a child')
		if child in self.children:
			raise IndexError('This object is already a child of this node')
		self.children.append(child)
		child.parent = self.parent
		if isinstance(child,Node):
			Bond(self.atom, child.atom)
		elif isinstance(child,Chain):
			Bond(self.atom, child.atomList[-1])
		elif isinstance(child, Loop):
			child.startNode = node0
			Bond(self.atom, child.atomList[0])
			Bond(self.atom, child.atomList[-1])
		if self.parent != None:
			self.parent.set_parent(child)


class Bond:
	def __init__(self, atom1, atom2):
		if not isinstance(atom1, Atom) or not isinstance(atom2, Atom):
			raise TypeError('Object: Bond constructor malformed')
		elif atom1 == atom2:
			raise ValueError('Object: Bond Atoms cannot be the same')

		self.atom1 = atom1
		self.atom2 = atom2
		#Define a string unique to each bond TYPE, ex. '1-1'
		#Want this string to be in ascending order of atom types, i.e. no '3-2', only '2-3'
		self.bondType = '%i-%i' % tuple(sorted( (atom1.atomType, atom2.atomType) )) #sorted returns a list
		atom1.bondList.append(self)
		atom2.bondList.append(self)

class Atom:
	currentCloneID = 0
	def __init__(self, parent, atomType):
		if (not isinstance(parent, Box) and not isinstance(parent, Chain) and not isinstance(parent, Node)) or not isinstance(atomType, int):
			raise TypeError('Object: Atom constructor malformed')
		self.parent = parent
		self.atomType = atomType
		self.pos = None#np.array([])
		self.bondList = []
		self.moleculeID = -1

	def clone(self):
		clone = Atom(self.parent, self.atomType)
		clone.pos = self.pos[:]
		clone.bondList = [] #empty on purpose
		self.cloneID = Atom.currentCloneID
		clone.cloneID = Atom.currentCloneID
		Atom.currentCloneID += 1
		#print Atom.currentCloneID
		return clone




if __name__ == '__main__':
	print 'Example polymer generation'
	print 'Make some arbitrary simple polymer'

	node0 = Node(1) #Get the first node object in the polymer 
	poly0 = Polymer(node0) #First node in polymer is an atom of type 1

	NArms = 4
	N = 3*[4]
	currAtomType = 1

	endChain = Chain(N[0]-1, currAtomType)
	node0.add_child(endChain)
	newNode = Node(currAtomType)
	endChain.add_child(newNode)

	for n in range(NArms-1):
		spineChain = Chain(N[1],currAtomType)
		armChain = Chain(N[2], currAtomType+1)
		newNode.add_child(spineChain)
		newNode.add_child(armChain)

		newNode = Node(currAtomType)
		spineChain.add_child(newNode)
		
	armChain = Chain(N[2], currAtomType+1)
	newNode.add_child(armChain)

	endChain = Chain(N[0], currAtomType)
	newNode.add_child(endChain)

	print 'Make a ring polymer'
	node0 = Node(1)
	poly1 = Polymer(node0)
	chain0 = Chain(10,2)

	loop = Loop( [chain0] )
	node0.add_child(loop)

	print 'Make a more complicated ring polymer'
	node0 = Node(1)
	poly2 = Polymer(node0)

	chain0 = Chain(10,2)
	node1 = Node(1)
	chain1 = Chain(10,2)
	node2 = Node(1)
	chain2 = Chain(10,2)
	loop = Loop([chain0,node1,chain1,node2,chain2])

	node0.add_child(loop)

	print 'Sawblade Polymer'
	node0 = Node(1)
	poly3 = Polymer(node0)

	loopList = []
	for i in range(5):
		loopList.append( Chain(10,i+1) )
		node = Node(1)
		loopList.append( node )
		node.add_child( Chain(5,3) )

	loopList.append( Chain(10,2) )
	node0.add_child( Chain(5,3) )
	node0.add_child( Loop(loopList) )

	box = Box([100,100,100], pointSpacing=2) #Initialize a 100x100x100 box
	box.add_polymer(poly0, startPos=[-35,0,0]) #create that polymer one time in the box, starting at the position (-35,-35,-35)
	#box.add_polymer(poly1, startPos=[ 35, 35, 35]) 
	box.add_polymer(poly2, startPos=[ 35,-35,-35])
	box.add_polymer(poly3, startPos=[ 35, 35,-35])
	#box.add_solvents(100, 3) #generate solvents in the box

	box.write_box('./polymer0.data')
	import subprocess as prcs
	simulate_str = prcs.check_output("./other-scripts/lmp_serial -sf omp -pk omp 4 -in ./other-scripts/polymer.in")