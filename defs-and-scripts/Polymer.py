#Still not sure how to handle closed loops with polymer generation, I guess
#For future: consider making atomType class so that mass, radius etc... can be specified for those
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
		self.occupiedPoints = np.array([])

	def add_polymer(self, polymer, startPos=[0,0,0], minDistance=None):
		'''Adds a list of atoms for this polymer, which have Position, Type and ID information
		startPos is the starting point for this polymer to be built
		boxDims is the three side lengths of the box
		This assumes that the origin of the box is in the center'''
		#Make sure polymer is ready to be written
		if minDistance == None:
			minDistance=self.pointSpacing/2.#this used to be an argument to this function and I'm too lazy to replace it with self.pointSpacing everywhere so I'm just redefining it here
		for i in polymer.chainList:
			if i.nodeParent == None:
				raise ValueError('Object: Box cannot accept a polymer with a chain that has no parent node')
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

		#Check that all loops are closed
		for i in polymer.loopList:
			i.loop.built = False#reset these to rebuild a polymer
			i.loop.currentAtomNum = 0
			if not i.loop.closed:
				raise ValueError('Attempted to write polymer containing unclosed loop')

		self.numPolymers += 1
		self.polmerAtomList = []
		if not isinstance(polymer, Polymer):
			raise TypeError('Object: Box add_polymer call malformed')
		print 'Adding polymer'

		startPos = np.array(startPos)
		# polymer.nodeList[0].atom.pos = startPos
			
		forwardVec = np.random.random(3) - 0.5
		forwardVec = forwardVec / LA.norm(forwardVec)
		self.generate(polymer.nodeList[0], startPos, startPos+forwardVec, minDistance)

		#clone atoms and bonds in this box so that changes to them in the polymer object will not affect the ones present in the box:
		atomList = []

		for atom in self.polymerAtomList:
			clone = atom.clone()
			clone.moleculeID = self.numPolymers
			clone.parent = self
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
					temp.append( Bond(atom1, atom2) )#adds this bond to both atom1 and atom2, which are clones, and these bonds are now children of this box
					bondList.append(bond) #since bond is known to two atoms by definition, we don't want it to get written twice
		self.bondList = temp

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
			for objChild in obj.attached:
				if isinstance(objChild, tuple):#this check is not really necessary
					if isinstance(objChild[0], Chain) and objChild[1] == 'myChild': # if the node is the parent of this chain
						self.generate(objChild[0], prevPos1[:], prevPos2[:], minDistance)
						#important that these lists are not passed by reference if there's ever two chains attached to one node
					if isinstance(objChild[0], LoopSegment) and objChild[1] == 'myChild': #elif (objChild.loop.nodeParentSuper == obj and not objChild.loop.built) or objChild.loop == loop: #Only begin to generate this loop if you are starting at the loop's beginning or in the process of generating it
						self.generate_closed_loop(objChild[0], prevPos1[:], minDistance)
				else:
					raise ValueError('Object: Node has a child that is not of the correct type')
		elif isinstance(obj, Chain):
			self.generate(obj.nodeChild, prevPos1[:], prevPos2[:], minDistance)

	def generate_closed_loop(self, loopseg, prevPos, minDistance):
		#choose random direction, this will point to circle center
		#loop.center will be this

		if loopseg.loop.currentAtomNum == 0:#if this is the start of the loop
			loopseg.loop.Position = []
			circleRadius = loopseg.loop.numAtoms*self.pointSpacing/2./math.pi
			randDirection = np.random.random(3)-0.5
			randDirection /= LA.norm(randDirection)

			A = prevPos / LA.norm(prevPos)
			C = np.cross(A,randDirection)#random vector perpendicular to A
			
			B = np.cross(A,C)
			B /= LA.norm(B)

			for i in range(loopseg.loop.numAtoms): 
				theta = 2.*math.pi*(i+1.)/(loopseg.loop.numAtoms)
				Pos = circleRadius*(math.cos(theta)-1)*A + circleRadius*math.sin(theta)*B +prevPos
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
				loopseg.loop.Positions.append(Pos)
				self.occupiedPoints = np.vstack([self.occupiedPoints, Pos[:]])

		for i in range(len(loopseg.atomList)):
			atom = loopseg.atomList[i]
			i += loopseg.loop.currentAtomNum
			self.polymerAtomList.append(atom)
			atom.pos = loopseg.loop.Positions[i][:]
		loopseg.loop.currentAtomNum += len(loopseg.atomList)+1 #+1 because of the node after this
		prevPos1 = loopseg.atomList[-1].pos
		prevPos2 = loopseg.atomList[-2].pos

		if loopseg.nodeChild == loopseg.loop.nodeParentSuper:
		#if this has reached back to the beginning of the loop
			return
		else:
			loopseg.loop.built = True
			self.generate(loopseg.nodeChild, prevPos1[:], prevPos2[:], minDistance, loop=loopseg.loop)



class Polymer:
	def __init__(self, initAtomType):
		if not isinstance(initAtomType, int):
			raise TypeError('Object: Polymer constructor malformed')
		self.currentChainID = 0 #ID of the most recently added chain
		self.currentNodeID = 0
		self.currentLoopID = 0
		self.chainList = []
		self.nodeList = [ Node(initAtomType, self.currentNodeID) ]
		self.loopList = []
		self.bondList = []
		self.currentNodeID += 1

	def add_chain(self, chainLen, atomType):
		self.chainList.append( Chain(chainLen, atomType, self.currentChainID) )
		self.currentChainID += 1
		return self.chainList[-1]

	def add_loop_segment(self, segLen, atomType, loop):
		self.loopList.append( LoopSegment(segLen, atomType, loop, self.currentLoopID) )
		self.currentLoopID += 1
		return self.loopList[-1]

	def add_node(self, atomType):
		self.nodeList.append( Node(atomType, self.currentNodeID) )
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

	def get_loop_by_id(self, loopSegmentID):
		for i in self.loopList:
			if i.loopSegmentID == loopSegmentID:
				return i
		else: #if for loop doesn't get broken out of
			raise NameError('Could not find loopSegment with ID %i' % loopSegmentID)


class Chain:
	def __init__(self, chainLen, atomType, chainID):
		if not isinstance(chainLen, int) or not isinstance(atomType, int) or not isinstance(chainID, int):
			raise TypeError('Object: Chain constructor malformed')
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
			Bond(self.atomList[i], self.atomList[i+1])
	
	def attach_parent_node(self, node):
		#Attach this chain to its parent node
		#Bond here too
		if self.nodeParent != None:
			raise TypeError('Object: A chain cannot detach a parent node')
		if not isinstance(node, Node):
			raise TypeError('Object: Chain expected an instance of Node')
		node.attach_chain( self, 'myChild' ) #this chain is a child of node
		self.nodeParent = node
		atom1 = node.atom
		atom2 = self.atomList[0]
		bond = Bond(atom1, atom2)

	def attach_child_node(self, node):
		#Attach this chain to its parent node
		#Bond here too
		if self.nodeChild != None:
			raise TypeError('Object: A chain cannot detach a child node')
		if not isinstance(node, Node):
			raise TypeError('Object: Chain expected an instance of Node')
		node.attach_chain( self, 'myParent' ) #this chain is the parent of node
		self.nodeChild = node
		atom1 = node.atom
		atom2 = self.atomList[-1]
		bond = Bond(atom1, atom2)


class Loop:
	def __init__(self, node):
		self.numAtoms = 0
		self.nodeParentSuper = node
		self.center = []
		self.closed = False
		self.built=False
		self.currentAtomNum = 0
		self.Positions = []

class LoopSegment: #basically no different than a Chain
	def __init__(self, segLen, atomType, loop, loopsegID):
		if loop.closed == True:
			raise IndexError('Attempted to add LoopSegment to loop that is already closed')
		self.nodeParent = None
		self.nodeChild = None
		self.segmentLen = segLen
		self.loop = loop
		self.loopSegmentID = loopsegID
		self.loop.numAtoms += self.segmentLen + 1
		#+1 from parent node, this number doesn't need to be exact b/c it's only for calculating the spatial radius when generating the loop
		self.atomList = []
		for i in range(segLen): #Add atoms
			self.atomList.append( Atom(self, atomType) )
		for i in range(segLen-1): #Add chainLen-1 bonds between atoms
			Bond(self.atomList[i], self.atomList[i+1])

	def attach_parent_node(self, node): #Bond here too
		if self.nodeParent != None:
			raise TypeError('Object: A LoopSegment cannot detach a parent node')
		if not isinstance(node, Node):
			raise TypeError('Object: LoopSegment expected an instance of Node')
		node.attach_loop_segment( self, 'myChild') #this chain is a child of node
		self.nodeParent = node
		atom1 = node.atom
		atom2 = self.atomList[0]
		bond = Bond(atom1, atom2)

	def attach_child_node(self, node): #Bond here too
		if self.nodeChild != None:
			raise TypeError('Object: A LoopSegment cannot detach a child node')
		if not isinstance(node, Node):
			raise TypeError('Object: LoopSgement expected an instance of Node')
		node.attach_loop_segment( self, 'myParent' ) #this chain is the parent of node
		self.nodeChild = node
		atom1 = node.atom
		atom2 = self.atomList[-1]
		bond = Bond(atom1, atom2)
		if node == self.loop.nodeParentSuper:
			self.loop.closed = True


class Node:
	def __init__(self, atomType, nodeID):
		if not isinstance(nodeID, int) or not isinstance(atomType, int):
			raise TypeError('Object: Node constructor malformed')
		self.attached = []
		self.nodeID = int(nodeID)
		self.atom = Atom(self, int(atomType))

	def attach_chain(self, chain, relationship):
		for child in self.attached:
			if chain in child:
				raise IndexError('This chain is already in a relationship with this node') 
		self.attached.append( (chain, relationship) )

	def attach_loop_segment(self, loopseg, relationship):
		for child in self.attached:
			if loopseg in self.attached[0]:
				return
		self.attached.append( (loopseg, relationship) )


class Bond:
	def __init__(self, atom1, atom2):
		if not isinstance(atom1, Atom) or not isinstance(atom2, Atom):
			raise TypeError('Object: Bond constructor malformed')
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
		if (not isinstance(parent, Box) and not isinstance(parent, Chain) and not isinstance(parent, Node) and not isinstance(parent, LoopSegment)) or not isinstance(atomType, int):
			raise TypeError('Object: Atom constructor malformed')
		self.parent = parent #parent is chain or node
		self.atomType = atomType
		self.pos = 'j'#np.array([])
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
	print 'Make some polymer with an arbitrary number of arms'
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

	chain4 = poly0.add_chain(10,5)
	chain4.attach_parent_node(node1)

	print 'Make a ring polymer'
	poly1 = Polymer(1)

	node0 = poly1.get_node_by_id(0)
	loop = Loop(node0)
	loopseg = poly1.add_loop_segment(10,2, loop)
	loopseg.attach_parent_node(node0)
	loopseg.attach_child_node(node0)

	print 'Make a more complicated ring polymer'
	poly2 = Polymer(1)

	node0 = poly2.get_node_by_id(0)
	loop0 = Loop(node0)
	loopseg0 = poly2.add_loop_segment(10,2, loop0)
	loopseg0.attach_parent_node(node0)

	node1 = poly2.add_node(1)
	loopseg0.attach_child_node(node1)
	loopseg1 = poly2.add_loop_segment(10,2, loop0)
	loopseg1.attach_parent_node(node1)
	loopseg1.attach_child_node(node0)

	node2 = poly2.add_node(1)
	loop1 = Loop(node1)
	loopseg2 = poly2.add_loop_segment(20,3, loop1)
	loopseg2.attach_parent_node(node1)
	loopseg2.attach_child_node(node2)
	loopseg3 = poly2.add_loop_segment(20,3, loop1)
	loopseg3.attach_child_node(node1)
	loopseg3.attach_parent_node(node2)

	print 'Sawblade Polymer'
	poly3 = Polymer(1)
	node0 = poly3.get_node_by_id(0)
	loop = Loop(node0)
	for i in range(5):
		chain = poly3.add_chain(5,2)
		chain.attach_parent_node(node0)
		loopseg0 = poly3.add_loop_segment(10,i+3, loop)
		loopseg0.attach_parent_node(node0)
		node0 = poly3.add_node(1)
		loopseg0.attach_child_node(node0)
	loopseg0 = poly3.add_loop_segment(10,i+4, loop)
	loopseg0.attach_parent_node(node0)
	chain = poly3.add_chain(5,2)
	chain.attach_parent_node(node0)
	loopseg0.attach_child_node(poly3.get_node_by_id(0))

	
	box = Box([100,100,100], pointSpacing=2) #Initialize a 100x100x100 box
	box.add_polymer(poly0, startPos=[-35,-35,-35]) #create that polymer one time in the box, starting at the position (-35,-35,-35)
	box.add_polymer(poly1, startPos=[ 35, 35, 35]) 
	box.add_polymer(poly2, startPos=[ 35,-35,-35])
	box.add_polymer(poly3, startPos=[ 35, 35,-35])
	box.add_solvents(100, 6) #generate solvents in the box

	box.write_box('./polymer0.data')
	import subprocess as prcs
	simulate_str = prcs.check_output("./lmp_serial -sf omp -pk omp 4 -in polymer.in")