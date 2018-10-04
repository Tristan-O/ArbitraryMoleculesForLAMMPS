#TODO: Add functionality for Angles, Dihedrals, Bond Coeffs
#Ways to improve:
#	- Position generation is done as you add polymers, but it could be done when you try to write it out, 
#	  and then the box would know how many polymers, how much "volume" they will take, could even possibly
#	  do something fancy with splitting the polymers into smaller boxes or something
#	- Look at scipy, they have a whole spatial module that could help with placing atoms. See KDTree
import random
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
		self.boxDims = np.copy(boxDims)
		for i in range(3):
			self.boxDims[i] = self.boxDims[i] / 2.
		self.currentAtomID = 1
		self.atomList = []
		self.polymerAtomList = []#this is a temporary sort of thing, all of the atoms in this get added to the atomList after polymer has been added
		self.bondList = []
		self.pointSpacing = float(pointSpacing) #must be a float because I sometimes scale by this factor
		self.numPolymers = 0
		self.polymerList = []

		numPointsX = int(2.*self.boxDims[0]/self.pointSpacing)
		numPointsY = int(2.*self.boxDims[1]/self.pointSpacing)
		numPointsZ = int(2.*self.boxDims[2]/self.pointSpacing)

		self.lat_points_x = np.linspace(-0.95*self.boxDims[0],0.95*self.boxDims[0],numPointsX)
		self.lat_points_y = np.linspace(-0.95*self.boxDims[1],0.95*self.boxDims[1],numPointsY)
		self.lat_points_z = np.linspace(-0.95*self.boxDims[2],0.95*self.boxDims[2],numPointsZ)
		#print len(self.lat_points_z), numPointsZ, self.pointSpacing, self.lat_points_z

		self.occupiedPoints = []

		self.atomTypes = {}

		self.numSpherePoints = 50

	def define_atom_type(self, atomType, mass=1., diameter=1., density=1.):
		if not isinstance(atomType, int): #in case the user gives a string or float or something
			atomType = int(atomType)
		self.atomTypes[atomType] = {'mass':mass, 'diameter':diameter, 'density':density}
		return	

	def write_box(self, outputFileLoc):
		'''Write everything in this box to LAMMPS data file'''
		print
		print 'If you have trouble writing polymers to a box, try increasing the box size and then using the LAMMPS command fix deform to shrink box to smaller size.'
		#sort polymers by largest average bead size first
		polymersGenerated = []

		for i,polymer in enumerate(self.polymerList):
			# assign positions to every atom
			polymersGenerated.append(polymer)
			self.generate(polymer, polymerNeighbors=polymersGenerated) #fills atomList, bondList and assigns positions to those atoms
			print 'Polymer {} successfully added'.format(i+1)

		print 'Writing box to file'
		print len(self.atomList),"atoms"
		with open(str(outputFileLoc), "wb") as f:
	    
		    #LAMMPS DESCRIPTION:
			f.write("This is the LAMMPS data file %s\n" % outputFileLoc)
		    
		    #header - box bounds **************************************************************
			f.write("\n%i atoms\n" % len(self.atomList))
			f.write("%i bonds\n" % len(self.bondList))
			f.write("0 angles\n")
			f.write("0 dihedrals\n")
			f.write("0 impropers\n")

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

			for i,atom in enumerate(self.atomList):
				atom.atomID = i+1
				f.write("%i %i " % (atom.atomID, atom.atomType))
				f.write("%3.5f %3.5f %3.5f " % (atom.pos[0], atom.pos[1], atom.pos[2])) # positions, using old string format because I dont want to get rid of it

				f.write("{} ".format(atom.moleculeID)) # molecular tag ***************######################################################################################
				f.write("{} ".format(self.atomTypes[atom.atomType]['diameter'])) # diameter
				f.write("{} ".format(self.atomTypes[atom.atomType]['density'])) # density
				f.write("0 0 0\n")

				#print len(atom.bondList)
				# for b in atom.bondList:
				# 	if b not in self.bondList:
				# 		self.bondList.append(b)

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
					#print len(self.bondList)
					#print bond.atom1.parent.parent
					f.write("%i " % i)
					f.write("%i " % bondDict[bond.bondType]) # the bond is hydrophobic-to-philic
					f.write("%i %i\n" % (bond.atom1.atomID, bond.atom2.atomID)) 

				print "Bond types and their corresponding atom type definitions:"
				print bondDict

	@staticmethod
	def dist(r1,r2):
		r=[]
		for i,_ in enumerate(r1):
			r.append(r1[i]-r2[i])
		return LA.norm(r)

	def sphere_of_points(self, center, num_pts, radius):
		#https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/26127012#26127012

		#num_pts = int(4.*math.pi*Rc**2.) # Number of area units on sphere
		indices = np.arange(0, num_pts, dtype=float) + 0.5

		phi = np.arccos(1 - 2*indices/num_pts)
		theta = math.pi * (1. + 5.**0.5) * indices

		points = np.transpose(np.array([(radius+0.01) * np.cos(theta) * np.sin(phi), (radius+0.01) * np.sin(theta) * np.sin(phi), (radius+0.01) * np.cos(phi)]))
		return center + points

	def outsideBox(self, pnt):
		for i in range(3):
			if pnt[i]>0.9*self.boxDims[i] or pnt[i]<-0.9*self.boxDims[i]:
				return True
		else:
			return False

	def generate(self, obj, prevPos1=None, prevradius=0, polymerNeighbors=[], loop=False, startPos=None):
		# recursively build polymer
		# assign positions to atoms of obj

		# do this by choosing N points the correct distance from the center of the previous atom,
		#	then  choose the point with the least overlap between the new atom and already occupied lattice points

		if obj == None:
			return

		# if prevPos1 == None:
		# 	print obj
		if isinstance(obj, Polymer):
			if obj.startPos == None:
				x = random.choice(self.lat_points_x)
				y = random.choice(self.lat_points_y)
				z = random.choice(self.lat_points_z)
				while not [x,y,z] not in self.occupiedPoints and not self.outsideBox([x,y,z]):#choose some random unoccupied spot
					x = random.choice(self.lat_points_x)
					y = random.choice(self.lat_points_y)
					z = random.choice(self.lat_points_z)
				prevPos1 = [x,y,z]
				obj.startPos = [x,y,z]
			else:
				prevPos1 = np.copy(obj.startPos)

			#find neighbor polymers
			actualNeighbors = []
			for possibleNeigh in polymerNeighbors:
				if Box.dist(possibleNeigh.startPos, obj.startPos) < math.sqrt(float(possibleNeigh.length))+math.sqrt(float(obj.length)):
					actualNeighbors.append(possibleNeigh)

			polymerNeighbors = actualNeighbors

			atomList=[]
			obj=obj.nodeList[0]
			obj.atom.pos = np.copy(prevPos1)
			if obj.atom.atomType not in self.atomTypes.keys():
				self.define_atom_type(obj.atom.atomType)


		elif isinstance(obj, Node):
			atomList = [ obj.atom ]
		elif isinstance(obj, Chain):
			atomList = obj.atomList
		else:
			atomList = []
		
		# minOverlapPnt = np.array([0,0,0])
		startPos_2 = np.copy(prevPos1)
		for atom in atomList:
			#find a free spot in lattice
			if atom.atomType not in self.atomTypes.keys():
				self.define_atom_type(atom.atomType)
			radius = self.atomTypes[atom.atomType]['diameter']/2.
			netRadius = radius + prevradius
			prevradius = radius

			prevPos1 = np.array(prevPos1)

			#then find the position with the minimum amount of overlap in occupied points
			minOverlapPnt = None
			minOverlap = None

			sphereOfPoints = self.sphere_of_points(prevPos1, self.numSpherePoints, netRadius)
			if not loop:
				np.random.shuffle(sphereOfPoints)
			else:
				sphereOfPoints = list(sphereOfPoints)
				sphereOfPoints.sort(key = lambda p: ((p[0]-startPos[0])**2+(p[1]-startPos[1])**2+(p[2]-startPos[2])**2))

			for i,spherePoint in enumerate(sphereOfPoints):
				# for each point on the sphere of possible points to choose, 
				# pick the one with the least overlap of occupied points on the lattice

				if self.outsideBox(spherePoint):
					# dont generate atoms outside of box
					# This can cause issues if every point on the sphere is outside of the box, 
					# so there is a check farther down to see if minOverlap was ever set
					continue

				tempCount = 0
				for poly in polymerNeighbors:
					for a in poly.atomList:
						try:
							if Box.dist(a.pos, spherePoint) < self.atomTypes[a.atomType]['diameter']/2.+radius:
								tempCount += 1#./Box.dist(spherePoint,startPos)
						except TypeError: #atom.pos possibly not specifed yet in current polymer
							continue

				# find minimum overlap point
				if minOverlap == None:
					minOverlapPnt = np.copy(spherePoint)
					minOverlap = tempCount
				elif tempCount < minOverlap:
					minOverlap=tempCount
					minOverlapPnt = np.copy(spherePoint)

			if minOverlap == None:
				# failed to find a good point to place the atom
				# try again with a larger radius
				# see the outside box check above
				print 'regenerating...'
				return self.generate(obj, startPos_2, netRadius, polymerNeighbors, loop=loop, startPos=startPos)
				

			#being safe with pass by reference bs, just np.copy()'ing everything
			atom.pos = np.copy(minOverlapPnt)
			minOverlapPnt = list(np.copy(minOverlapPnt))

			#set affected lattice points to be occupied
			sublatticeX = []
			sublatticeY = []
			sublatticeZ = []
			for tempx in self.lat_points_x:
				if tempx>minOverlapPnt[0]-radius and tempx<minOverlapPnt[0]+radius:
					sublatticeX.append(tempx)
			for tempy in self.lat_points_y:
				if tempy>minOverlapPnt[1]-radius and tempy<minOverlapPnt[1]+radius:
					sublatticeY.append(tempy)
			for tempz in self.lat_points_z:
				if tempz>minOverlapPnt[2]-radius and tempz<minOverlapPnt[2]+radius:
					sublatticeZ.append(tempz)
			for x in sublatticeX:
					for y in sublatticeY:
						for z in sublatticeZ:
							if Box.dist([x,y,z],prevPos1) <radius:
								self.occupiedPoints.append([x,y,z])

			prevPos1 = np.copy(np.array(minOverlapPnt)) #last location an atom was placed

		#test if there are children remaining in the tree (polymer is like a tree data structure)
		
		if isinstance(obj, Node):#if there are, call this function again on each of those children
			for child in obj.children:
				if isinstance(child, Chain) or isinstance(child, Node):
					self.generate(child, np.copy(prevPos1[:]), prevradius, polymerNeighbors)
					#important that these lists are not passed by reference if there's ever two chains attached to one node
				elif isinstance(child, Loop): #elif (objChild.loop.nodeParentSuper == obj and not objChild.loop.built) or objChild.loop == loop: #Only begin to generate this loop if you are starting at the loop's beginning or in the process of generating it
					# self.generate_closed_loop(child, np.copy(prevPos1[:]),prevradius, polymerNeighbors)
					loopStart= np.copy(prevPos1)
					for loopChild in child.loopElems:
						prevPos1 = self.generate(loopChild, prevPos1, prevradius, polymerNeighbors, loop=True, startPos=loopStart) #loopStart is passed by reference, changes every iteration
				else:
					raise ValueError('Object: Node has a child that is not of the correct type')
		elif isinstance(obj, Chain):
			self.generate(obj.child, prevPos1, prevradius, polymerNeighbors)


		return prevPos1

	def add_polymer(self, polymer, startPos=None):
		'''Adds a list of atoms for this polymer, which have Position, Type and ID information
		startPos is the starting point for this polymer to be built
		boxDims is the three side lengths of the box
		This assumes that the origin of the box is in the center'''
		#Objects are pass by reference, make clone

		polymerClone = polymer.clone()
		polymerClone.startPos=startPos
		self.polymerList.append(polymerClone)

		polymerClone.atomList = []

		children = [polymerClone.nodeList[0]]
		while children:
			child = children[0]
			del children[0]
			if isinstance(child,Node):
				#print child.atom.bondList
				children += child.children
				polymerClone.atomList.append(child.atom)

			elif isinstance(child, Chain):
				polymerClone.atomList += child.atomList
				if child.child != None:
					children += [child.child]

			elif isinstance(child, Loop):
				children += child.loopElems
				continue

			else:
				continue

		self.atomList += polymerClone.atomList
		tempBondList = [] #for speed, bonds are never going to appear the
		# same between polymers so I will just check for this polymer's bonds
		for atom in polymerClone.atomList:
			for bond in atom.bondList:
				if bond not in tempBondList:
					tempBondList.append(bond)
		self.bondList+=tempBondList

		# assign molecule IDs (not sure what these are used for but whatever)
		for atom in polymerClone.atomList:
			atom.moleculeID = len(self.polymerList)

class Polymer:
	def __init__(self, node0):
		if not isinstance(node0, Node):
			raise TypeError('Object: Polymer constructor malformed')
		self.chainList = []
		self.nodeList = []
		self.loopList = []
		self.bondList = []
		self.length = 0 #1 added at set_parent below
		# self.add_node(node0)
		self.set_parent(node0)


	# def add_chain(self, chain):
	# 	self.chainList.append( chain )

	# def add_node(self, node):
	# 	self.nodeList.append( node )
 
	# def add_loop(self, loop):
	# 	self.loopList.append( loop )

	def set_parent(self, obj):
		if obj == None:
			return
		if obj.parent == self:
			return
		elif isinstance(obj, Chain):
			children = [obj.child]
			self.chainList.append(obj)
			self.length += obj.chainLen
		elif isinstance(obj, Node):
			children = obj.children
			self.nodeList.append(obj)
			self.length += 1
		elif isinstance(obj, Loop):
			children = obj.loopElems
			self.loopList.append(obj)
			#because of the way children are set, don't need to add to length here
		else:
			return

		obj.parent = self
		for child in children:
			self.set_parent(child)

	def get_chain_by_id(self, chainID):
		return self.chainList[chainID]

	def get_node_by_id(self, nodeID):
		return self.nodeList[nodeID]

	def clone(self):
		node0 = self.nodeList[0].clone()
		return Polymer(node0)


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

	def add_bond(self,obj):
		if isinstance(obj,Node):
			Bond(self.atomList[-1], obj.atom)
		elif isinstance(obj,Chain):
			Bond(self.atomList[-1], obj.atomList[0])

	def clone(self):
		myClone = Chain(self.chainLen, self.atomType)
		if self.child != None:
			myClone.add_child(self.child.clone())
		return myClone


class Loop:
	def __init__(self, loopElemList):
		self.numAtoms = 0
		self.loopElems = loopElemList #loop below requires two elements
		self.atomList = [ ]
		for i in range(len(loopElemList)):
			#self.loopElems.append(loopElemList[i])
			if len(loopElemList) > 1 and i != len(loopElemList)-1:
				self.loopElems[i].add_bond(self.loopElems[i+1])

			if isinstance(loopElemList[i], Chain):
				self.numAtoms += loopElemList[i].chainLen
				for j in loopElemList[i].atomList:
					self.atomList.append(j)
			elif isinstance(loopElemList[i], Node):
				self.numAtoms += 1
				self.atomList.append(loopElemList[i].atom)
		self.Positions = []

	def clone(self):
		cloneList = []
		for elem in self.loopElems:
			cloneList.append(elem.clone())
		return Loop(cloneList)


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
			Bond(self.atom, child.atomList[0])
		elif isinstance(child, Loop):
			# child.startNode = node0
			Bond(self.atom, child.atomList[0])
			Bond(self.atom, child.atomList[-1])
		if self.parent != None:
			self.parent.set_parent(child)

	def add_bond(self, obj):
		if isinstance(obj,Node):
			Bond(self.atom, obj.atom)
		elif isinstance(obj,Chain):
			Bond(self.atom, obj.atomList[0])


	def clone(self):
		myClone = Node(self.atom.atomType)
		for child in self.children:
			childClone = child.clone()
			myClone.add_child(childClone)
		return myClone


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
	# currentCloneID = 0
	def __init__(self, parent, atomType):
		if (not isinstance(parent, Box) and not isinstance(parent, Chain) and not isinstance(parent, Node)) or not isinstance(atomType, int):
			raise TypeError('Object: Atom constructor malformed')
		self.parent = parent
		self.atomType = atomType
		self.pos = None
		self.bondList = []
		self.moleculeID = -1

	# def clone(self):
	# This is the old way of copying polymers over
	# 	clone = Atom(self.parent, self.atomType)
	# 	clone.pos = self.pos[:]
	# 	clone.bondList = [] #empty on purpose
	# 	self.cloneID = Atom.currentCloneID
	# 	clone.cloneID = Atom.currentCloneID
	# 	Atom.currentCloneID += 1
	# 	#print Atom.currentCloneID
	# 	return clone




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

	node0 = Node(1)
	poly1 = Polymer(node0)
	chain = Chain(5,2)
	node0.add_child( chain )
	node0.add_child( Chain(5,3))
	node1 = Node(1)
	chain.add_child( node1 )
	node1.add_child( Chain(5,4))

	node0 = Node(1)
	poly0 = Polymer(node0)
	node1 = Node(3)
	node0.add_child( Loop([Chain(5,2),node1,Chain(5,2)]) )


	box = Box([100,100,100]) #Initialize a 100x100x100 box
	box.add_polymer(poly0, startPos=[0,0,0]) #create that polymer one time in the box, starting at the position (-35,-35,-35)
	box.add_polymer(poly1, startPos=[ 35, 35, 35]) 
	box.add_polymer(poly2, startPos=[ 35,-35,-35])
	box.add_polymer(poly3, startPos=[ 35, 35,-35])

	box.write_box('./other-scripts/polymer0.data')
	import subprocess as prcs
	simulate_str = prcs.check_output("./other-scripts/lmp_serial.exe -in ./other-scripts/polymer.in")