
class Polymer:
	def __init__(self, initChainLen, initAtomType):
		self.currentAtomID = 1 #ID of the most recently built atom
		self.currentChainID = 0 #ID of the most recently added chain
		self.currentNodeID = 0
		self.chainList = [ chain(self, initChainLen, initAtomType, self.currentChainID) ]
		self.nodeList = []
		self.currentchainID += 1
		pass

	def build(self, Pos, aveDist, boxDims):
	'''Returns a list of atoms for this polymer, which have Position, Type and ID information
	Pos is the starting point for this polymer to be built
	This assumes that the origin of the box is in the center'''
		pass

	def add_chain(self, chainLen, atomType):
		self.chainList.append( Chain(self, chainLen, atomType, self.currentchainID) )
		self.currentChainID += 1
		return self.chainList[-1]

	def add_node(self, atomType):
		self.nodeList.apppend( Node(self, self.currentNodeID) )
		self.currentNodeID += 1
		return self.nodeList[-1]

	def get_chain(self, chainID):
		for i in self.chainList:
			if i.chainID == chainID:
				return i
		else: #if for loop doesn't get broken out of
			raise NameError('Could not find chain with ID %i' % chainID)

	def get_node(self, nodeID):
		for i in self.nodeList:
			if i.nodeID == nodeID:
				return i
		else: #if for loop doesn't get broken out of
			raise NameError('Could not find node with ID %i' % nodeID)


class Chain:
	def __init__(self, parentPolymer, chainLen, atomType, chainID):
		self.parent = parentPolymer
		self.chainLen = int(chainLen)
		self.atomType = int(atomType)
		self.chainID = int(chainID)
		self.nodeList = []
	
	def attach(self, node, end):
		if self.parent.chainList[0].chainID == self.chainID:
			raise 
		node.attach( self, end )
		self.nodeList.append( (node, end) )


class Node:
	def __init__(self, parentPolymer, nodeID, atomType):
		self.parent = parentPolymer
		self.attachedChains = []
		self.nodeID = int(nodeID)
		self.atomType = int(atomType)

	def attach(self, chain, end):
		self.attachedChains.append( (chain, end) )


class Bond:
	def __init__(self):
		pass


class Atom:
	def __init__(self):
		pass
