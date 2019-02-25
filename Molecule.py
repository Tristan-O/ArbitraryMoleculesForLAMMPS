import numpy as np
from numpy import linalg as LA
from sklearn import manifold

try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
	from matplotlib import cm
	from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
	import matplotlib
except Exception as e:
	print 'Could not import matplotlib (using Tkinter). Do not use debug argument in \'Box\' Initialization or this will fail'
	print e
	print 'Continuing...'

class plot3dClass( object ):

	def __init__( self, dims ):
		# self.plot = self.ax.scatter( 
		#     [0],[0],[0], 
		#     cmap=cm.jet, linewidth=0, antialiased=False )
		# # plt.draw() maybe you want to see this frame?
		matplotlib.interactive(True)
		self.points = []
		self.cube_drawn = False
		self.points_drawn = False
		self.dims = dims

	def draw_now( self, point=None):
		if self.points_drawn == False:
			self.fig = plt.figure()
			self.ax = self.fig.add_subplot( 111, projection='3d' )
			self.ax.set_xlim3d( -self.dims[0], self.dims[0] )
			self.ax.set_ylim3d( -self.dims[1], self.dims[1] )
			self.ax.set_zlim3d( -self.dims[2], self.dims[2] )
			self.points_drawn = True
		if not (point is None):
			self.points.append(point)
		# try:
			# self.plot.remove()
		# except AttributeError:
			# pass
		self.plot = self.ax.scatter( 
			np.transpose(self.points)[0], np.transpose(self.points)[1], np.transpose(self.points)[2],
			c=[0,0,0], linewidth=0, antialiased=False )
		plt.draw()  
		# self.fig.canvas.flush_events() # redraw the canvas
		plt.pause(.0001)

	def plot_cube(self, cube_definition):
		if self.cube_drawn == False:
			self.fig_cube = plt.figure()
			self.ax_cube = self.fig_cube.add_subplot(111, projection='3d')
			self.cube_drawn = True
		cube_definition_array = [
			np.array(list(item))
			for item in cube_definition
		]

		points = []
		points += cube_definition_array
		vectors = [
		    cube_definition_array[1] - cube_definition_array[0],
		    cube_definition_array[2] - cube_definition_array[0],
		    cube_definition_array[3] - cube_definition_array[0]
		]

		points += [cube_definition_array[0] + vectors[0] + vectors[1]]
		points += [cube_definition_array[0] + vectors[0] + vectors[2]]
		points += [cube_definition_array[0] + vectors[1] + vectors[2]]
		points += [cube_definition_array[0] + vectors[0] + vectors[1] + vectors[2]]

		points = np.array(points)

		edges = [
		    [points[0], points[3], points[5], points[1]],
		    [points[1], points[5], points[7], points[4]],
		    [points[4], points[2], points[6], points[7]],
		    [points[2], points[6], points[3], points[0]],
		    [points[0], points[2], points[4], points[1]],
		    [points[3], points[6], points[7], points[5]]
		]

		faces = Poly3DCollection(edges, linewidths=1, edgecolors='k')
		faces.set_facecolor((0,0,1,0.1))

		self.ax_cube.add_collection3d(faces)

		# Plot the points themselves to force the scaling of the axes
		self.ax_cube.scatter(points[:,0], points[:,1], points[:,2], s=5)

		self.ax_cube.set_aspect('equal')


class Box:
	#This can keep track of polymers, atomIDs, etc.
	def __init__(self, boxDims, debug = False):
		# self.boxDims = np.copy(boxDims)
		# for i in range(3):
		self.boxDims = np.array(boxDims) / 2.

		self.moleculeList = []
		self.otherSections = {}
		self.atomTypes = {}
		# self.bondTypes = {}
		# self.angleTypes = {}
		# self.dihedralTypes = {}
		# self.improperTypes = {}

		self.debug = debug
		if self.debug:
			self.interactivePlotter = plot3dClass(np.copy(self.boxDims))

	def define_atom_type(self, atomType, mass=1., diameter=1., density=1.):
		atomType = int(atomType)
		self.atomTypes[atomType] = {'mass':mass, 'diameter':diameter, 'density':density}
		return	

	# def define_bond_type(self, bondType, params):
	# 	bondType = int(bondType)
	# 	self.bondTypes[bondType] = params[:]
	# 	return

	# def define_angle_type(self, angleType, params):
	# 	angleType = int(angleType)
	# 	self.angleTypes[angleType] = params[:]
	# 	return

	# def define_dihedral_type(self, dihedralType, params):
	# 	dihedralType = int(dihedralType)
	# 	self.dihedralTypes[dihedralType] = params[:]
	# 	return	
		
	# def define_bond_type(self, improperType, params):
	# 	improperType = int(improperType)
	# 	self.improperTypes[improperType] = params[:]
	# 	return	

	def define_other_section(self, header, line):
		''' '''
		if header not in self.otherSections.keys():
			self.otherSections[header] = []
		self.otherSections[header].append(line)
		return

	def add_molecule(self, m):
		'''Add a molecule to a box. This clones the molecule object, so one molecule object can be added several times.'''
		#Objects are pass by reference, make clone
		mClone = m.clone()
		self.moleculeList.append(mClone)

	def write_box(self, outputFileLoc):
		'''Write everything in this box to LAMMPS data file'''
		# assign positions in this box to every atom
		self._generate()

		#data setup
		atomList = []
		bondList = []
		angleList = []
		dihedralList = []
		improperList = []

		bondTypes = []
		angleTypes = []
		dihedralTypes = []
		improperTypes = []

		for i_m,m in enumerate(self.moleculeList):

			for atom in m.atomList:
				atom.moleculeID = i_m+1

			atomList += m.atomList
			bondList += m.bondList
			angleList += m.angleList
			dihedralList += m.dihedralList
			improperList += m.improperList

			for b in m.bondList:
				if b.bondType not in bondTypes:
					bondTypes.append( b.bondType )

			for ang in m.angleList:
				if ang.angleType not in angleTypes:
					angleTypes.append( ang.angleType )

			for di in m.dihedralList:
				if di.dihedralType not in dihedralTypes:
					dihedralTypes.append( di.dihedralType )

			for imp in m.improperList:
				if imp.improperType not in improperTypes:
					improperTypes.append( imp.improperType )

		with open(str(outputFileLoc), "wb") as f:
	    
		    #LAMMPS DESCRIPTION:
			f.write("This is the LAMMPS data file %s\n" % outputFileLoc)
		    
		    # HEADER *******************************************************************************
			f.write("\n{} atoms\n".format(len(atomList)))
			f.write("{} bonds\n".format(len(bondList)))
			f.write("{} angles\n".format(len(angleList)))
			f.write("{} dihedrals\n".format(len(dihedralList)))
			f.write("{} impropers\n".format(len(improperList)))		

			f.write("\n{} atom types\n".format(len(self.atomTypes)))
			f.write("{} bond types\n".format(len(bondTypes)))
			f.write("{} angle types\n".format(len(angleTypes)))
			f.write("{} dihedral types\n".format(len(dihedralTypes)))
			f.write("{} improper types\n".format(len(improperTypes)))
		    
		    # BOX SIZE: ****************************************************************************
			Lx = self.boxDims[0]
			Ly = self.boxDims[1]
			Lz = self.boxDims[2]
			f.write("\n%2.3f %2.3f xlo xhi"   % (-1.*Lx , 1.*Lx))
			f.write("\n%2.3f %2.3f ylo yhi"   % (-1.*Ly , 1.*Ly))
			f.write("\n%2.3f %2.3f zlo zhi\n" % (-1.*Lz , 1.*Lz))
		    
		    # MASSES: ******************************************************************************
			f.write("\nMasses\n")
			f.write("\n")
			for i in sorted(self.atomTypes.keys()):
				f.write("{} {}\n".format(i, self.atomTypes[i]['mass']))
		    
		    # ATOMS: *******************************************************************************
			f.write("\nAtoms\n")
			f.write("\n")

			for i,atom in enumerate(atomList):
				atom.atomID = i+1
				f.write("%i %i " % (atom.atomID, atom.atomType))
				f.write("%3.5f %3.5f %3.5f " % (atom.pos[0], atom.pos[1], atom.pos[2])) # positions, using old string format because I am lazy

				f.write("{} ".format(atom.moleculeID)) # molecular tag
				f.write("{} ".format(self.atomTypes[atom.atomType]['diameter'])) # diameter
				f.write("{} ".format(self.atomTypes[atom.atomType]['density'])) # density
				f.write("0 0 0\n")

		    # BONDS: *******************************************************************************
			if len(bondList)>0:
				f.write("\nBonds\n")
				f.write("\n")

				for i,bond in enumerate(bondList):
					f.write("%i " % (i+1)) # bond ID
					f.write("%i " % bond.bondType) # bond type
					f.write("%i %i\n" % tuple( sorted([bond.atoms[0].atomID, bond.atoms[1].atomID] )) ) # atom IDs 

			# ANGLES: ******************************************************************************

			if len(angleList)>0:
				f.write("\nAngles\n")
				f.write("\n")

				for i,ang in enumerate(angleList):
					f.write("%i " % (i+1))
					f.write("%i " % ang.angleType)
					f.write( "%i %i %i\n" % (ang.atoms[0].atomID, ang.atoms[1].atomID, ang.atoms[2].atomID ) ) 

			# DIHEDRALS: ***************************************************************************

			if len(dihedralList)>0:
				f.write("\nDihedrals\n")
				f.write("\n")

				for i,di in enumerate(dihedralList):
					f.write("%i " % (i+1))
					f.write("%i " % di.dihedralType)
					f.write( "%i %i %i %i\n" % (di.atoms[0].atomID, di.atoms[1].atomID, di.atoms[2].atomID, di.atoms[3].atomID) )

			# IMPROPERS: ***************************************************************************

			if len(improperList)>0:
				f.write("\nImpropers\n")
				f.write("\n")

				for i,imp in enumerate(improperList):
					f.write("%i " % (i+1))
					f.write("%i " % imp.improperType)
					f.write( "%i %i %i %i\n" % ( imp.atoms[0].atomID, imp.atoms[1].atomID, imp.atoms[2].atomID, imp.atoms[3].atomID) )  

			# OTHER SECTIONS **************************************************************************
			
			for header in self.otherSections.keys():
				f.write( "\n{}\n".format( header ) )
				f.write( "\n" )
				for line in self.otherSections[header]:
					for item in line:
						#if an object is passed in, find all combinations
						# if isinstance(item, Molecule.Atom) or isinstance(item, Molecule.Bond) or isinstance(item, Molecule.Angle) or isinstance(item, Molecule.Dihedral) or isinstance(item, Molecule.Improper):
						# 	for clone in item.clones:	
						f.write( str(item) )
						f.write( ' ' )
					f.write('\n')

	def _partition(self):
		class _volume:
			def __init__(self,molList, box, corner1, corner2, depth):
				self.molList = molList
				self.parentbox = box
				self.corner1 = np.array(corner1)
				self.corner2 = np.array(corner2)
				if self.parentbox.debug:
					corner3 = (self.corner2-self.corner1)*np.array([1,0,0]) + self.corner1
					corner4 = (self.corner2-self.corner1)*np.array([0,1,0]) + self.corner1
					corner5 = (self.corner2-self.corner1)*np.array([0,0,1]) + self.corner1
					self.parentbox.interactivePlotter.plot_cube([corner1,corner3,corner4,corner5])
					# plt.pause(.05)
				
				self.depth = depth
				volumes = np.array(list((x.effective_radius(self.parentbox.atomTypes) for x in self.molList)))
				volumes = (2**3)*volumes*volumes*volumes
				#have the "box volumes" of polymers, scale to fit the box exactly:
				#"normalize" volumes to sum to 1, working with unit volume is easier, stretch after
				if len(self.molList) == 1:
					#gone as far as I can with partitioning
					self.molList[0].corner1 = self.corner1*(2*0.9*self.parentbox.boxDims) - 0.9*self.parentbox.boxDims
					self.molList[0].corner2 = self.corner2*(2*0.9*self.parentbox.boxDims) - 0.9*self.parentbox.boxDims
					return
				depth = self.depth + 1
				for i,part in enumerate(self._partition_list()):
					corner1 = self.corner1 + 0.5**depth*np.array([i%2, (i//2)%2, (i//4)%2])
					corner2 = self.corner2 + 0.5*(self.corner1-self.corner2) + 0.5**depth*np.array([i%2, (i//2)%2, (i//4)%2])

					_volume( part, self.parentbox, corner1, corner2, depth)
			def _partition_list(self):
				#partition list into 8 parts of near equal volumes
				if len(self.molList)<=8:
					return ([x] for x in self.molList)

				molList = []
				for i,m in enumerate(self.molList):
					molList.append( (self.volumes[i],m) )
				molList.sort()

				partitions = [{'volume':0, 'a':[]},{'volume':0, 'a':[]},{'volume':0, 'a':[]},{'volume':0, 'a':[]},{'volume':0, 'a':[]},{'volume':0, 'a':[]},{'volume':0, 'a':[]},{'volume':0, 'a':[]}]
				while molList:
					#sort partitions by volume, smallest will be at front
					partitions.sort(key=lambda a: a['volume'])
					partitions[0]['a'].append(molList[0][1])
					#ugly but should work, first element in partitions is min, 'a' to ref list of polys, [-1] to ref most recently added (vol,poly) tuple, [0] to get vol
					partitions[0]['volume'] += molList[0][0]
					polyList = polyList[1:]

				temp = []
				for d in partitions:
					temp.append( d['a'] )
				return temp

		self.moleculeList.sort(key=lambda x: x.effective_radius(self.atomTypes))
		_volume( self.moleculeList, self, [0,0,0],[1,1,1], 0)

	def _generate(self):
		self._partition()
		# build polymer using MDS
		for m in self.moleculeList:
			if len(m.atomList) == 0:
				continue
			elif len(m.atomList) == 1:
				m.atomList[0].pos = m.corner1 + (m.corner2 - m.corner1) / 2.
				if self.debug:
					self.interactivePlotter.draw_now(np.copy(m.atomList[0].pos))
			else:
				m.determine_dissimilarity(self.atomTypes)
				mds = manifold.MDS(n_components=3, max_iter=30000, eps=1e-9, dissimilarity="precomputed", n_jobs=1)
				pos = mds.fit(m.dissimilarity).embedding_

				# find box enclosing all points
				xMax = max(pnt[0] for pnt in pos)
				yMax = max(pnt[1] for pnt in pos)
				zMax = max(pnt[2] for pnt in pos)
				xMin = min(pnt[0] for pnt in pos)
				yMin = min(pnt[1] for pnt in pos)
				zMin = min(pnt[2] for pnt in pos)
				corner1 = np.array([xMin,yMin,zMin])
				corner2 = np.array([xMax,yMax,zMax])
				
				for i,pnt in enumerate(pos):
					#move points to be within unit cube [0,1]^3
					pnt -= corner1
					pnt /= (corner2-corner1)

					#translate point to partitioned volume
					pnt *= (m.corner2-m.corner1)
					pnt += m.corner1

					m.atomList[i].pos = np.copy(pnt)
					if self.debug:
						self.interactivePlotter.draw_now(np.copy(pnt))


class Molecule:

	class Atom:
		def __init__(self, atomType, atomNum):
			self.atomType = atomType
			self.bondedAtoms = []

			self.atomNum = atomNum
			return
		def clone(self, mClone):
			mClone.add_atom( self.atomType )
		def __str__(self):
			return str( self.atomID )

	class Bond:
		def __init__(self, bondType, atom1, atom2):
			if atom1 in atom2.bondedAtoms or atom2 in atom1.bondedAtoms:
				raise ValueError('Bond already exists')
			self.bondType = bondType
			self.atoms = (atom1, atom2)
			atom1.bondedAtoms.append(atom2)
			atom2.bondedAtoms.append(atom1)
			return
		def clone(self, mClone):
			atom1 = mClone.atomList[self.atoms[0].atomNum]
			atom2 = mClone.atomList[self.atoms[1].atomNum]
			mClone.bond_atoms(self.bondType, atom1, atom2)

	class Angle:
		def __init__(self, angleType, atom1, atom2, atom3):
			self.angleType = angleType
			self.atoms = [atom1, atom2, atom3]
		def clone(self, mClone):
			atom1 = mClone.atomList[self.atoms[0].atomNum]
			atom2 = mClone.atomList[self.atoms[1].atomNum]
			atom3 = mClone.atomList[self.atoms[2].atomNum]
			mClone.angle_atoms(self.angleType, atom1, atom2, atom3)

	class Dihedral:
		def __init__(self, dihedralType, atom1, atom2, atom3, atom4):
			self.dihedralType = dihedralType
			self.atoms = [atom1, atom2, atom3, atom4]
		def clone(self, mClone):
			atom1 = mClone.atomList[self.atoms[0].atomNum]
			atom2 = mClone.atomList[self.atoms[1].atomNum]
			atom3 = mClone.atomList[self.atoms[2].atomNum]
			atom4 = mClone.atomList[self.atoms[3].atomNum]
			mClone.dihedral_atoms(self.dihedralType, atom1, atom2, atom3, atom4)

	class Improper:
		def __init__(self, improperType, atom1, atom2, atom3, atom4):
			self.improperType = improperType
			self.atoms = [atom1, atom2, atom3, atom4]
		def clone(self, mClone):
			atom1 = mClone.atomList[self.atoms[0].atomNum]
			atom2 = mClone.atomList[self.atoms[1].atomNum]
			atom3 = mClone.atomList[self.atoms[2].atomNum]
			atom4 = mClone.atomList[self.atoms[3].atomNum]
			mClone.improper_atoms(self.improperType, atom1, atom2, atom3, atom4)

	def __init__(self):
		self.atomList = []
		self.bondList = []
		self.angleList = []
		self.dihedralList = []
		self.improperList = []
		self.netradius = None
		return

	def add_atom(self, atomType):
		atom = self.Atom(atomType, len(self.atomList) ) # do not need -1 in last argument
		self.atomList.append( atom )
		return atom

	def bond_atoms(self, bondType, atom1, atom2):
		if not (atom1 in self.atomList and atom2 in self.atomList):
			raise ValueError('Attempted to bond atoms not defined within this molecule')
		elif atom1 is atom2:
			raise ValueError('Attemped to bond an atom to itself')
		else:
			self.bondList.append( self.Bond(bondType, atom1,atom2) )

	def angle_atoms(self, angleType, atom1, atom2, atom3):#order matters!
		self.angleList.append( self.Angle(angleType, atom1, atom2, atom3) )
		return

	def dihedral_atoms(self, dihedralType, atom1, atom2, atom3, atom4):#order matters!
		self.dihedralList.append( self.Dihedral(dihedralType, atom1, atom2, atom3, atom4) )
		return

	def improper_atoms(self, improperType, atom1, atom2, atom3, atom4):#order matters!
		self.improperList.append( self.Improper(improperType, atom1, atom2, atom3, atom4) )
		return

	def define_atoms(self,atomTypeInfo):
		#expect atomTypeInfo to be dict with keys atomType, of dicts with keys 'diameter', 'mass', 'density'
		for atom in self.atomList:
			atom.radius = atomTypeInfo[atom.atomType]['diameter']/2.

	def determine_dissimilarity(self,atomTypeInfo):
		self.define_atoms(atomTypeInfo)
		self.dissimilarity = np.zeros( (len(self.atomList), len(self.atomList)) )
		for i,a1 in enumerate(self.atomList):
			for j,a2 in enumerate(self.atomList):
				self.dissimilarity[i,j] = self.find_shortest_path(a1,a2)
				self.dissimilarity[j,i] = self.dissimilarity[i,j]

	def find_shortest_path(self,a1,a2):
		#take shortest path and return length, length is adjusted for atom radii
		paths = self.find_all_paths(a1,a2)
		minLen = -1
		for p in paths:
			pathLength = sum(2*x.radius for x in p) - p[0].radius - p[-1].radius
			if pathLength < minLen or minLen < 0:
				minLen = pathLength
		return minLen

	def find_all_paths(self, a1,a2, path = []):
		#find all paths from a1 to a2
		path = path + [a1]
		if a1 == a2:
		    return [path]
		if a1 not in self.atomList:
		    return []
		paths = []
		for a_mid in a1.bondedAtoms:
		    if a_mid not in path:
		        extended_paths = self.find_all_paths(a_mid, a2, path)
		        for p in extended_paths: 
		            paths.append(p)
		return paths

	def clone(self):
		m = Molecule()
		for x in self.atomList + self.bondList + self.angleList + self.dihedralList + self.improperList:
			x.clone( m )
		return m

	def effective_radius(self, atomTypeInfo):
		self.define_atoms(atomTypeInfo)
		if self.netradius is None:
			volume = 0
			for atom in self.atomList:
				volume += (atom.radius)**3.
			self.netradius = volume ** (1./3.) / 0.64
			# factor of 1/0.64 comes from the density of randomly packing spheres
			# https://en.wikipedia.org/wiki/Sphere_packing
		return self.netradius


if __name__ == '__main__':

	m = Molecule()
	a_1 = m.add_atom(1)
	a_previous = a_1
	for i in range(20):
		a_current = m.add_atom(2)
		m.bond_atoms(a_previous,a_current)
		a_previous = a_current

	m.bond_atoms(a_previous, a_1)

	box = Box([90,90,90],debug=False)

	for i in range(5):
		box.add_molecule(m)

	m = Molecule()
	a1 = m.add_atom(1)
	a2 = m.add_atom(2)
	a3 = m.add_atom(3)
	m.angle_atoms(1, a1,a2,a3)
	box.add_molecule(m)

	m = Molecule()
	a1 = m.add_atom(1)
	a2 = m.add_atom(2)
	a3 = m.add_atom(3)
	a4 = m.add_atom(4)
	m.dihedral_atoms(1, a1,a2,a3,a4)
	box.add_molecule(m)

	m = Molecule()
	a1 = m.add_atom(1)
	a2 = m.add_atom(2)
	a3 = m.add_atom(3)
	a4 = m.add_atom(4)
	m.improper_atoms(1, a1,a2,a3,a4)
	box.add_molecule(m)

	box.define_atom_type(1,diameter=1)
	box.define_atom_type(2,diameter=1)
	box.define_atom_type(3,diameter=10)
	box.define_atom_type(4)
	box.write_box('LAMMPS_sim.data')