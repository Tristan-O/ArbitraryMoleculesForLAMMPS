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
	print('Could not import matplotlib (using Tkinter). Do not use debug argument in \'Box\' Initialization or this will fail')
	print(e)
	print('Continuing...')

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
	def __init__(self, box_dims, num_threads = 1, verbose=True, debug = False):
		self.boxDims = np.array(box_dims) / 2.

		self.moleculeList = []
		self.otherSections = {}
		self.atomTypes = {}

		self.debug = debug
		self.verbose = verbose
		if self.debug:
			self.interactivePlotter = plot3dClass(np.copy(self.boxDims))

				#arbitrary definition of randomness assignment to positions:
		self.random_offset = 1#self.boxDims/sum(self.boxDims)
		self.num_threads = num_threads

	def define_atom_type(self, atomType, mass=1., diameter=1., density=1.):
		atomType = int(atomType)
		self.atomTypes[atomType] = {'mass':mass, 'diameter':diameter, 'density':density}
		return	

	def define_other_section(self, header, line):
		''' '''
		if header not in self.otherSections.keys():
			self.otherSections[header] = []
		self.otherSections[header].append(line)
		return

	def add_molecule(self, m, max_len = -1):
		'''Add a molecule to a box. This clones the molecule object, so one molecule object can be added several times.'''
		#Objects are pass by reference, make clone
		mClone = m.clone(self.atomTypes)
		mClone.max_len = max_len
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

	def _generate(self):
		# build polymer using MDS
		for m_i,m in enumerate(self.moleculeList):
			m.center = np.array([0.,0.,0.])
			m.radius = 0

			if len(m.atomList) == 0:
				continue
			elif len(m.atomList) == 1:
				m.atomList[0].pos = np.array([0.,0.,0.])
				m.radius = self.atomTypes[m.atomList[0].atomType]['diameter']/2.
				
			else:
				m.determine_dissimilarity(self.atomTypes)
				mds = manifold.MDS(n_components=3, max_iter=30000, eps=1e-9, dissimilarity="precomputed", n_jobs=self.num_threads)
				pos = mds.fit(m.dissimilarity).embedding_

				# find box enclosing all points
				
				for i,pnt in enumerate(pos):
					#move points to be within unit cube [0,1]^3
					m.atomList[i].pos = np.copy(pnt)
					# if self.debug:
						# self.interactivePlotter.draw_now(np.copy(pnt))
					m.center += np.copy(pnt)/len(m.atomList)
				for i,pnt in enumerate(pos):
					if m.radius < LA.norm(m.center-pnt)+self.atomTypes[m.atomList[i].atomType]['diameter']/2.:
						m.radius = LA.norm(m.center-pnt)+self.atomTypes[m.atomList[i].atomType]['diameter']/2.
					m.atomList[i].pos -= m.center
			m.center = np.array([0.,0.,0.])
			if self.verbose:
				print('Molecule {} atoms arranged via MDS'.format(m_i))
		#molecules have been built, but not placed into box
		# Strategy: Warped Lattice Technique
		#   Build linear lattice and associate molecules with points on it
		#   such that there are approximately the same number of lattice
		#   points as there are molecules. Then use MDS again, taking distance
		#   between two molecules to be the either the lattice distance or the
		#   minimum distance (the sum of the molecule's radii)

		#building lattice:
		i = 0
		box_vol = 8.*self.boxDims[0]*self.boxDims[1]*self.boxDims[2]
		vol_per_point = box_vol / len(self.moleculeList)
		linear_space_per_point = vol_per_point**(1./3.)

		np.random.shuffle(self.moleculeList)
		for x in np.arange(-self.boxDims[0], self.boxDims[0],linear_space_per_point):
			for y in np.arange(-self.boxDims[1], self.boxDims[1], linear_space_per_point):
				for z in np.arange(-self.boxDims[2], self.boxDims[2], linear_space_per_point):
					if i == len(self.moleculeList):
						break
					self.moleculeList[i].latticePosition = np.array([x,y,z])
					i += 1
		dissimilarity = np.zeros( (len(self.moleculeList), len(self.moleculeList)) )

		#fit molecules into box:
		for i1,m1 in enumerate(self.moleculeList):
			for i2,m2 in enumerate(self.moleculeList):
				if m1 is m2:
					dissimilarity[i1,i2] = 0.
				else:
					dissimilarity[i1,i2] = max( m1.radius+m2.radius, LA.norm(m1.latticePosition-m2.latticePosition) )
		mds = manifold.MDS(n_components=3, max_iter=30000, eps=1e-9, dissimilarity="precomputed", n_jobs=1)
		pos = mds.fit(dissimilarity).embedding_
		for i,m in enumerate(self.moleculeList):
			for a in m.atomList:
				m.center = pos[i]
				a.pos += pos[i] + self.random_offset*np.random.rand(3)

				for j,_ in enumerate(a.pos):
					while not (-self.boxDims[j] < a.pos[j] < self.boxDims[j]):
						if a.pos[j] < -self.boxDims[j]:
							a.pos[j] += 2*self.boxDims[j]
						elif a.pos[j] > self.boxDims[j]:
							a.pos[j] -= 2*self.boxDims[j]
					 
				if self.debug:
					self.interactivePlotter.draw_now( np.copy(a.pos) )
		if self.verbose:
			print('Molecules placed in box.')

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
		self.dissimilarity = None
		self.minLenMax = 3
		return

	def add_atom(self, atomType):
		atom = self.Atom(atomType, len(self.atomList) ) # do not need -1 in last argument
		self.atomList.append( atom )
		if self.dissimilarity is not None:
			self.dissimilarity = None
		return atom

	def bond_atoms(self, bondType, atom1, atom2):
		if not (atom1 in self.atomList and atom2 in self.atomList):
			raise ValueError('Attempted to bond atoms not defined within this molecule')
		elif atom1 is atom2:
			raise ValueError('Attemped to bond an atom to itself')
		else:
			self.bondList.append( self.Bond(bondType, atom1,atom2) )
		if self.dissimilarity is not None:
			self.dissimilarity = None

	def angle_atoms(self, angleType, atom1, atom2, atom3):#order matters!
		self.angleList.append( self.Angle(angleType, atom1, atom2, atom3) )
		if self.dissimilarity is not None:
			self.dissimilarity = None

	def dihedral_atoms(self, dihedralType, atom1, atom2, atom3, atom4):#order matters!
		self.dihedralList.append( self.Dihedral(dihedralType, atom1, atom2, atom3, atom4) )
		if self.dissimilarity is not None:
			self.dissimilarity = None

	def improper_atoms(self, improperType, atom1, atom2, atom3, atom4):#order matters!
		self.improperList.append( self.Improper(improperType, atom1, atom2, atom3, atom4) )
		if self.dissimilarity is not None:
			self.dissimilarity = None

	def define_atoms(self,atomTypeInfo):
		#expect atomTypeInfo to be dict with keys atomType, of dicts with keys 'diameter', 'mass', 'density'
		for atom in self.atomList:
			atom.radius = atomTypeInfo[atom.atomType]['diameter']/2.

	def determine_dissimilarity(self,atomTypeInfo):
		self.define_atoms(atomTypeInfo)
		self.dissimilarity = np.zeros( (len(self.atomList), len(self.atomList)) )
		for i,a1 in enumerate(self.atomList):
			for j,a2 in enumerate(self.atomList[i:]):
				j+=i
				self.dissimilarity[i,j] = self.find_shortest_path(a1,a2)
				self.dissimilarity[j,i] = self.dissimilarity[i,j]
		for i,row in enumerate(self.dissimilarity):
			for j,val in enumerate(row):
				if val < 0:
					val = np.max(self.dissimilarity) + (self.atomList[i].radius + self.atomList[j].radius)*(np.random.rand(1)*0.5 + 1)
	def find_shortest_path(self,a1,a2):
		#take shortest path and return length, length is adjusted for atom radii
		paths = self.find_all_paths(a1,a2)
		minLen = -1
		for p in paths:
			if paths[-1] is a2: #if this path actually goes all the way to a2, because of max_len this may not be satisfied 
				pathLength = sum(2*x.radius for x in p) - p[0].radius - p[-1].radius
				if pathLength < minLen or minLen < 0:
					minLen = pathLength
		return minLen

	def find_all_paths(self, a1,a2, path = []):
		#find all paths from a1 to a2
		path = path + [a1]
		if a1 == a2 or ( len(path) >= self.max_len and self.max_len > 0 ):
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

	def clone(self, atomTypeInfo):
		m = Molecule()
		for x in self.atomList + self.bondList + self.angleList + self.dihedralList + self.improperList:
			x.clone( m )
		return m


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