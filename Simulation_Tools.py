# This file provides the definition of the class "event", which is a class that allows one to store various
# parameters of simulations and then simulate under those conditions
# It also contains an algorithm that serves as an example usage, 
import matplotlib.pyplot as plt
import numpy as np
import scipy.misc as spmisc
from scipy.signal import savgol_filter
#import os
import math
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize': 8}) #makes endcaps of errorbars visible
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
#import sys
import subprocess as prcs
import Polymer as poly
#import pymbar
from decimal import * # values are often too large or too small (truncated to zero) with the exp and ln that happen below
# sys.path.insert(0, '../../defs-and-scripts')

def make_dat_file(Lx, Ly, Lz, n_p, n_s=None, N=None, NArms=None, NSegments=None, style=None, outputFile='./polymer0.data'):
	'''Makes a dat file called polymer0.dat with the parameters that are necessary'''
	node0 = Node(1) #Make a node object to be the first atom in the polymer
	poly0 = poly.Polymer(node0) #Make a polymer from that node
	currAtomType = node0.atom.atomType

	style = style.lower()
	if style == 'beadspring':
		if not isinstance(N,int):
			raise ValueError('Style beadspring requires argument N to be specified as int')

		chain = poly.Chain(N-2, currAtomType+1)
		node0.add_child(chain)

		node1 = poly.Node(currAtomType)
		chain.add_child(node1)

	elif style == 'multiblock':
		if not (isinstance(N, list) or isinstance(N, tuple)):
			raise ValueError('Style multiblock requires argument N to be a list or a tuple')
		if isinstance(N, int):
			N = [N]
		lastNode = Node(currAtomType)
		for n in N:
			currAtomType += 1
			chain = Chain(n-1, currAtomType)
			node0.add_child(chain)
			node0 = Node(currAtomType)
			chain.add_child(node0)
		node0.add_child(lastNode)

	elif style == 'star':
		if N==None or (isinstance(N,int) and not isinstance(NArms,int)) and not (isinstance(N,list) or isinstance(N,tuple)):
			raise ValueError('Style star requires either N,NArms to be specified as int or N to be specified as list')
		if isinstance(N, int):
			N=[N]*NArms
		else:
			if NArms != None and len(N)!=NArms:
				print 'WARNING: Ignoring NArms because you passed in a list of a different size'

		for n in N:
			currAtomType += 1
			chain = Chain(n, currAtomType)
			node0.add_child(chain) #all arms extend from this node

	elif style == 'spine':
		if not (isinstance(N,int) or isinstance(N,list) or isinstance(N, tuple)) and not isinstance(NArms,int):
			raise ValueError('Style star requires either N,NArms to be specified as int or N to be specified as list')
		if isinstance(N, list) or isinstance(N, tuple):
			if len(N) != 3:
				raise ValueError('Style spine requires N to be a list of 3 ints - end-length and length-between-arms, arm-length')
		else:
			N =[N]*3

		endChain = poly.Chain(N[0]-1, currAtomType)
		node0.add_child(endChain)
		newNode = poly.Node(currAtomType)
		endChain.add_child(newNode)

		for n in range(NArms-1):
			spineChain = poly.Chain(N[1],currAtomType)
			armChain = poly.Chain(N[2], currAtomType+1)
			newNode.add_child(spineChain)
			newNode.add_child(armChain)

			newNode = poly.Node(currAtomType)
			spineChain.add_child(newNode)
			
		armChain = poly.Chain(N[2], currAtomType+1)
		newNode.add_child(armChain)

		endChain = poly.Chain(N[0], currAtomType)
		newNode.add_child(endChain)

	elif style == 'ring':
		if not (isinstance(NSegments, int) and isinstance(N, int)):
			raise ValueError('Style ring requires N,NSegments to be specified as int')
		
		loopList = []
		for i in range(NSegments-1):
			node1 = poly.Node(currAtomType)
			chain1 = poly.Chain(currAtomType)
			loopList.append(node1)
			loopList.append(chain1)
			currAtomType += 1
		
		node0.add_child( poly.Loop(loopList) )
	else:
		raise ValueError('Invalid Style: ' + str(style))

	box = poly.Box( [Lx,Ly,Lz], pointSpacing=2 ) #Initialize a box

	for i in range(n_p):
		#start each one in one of the upper edges of the box
		rx = 0.8*(-Lx/2. + Lx*i/float(n_p))
		box.add_polymer(poly0, startPos=[rx, 0,0])

	if n_s != None:
		box.add_solvents(n_s, 3) #generate solvents in the box
	box.write_box(outputFile)

	return

def clean_log(infileloc, outfileloc, header='Steps Temp KinE PotE TotE Press Vol'):
	'''Cleans the log so that it can be put into stats.py, header must have same number of columns'''
	#Look for 10 lines that are just numbers
	import re # regex module

	with open(infileloc, 'r') as f:
			simulate_str = f.read()

	prog = re.compile("^[0-9. ]+$", re.MULTILINE) # match regex thermo output lines 
	result = prog.findall(simulate_str) # (might be imperfect but probably fine, can fix if it starts incorrectly matching)

	#write in csv file format if the user specifies .csv
	if outfileloc[-4:-1] == '.csv':
		delim = ','
	else:
		delim = ' '

	with open(outfileloc, "wb") as f:
		f.write('# ' + header + '\n')
		f.writelines("%s\n" % re.sub("\s+", delim, res.strip()) for res in result)
	
	for i in range(len(result)):
		result[i] = result[i].split()

	return np.array(result, float)

def do_stats(params_list, infileloc, outfileloc=None):
	'''Calls stats.py and obtains mean and error on each quantity obtained from lammps simulation, returns list of mean and std dev'''
	if type(params_list) == type(' '):
		params_list = params_list.split()

	print 'May get "divide error". This is probably from trying to do stats on volume,\nwhich is always constant and thus has an uncertainty of 0.'
	stats_data = []
	for param_i in params_list:
		print 'Doing stats on ' + str(param_i)
		tmp_str = prcs.check_output("python ../../defs-and-scripts/stats.py -o %s -f %s"  % (param_i, infileloc))
		for x in tmp_str.split('\n'): #go line by line in the output of stats.py
			if 'Mean' in x: #take third entry in line containing 'Mean' to be mean
				# Make stats_data a list of lists
				# Each sublist is len 3 and has the parameter string, mean and std dev
				stats_data.append([param_i, float(x.split()[3]), float(x.split()[5])]) #split[4] is '+/-'	

	if outfileloc != None:
		with open(outfileloc, 'wb') as f:
			f.writelines("%s\n" % stat for stat in stats_data)

	return stats_data

def plot_params(outfile, data, stats_data='Steps Temp KinE PotE TotE Press Vol', numcols=None, numrows=None, t_lo=None, t_hi=None, txt = ''):
	'''Makes several plots vs time of all of the params output by lammps
	This is written to exclude the plot of timestep vs. timestep, assumed
	to be the first column of the data list.
	- The list 'data' is assumed to be in the same format as that of the 
	output of the clean_log() function
	- stats_data can be either the output of do_stats(), a list of the parameters
	in string format, or one string of the parameters separated by spaces'''

	if type(stats_data) == type(''):
		stats_data = stats_data.split()
	elif type(stats_data[0]) == type([]):
		temp = []
		for i in stats_data:
			temp.append(i[0][:]) # [:] - pass by value
		stats_data = temp

	if t_lo != None and t_hi != None:
		if t_lo > t_hi:
			print "Error: t_lo > t_hi"
			return

	data = data[:] #so that changes to data in this function do not affect the original argument
	# Only plot in range t_lo to t_hi
	# This is not the best way, it is just deleting every row whose 0th value is greater(less) than t_hi(t_lo)
	if t_hi != None:
		i = len(data)-1
		while i > 0: # timestep column
			if data[i][0] < t_hi:
				# as soon as this hits a step that is less than t_hi
				data = data[0:i]
				break
			i -= 1

	if t_lo != None:
		i = 0
		while i < len(data): # timestep column
			if data[i][0] > t_lo:
				# delete row
				data = data[i:len(data)]
				break
			i += 1

	if numcols == None:
		numcols = int(math.ceil(math.sqrt(len(data[0]))))
	if numrows == None:
		numrows = numcols -1

	fig = plt.figure(figsize=(25, 15))
	for i in range(1, len(data[0])):
		plt.subplot(int('%i%i%i' % (numrows, numcols, i)))
		plt.plot(data[:,0], data[:,i])
		fontsize = 18
		plt.title('%s vs Time' % stats_data[i], fontsize=fontsize)
		plt.ylabel(stats_data[i], fontsize=fontsize)
		plt.xlabel('Timestep', fontsize= fontsize)

		ax = plt.gca()

		# For the minor ticks
		ax.xaxis.set_minor_locator(AutoMinorLocator())
		ax.yaxis.set_minor_locator(AutoMinorLocator())
		ax.tick_params(which='both', width=1)
		ax.tick_params(which='major', length=10)
		ax.tick_params(which='minor', length=6)


	if txt == '':
		txt = outfile
	fig.text(0.5, 0.05, txt, ha='center', fontsize=fontsize)
	fig.savefig(outfile)#, bbox_inches='tight')

def analyze_dump(infile, style='beadspring', title='', POLY_END_TYPE = 1, POLY_MID_TYPES = [2], COLLOID_TYPE = 4, TYPE_IGNORES = [3], id1 = None, id2 = None, nbins=100, expectedPotentialX=None, expectedPotentialY=None, minE=None, infE = None, show=True, save=True, noVolume=False, startStep=0): #It would make a lot more sense to just input which bond distances you want by atom type...
	'''Analyze an atom-type dump file
	INPUT:
		infile: The name of an atom-type dump file generated by LAMMPS
		style: 'polymer' or 'colloid', style of the contents of the dump file
		POLY_END_TYPE and POLY_MID_TYPES: only used if style is 'polymer',
			these are the respective types of the atoms that make up the polymer.
			POLY_MID_TYPES supports multiple types in list form.'''
	import pickle
	if 'pickle' in infile:
		with open(infile, 'rb') as f:
			dists = pickle.load(f)
		print 'Number of data points: ', len(dists)
	else:		
		if type(POLY_MID_TYPES) == int:
			POLY_MID_TYPES = [POLY_MID_TYPES]
		if POLY_END_TYPE in POLY_MID_TYPES:
			print 'ERROR: You specified that the end type of the polymer was in POLY_MID_TYPES.'
			raise ValueError
		if type(TYPE_IGNORES) == int:
			TYPE_IGNORES = [TYPE_IGNORES]

		# Create a class of atoms so that I can 
		class Atom:
			def __init__(self, id_, type_):
				# the _ is there because id and type are built in to python and can't be overridden
				self.Pos = []
				self.Box = []
				self.neighList = []
				self.atomType = int(type_)
				self.atomID = int(id_)

			def addCoord(self, pos, box):
				'''Add a coordinate to this atom
				The expected box is of the form [(xlo,xhi),(ylo,yhi),(zlo,zhi)]'''
				for i in range(3):
					pos[i] = float(pos[i])
					for j in range(2):
						box[i][j] = float(box[i][j])	
				self.Pos.append(pos)
				self.Box.append(box)

			def addNeighbor(self, neighbor):
				'''Specify an atom that is bonded to this one'''
				self.neighList.append(neighbor)


		#First must generate a list of atoms with their respective IDs and atom types
		print 'Creating atom list'
		with open(infile, 'rb') as f:
			dumplist = f.read().split('ITEM: TIMESTEP')

		atoms = []
		del dumplist[0] # this is an empty string


		for i in dumplist[0].split('ITEM: ATOMS id type xs ys zs')[1].split('\n'): 
		# This is a terrible way of looping through the lines that have the initial position info I need
			#print repr(i)
			line = i.split()
			if len(line) < 1:
				continue
			id_ = line[0]
			type_ = line[1]
			atoms.append(Atom(id_, type_))
		atoms.sort(key=lambda atom: atom.atomID) #Sort atoms by ID


		#Fill atoms with position data
		print 'Filling position values'
		skipCount = 0
		for indx,timestepData in enumerate(dumplist):
			timestep = dumplist[indx].split('\n')[1]
			if int(timestep)<int(startStep):
				skipCount+=1
				continue
			temp = timestepData.split('ITEM: ATOMS id type xs ys zs')[1].split('\n')
			temp2 = timestepData.split('ITEM: ATOMS id type xs ys zs')[0].split('ITEM: BOX BOUNDS pp pp pp')[1].split('\n')
			box = []
			for i in temp2:
				if i != '':
					box.append(i.split())
					for i in range(len(box)):
						for j in range(len(box[i])):
							box[i][j] = float(box[i][j])

			if [] in box:
				del box[0]

			for atom_data in temp:
				# print repr(atom_data)
				atom_data = atom_data.split()
				if len(atom_data) == 0:
					continue
				id_ = int(atom_data[0]) # id
				Pos = [float(atom_data[2]), float(atom_data[3]), float(atom_data[4])]

				for atom in atoms:
					if atom.atomID == id_:
						for i in range(3):
							Pos[i] = (box[i][1] - box[i][0])*Pos[i] + box[i][0]
						atom.addCoord(Pos, box)
						break
				else:
					print "ID not found ", id_
					print timestepData.split('ITEM: ATOMS id type xs ys zs')[0]
					raise ValueError

		print 'Number of data points per atom (timesteps recorded):', len(atoms[0].Pos)
		print 'Number of points skipped:', skipCount
		#list of atoms has been filled with ID, type and position data
		#Now to add neighbors
		print 'Adding neighbors'
		if style == 'beadspring':
			polyEnd = False
			for i in range(len(atoms)):
				if atoms[i].atomType == POLY_END_TYPE:
					if not polyEnd:
					# This atom is on the end of a polymer, and it is the beginning of a polymer
						atoms[i].addNeighbor(atoms[i+1])
						polyEnd = True
					else:
						atoms[i].addNeighbor(atoms[i-1])
						polyEnd = False
				elif atoms[i].atomType in POLY_MID_TYPES:
					atoms[i].addNeighbor(atoms[i+1])
					# atoms[i].addNeighbor(atoms[i-1])
				elif atoms[i].atomType not in TYPE_IGNORES:
					print "WARNING: Atom of unknown type encountered."
					print "atom ID ", atoms[i].atomID

		elif style == 'colloid':
			#distance between atoms of one type
			colloids = []
			for i in range(len(atoms)):
				if atoms[i].atomType != COLLOID_TYPE and atoms[i].atomType not in TYPE_IGNORES:
					print "WARNING: Atom of unknown type encountered."
					print "Atom ID ", atoms[i].atomID
				else:
					colloids.append(atoms[i])
			for i in range(len(colloids)):
				for j in range(i+1, len(colloids)):
					atoms[i].addNeighbor(atoms[j])
			atoms = colloids[:]

		elif style == 'id':
			#distance between the atoms with specified ids
			temp = []
			for atom in atoms:
				if atom.atomID == id1:
					temp.append( atom )
					id1 = None
				elif atom.atomID == id2:
					temp.append( atom )
					id2 == None
			atoms = temp

			for i in range(len(atoms)):
				for j in range(i+1, len(atoms)):
					atoms[i].addNeighbor(atoms[j])


		print 'Calculating distance data of neighbors'
		#generate distance data of neighbors
		dists = []

		def nint(x):
			for i,j in enumerate(x):
				x[i] = round(j)
			return x

		for i in atoms:
			for j in i.neighList:
				for k, x1 in enumerate(i.Pos):
					#Find minimum distance - because periodic boundaries are a thing the bond might be across two walls
					#print math.sqrt(sum(min( (i.Pos[k][l]-j.Pos[k][l])**2. , (i.Pos[k][l]+j.Box[k][l][0]-j.Pos[k][l]-j.Box[k][l][1])**2., (i.Pos[k][l]+j.Box[k][l][1]-j.Pos[k][l]-j.Box[k][l][0])**2. ) for l in range(3)))
					# Pretty sure the above formula is right, it's the distance between atom1 and the wall + the distance between atom2 and the opposite wall
					x2 = j.Pos[k]
					L = np.array([0.,0.,0.])
					for box_indx,box_boundary in enumerate(i.Box[k]):
						L[box_indx] = box_boundary[1]-box_boundary[0]
					r_ij = np.array(x2)-np.array(x1)
					r0_ij = r_ij-L*(nint(r_ij/L))
					dist = np.linalg.norm(r0_ij)
					if dist < min(L)/2.:
						dists.append(dist)
		if save:
			with open('dump{}.pickle'.format(title), 'wb') as f:
				pickle.dump(dists, f)
				print 'Distance data saved into dump{}.pickle'.format(title)
	#Finally, generate some histograms
	#I think I want all distances involved, as well as some individual bond distances vs. time although unless the timesteps are small for this, they are going to be useless
	print 'Generating plots'

	fontsize = 18
	counts, bins, bars = plt.hist(dists, nbins)
	U = np.copy(counts)
	plt.title('Histogram of Distances\n{}'.format(infile), fontsize=fontsize)

	fig = plt.gcf()
	if save:
		fig.savefig('Hist{}-{}bins.png'.format(title, nbins))
	if show:
		plt.show()
	plt.close()

	# for i in counts: #printing to see if uncertainty on histograms cannot be just sqrt(N)
	# 	if float(i)/sum(counts) > 0.01:
	# 		print float(i)/sum(counts)
	if show:
		plt.subplot(121)
		plt.plot(bins[:-1], counts,'b-')
	if not noVolume:
		for i,x in enumerate(U):
			U[i] /= bins[i+1]**3. - bins[i]**3.

	if show:
		plt.subplot(122)
		plt.plot(bins[:-1], counts,'r-')
		plt.title('Shape of histogram before and after weighting\nwith radial volume correction')

		plt.show()
	plt.close()

	while 0. in counts:
		for i,x in enumerate(counts):
			if x == 0:
				print "0 counts in bin {}, removing...".format(i)
				counts = np.delete(counts, i)
				U = np.delete(U, i)
				bins = np.delete(bins, i)
				break

	U = -np.log(U)

	plt.title('Potential Energy\n{} - {} bins'.format(infile,nbins), fontsize=fontsize)
	x = []
	for i,b in enumerate(bins):
		if i < len(bins)-1:
			x.append( (bins[i] + bins[i+1])/2. )
	if isinstance(expectedPotentialY, np.ndarray) and isinstance(expectedPotentialX, np.ndarray):
		U -= min(U)-min(expectedPotentialY)
		plt.plot(x, U, 'r.', expectedPotentialX, expectedPotentialY, 'b-')
	elif minE != None:
		U -= min(U)-minE #shift minimum potential energy to min E
		plt.plot(x, U)
	elif infE != None: #behavior as r goes to inf
		U -= U[-1]-infE
		plt.plot(x,U)

	fig = plt.gcf()
	if save:
		fig.savefig('Potential{}-{}bins.png'.format(title,nbins))
	if show:
		plt.show()
	plt.close()

	plt.plot(dists)
	plt.title('Time evolution of distances\n{}'.format(infile), fontsize=fontsize)

	fig = plt.gcf()
	if save:
		fig.savefig('TimeEv{}.png'.format(title))
	plt.close()
	return bins[1:],U,dists,x

def find_density(infile, radius, centerAtomID, surroundAtomTypes, scaledCoords = True):
	'''Finds the density of certain atom types around another atom type using an atom-style dump file'''
	class Atom:
		def __init__(self, id_, type_):
			# the _ is there because id and type are built in to python and can't be overridden
			self.Pos = []
			self.Box = []
			self.neighList = []
			self.atomType = int(type_)
			self.atomID = int(id_)

		def addCoord(self, pos, box):
			'''Add a coordinate to this atom
			The expected box is of the form [(xlo,xhi),(ylo,yhi),(zlo,zhi)]'''
			for i in range(3):
				pos[i] = float(pos[i])
				for j in range(2):
					box[i][j] = float(box[i][j])	
			self.Pos.append(pos)
			self.Box.append(box)
	if not (isinstance(surroundAtomTypes, list) or isinstance(surroundAtomTypes, tuple)):
		surroundAtomTypes = [surroundAtomTypes]

	#First must generate a list of atoms with their respective IDs and atom types
	print 'Creating atom list'
	with open(infile, 'rb') as f:
		dumplist = f.read().split('ITEM: TIMESTEP')

	atoms = []
	del dumplist[0] # this is an empty string
	for i in dumplist[0].split('ITEM: ATOMS id type xs ys zs')[1].split('\n'): 
	# This is a terrible way of looping through the lines that have the initial position info I need
		#print repr(i)
		if i == '':
			continue
		line = i.split()
		id_ = line[0]
		type_ = line[1]
		atoms.append(Atom(id_, type_))
	atoms.sort(key=lambda atom: atom.atomID) #Sort atoms by ID


	#Fill atoms with position data
	print 'Filling position values'
	for timestepData in dumplist:
		temp = timestepData.split('ITEM: ATOMS id type xs ys zs')[1].split('\n')
		temp2 = timestepData.split('ITEM: ATOMS id type xs ys zs')[0].split('ITEM: BOX BOUNDS pp pp pp')[1].split('\n')
		box = []
		for i in temp2:
			if i != '' and not scaledCoords:
				box.append(i.split())	
			elif i != '' and scaledCoords: #coords runs from 0-1 in every direction
				box.append( [0,1] )

		for atom_data in temp:
			# print repr(atom_data)
			if atom_data == '':
				continue
			atom_data = atom_data.split()
			id_ = int(atom_data[0]) # id
			Pos = [float(atom_data[2]), float(atom_data[3]), float(atom_data[4])]

			for atom in atoms:
				if atom.atomID == id_:
					atom.addCoord(Pos, box)
					break
			else:
				print "ID not found ", id_
				print timestepData.split('ITEM: ATOMS id type xs ys zs')[0]
				raise ValueError

	print 'Finding Density at each timestep'
	vol = 4.*math.pi* radius**3. /3.
	atomCount = []
	for i in atoms[0].Pos:
		atomCount.append(0)
	for atom in atoms:
		if atom.atomType in surroundAtomTypes:
			for i in range(len(atom.Pos)):
				Pos = atom.Pos[i]
				center = atoms[centerAtomID].Pos[i]
				Box = atom.Box[i]
				dist = math.sqrt(sum(min( (Pos[l]-center[l])**2. , (Pos[l]-Box[l][0]-center[l]+Box[l][1])**2., (-Pos[l]+Box[l][1]+center[l]-Box[l][0])**2. ) for l in range(3)))
				if dist <= radius:
					atomCount[i] += 1

	print 'Average Density'
	aveDensity = 0
	for i in atomCount:
		i/=vol
		aveDensity += i

	aveDensity /= len(atomCount)
	print aveDensity

def analyze_dumps_umbrella_sampled(infiles, k, r0, nbins, id1,id2, tolerance=0.00001, startStep=0):
	'''k here is a tuple and assumes the potential is of the form U=k(r-r0)^2
	notice that the factor of 0.5 is missing (consistent with LAMMPS notation)'''
	if len(infiles) != len(k) or len(infiles) != len(r0):
		raise ValueError('Arguments must have the same size')

	# x = []
	# y = []
	allDists = []
	temp = []
	for i,f in enumerate(infiles): #because analyze_dump calls plt.close, we cannot plot inside of this function
		print 'Umbrella ',i
		if 'pickle' in f:
			import pickle
			with open(f, 'rb') as fi:
				dists = pickle.load(fi)
			print 'Number of data points: ', len(dists)
		else:
			_, _, dists, _ = analyze_dump(f, title='k{}_r0{}'.format(k[i],r0[i]), nbins=nbins, style='id', id1=id1, id2=id2, noVolume=False, show=False, save=False, startStep=startStep)
		allDists += dists
		temp.append(dists)

	print 'Combined'


	# WHAM Method Implementation:
	# See Scott Shell's lecture 210D notes

	Beta = 1.0 #1/kT
	counts, bins, bars = plt.hist(allDists, nbins)
	plt.close()
	numTot = len(dists)
	BF_z = [] # F as function of z
	z_z = [] # list of z values
	
	getcontext().prec = 200 #number of decimal places to track, arbitrarily large, part of decimal package
	bins = list(bins)

	for i,_ in enumerate(bins):
		bins[i] = bins[i]
	for i,_ in enumerate(bins[:-1]):
		z_z.append((bins[i+1]+bins[i])/2) # bins are the left edges, last element is the right edge, but z should be the center, also bins may have irregular width ()

	# del bins[-1]

	Aj = []
	prevAj = []
	for i,_ in enumerate(infiles): # loop "number of umbrellas" times
		Aj.append(0.)
		prevAj.append( tolerance + 1 )

	debug = 0
	allAj = []
	try:
		while sum(abs(prevAj[i]-Aj[i]) > tolerance for i,_ in enumerate(Aj)):
			debug +=1
			if debug > 400:
				print 'stuck at z=', z
				print BF_z
				print
				print prevAj
				print
				print Aj
				break
			if debug%200 == 0:
				print
				print
				print debug, 'th iteration'
				for i,_ in enumerate(Aj):
					print abs(prevAj[i]-Aj[i])
					print
				#raise KeyboardInterrupt
			BF_z = []
			for i,z in enumerate(z_z):
				shiftTerm = []
				for j,_ in enumerate(infiles):
					eta_jz = -1 * Beta * k[j] * (( z-r0[j] )**2)
					shiftTerm.append(eta_jz + Beta*Aj[j])
				try:
					BF_z.append( -1*math.log( float(counts[i])/float(numTot) ) + spmisc.logsumexp(shiftTerm))
				except ValueError:
					print counts[i]
					print shiftTerm
					raise ValueError
			prevAj = []
			allAj.append([])
			for x in Aj:
				prevAj.append(x) #want new object so call constructor
				allAj[-1].append(x)

			Aj = [0.]
			for j,_ in enumerate(infiles):
				if j != 0:
					Aj.append( [] )
					for indx,z in enumerate(z_z):
						Aj[-1].append(-1*BF_z[indx] - Beta*k[j]*( z-r0[j] )**2 ) #1e6 to avoid overflow errors
					Aj[-1] = -1* spmisc.logsumexp( Aj[-1] )

			# print Aj

	except KeyboardInterrupt:
		pass

	print debug, ' iterations'
	plt.title('{} bins'.format(nbins)) #With ( +ln(v) )

	for i,_ in enumerate(z_z):
		BF_z[i] += 2/Beta*math.log(z_z[i]) # entropy correction
	plt.plot(z_z, BF_z)
	plt.show()
	plt.close()

	allAj = np.array(allAj)
	allAj = np.transpose(allAj)
	print np.shape(allAj)
	for i in allAj:
		plt.plot(i)
	plt.show()
	plt.close()
		# bins = list(bins)
		# del bins[0]
		# del bins[0]
		# del bins[-1]
		# del bins[-1]
		# bins = np.array(bins)
	# 	try:
	# 		bins_kn = np.vstack((bins_kn,bins))
	# 	except UnboundLocalError:
	# 		bins_kn = np.copy(bins)
	# 	x.append(np.copy(r))
	# 	y.append(np.copy(BU))

	# 	BU -= k[i]*( r-r0[i] )**2.

	# 	x.append(np.copy(r))
	# 	y.append(np.copy(BU))

	plt.title('Briefly all of the samples on one plot')
	x_i = []
	y_i = []
	for i,x in enumerate(temp):
		counts,bins,_ = plt.hist(x, nbins)
		# plt.show()
		# plt.close()
		counts = -np.log(counts)
		plt.plot(bins[:-1],counts)
		plt.savefig('{}.png'.format(i))
		# plt.show()
		# plt.close()
		y_i.append(np.copy(counts))
		x_i.append(np.copy(bins[:-1]))

	for i,_ in enumerate(x_i):
		y_i[i] = y_i[i] - float(k[i])*(x_i[i]-float(r0[i]))**2 + i
		plt.plot(x_i[i], y_i[i])
	plt.show()
	plt.close()




	# plt.xlabel('r')
	# plt.ylabel(r'$\beta$U(r)')
	# #plt.show()
	# #plt.close()

	# for i,_ in enumerate(x):
	# 	if i%2 == 0:
	# 		continue

	# 	Z = [j for _,j in sorted(zip(x[i],y[i]))]
	# 	y1=savgol_filter(Z, 51,7)
	# 	plt.subplot(313)
	# 	plt.plot(x[i],y1, label='Savitsky-Golay smoothed curve')

		#plt.plot(x[i], y[i], 'k-', label='original potential reconstructed from histogram')
		#plt.legend()
	#plt.show()
	# plt.close()

	# #*********************************************************************************************************************
	# # Allocate storage for simulation data

	# N_k = np.zeros(len(allDists), np.int32) # N_k[k] is the number of snapshots from umbrella simulation k
	# N_max = -1
	# for i,x in enumerate(allDists):
	# 	N_k[i] = len(x)
	# 	N_max = max(N_max,len(x))

	# K_k = np.array(k, np.float64) # K_k[k] is the spring constant (in kJ/mol/nm**2) for umbrella simulation k
	# dist0_k = np.array(r0, np.float64) # dist0_k[k] is the spring center location (in nm) for umbrella simulation k
	# dist_kn = np.array(allDists, np.float64) # dist_kn[k,n] is the distance set point (in nm) for snapshot n from umbrella simulation k
	# u_kn = np.zeros([len(k),N_max], np.float64) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
	# # this can just be zeros because it cancels out... Or something
	# u_kln = np.zeros([len(k), len(k),N_max], np.float64)

	# for k_i,_ in enumerate(k):
	# 	for n in range(N_k[k_i]):
	# 		# Compute minimum-image distance deviation from umbrella center l
	# 		ddist = dist_kn[k_i,n] - dist0_k

	# 		# Compute energy of snapshot n from simulation k in umbrella potential l
	# 		u_kn[k_i,:,n] = (K_k[k_i]) * ddist**2

	# Initialize MBAR.

	# print "Running MBAR..."
	# mbar = pymbar.MBAR(u_kln, N_k, verbose = True)

	# allDists = np.array(allDists)
	# x_n_sorted = np.sort(allDists) # unroll to n-indices
	# bins = np.append(x_n_sorted[0::int(N_max/nbinsPer)], x_n_sorted.max()+0.1)

	# bin_n = np.zeros(allDists.shape, np.int64)
	# bin_n = np.digitize(allDists, bins) - 1
	# # Compute PMF in unbiased potential (in units of kT).
	# (f_i, df_i) = mbar.computePMF(u_kn, bin_n, nbinsPer)
	# plt.plot(df_i)#, bins)
	# plt.show()
	# plt.close()

	return


if __name__ == "__main__":
# 	#simulate_str = prcs.check_output("../../defs-and-scripts/lmp_serial -sf omp -pk omp 4 -in polymer.in")
	infile = '../R12warmup.coords.dat'
	surroundAtomTypes = [3]
	centerAtomID = 44624
	radius = 2
	find_density(infile, radius, centerAtomID, surroundAtomTypes)