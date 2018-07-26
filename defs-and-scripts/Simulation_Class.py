# This file provides the definition of the class "event", which is a class that allows one to store various
# parameters of simulations and then simulate under those conditions
# It also contains an algorithm that serves as an example usage, 
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import os
import math
import subprocess as prcs
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize': 8}) #makes endcaps of errorbars visible
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import sys
import time

class Simulation: # a simulation run for one set of parameters, defined in __init__
	L_default = 40.
	N_default = 100
	T_default = 1.0
	k_default = 4.0
	r0_default = 0.0
	A_default = -22.44839027 #
	B_default = 0.25
	ts_default = 0.005
	runs_default = 2000
	exampleInFileName_default = 'polymer.in.example'
	def __init__(self, n_p, n_s, Lx=L_default, Ly=L_default, Lz=L_default, N=N_default, T=T_default, k=k_default, r0=r0_default, A=A_default, B=B_default, timestep=ts_default, runs=runs_default, exampleInFileName=exampleInFileName_default):
		self.Lx = float(Lx)
		self.Ly = float(Ly)
		self.Lz = float(Lz)
		self.N = int(N)
		self.T = float(T)
		self.n_p = int(n_p)
		self.n_s = int(n_s)
		self.k = float(k)
		self.r0 = float(r0)
		self.A = float(A)
		self.B = float(B)
		self.ts = float(timestep)
		self.runs = int(runs)
		self.exampleInFileName = exampleInFileName

		self.dat_made = False
		self.log_cleaned = False

		self.V = Lx*Ly*Lz
		self.timeCreated = time.strftime("%m/%d/%Y %H:%M:%S")


	def make_in_file(self, outdir='./'):
		'''Makes a polymer.in file with the parameters that are necessary'''
		# read from example file, fill in necessary params, output
		# example file has %i, %f written in where applicable, and so just reading the whole string in I can change stuff
		# order of parameters: k, r0, A, B, timestep, runs
		if outdir[-1] != '/': #make outdir a proper directory extension
			outdir += '/'

		if not os.path.exists(self.exampleInFileName):
			print 'Error: example input file does not exist!'
			sys.exit()
		filestr = ''
		with open('polymer.in.example', 'r') as infile:
			filestr = infile.read()
		with open('polymer.in', 'wb') as outfile:
			# k, f, A, B, timestep, runs
			outfile.write(filestr % (self.k, self.r0, self.A, self.B, outdir,self.ts, self.runs, outdir, outdir, outdir, outdir))
		return


	def make_dat_file(self):
		'''Makes a dat file called polymer0.dat with the parameters that are necessary'''
		self.dat_made = True
		if os.path.exists("polymer0.data"):
			try:
				os.remove("polymer0.data")
				os.remove("anim0.pdb") # this is just for my convenience
			except WindowsError:
				pass # never do this kids

		self.dat_str = '*********\nn_p = ' + str(self.n_p) + '\n'
		self.dat_str += prcs.check_output("python ../../defs-and-scripts/BSpolymerExpSolventInit.py -Lx %i -Ly %i -Lz %i -N %i -np %i -ns %i" % (self.Lx, self.Ly, self.Lz, self.N, self.n_p, self.n_s))
		return self.dat_str

	def simulate(self, outdir,simulation_log = None):
		'''Runs a lammps simulation for the given params.
	    Calls do_stats on the results.
		Returns a tuple of the output string of lammps and the output string of stats.py
		- Can pass simulation_log arg a file and skip the simluation part to just load
		results from that file'''
		if simulation_log == None:
			if not self.dat_made:
				print 'Error: Never made polymer0.data file for %i' % self.n_p
				sys.exit()
			if self.runs <= 100:
				print "Warning, not enough runs, this probably can't analyze..."
			self.simulate_str = prcs.check_output("../../defs-and-scripts/lmp_serial -sf omp -pk omp 4 -i polymer.in")
		else: 
			print 'Don\'t forget to include a load in of parameters somehow, N, np etc'
			with open(simulation_log) as f:
				self.simulate_str = f.read()
				# for row in f:
				# 	self.simulate_str += row
		#print self.simulate_str
		self.cleanLog(outdir)
		#print self.clean_str
		self.params_list = self.clean_str.split('\n')[0].split()
		del self.params_list[0] # the first entry is the header marker for stats.py, not important to anything

		self.do_stats(outdir)
		return (self.simulate_str, self.stats_str)


	def cleanLog(self, outdir):
		'''Cleans the log so that it can be put into stats.py'''
		#Look for 10 lines that are just numbers
		print 'Cleaning'
		self.clean_str = ''
		outfile = open(outdir + "/temp.txt", "w")
		validRowCount = 0

		lineStrList_ten = []
		ten_valid = False
		for line in self.simulate_str.split('\n'): # slow but whatever
			line += '\n' # so line includes newline afterward that split('\n') gets rid of
			if not ten_valid:
				lineStrList_ten.append(line)
				if self.is_valid_line(line):
					validRowCount += 1
				else:
					validRowCount = 0

				if len(lineStrList_ten) == 12:
					del lineStrList_ten[0] # keep only previous 11 lines

				if validRowCount >= 10:
					ten_valid = True
					self.clean_str += '# '
					for x in lineStrList_ten:
						self.clean_str += x

			else:
				if self.is_valid_line(line):
					self.clean_str += line

		outfile.write(self.clean_str)
		outfile.close()
		self.log_cleaned = True
		print 'Done cleaning'

	def is_valid_line(self, line):
		'''test for valid line'''
		temp = []
		for cell in line.split(): 
			try:
				temp.append(float(cell))
			except ValueError:
				return False #if one of the entries in this list was not a number
		for i in temp:
			if i != 0.:
				return True #if entire row is zeros, then it is not a valid line
		return False


	def do_stats(self, fileloc):
		'''Calls stats.py and obtains mean and error on each quantity obtained from lammps simulation'''
		if not self.log_cleaned:
			print "Logfile never cleaned for n_p %i" % self.n_p
			sys.exit()

		print 'May get "divide error". This is probably from trying to do stats on volume,\nwhich is always constant and thus has an uncertainty of 0.'
		self.stats_str = ''
		self.stats_data = []
		for param_i in self.params_list:
			print 'Doing stats on ' + str(param_i)
			tmp_str = prcs.check_output(str("python ../../defs-and-scripts/stats.py -o %s -f " + fileloc + "/temp.txt") % param_i)
			for x in tmp_str.split('\n'): #go line by line in the output of stats.py
				if 'Mean' in x: #take third entry in line containing 'Mean' to be mean
					self.stats_data.append([param_i, float(x.split()[3]), float(x.split()[5])]) #split[4] is '+/-'
			self.stats_str += tmp_str
			self.stats_str += '\n\n\n'


	def save(self, outfilestats, outfileraw):
		'''Save mean and std. dev, as well as bond params and everything else relevant to a csv file
		outfilestats is a summary of the data obtained, outfileraw contains everything'''

		infotuple = (self.n_p,self.N,self.n_s)
		if not os.path.exists(outfilestats):
			with open(outfilestats, 'w') as f:
				f.write('np,N,ns,') #write basic header info
				for i in self.stats_data:
					f.write('%s,d(%s),' % (i[0],i[0])) #write stats header info
				f.write('\n')
						#np, N, L, T, k,r0, A, B,ts,runs,
				f.write('%i,%i,%i,' % infotuple) #write first basic data
				for i in self.stats_data: #write data from doing stats
					f.write('%s,%s,' % (i[1],i[2]))
		else:
			with open(outfilestats, 'a') as f: #if file already exists, skip writing header info and just write data
				f.write('\n')
						#np, N, L, T, k,r0, A, B,ts,runs,
				f.write('%i,%i,%i,' % infotuple) #write first basic data
				for i in self.stats_data: #write data from doing stats
					f.write('%s,%s,' % (i[1],i[2]))

		with open(outfileraw, 'a') as f:
			f.write('np,N,ns\n') #write basic header info

					#np, N, L, T, k,r0, A, B,ts,runs,
			f.write('%i,%i,%i\n' % infotuple) # write basic data
			
			for line in self.clean_str.split('\n'): #write raw data, no stats
				for cell in line.split()[1:]:
					f.write('%s,' % cell)
				f.write('\n')

	# def load(self, infilestats=None, infileraw=None):
	# 	'''Loads values from stats data file for computed values, and raw data file for raw simulation data.
	# 	These are expected to be in csv format
	# 	THIS FUNCTION IS UNFINISHED! DO NOT USE!'''


	# 	if infilestats != None:
	# 		if infilestats[-4:0] != '.csv':
	# 			print infilestats[-4:0]
	# 			print 'Warning: File\t' +infilestats+ '\n is not csv. Attempting read anyway...'
	# 		if not os.path.exists(infilestats):
	# 			print 'Error: File\t' +infilestats+ '\ndoes not exist! No data loaded'
	# 		with open(infilestats, 'r') as f:
	# 			list_f = list(f)
	# 			self.params_list = list_f[0].split(',')
	# 			del self.params_list[-1] # newline
	# 			wrong_params_list = ['np','N','L','T','k','r0','A','B','ts','runs','time created']
	# 			for i in range(len(wrong_params_list)):
	# 				indx = self.params_list.index(wrong_params_list[i])
	# 				del self.params_list[indx]
	# 			templist = list(self.params_list) # because deleting an element from a list while iterating causes you to skip an element
	# 			for i in range(len(templist)-1):
	# 				if templist[i+1] == 'd(%s)' % templist[i]:
	# 					indx = self.params_list.index(templist[i+1])
	# 					del self.params_list[indx]

	# 			print "Params list after cleanup:"
	# 			print self.params_list
	# 			params_indices = [] #used later
	# 			# this could technically pick up the d() signifying the uncertainty part if something got changed in save(), hope that doesn't happen
	# 			for i in self.params_list:
	# 				params_indices.append(list_f[0].split(',').index(i))

	# 			for row in list_f:
	# 				cells = row.split(',')
	# 				expected_vals = [self.n_p,self.N,self.L,self.T,self.k,self.r0,self.A,self.B,self.ts,self.runs]
	# 				foundMyRow = True
	# 				for i in range(len(expected_vals)):
	# 					if expected_vals[i] != float(cells[i]):
	# 						print expected_vals[i]
	# 						print cells[i]
	# 						foundMyRow = False

	# 				if foundMyRow:
	# 					self.stats_data = []
	# 					for i in params_indices:
	# 						self.stats_data.append([self.params_list[params_indices.index(i)], float(cells[i]), float(cells[i+1])])
	# 					break


	# 	if infileraw != None:
	# 		if infileraw[-4:0] != '.csv':
	# 			print 'Warning: File\t' +infileraw+ '\n is not csv. Attempting read anyway...'
	# 		if not os.path.exists(infileraw):
	# 			print 'Error: File\t' +infileraw+ '\ndoes not exist! No data loaded'
	# 		with open(infileraw, 'r') as f:
	# 			list_f = list(f)
	# 			header_row = False
	# 			indx = 0
	# 			for row in f:
	# 				cells = row.split(',')

	# 				# If the previous row was a header row, check if the following row has the same values as self does
	# 				my_entry_found = True
	# 				if header_row == True:
	# 					self_vals = [self.n_p,self.N,self.L,self.T,self.k,self.r0,self.A,self.B,self.ts,self.runs,self.timeCreated]
	# 					if len(cells) == len(self_vals):
	# 						for i in range(len(cells)):
	# 							if str(self_vals[i]) != cells[i]:
	# 								# if self value is not the same as what is in the cell... 
	# 								# however I'm not sure how well this works because decimal points could be off
	# 								my_entry_found = False

	# 				# If entry was found (still true after check above)
	# 				if my_entry_found == True:
	# 					indx = list_f.index(row)
	# 					break
									

	# 				# Check if this row is a header row
	# 				expected_vals = ['np','N','L','T','k','r0','A','B','ts','runs','time created']
	# 				header_row = True
	# 				if len(cells) == len(expected_vals):
	# 					for i in range(len(cells)):
	# 						if expected_vals[i] != cells[i]:
	# 							header_row = False
	# 			self.clean_str = ''
	# 			for el in list_f[indx:]:
	# 				# Check if this row is a header row
	# 				cells = el.split(',')
	# 				expected_vals = ['np','N','L','T','k','r0','A','B','ts','runs','time created']
	# 				header_row = True
	# 				if len(cells) == len(expected_vals):
	# 					for i in range(len(cells)):
	# 						if expected_vals[i] != cells[i]:
	# 							header_row = False
	# 				if header_row == True:
	# 					break
	# 				else:
	# 					self.clean_str += el




	def plot_params(self, outfile, numcols=None, numrows=None, t_lo=None, t_hi=None):
		'''Makes several plots vs time of all of the params output by lammps
		This is written to exclude the plot of timestep vs. timestep, assumed
		to be the first column of the clean_str'''
		if not self.log_cleaned:
			print "Error: Attempted to plot before log was cleaned"
			return
		elif t_lo > t_hi:
			print "Error: t_lo > t_hi"
			return

		rows = self.clean_str.split('\n')[1:] # skip zeroth row because it is the header row
		del rows[-1] #last row is empty line
		rawdat = np.zeros((len(rows)-1, len(self.params_list))) # runs x num params

		for rowNum in range(len(rows)):
			cells = rows[rowNum].split()
			for colNum in range(len(cells)):
				# print cells
				rawdat[rowNum-1, colNum] = float(cells[colNum])


		# Only plot in range t_lo to t_hi
		# Do t_hi first to avoid the indices getting messed up
		# This is not the best way of coding, it technically is just deleting every row whose 0th value is greater(less) than t_hi(t_lo)
		if t_hi != None:
			temp = np.zeros((1, len(self.params_list)))
			for i in range(np.shape(rawdat)[0]): # timestep column
				if rawdat[i,0] < t_hi:
					# delete row
					temp = np.vstack((temp, rawdat[i,:]))
			#delete first row (all zeros from definition necessary for np.append())
			rawdat = np.delete(temp, 0, axis=0)

		if t_lo != None:
			temp = np.zeros((1, len(self.params_list)))
			for i in range(np.shape(rawdat)[0]): # timestep column
				if rawdat[i,0] > t_lo:
					# delete row
					temp = np.vstack((temp, rawdat[i,:]))
			#delete first row (all zeros from definition necessary for np.append())
			rawdat = np.delete(temp, 0, axis=0)

		if numcols == None:
			numcols = int(math.ceil(math.sqrt(len(self.params_list))))
		if numrows == None:
			numrows = numcols -1

		fig = plt.figure(figsize=(25, 15))
		for i in range(1, len(self.params_list)):
			plt.subplot(int('%i%i%i' % (numrows, numcols, i)))
			plt.plot(rawdat[:,0], rawdat[:,i])
			fontsize = 18
			plt.title('%s vs Time' % self.params_list[i], fontsize=fontsize)
			plt.ylabel(self.params_list[i], fontsize=fontsize)
			plt.xlabel('Time', fontsize= fontsize)

			ax = plt.gca()

			# For the minor ticks
			ax.xaxis.set_minor_locator(AutoMinorLocator())
			ax.yaxis.set_minor_locator(AutoMinorLocator())
			ax.tick_params(which='both', width=1)
			ax.tick_params(which='major', length=10)
			ax.tick_params(which='minor', length=6)
		txt = "n_p = %i, N = %i, T = %f" % (self.n_p, self.N, self.T)
		fig.text(0.5, 0.05, txt, ha='center', fontsize=fontsize)
		fig.savefig(outfile)#, bbox_inches='tight')
						

# ************************************** END CLASS DEF ***************************************
# ***************************** USER FILL SIMULATIONS LIST HERE ******************************
# simulations = [event(1, runs = 10000000), event(2, runs=10000000), event(10,runs=1000000), event(50)]
# Note: Can change defaults directly in class definition

if __name__ == "__main__":
	simulations = [Simulation(128, 115200, N=100, runs=20000)]#, Sim(10, 10, runs=2000)]

	# ****************************************** MAIN ********************************************




	#logfile = open(outdir+ '/' + outdir + '.LOG', 'w')

	# if file already exists, then don't redo simulation


	# if os.path.exists(outdir):
		# print "File: \t"+outdir+"\nDoes not exist. Must run simulations to generate data.\nThis will take a while..."

	for sim in simulations:
		outdir = "Mean-vs-np_%s" % (time.strftime("%m-%d-%Y_%H-%M-%S"))
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		print "Simulating n_p = %i..." % sim.n_p

		sim.make_in_file(dumpfile=outdir+'/dump.coords.dat')
		sim.make_dat_file()
		#sim.simulate(outdir, simulation_log=)
		#output of these commands are saved to string s, disregarded
		#sim.save(outdir + '/stats.csv', outdir + '/raw.csv')
		#sim.plot_params(outdir + '/plots' + str(sim.n_p) + '.png')

	# plt.errorbar(n_pList, meanList, yerr=sdList, color='r', marker='.', markersize=2, linestyle='None', label="Simulation")
	# fontsize = 18
	# plt.title('Pressure vs Number of polymers\n Each with %i identical atoms in a box of length %i' % (N,L), fontsize=fontsize)
	# plt.ylabel('Pressure', fontsize=fontsize)
	# plt.xlabel('Num Polymers', fontsize= fontsize)
	# #prms.update({'font.size': 22})

	# from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
	# ax = plt.gca()

	# # For the minor ticks
	# ax.xaxis.set_minor_locator(AutoMinorLocator())
	# ax.yaxis.set_minor_locator(AutoMinorLocator())
	# ax.tick_params(which='both', width=1)
	# ax.tick_params(which='major', length=10)
	# ax.tick_params(which='minor', length=6)

	# # The actual curve fitting happens here
	# # get list of n_p's, means before fitting
	# optimizedParameters, pcov = opt.curve_fit(func, n_pList, meanList);

	# x_fit = np.array([])
	# for x in range(0, n_pHigh):
	# 	for i in range(0,100):
	# 		x_fit = np.append(x_fit, x+i/100.)

	# perr = np.sqrt(np.diag(pcov))
	# # Use the optimized parameters to plot the best fit
	# plt.plot(x_fit, func(x_fit, *optimizedParameters), label="Fit, B(T*=1) = %f+/-%f" % (optimizedParameters[0], perr));
	# plt.plot(n_pList, idealP, color='k', marker='s', markerfacecolor='None', markersize=10, linestyle='None', label="Ideal Gas")
	# plt.plot(n_pList, virial, color='g', marker='o', markerfacecolor='None', markersize=12, linestyle='None', label="Virial")

	# print "\nOptimized Params:"
	# print optimizedParameters
	# print perr
	# print ""

	# plt.legend(fontsize=fontsize)

	# plt.show()