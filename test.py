
import numpy as np
import time
import os
import sys

sys.path.insert(0, '../../defs-and-scripts')
import Simulation_Class as sim # class definition of simulation with useful functions for analysis


simulations = [sim.Simulation(10, 1, N=100, Lx=100,Ly=60,Lz=60, runs=10)]
print "Have you edited the prefactors in beadwall?"
for sim in simulations:
	#outdir = "Mean-vs-np_%s" % (time.strftime("%m-%d-%Y_%H-%M-%S"))
	#if os.path.exists(outdir):
	#s.system('rm -rf ./' + outdir)
	#os.makedirs(outdir)

	print "Simulating N atoms = %i..." % (sim.n_p*sim.N + sim.n_s)
	#Run *********
	#sim.make_in_file()
	sim.make_dat_file()
	print sim.dat_str
	
	#output of these commands are saved to string s, disregarded
	#sim.save(outdir + '/stats.csv', outdir + '/raw.csv')

	#Analysis ****
	sim.simulate('.')
