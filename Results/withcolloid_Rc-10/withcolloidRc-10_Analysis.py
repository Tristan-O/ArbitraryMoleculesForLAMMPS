
import numpy as np
import time
import os
import sys
import subprocess as prcs
sys.path.insert(0, '../../defs-and-scripts')
import Simulation_Tools as sim # class definition of simulation with useful functions for analysis

Lx = 100*2
Ly = 60*2
Lz = 60*2
n_p = 26
n_s = 23251
N = 100
Rc = 10
try:
	os.remove('polymer0.data')
except Exception as e:
	pass
	
print sim.make_dat_file(Lx, Ly, Lz, n_p, N, n_s, Rc=Rc)

simulate_str = prcs.check_output("../../defs-and-scripts/lmp_serial -sf omp -pk omp 4 -in input/polymer.posterity.in")