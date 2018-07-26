
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
n_s = 1
N1 = 10
N2 = 10
Rc = 0
try:
	os.remove('polymer0.data')
except Exception as e:
	pass

#sim.make_dat_file(Lx, Ly, Lz, n_p, n_s, N1=N1, N2=N2, Rc=Rc, style='diblock')
sim.analyzeDump('dump.coords.dat', style='beadspring', POLY_END_TYPE = 1, POLY_MID_TYPES = [2], TYPE_IGNORES = [3,4])
#simulate_str = prcs.check_output("../../defs-and-scripts/lmp_serial -sf omp -pk omp 4 -in input/polymer.posterity.in")