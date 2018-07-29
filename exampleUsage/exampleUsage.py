
import numpy as np
import time
import os
import sys
import subprocess as prcs
sys.path.insert(0, '../defs-and-scripts')
import Simulation_Tools as sim # class definition of simulation with useful functions for analysis
import Polymer as poly

Lx = 100*2.
Ly = 60*2.
Lz = 60*2.
n_p = 26
n_s = 25000
N=100#[5,10,5]
Rc = 0
NArms=5
try:
	os.remove('polymer0.data')
except Exception as e:
	pass

sim.make_dat_file(Lx, Ly, Lz, n_p, n_s=n_s, N=N, NSegments=NArms, NArms=NArms, style='beadspring')
#sim.analyzeDump('dumpcolloids.coords.dat', style='colloid', POLY_END_TYPE = 1, POLY_MID_TYPES = [2], COLLOID_TYPE=4)#, TYPE_IGNORES = [3,4])

simulate_str = prcs.check_output("../defs-and-scripts/lmp_serial -sf omp -pk omp 4 -in ../defs-and-scripts/polymer.in")