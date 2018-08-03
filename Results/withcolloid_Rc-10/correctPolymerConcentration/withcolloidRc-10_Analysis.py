
import numpy as np
import time
import os
import sys
import subprocess as prcs
sys.path.insert(0, '../../../defs-and-scripts')
import Simulation_Tools as sim # class definition of simulation with useful functions for analysis
import BeadSpringInit_WithColloid as bs

Lx=100*2.
Ly=60*2.
Lz=60*2.
n_p = 57
N = 100
n_s = 20134
Rc = 10

bs.AtomArange(Lx, Ly, Lz, n_p, N, n_s, Rc)
# print 'Cleaning input file'
# data = sim.clean_log('LAMMPS_WCol.o21933', 'cleaned.txt', header='Steps Temp KinE PotE TotE Press Vol')
# print 'Plotting'
# sim.plot_params('post-warmup.png', data, stats_data='Steps Temp KinE PotE TotE Press Vol', numcols=None, numrows=None, t_lo=3e6, t_hi=None, txt = '')