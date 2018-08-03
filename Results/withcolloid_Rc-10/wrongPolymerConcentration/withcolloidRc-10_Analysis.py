
import numpy as np
import time
import os
import sys
import subprocess as prcs
sys.path.insert(0, '../../../defs-and-scripts')
import Simulation_Tools as sim # class definition of simulation with useful functions for analysis

# print 'Cleaning input file'
# data = sim.clean_log('LAMMPS_WCol.o21933', 'cleaned.txt', header='Steps Temp KinE PotE TotE Press Vol')
# print 'Plotting Thermo'
# sim.plot_params('post-warmup.png', data, stats_data='Steps Temp KinE PotE TotE Press Vol', numcols=None, numrows=None, t_lo=3e6, t_hi=None, txt = '')
print 'Plotting Bond Lengths'
sim.analyze_dump('dumpcolloids.coords.dat', style='colloid', COLLOID_TYPE = 4)
