
import numpy as np
import time
import os
import sys

sys.path.insert(0, '../../defs-and-scripts')
import Simulation_Tools as sim # class definition of simulation with useful functions for analysis

logfileLoc = 'nocolloid.out.log'
cleanedfileLoc = 'temp.csv'
header = 'Steps Temp KinE PotE TotE Press Vol'

data = sim.cleanLog(logfileLoc, cleanedfileLoc, header)
print 'done cleaning'

stats_data = header
# stats_data = sim.do_stats(header, cleanedfileLoc, 'stats.data')
# print 'done stats'

sim.plot_params('fullrangeplots.png', data, stats_data)
print 'done plotting 1'

sim.plot_params('aftermixingplots.png', data, stats_data, t_lo=1000000)
print 'done plotting 2'

# sim.plot_params('t.png', data, stats_data, t_lo=0.96125e7, t_hi=0.9615e7)
# print 'done plotting 3'

print 'Done'