#!/usr/bin/env python
import numpy as np
import os

# Folder to store input parameters and initial conditions, output data and
# generated plots:
outfolder = "./results/"
pfiles = ['benchmark_params.txt',
          'awkwark_front_fork.txt',
          'BatavusBrowserPar.txt',
          'test_params.txt']
parameters = "../parameters/" + pfiles[0]

# options for eigenvalue generation
vi = 0.0
vf = 20.0
N = 300
evalfile = outfolder + "eigenvalues.dat"

os.system('./src/whippleeig ' +
          ' -m ' + parameters +
          ' -i ' + str(vi) +
          ' -f ' + str(vf) +
          ' -n ' + str(N) + 
          ' -o ' + evalfile)

# options for simulation generation
sfiles = ['benchmark_ic.txt',
          'basumandal_states.txt',
          'test.txt',
          'upright.txt',
          'uprightsteady_4.6.txt']
initial_state = '../state/' + sfiles[0]
tf = 5.0
fps = 60.0

os.system('./src/whipplesim' +
          ' -m ' + parameters +
          ' -s ' + initial_state +
          ' -t ' + str(tf) +
          ' -f ' + str(fps) +
          ' -o ' + outfolder)

os.system("cp " + parameters + " " + outfolder + "bike_parameters.txt")
os.system("cp " + initial_state + " " + outfolder + "initial_conditions.txt")
os.system("cp plotfunctions.py " + outfolder)

makeplots=\
"""import numpy as np
import plotfunctions as pf
import matplotlib.pyplot as plt

from record import record_dt, eval_dt

# Get the data from file and put into a custom data type -- examine
# ./simulation.data for details on all the data fields.
sim_data = np.fromfile("simulation.dat", dtype=record_dt)
eval_data = np.fromfile("eigenvalues.dat", dtype=eval_dt)

plot_dict = {'evals': True,
             'orientation': True,
             'contact': True,
             'energy' : True,
             'constraintforces': True,
             'constraints': True}

pf.makeplots(plot_dict, eval_data, sim_data, '')

plt.show()
"""

fp = open(outfolder + "plotter.py", "w")
fp.write(makeplots)
print("Plotter code written to " + fp.name);
fp.close()

import time
ts = time.localtime()
t = (str(ts[0]) + '_'  + str(ts[1]) + '_' + 
     str(ts[2]) + '_' + str(ts[3]) + '_' +
     str(ts[4]) + '_' + str(ts[5]))
tarball_name = "results" + t + ".tar.bz2"

os.system('cd results && python plotter.py')
os.system('tar cjf ' + outfolder + tarball_name + " " + 
           outfolder + '*.pdf ' +
           outfolder + '*.txt ' +
           outfolder + '*.dat ' +
           outfolder + '*.py')
print("tarball written to " + outfolder + tarball_name)

os.system('rm -rf ' + 
           outfolder + '*.pdf ' +
           outfolder + '*.txt ' +
           outfolder + '*.dat ' +
           outfolder + '*.py*')
print("results folder cleanup performed.")
