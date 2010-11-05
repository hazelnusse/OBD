#!/usr/bin/env python
import numpy as np
import os

# Folder to store all simulation input and output for each run
outfolder = "./results/"

# Parameter files: TODO change where these are located
pfiles = ['benchmark_params.txt',
          'awkwark_front_fork.txt',
          'BatavusBrowserPar.txt',
          'test_params.txt']
parameters = "../parameters/" + pfiles[0]

# TODO add this information to a text file of some sort that is stored along
# with the eigenvalue output data
# options for eigenvalue generation
vi = 0.0
vf = 20.0
N = 300

os.system('./src/whippleeig ' +
          ' -m ' + parameters +
          ' -i ' + str(vi) +
          ' -f ' + str(vf) +
          ' -n ' + str(N) + 
          ' -o ' + outfolder)

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

# Copy the parameters and initial state that were used as inputs to the
# eigenanalysis / simulation
os.system("cp " + parameters + " " + outfolder + "bike_parameters.txt")
os.system("cp " + initial_state + " " + outfolder + "initial_conditions.txt")
os.system("cp plotfunctions.py " + outfolder)

# Code that is packaged along with data so that people can regenerate the
# plots.
makeplots=\
"""import numpy as np
import plotfunctions as pf
import matplotlib.pyplot as plt
from sim_record import sim_dt
from eval_record import eval_dt

# Get the data from file and put into a custom data type -- examine
# ./record.py for details on all the data fields.
sim_data = np.fromfile("simulation.dat", dtype=sim_dt)
eval_data = np.fromfile("eigenvalues.dat", dtype=eval_dt)

# Choose which plots to generate here
plot_dict = {'evals': True,
             'orientation': True,
             'contact': True,
             'energy' : True,
             'constraintforces': True,
             'constraints': True}

# Make the plots and save them to file
pf.makeplots(plot_dict, eval_data, sim_data, folder='')

# Display the plots on screen
# plt.show()
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
