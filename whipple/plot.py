#!/usr/bin/env python
import numpy as np
import os

#TODO Need to fix all occurances of './' and '/' so that everything runs under
# cygwin.
# Folder to store all simulation input and output for each run
outfolder = "results/"

# Parameter files:
pfiles = ['benchmark.txt',
          'benchmark_toroidal.txt',
          'test.txt']
parameters = "bikeparameters/" + pfiles[0]

# TODO add this information to a text file of some sort that is stored along
# with the eigenvalue output data
# options for eigenvalue generation
vi = 0.0
vf = 10.0
N = 301

os.system('whippleeig ' +
          ' -m ' + parameters +
          ' -i ' + str(vi) +
          ' -f ' + str(vf) +
          ' -n ' + str(N) + 
          ' -o ' + outfolder)

# options for simulation generation
sfiles = ['benchmark_ic.txt',
          'upright_static.txt',
          'test.txt']
ic = 'initialconditions/' + sfiles[0]
tf = 5.0
fps = 60.0

os.system('whipplesim' +
          ' -m ' + parameters +
          ' -s ' + ic +
          ' -t ' + str(tf) +
          ' -f ' + str(fps) +
          ' -o ' + outfolder)

# options for steady turning plots
os.system('whipplesteady' +
          ' -o ' + outfolder)


# Copy the parameters and initial state that were used as inputs to the
# eigenanalysis / simulation
os.system("cp " + parameters + " " + outfolder + "bike_parameters.txt")
os.system("cp " + ic + " " + outfolder + "initial_conditions.txt")
os.system("cp plotfunctions.py " + outfolder)

# Code that is packaged along with data so that people can regenerate the
# plots.
makeplots=\
"""import numpy as np
import plotfunctions as pf
import matplotlib.pyplot as plt
from sim_record import sim_dt
from eval_record import eval_dt
from boundary_record import boundary_dt

# Get the data from file and put into a custom data type -- examine
# sim_record.py, eval_record.py for details on data fields.
sim_data = np.fromfile("simulation.dat", dtype=sim_dt)
eval_data = np.fromfile("eigenvalues.dat", dtype=eval_dt)
boundary_data = np.fromfile("boundary.dat", dtype=boundary_dt)

# Choose which plots to generate here
plot_dict = {'evals': True,
             'orientation': True,
             'contact': True,
             'energy' : True,
             'constraintforces': True,
             'constraints': True,
             'cfglim' : True,
             'feasibleboundary' : True}

# Make the plots and save them to file
pf.plotcontroller(plot_dict,
                  sim_data,
                  [eval_data],
                  steady_data=[boundary_data],
                  folder="")

# Display the plots on screen
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

os.system('cd ' + outfolder + ' && python plotter.py && ' +
          'tar cjf ' + tarball_name + ' ' + 
           '*.pdf *.txt *.dat *.py')
print("tarball written to " + outfolder + tarball_name)

os.system('rm -rf ' + 
           outfolder + '*.pdf ' +
           outfolder + '*.txt ' +
           outfolder + '*.dat ' +
           outfolder + '*.py*')
print("results folder cleanup performed.")
