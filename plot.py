#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import plotfunctions as pf

eval_dt = np.dtype([('v', np.float64),
                    ('lambda1', np.float64),
                    ('lambda2', np.float64),
                    ('lambda3', np.float64),
                    ('lambda4', np.float64)])

os.system('rm -rf *.dat')

pfiles = ['benchmark_params.txt',
          'awkwark_front_fork.txt',
          'BatavusBrowserPar.txt',
          'test_params.txt']
parameters = "../parameters/" + pfiles[0]

# options for eigenvalue generation
vi = 0.0
vf = 20.0
N = 300
evalfile1 = "eigenvalues1.dat"
evalfile2 = "eigenvalues2.dat"

os.system('./src/whippleeig ' +
          ' -m ' + parameters +
          ' -i ' + str(vi) +
          ' -f ' + str(vf) +
          ' -n ' + str(N) + 
          ' -o ' + evalfile1)

os.system('./src/whippleeig_gyro ' +
          ' -m ' + parameters +
          ' -i ' + str(vi) +
          ' -f ' + str(vf) +
          ' -n ' + str(N) + 
          ' -o ' + evalfile2)

eval_data1 = np.fromfile(evalfile1, eval_dt)
eval_data2 = np.fromfile(evalfile2, eval_dt)

# options for simulation generation
sfiles = ['benchmark_ic.txt',
          'basumandal_states.txt',
          'test.txt',
          'upright.txt',
          'uprightsteady_4.6.txt']
initial_state = '../state/' + sfiles[0]
tf = 5.0
fps = 60.0
sim_file1 = 'simulation1.dat'
sim_file2 = 'simulation2.dat'

os.system('./src/whipplesim' +
          ' -m ' + parameters +
          ' -s ' + initial_state +
          ' -t ' + str(tf) +
          ' -f ' + str(fps) +
          ' -o ' + sim_file1)

os.system('./src/whipplesim_gyro' +
          ' -m ' + parameters +
          ' -s ' + initial_state +
          ' -t ' + str(tf) +
          ' -f ' + str(fps) +
          ' -o ' + sim_file2 +
          ' -v')

from record import record_dt

# Get the data from file and put into a custom data type -- examine
# ./simulation.data for details on all the data fields.
sim_data1 = np.fromfile(sim_file1, dtype=record_dt)
sim_data2 = np.fromfile(sim_file2, dtype=record_dt)

plot_dict = {'evals': True,
             'orientation': True,
             'contact': True,
             'energy' : True,
             'constraintforces': True,
             'constraints': True}

pf.makeplots(plot_dict, eval_data1, sim_data1, '')
pf.makeplots(plot_dict, eval_data2, sim_data2, '')

plt.show()
