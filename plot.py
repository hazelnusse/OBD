#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os

os.system('rm -rf simulation.dat')
os.system('./src/whipplesim')
from record import record_dt

eval_dt = np.dtype([('v', np.float64),
                    ('lambda1', np.float64),
                    ('lambda2', np.float64),
                    ('lambda3', np.float64),
                    ('lambda4', np.float64)])

os.system('rm -rf eigenvalues.dat')
os.system('./src/whippleeig')
eval_data = np.fromfile('./eigenvalues.dat', eval_dt)

# Eigenvalue plot
plt.figure()
plt.plot(eval_data[:]['v'], eval_data[:]['lambda1'], 'k,')
plt.plot(eval_data[:]['v'], eval_data[:]['lambda2'], 'k,')
plt.plot(eval_data[:]['v'], eval_data[:]['lambda3'], 'k,')
plt.plot(eval_data[:]['v'], eval_data[:]['lambda4'], 'k,')
plt.axis([0, 10, -10, 10])

# Get the data from file and put into a custom data type -- examine
# ./simulation.data for details on all the data fields.
data = np.fromfile('./simulation.dat', dtype=record_dt)

# Configuration variable plots
f1, (f1a1, f1a2, f1a3) = plt.subplots(3, sharex=True, sharey=False)
f1a1.plot(data[:]['t'], data[:]['q0'], label='Frame Yaw')
f1a1.plot(data[:]['t'], data[:]['q1'], label='Frame Lean')
f1a1.plot(data[:]['t'], data[:]['q2'], label='Frame Pitch')
f1a1.legend(loc=0)
f1a1.set_title('Bicycle orientation angles')
f1a1.set_yticks(f1a1.get_yticks()[1:])
f1a2.plot(data[:]['t'], data[:]['q3'], label='Steer')
f1a2.legend(loc=0)
f1a2.set_yticks(f1a2.get_yticks()[1:-1])
f1a3.plot(data[:]['t'], data[:]['fa_yaw'], label='Fork Yaw')
f1a3.plot(data[:]['t'], data[:]['fa_lean'], label='Fork Lean')
f1a3.plot(data[:]['t'], data[:]['fa_pitch'], label='Fork Pitch')
f1a3.legend(loc=0)
f1a3.set_yticks(f1a3.get_yticks()[:-1])
f1.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f1.axes[:-1]], visible=False)
f1a3.axes.set_xlabel('seconds')
plt.savefig('orientation.eps')

# Contact point plots
plt.figure()
plt.plot(data[:]['q6'], data[:]['q7'], 'r-', label='Rear contact')
plt.plot(data[:]['fnx'], data[:]['fny'], 'g-', label='Front contact')
plt.xlabel('x [meters]')
plt.ylabel('y [meters]')
plt.axis('tight')
plt.title('Contact point locations in x-y plane')
plt.legend(loc=0)
plt.savefig('xy.eps')

# Energy plots
plt.figure()
plt.plot(data[:]['t'], data[:]['ke'], 'r-', label='KE')
plt.plot(data[:]['t'], data[:]['pe'], 'g-', label='PE')
plt.plot(data[:]['t'], data[:]['pe'] + data[:]['ke'], 'b-', label='TE')
plt.xlabel('seconds')
plt.ylabel('Joules')
plt.title('Energy')
plt.legend(loc=0)
plt.savefig('energy.eps')

# Constraint force plots
f2, (f2a1, f2a2, f2a3) = plt.subplots(3, sharex=True, sharey=False)
f2a1.plot(data[:]['t'], data[:]['Rx'], 'r-', label='$R_x$')
f2a1.plot(data[:]['t'], data[:]['Fx'], 'g-', label='$F_x$')
f2a1.set_title('Wheel contact forces')
f2a1.legend(loc=0)
f2a1.set_yticks(f2a1.get_yticks()[1:])
f2a2.plot(data[:]['t'], data[:]['Ry'], 'r-', label='$R_y$')
f2a2.plot(data[:]['t'], data[:]['Fy'], 'g-', label='$F_y$')
f2a2.legend(loc=0)
f2a2.set_yticks(f2a2.get_yticks()[1:-1])
f2a3.plot(data[:]['t'], data[:]['Rz'], 'r-', label='$R_z$')
f2a3.plot(data[:]['t'], data[:]['Fz'], 'g-', label='$F_z$')
f2a3.legend(loc=0)
f2a3.set_yticks(f2a3.get_yticks()[:-1])
f2a3.axes.set_xlabel('seconds')
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f2.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f2.axes[:-1]], visible=False)
plt.savefig('constraintforces.eps')

# Constraint
f3, (f3a1, f3a2, f3a3, f3a4, f3a5) = plt.subplots(5, sharex=True, sharey=False)
f3a1.plot(data[:]['t'], data[:]['fnz'], label='Holonomic constraint')
f3a1.set_title('Numerical satisfaction of constraints')
f3a1.legend(loc=0)
f3a1.set_yticks(f3a1.get_yticks()[1:])
f3a2.plot(data[:]['t'], data[:]['nh1'], label='Non-holonomic constraint 1')
f3a2.legend(loc=0)
f3a2.set_yticks(f3a2.get_yticks()[1:-1])
f3a3.plot(data[:]['t'], data[:]['nh2'], label='Non-holonomic constraint 2')
f3a3.legend(loc=0)
f3a3.set_yticks(f3a3.get_yticks()[1:-1])
f3a4.plot(data[:]['t'], data[:]['nh3'], label='Non-holonomic constraint 3')
f3a4.legend(loc=0)
f3a4.set_yticks(f3a4.get_yticks()[1:-1])
f3a5.plot(data[:]['t'], data[:]['ke'] + data[:]['pe'] -
                       (data[0]['ke'] + data[0]['pe']), label='$\Delta(KE+PE)$')
f3a5.legend(loc=0)
f3a5.set_yticks(f3a5.get_yticks()[:-1])
f3.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f3.axes[:-1]], visible=False)
f3a5.axes.set_xlabel('seconds')
plt.savefig('constraints.eps')
plt.show()
