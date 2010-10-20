#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os

os.system('rm -rf simulation.dat')
os.system('./src/simulate')
from record import record_dt

# Get the data from file and put into a custom data type -- examine
# ./simulation.data for details on all the data fields.
data = np.fromfile('./simulation.dat', dtype=record_dt)
plt.figure()
plt.plot(data[:]['t'], data[:]['q0'], label='Yaw')
plt.plot(data[:]['t'], data[:]['q1'], label='Lean')
plt.plot(data[:]['t'], data[:]['q2'], label='Pitch')
plt.plot(data[:]['t'], data[:]['q3'], label='Steer')
plt.xlabel('seconds')
plt.xlabel('radians')
plt.title('Configuration Variables')
plt.legend()


plt.figure()
plt.plot(data[:]['q6'], data[:]['q7'], 'r-', label='Rear contact')
plt.plot(data[:]['fnx'], data[:]['fny'], 'g-', label='Front contact')
plt.xlabel('meters')
plt.xlabel('meters')
plt.title('Contact point locations in x-y plane')
plt.legend()

plt.figure()
plt.plot(data[:]['t'], data[:]['ke'], 'r-', label='KE')
plt.plot(data[:]['t'], data[:]['pe'], 'g-', label='PE')
plt.plot(data[:]['t'], data[:]['pe'] + data[:]['ke'], 'b-', label='TE')
plt.xlabel('seconds')
plt.xlabel('Joules')
plt.title('Energy')
plt.legend()

plt.figure()
plt.subplot(311)
plt.plot(data[:]['t'], data[:]['Rx'], 'r-', label='$R_x$')
plt.plot(data[:]['t'], data[:]['Fx'], 'g-', label='$F_x$')
plt.title('Wheel contact forces')
plt.legend()
plt.subplot(312)
plt.plot(data[:]['t'], data[:]['Ry'], 'r-', label='$R_y$')
plt.plot(data[:]['t'], data[:]['Fy'], 'g-', label='$F_y$')
plt.legend()
plt.subplot(313)
plt.plot(data[:]['t'], data[:]['Rz'], 'r-', label='$R_z$')
plt.plot(data[:]['t'], data[:]['Fz'], 'g-', label='$F_z$')
plt.legend()
plt.xlabel('seconds')

plt.show()
