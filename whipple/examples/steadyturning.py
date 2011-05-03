#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import os

steady_data = np.dtype([('lean', np.float64),
                        ('steer', np.float64),
                        ('v', np.float64),
                        ('Ts', np.float64),
                        ('ke', np.float64),
                        ('pe', np.float64),
                        ('te', np.float64),
                        ])

leansteer_data = np.dtype([('lean', np.float64),
                           ('steer', np.float64)])

N_phi = 1001
N_delta = 1001
phi_min = 0.001
delta_min = 0.001
phi_max = 17.0*np.pi/180.
delta_max = 6.0*np.pi/180.
p = os.system('./whipplesteady_origin ' +
              '{0} {1} {2} {3} {4} {5}'.format(N_phi,
                                               N_delta,
                                               phi_min,
                                               delta_min,
                                               phi_max,
                                               delta_max))

data = np.fromfile('steadyturning.dat', dtype=steady_data)
stable_i = np.fromfile('stableimag.dat', dtype=leansteer_data)
stable_r = np.fromfile('stablereal.dat', dtype=leansteer_data)

L = data[:]['lean'].reshape((N_phi,N_delta))
S = data[:]['steer'].reshape((N_phi,N_delta))
Z_v = np.abs(data[:]['v']).reshape((N_phi, N_delta))
Z_Ts = data[:]['Ts'].reshape((N_phi, N_delta))
Z_ke = data[:]['ke'].reshape((N_phi, N_delta))
Z_pe = data[:]['pe'].reshape((N_phi, N_delta))
Z_te = data[:]['te'].reshape((N_phi, N_delta))

#width = 6.0
#height = 6.0
xll = .1
yll = .1
xur = .8
yur = .8
fig = plt.figure()#figsize=(width, height))
ax = fig.add_axes((xll, yll, xur, yur))
ax.grid(b=True)
ax.hold(True)
ax.set_xlabel(r'$\phi$ [degrees]')
ax.set_ylabel(r'$\delta$ [degrees]')
ax.set_xlim(phi_min, phi_max, auto=True)
ax.set_ylim(delta_min, delta_max, auto=True)
ax.set_aspect('equal')

CS_Ts = ax.contour(L, S, Z_Ts,
        levels=[-2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -.75, -0.5, -0.25, 0,
            .25, .5],
        antialiased=True,
        colors='k')

CS_te = ax.contour(L, S, Z_te,
        levels=[1693.41, 2565.5],
        antialiased=True,
        colors='k')

plt.clabel(CS_Ts, inline=True, fontsize=6, fmt='$T_\delta$ = %.2f', manual=True)
plt.clabel(CS_te, inline=True, fontsize=6, fmt='ke + pe = %.2f', manual=True)
real = plt.plot(stable_r['lean']*180./np.pi,
                stable_r['steer']*180./np.pi, 'r,',
                alpha=0.05)
imag = plt.plot(stable_i['lean']*180./np.pi,
                stable_i['steer']*180./np.pi, 'g,',
                alpha=0.05)
"""
CS_v = ax.contour(L, S, Z_v, levels=[0, 1.0, 2.0, 4.0, 8.0])
plt.clabel(CS_v, inline=True, fontsize=6, fmt='%0.1f')
tri_real = tri.Triangulation(stable_r['lean']*180./np.pi,
                             stable_r['steer']*180.0/np.pi)
tri_imag = tri.Triangulation(stable_i['lean']*180./np.pi,
                             stable_i['steer']*180.0/np.pi)

# Want to mask all triangles except those who have a neighbor of -1, which
# indicates they are a boundary triangle.

mask_real = np.where(np.where(tri_real.neighbors < 0, 0, 1).all(axis=1), 1, 0)
mask_imag = np.where(np.where(tri_imag.neighbors < 0, 0, 1).all(axis=1), 1, 0)
tri_real.set_mask(mask_real)
tri_imag.set_mask(mask_imag)

plt.plot(tri_real.x[tri_real.edges], tri_real.y[tri_real.edges],  'r-')
plt.plot(tri_imag.x[tri_imag.edges], tri_imag.y[tri_imag.edges], 'g-')
#plt.triplot(tri_imag, 'go-')
"""
plt.title('Steer Torque and Stability Level Curves')
#plt.savefig('steady_benchmark_tau.pdf')
plt.savefig('steady_benchmark_tau.png', dpi=300)
