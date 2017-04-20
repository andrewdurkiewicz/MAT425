#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.mlab as mlab

error = []
    ## Initial conditions
T = 1
x0, x1 = 0, 1

dt = .2
dx = .1

c = 1

def u_t0(x):
    # Initial positions
    return np.sin(2 * np.pi * x)
def v_t0(x):
    # Initial velocities
    return 2 * np.pi * np.sin(2 * np.pi * x)

## Setup Parameters
intX = int(1 + round((x1-x0) / dx))
intT = 1 + round(T / dt)
x_range = np.linspace(x0, x1, intX)
t_range = np.linspace(0, T, intT)

#use our t values to make our u matrix
u_prev = u_t0(x_range)  # t_i=0
u_mat = u_prev + dt * v_t0(x_range)  # t_i=1

# u_{n+1} = A (dot) u_{n}
A = np.zeros((intX, intX))
range1 = np.arange(intX)
range2 = np.arange(intX-1)
ld = c**2 * dt**2 / dx**2
A[range1, range1] = 2 - 2*ld
A[range2, range2+1] = ld
A[range2+1, range2] = ld
A[0, -1] = ld  # Coefficients for periodic boundary condition
A[-1, 0] = ld

plot_data = np.zeros((intT, intX))
plot_data[0, :] = u_prev

def solution(x, t):
    return np.sin(2*np.pi*x) * (np.cos(2*np.pi*t) + np.sin(2*np.pi*t))
plot_our_solution = np.zeros_like(plot_data)

## Solve - iterate over time
for t_i in range(intT-1):
    u_prev[:], u_mat = u_mat, A.dot(u_mat) - u_prev

    plot_data[t_i, :] = u_mat
    plot_our_solution[t_i, :] = solution(x_range, t_range[t_i])

error.append(np.max(np.abs(plot_data - plot_our_solution)))

# Setup plot
X, Y = np.meshgrid(range(plot_our_solution.shape[0]), range(plot_our_solution.shape[1]))
Z = plot_our_solution[X, Y]
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot_wireframe(X*dt, Y*dx, Z)
ax.set_xlabel(r"time $t$")
ax.set_ylabel(r"position $x$")
ax.set_title(r'$\Delta x\/=\/10^{-1}\/,\/t\/=\/2 \/ \cdot \/10^{-1}$')

plt.savefig('1b_final_true.png')
plt.show()
