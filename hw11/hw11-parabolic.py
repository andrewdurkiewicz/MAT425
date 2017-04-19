#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
a = .5**2
l = 10
T = 30
dx = .4
dt = .2
ld = a/dx**2*dt #lambda
Nx = int(1+round(10/dx))
Nt = 1 + round(T/dt)
x_range = np.linspace(0,l,Nx)
t_range = np.linspace(0,T,Nt)


def u0(x):
    #initial conditions for when t = 0
    return np.sin(np.pi * x / l)

#initialize u matrix:
u = u0(x_range)

#U_{n+1} = A.u_n
A = np.zeros((Nx,Nx))
#using the matrix form of our equation:
i = np.arange(Nx)
j = np.arange(Nx-1)

A[i,i] = 1-2*ld
A[j,j+1] = ld
A[j+1,j] = ld
A[0,-1] = ld
A[-1,0] = ld

step = 20  # Plot every 10 time steps
plot_data = np.zeros((Nt//step+2, Nx))
plt_i = 0
masses = np.zeros(Nt)


for t_i in range(Nt):
    u = A.dot(u)
    masses[t_i] = np.sum(u)  # Probably off by a scalar factor
    if t_i % step == 0 or t_i == Nt-1:
        plot_data[plt_i, :] = u
        plt_i += 1
graph = 0
while graph < 3:
    if graph == 0:
        # Setup plot
        plt.plot(t_range, masses)
        plt.title(r'$System\/Mass\/At\/Time\/t$')
        plt.xlabel(r'$Time$')
        plt.ylabel(r'$Total\/mass$')
        plt.savefig('mass_at_time_t.png')
    if graph == 2:
        X, Y = np.meshgrid(range(plot_data.shape[0]), range(plot_data.shape[1]))
        Z = plot_data[X, Y]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot_wireframe(X*dt, Y*dx, Z)
        ax.set_xlabel(r"time $t$")
        ax.set_ylabel(r"position $x$")
        ax.set_zlabel(r'$u(x,t)$')


        plt.savefig("3d.png")
    #lets find how much the mass deviates from the initial mass:
    if graph == 1:
        def mass_change(m0,m):
            return abs(m-m0)
        change_mass = np.array([])
        for i in range(len(masses)):
            change_mass=np.append(change_mass,(mass_change(masses[0],masses[i])))
        plt.plot(t_range,change_mass)
        plt.ylim(-1,1)
        plt.xlabel(r'$Time(s)$')
        plt.ylabel(r'$\Delta m$')
        plt.title(r'$Mass\/Evolution\/vs\/Time$')
        plt.savefig('Mass_Evolution_vs_time.png')
    graph+=1