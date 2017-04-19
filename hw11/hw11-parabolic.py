#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

'''
 Solve the elliptic PDE
    ∂_t u = D ∂_x^2 u      ,  0<x<l, t>0
     u(x,0) = u0(x)        ,  0<x<l, t=0
     u(0,t) = u(l,t) = 0   ,     t≥0
'''
import numpy as np
import matplotlib.pyplot as plt
import scipy

##-- parameters
D = .5**2
l = 10
T = 1
def u0(x):
    # the initial condition
    return np.sin(np.pi * x / l)

##-- numerical paramters
dx = .2
dt = .005
ld = D/dx**2*dt
mass_total = np.array([])
## init
nTime = int(T/dt+.5) + 1
nX    = int(l/dx+.5)
## IC
intX  = np.linspace(0,l,nX)
intT = np.linspace(0,T,nTime)
u     = u0(intX)

# init plot
#plt.ion()
#plt.clf()
for index in mass_total:
    print(index)
pltU = plt.plot(intX, u)[0]
# deco
plt.axis([0, l, -.1, 1.2])
plt.xlabel(r'position $x$', fontsize=18)
plt.ylabel('density', fontsize=18)
theTitle = plt.title('Solution at T=' + '{:04.2f}'.format(0*dt),
                     horizontalalignment='center', fontsize=20)

##----------------------------------##
##-------      loop         --------##
##----------------------------------##
for n in range(nTime):
    ## new time step
    u = (1-2*ld)*u + ld*(np.append(u[1:], 0) + np.append(0, u[:-1]))
    ##trapz(intX,u)
    ## plot
    if (abs(u[0] - u[-1])) > 1e-14:
        #this allows us to check whether the boundry conditions are holding
        raise ValueError('u(0,t) =/= u(l,t) for t>=0')
    pltU.set_ydata(u)
    mass = 0
    for i in u:
        mass += i
    mass_total = np.append(mass_total,mass)
    theTitle.set_text('Time t=' + '{:04.2f}'.format(n*dt))
    plt.pause(0.0001)
plt.plot(intT,mass_total)
plt.ylabel(r'$M(t)$')
plt.xlabel(r'$Time$')
plt.title(r'$Mass\/Evolution\/Over\/Time$')
plt.savefig('mass_vs_time.png')
