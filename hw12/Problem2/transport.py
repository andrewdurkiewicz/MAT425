#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
'''
 Solve the transport eq. with periodic BC
    ∂_t ρ + ∂_x x(1-x)ρ  ,  0<x<l, t>0
     ρ(x,0) = ρ0(x)  ,  0<x<l, t=0
     ρ(0,t) = ρ(l,t)   ,     t≥0
'''
import numpy as np
import matplotlib.pyplot as plt

##-- parameters
l = 2
T = 6
def c_of_x(x):
    return x*(1-x)
def rho0(x):
    return 1*(0.25<x)*(x<1.75)
##-- numerical paramters
dx = .1
dt = .05
nTime = int(T/dt+.5)
nX    = int(l/dx+.5)+1
## IC
intX    = np.linspace(0,l,nX)
intX_12 = intX[:-1] + dx/2 # dx/2 3dx/2 ... l-dx/2
rho     = np.array(rho0(intX))
## trick
intL = range(nX-1)
intR = range(1,nX)
# init plot
plt.ion()
plt.clf()
pltU = plt.plot(intX, rho)[0]
# deco
plt.axis([0, l, -.1, 3.2])
plt.xlabel(r'position $x$', fontsize=18)
plt.ylabel('density', fontsize=18)
theTitle = plt.title('Solution at T=' + '{:04.2f}'.format(0*dt),
                     horizontalalignment='center', fontsize=20)

##----------------------------------##
##-------      loop         --------##
##----------------------------------##
for k in range(nTime):
##-- estimation flux
  c_i     = np.array(c_of_x(intX))
  c_i_12  = np.array(c_of_x(intX_12))
  for index in c_i_12:
    if index < 0:
      a = 0
      b = 1
    else:
      a = 1
      b = 0
  F_inter = a*c_i[:-1]*rho[:-1] + b*c_i[1:]*rho[1:]
  F_inter = np.concatenate(([F_inter[-1]], F_inter, [F_inter[0]]))
  ## update: ρ = ρ - dt/dx⋅(F_+1/2 - F_-1/2)
  rho = rho - dt/dx*(F_inter[1:] - F_inter[:-1])
  ## plot
  pltU.set_ydata(rho)
  theTitle.set_text('Time t=' + '{:04.2f}'.format(k*dt))
  plt.pause(0.006)
  if k*dt == 4:
    plt.savefig("cNEW.png")