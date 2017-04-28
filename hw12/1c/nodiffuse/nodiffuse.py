#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

##-- parameters
c = 2
l = 20
T = 100
def rho0(x):
    return np.exp(-(x-10)**(2))

##-- numerical paramters
dx = .02
dt = .01
ld = c*dt/dx

nTime = int(T/dt+.5)
nX    = int(l/dx+.5)+1
## IC
intX  = np.linspace(0,l,nX)
rho   = rho0(intX)
## Warning
CFL = c*dt/dx
if (CFL>1):
    print("!!!! Danger !!!! CFL>1\n")
    print("CFL = {:04.2f}".format(CFL))
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
    ## update
    rho = (1-CFL)*rho + CFL*(np.append(rho[1:],0))
    ## plot
    pltU.set_ydata(rho)
    theTitle.set_text('Time t=' + '{:04.2f}'.format(k*dt))
    plt.pause(0.05)
    if k*dt == 1 or k*dt == 4:
        plt.savefig("rho_time%i_dt-%d_T-%i_dx-%d.png"%(k*dt,dt,T,dx))
    if k*dt > 5:
        break