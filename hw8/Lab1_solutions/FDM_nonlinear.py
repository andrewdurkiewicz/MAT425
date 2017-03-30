# -*- coding: utf-8 -*-
''' Solve the BVP with the finite difference method
     y'' = f(y)
     y(a) = α, y(b) = β
'''

import numpy as np
import matplotlib.pyplot as plt

## Boundary value problem
def F(y):
    return y**2
def DF(y):
    return np.diag(2*y)
a = 0
b = 1
alpha = 1
beta  = 2

## Parameters
N  = 100-1
dx = 1.0/(N+1)
x  = a+dx*np.arange(1,N+1)
# explanation
#------------
# 0      1      2                          N     N+1
# |------|------|------| ... |------|------|------|
# a     a+dx   a+2dx                              b

## numerical parameter
nbIterations = 5
stockY = np.zeros((nbIterations+1,N))

## Initialization
Y = 1+x			# choice first iteration
L = 2*np.diag(np.ones(N)) - np.diag(np.ones(N-1),1) - np.diag(np.ones(N-1),-1)
vecBC = np.zeros(N)
vecBC[0] = alpha
vecBC[-1] = beta
stockY[0,:] = Y


##-----       loop        ---##
##---------------------------##
for k in range(nbIterations):
    ## Newton method
    DJ = L + dx**2*DF(Y)
    tp = np.linalg.lstsq(DJ,L.dot(Y) + dx**2*F(Y)- vecBC)
    Y  = Y - tp[0]
    ## save  
    stockY[k+1,:] = Y
##---------------------------##

## plot
plt.ion()
plt.figure(1)
plt.plot(x,Y)
plt.xlabel(r'$x$',fontsize=25)
plt.ylabel(r'$y$',fontsize=25)
plt.title(r'$y^{\prime\prime}=y^2,\, y(0)=1,y(1)=2$',fontsize=20)
plt.show()

plt.figure(2)
for k in range(nbIterations+1):
    plt.plot(x,stockY[k,:])
plt.xlabel(r'$x$',fontsize=25)
plt.ylabel(r'$y$',fontsize=25)
plt.legend(['k = 0','k = 1','k = 2','k = 3','k = 4'],loc=2)
plt.title("Iterations Newton method",fontsize=20)
plt.show()
