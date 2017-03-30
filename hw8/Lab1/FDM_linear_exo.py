# -*- coding: utf-8 -*-
''' Solve the BVP with the finite difference method
   y'' = p(x) y' + q(x) y + r(x)
   y(0)= α, y(b)= β
'''
import numpy as np
import matplotlib.pyplot as plt

## Boundary value problem
a     = 0
b     = 1
alpha = 2
beta  = 2

def p(x):
    return -1./(x+1)
def q(x):
    return 2 + 0*x
def r(x):
    return (1-x**2)*np.exp(-x)

## Parameters
N  = 10-1
dx = 1.0/(N+1)
x  = a+dx*np.arange(1,N+1)
# explanation
#------------
# 0      1      2                          N     N+1
# |------|------|------| ... |------|------|------|
# a     a+dx   a+2dx                              b


## Initialization
A = ...
vecB = -dx**2*r(x)
## boundary condition
vecB[0] = ...
vecB[N-1] = ...

##----------------------------##
##-----       solve        ---##
##----------------------------##

Y = np.linalg.lstsq(A,vecB)
## boundary condition
y_sol = np.concatenate(([alpha], Y[0].reshape(-1), [beta]))


##----- plot
x2  = a+dx*np.arange(0,N+2)
plt.plot(x2,y_sol)
plt.xlabel(r'$x$',fontsize=25)
plt.ylabel(r'$y$',fontsize=25)
plt.title(r'$y^{\prime\prime}=p y + q y^\prime + r,\, y(a)=\alpha,y(b)=\beta$',fontsize=20)
plt.show()
