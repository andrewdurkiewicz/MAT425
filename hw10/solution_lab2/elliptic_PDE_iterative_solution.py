# -*- coding: utf-8 -*-
''' Solve the elliptic PDE
    ∂_x^2 u + ∂_y^2 u  =  f    on Ω=[a,b]×[c,d]
          u = g                on ∂Ω
'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
## the functions
def f(x,y):
    return 4
def g(x,y):
    return (x-y)**2
## domain Ω=[a,b]×[c,d]
a = 0.0
b = 1.0
c = 0.0
d = 2.0
## numerical parameters
n = 50
m = 40
dx = (b-a)/n
dy = (d-c)/m
nIterMax = 1000
TOL = 10**(-10)
## init
ld = dx**2/dy**2
error = 10*TOL
nIter = 0

##----------- initialization
##-- matB
matB = np.zeros((n-1,m-1))
for i in range(n-1):
  for j in range(m-1):
      ## init
      xi   = a + (i+1)*dx
      yj   = c + (j+1)*dy
      ## the vector b
      matB[i,j] = dx**2*f(xi,yj)
##-- u
u     = np.zeros((n-1,m-1))
uBC   = np.zeros((n+1,m+1))
## boundary condition
##-- i=1 (x=a) and i=n+1 (x=b)
for j in range(m+1):
    yj = c + j*dy
    uBC[0,j] = g(a,yj)
    uBC[n,j] = g(b,yj)
##-- j=1 (y=c) and j=m+1 (y=d)
for i in range(n+1):
    xi = a + i*dx
    uBC[i,0] = g(xi,c)
    uBC[i,m] = g(xi,d)

##-----------------------------------------##
##---     Solution iterative method     ---##
##-----------------------------------------##

while (error>TOL) and (nIter<nIterMax):
    u_new = -1/(2+2*ld)*(matB - uBC[2:(n+1),1:m] - uBC[:-2,1:m] - ld*uBC[1:n,2:(m+1)] - ld*uBC[1:n,:-2])
    ## update
    error = np.max(abs(u_new-u))
    u     = np.copy(u_new)
    nIter += 1
    uBC[1:n,1:m] = u 

##------------------------##
##---     Solution     ---##
##------------------------##
##-- plot
x = np.linspace(a,b,n+1)
y = np.linspace(c,d,m+1)
X,Y = np.meshgrid(x,y)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, np.transpose(uBC))
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()
