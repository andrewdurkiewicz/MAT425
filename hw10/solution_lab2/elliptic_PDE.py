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
##--------------------------------##
##-----         FDM          -----##
##--------------------------------##
lMax = (n-1)*(m-1)
A    = np.zeros((lMax,lMax))
vecB = np.zeros(lMax)
def ij_to_l(i,j):
    return i + (m-1-j)*(n-1)
for i in range(1,n):
    for j in range(1,m):
        ## init
        l    = ij_to_l(i,j)
        xi   = a+i*dx
        yj   = c+j*dy
        l_ip = ij_to_l(i+1,j)
        l_im = ij_to_l(i-1,j)
        l_jp = ij_to_l(i,j+1)
        l_jm = ij_to_l(i,j-1)
        ## the matrix A: stencil with 5 points
        A[l-1,l-1]  = -(2+2*dx**2/dy**2)
        vecB[l-1] = dx**2*f(xi,yj)
        ## extra diago: !! the border !!
        if (i!=1):
            A[l-1,l_im-1] = 1
        else:
            ## x_i-1=a (left)
            vecB[l-1] -= g(a,yj)
        if (i!=(n-1)):
            A[l-1,l_ip-1] = 1
        else:
            ## x_i+1=b (right)
            vecB[l-1] -= g(b,yj)
        if (j!=1):
            A[l-1,l_jm-1] = dx**2/dy**2
        else:
            ## y_i-1=c (down)
            vecB[l-1] -= dx**2/dy**2*g(xi,c)
        if (j!=(m-1)):
            A[l-1,l_jp-1] = dx**2/dy**2
        else:
            ## y_i+1=d (up)
            vecB[l-1] -= dx**2/dy**2*g(xi,d)
##------------------------##
##---     Solution     ---##
##------------------------##
vecU_tp = np.linalg.lstsq(A,vecB)
vecU = vecU_tp[0].reshape(-1)
##--reconstruction matrix solution
matU = np.zeros((n+1,m+1))
for i in range(n+1):
  for j in range(m+1):
    if ((i==0) or (i==n) or (j==0) or (j==m)):
      ## at the border
      xi = a+i*dx
      yj = c+j*dy
      matU[i,j] = g(xi,yj)
    else:
      ## inside
      l = ij_to_l(i,j)
      matU[i,j] = vecU[l-1]
##-- plot
x = np.linspace(a,b,n+1)
y = np.linspace(c,d,m+1)
X,Y = np.meshgrid(x,y)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, np.transpose(matU))
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()
