#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt


def derivative_matrix(N):
    M = np.zeros((N, N))
    for i in range(N):
        line = np.array([0] * (i-1) +
                        [-1] * (i >= 1) +
                        [2] +
                        [-1] * (i < N-1) +
                        [0] * (N - i - 2))
        M[i] = line
    return M

def derivative(Y, dx):
    dY = np.zeros(Y.shape)
    dY[0] = (Y[1] - Y[0]) / dx
    for i in range(1, len(Y)-1):
        dY[i] = (Y[i-1] - Y[i+1]) / (2*dx)
    dY[-1] = (Y[-1] - Y[-2]) / dx
    return dY

def F(Y, dx):
    return Y**3 - Y * derivative(Y, dx)

def DF(Y, dx):
    N = len(Y)
    DF = np.zeros((N, N))
    for i, y_i in enumerate(Y):
        if i > 0:
            DF[i, i-1] = y_i / (2*dx)

        if i == 0:
            DF[i, i] = 3 * y_i**2 + (2*y_i - Y[i+1]) / dx
        elif i == N-1:
            DF[i, i] = 3 * y_i**2 + (Y[i-1] - 2*y_i) / dx
        else:
            DF[i, i] = 3 * y_i**2 - (Y[i+1] - Y[i-1]) / (2*dx)

        if i < N-1:
            DF[i, i+1] = -1 * y_i / (2*dx)
    return DF
plt.figure()
for k, N in enumerate((1001,)):
    x0 = 1
    y0 = 1/2
    xb = 2
    yb = 1/3
    dx = (xb - x0) / (N-1)
    x_range = np.linspace(x0, xb, N)
    Y0 = np.zeros(N)
    Y0[0] = y0
    Y0[-1] = yb
    L = derivative_matrix(N)
    Y = np.ones(N)  
    for i in range(5):
        J = L.dot(Y) + dx**2 * F(Y, dx) - Y0
        DJ = L + DF(Y, dx) * dx**2
        dY = np.linalg.solve(DJ, J)
        Y -= dY
plt.plot(x_range, Y, label='$\Delta x=%.3f$' % dx)
plt.title(r'$Solution\/ for\/ y^{\prime\prime} \/= \/y ^ {3}-yy^{\prime}$')
plt.xlabel(r'$x$')
plt.ylabel(r'$y(x)$')
plt.legend()
plt.savefig('Ex3-solution.png')

