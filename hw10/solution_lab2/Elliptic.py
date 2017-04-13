#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
import time

# Initial conditions
x0, x1 = 1, 3
y0, y1 = 0, 1
Ny = 9
error = []
x_range = []
for intX in range(5, 39):
    dx = (x1 - x0) / (intX + 1)
    dy = (y1 - y0) / (Ny + 1)
    l_max = (intX-1)*(Ny-1)
    def f(x, y):
        # u_xx + u_yy = 0
        return 0
    #conditions given in the problem
    def g_x0(y):
        return 1/(1 + y**2)
    def g_x1(y):
        return 3/(9 + y**2)
    def g_y0(x):
        return 1/x
    def g_y1(x):
        return x/(x**2 + 1)
    def ij_l(i_x, i_y):
        #making our matrix for l
        return i_x + (Ny - i_y - 1) * (intX - 1) - 1
    #Setting the initial conditions for our matrix
    A = np.zeros((l_max, l_max))
    B = np.zeros(l_max)

    start_time = time.time()
    for i in range(1, intX):
        for j in range(1, Ny):
            l = ij_l(i, j)
            x_i = x0 + i*dx
            y_i = y0 + j*dy
            l_ip = ij_l(i+1, j)
            l_im = ij_l(i-1, j)
            l_jp = ij_l(i, j+1)
            l_jm = ij_l(i, j-1)
            A[l, l] = -2 * (1 + dx**2 / dy**2)
            B[l] = dx**2 * f(x_i, y_i)
            #enforce BC
            if i != 1:
                A[l_im, l] = 1
            else:
                B[l] -= g_x0(y_i)
            if i != intX-1:
                A[l_ip, l] = 1
            else:
                B[l] -= g_x1(y_i)
            if j != 1:
                A[l, l_jm] = dx**2 / dy**2
            else:
                B[l] -= g_y0(x_i)
            if j != Ny-1:
                A[l, l_jp] = dx**2 / dy**2
            else:
                B[l] -= g_y1(x_i)

    u_vec = np.linalg.solve(A, B)
    u_solution = np.zeros((intX+1, Ny+1))
    # Turn u_solution back into a 2d matrix
    for i in range(intX+1):
        for j in range(Ny+1):
            # At border. These values are not in u_vec
            if i==0 or i==intX or j==0 or j==Ny:
                x_i = x0 + i*dx
                y_i = y0 + j*dy
                if i==0:
                    u_solution[i, j] = g_x0(y_i)
                elif i==intX:
                    u_solution[i, j] = g_x1(y_i)
                if j==0:
                    u_solution[i, j] = g_y0(x_i)
                elif j==Ny:
                    u_solution[i, j] = g_y1(x_i)
                continue
            l = ij_l(i, j)
            u_solution[i, j] = u_vec[l]

    def solution(x, y):
        return x / (x**2 + y**2)

    u_exact_soln = np.zeros((intX+1, Ny+1))
    for i, x in enumerate(np.linspace(x0, x1, intX+1)):
        for j, y in enumerate(np.linspace(y0, y1, Ny+1)):
            u_exact_soln[i, j] = solution(x, y)

    if round(dx, 3) == .1:
        #FDM Solution
        cm = plt.get_cmap('binary')
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        ax1.imshow(u_solution, extent=[x0, x1, y1, y0], cmap=cm)
        ax1.set_title('FDM Solution')
        ax1.set_ylabel('y')
        #Exact Solution Figure
        im2 = ax2.imshow(u_exact_soln, extent=[x0, x1, y1, y0], cmap=cm)
        ax2.set_title('Exact Solution')
        ax2.set_xlabel('x')
        #Our graph locations and plots
        fig.subplots_adjust(right=.8)
        ax3 = fig.add_axes([.85, .15, .05, .7])
        fig.colorbar(im2, cax=ax3)
        plt.savefig('FDMvsExact.png')
    #appending our error so we can graph it later on another graph    
    error.append(np.max(np.abs(u_exact_soln - u_solution)))
    x_range.append(dx)

plt.figure()
plt.plot(x_range, error)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$Error$')
plt.title(r'$Accuracy\/Of\/Our\/Method\/for\/(\Delta y\/=\/.1)$')
plt.savefig('Error_Plot_Problem_1.png')