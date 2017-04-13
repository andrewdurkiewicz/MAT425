#!/usr/bin/env python3

#y'' = 4(y-x)
#y(0) = 0
#y(1) = 1

#exact solution y(x) = (e^2/(e^4 - 1))*(e^(2x) - e^(-2x)) + x

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def p(x):
    return 0+0*x
def q(x):
    return 4 + 0*x
def r(x):
    return -4*x
def bvp(a,b,alpha,beta,N):
	dx = (b-a)/(N)
	x  = a+np.arange(a,b,(dx))
	A = -np.diag(1+dx/2*p(x[1:N]),-1)  + np.diag(2+dx**2*q(x)) - np.diag(1-dx/2*p(x[0:(N-1)]),1)
	vecB = -dx**2*r(x)
	vecB[0] = vecB[0] + (1+dx/2*p(x[0]))*alpha
	vecB[N-1] = vecB[N-1] + (1-dx/2*p(x[N-1]))*beta
	Y = np.linalg.lstsq(A,vecB)
	y_sol = np.concatenate(([alpha], Y[0].reshape(-1), [beta]))
	x2  = a+dx*np.arange(0,N+2)
	return y_sol,x2
def real(x):
	return (np.exp(2)/(np.exp(4) - 1))*(np.exp(2*x) - np.exp(-2*x)) + x
def find_accuracy(y_sol,y_real):
	return np.abs(y_real - y_sol)
def plot1(x2,y_sol,y_real):
	plt.plot(x2,y_sol, label = 'Approx')
	plt.plot(x2,y_real, label = 'Exact')
	plt.legend()
	plt.title(r'$Finite \/Difference \/Method \/Approximation \/for \/y^{\prime\prime} = 4(y-x)$')
	plt.xlabel(r'$x$')
	plt.ylabel(r'$y(x)$')
	plt.savefig("Problem_1_Approx_vs_Real.png")

def plot2(x2,error):
	c = plt.plot(x2,error, label = 'Error')
	plt.legend()
	plt.title(r'$Finite\/ Difference\/ Method\/ Error\/ for \/y^{\prime\prime} = 4(y-x)$')
	plt.xlabel(r'$x$')
	plt.ylabel(r'$y(x)$')
	plt.ylim(0,2)
	plt.savefig("Problem_1_Error_Autofit.png")
y_sol, x2 = bvp(0,1,0,2,100)
y_real = real(x2)
error = find_accuracy(y_sol,y_real)
#plot1(x2,y_sol,y_real)
plot2(x2,error)





