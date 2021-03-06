#!/usr/bin/env/ python3
import numpy as np
import matplotlib.pyplot as plt
def F(y):
    return np.cos(y)
def DF(y):
    return np.diag(-np.sin(y))
def calculate_solution(a,b,alpha,beta,N):
	N = N-1
	dx = (b-a)/(N+1)
	x  = a+dx*np.arange(1,N+1)
	iterations = 5
	yaxis = np.diag(np.ones(N))- np.diag(np.ones(N-1),1) - np.diag(np.ones(N-1),-1)
	Y = 1+x			
	L = 2*np.diag(np.ones(N)) - np.diag(np.ones(N-1),1) - np.diag(np.ones(N-1),-1)
	vecBC = np.zeros(N)
	vecBC[0] = alpha
	vecBC[-1] = beta
	yaxis[0,:] = Y
	for k in range(iterations):
	    DJ = L + dx**2*DF(Y)
	    tp = np.linalg.lstsq(DJ,L.dot(Y) + dx**2*F(Y)- vecBC)
	    Y  = Y - tp[0]
	    yaxis[k+1,:] = Y
	return x,Y
 
a = 0
b = 1
alpha = 0
beta = 0
exponent = np.arange(1,4)
N = []
for i in exponent:
	N.append(10**(i))
ysolutions = [[]]*10
xsolutions = [[]]*10
tick = 0
while tick < len(N):
	x,y = calculate_solution(a,b,alpha,beta,N[tick])
	yvalues = np.array(y,dtype = object)
	xvalues = np.array(x,dtype = object)

	ysolutions[tick] = yvalues
	xsolutions[tick] = xvalues
	tick+=1
plotnumber = np.arange(0,len(N))

plt.plot(xsolutions[0],ysolutions[0],label = '$\Delta{x} = 10^{-%d}$' %exponent[0])
plt.plot(xsolutions[1],ysolutions[1],label = '$\Delta{x} = 10^{-%d}$' %exponent[1])
plt.plot(xsolutions[2],ysolutions[2],label = '$\Delta{x} = 10^{-%d}$' %exponent[2])

plt.title(r'$Finite \/Difference \/Method\/ and \/Newton \/Method\/Approximation \/for \/y^{\prime\prime} = \cos(y)$')
plt.xlabel(r'$x$')
plt.ylabel(r'$y(x)$')
plt.legend()
plt.savefig("Problem_2_Solution_vs_dx.png")