#using the labsolutions as my barebones code. Otherwise I have no clue how to make the graph move
#in time
import numpy as np
import matplotlib.pyplot as plt

D = 1
l = 10
T = 1
def u0(x):
    return np.sin(np.pi * x / l)


dx = .127
dt = 10**(-2)
ld = D/dx**2*dt

nTime = int(T/dt+.5) + 1
nX    = int(l/dx+.5)

intX  = np.linspace(0,l,nX)
intT = np.linspace(dt,T)
u     = u0(intX)
utrue = u0(intX)
error = np.abs(utrue-u)
plt.ion()
plt.clf()
pltU = plt.plot(intX, u,label = r'$FD$')[0]
pltUtrue = plt.plot(intX,u,label = r'$Exact\/Solution$')[0]
plterror = plt.plot(intX,error,label = r'$Error\/For\/dx=\/ $''%-.2f'%dx)[0]
plt.legend(loc = 1)

plt.axis([0, l, -.1, 1.2])
plt.xlabel(r'position $x$', fontsize=18)
plt.ylabel('density', fontsize=18)
theTitle = plt.title('Solution at T=' + '{:04.2f}, dx = %-.2f'.format(0*dt),
                     horizontalalignment='center', fontsize=20)

for n in range(nTime):
    u_true = np.exp(-1*n*dt*np.pi**(2)/l**(2))*utrue
    u = (1-2*ld)*u + ld*(np.append(u[1:], 0) + np.append(0, u[:-1]))
    u[0] = 0
    u[-1] = 0
    plterror.set_ydata(np.abs(u_true - u))
    pltU.set_ydata(u)
    pltUtrue.set_ydata(u_true)
    theTitle.set_text('Time t=' + '{:04.2f}'.format(n*dt))
    plt.pause(.001)
    if n ==100:
        plt.savefig("FD_dx_%-.2f_t_%-.1f.png"%(dx,n*dt))




