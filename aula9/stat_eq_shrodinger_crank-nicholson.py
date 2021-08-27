# Apply Crank-Nicholson method on the one dimensional Schrodinger Equation (null potential)
import copy
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt

# Function to generate gaussian initial condition
def init():
    L=100   # Size of the array
    dt=0.5    # Time increment
    dx=1.   # Spatial increment
    t0=0    # Initial time
    tmax=50.   # Maximum time
    k=10
    a=dt/4.j
    b=1+2*a
    z = 1.j
    L0=L/3.     
    S = L/10.
    x = np.arange(0,L,1)
    f = np.mat(np.exp(((-(x-L0)**2)/(2*S*S))+k*z*x)/(np.sqrt(S*np.sqrt(np.pi)))).T
    E = np.identity(L)*b
    for i in range(1,L):
        E[i-1,i]=-a
        E[i,i-1]=-a

    E[0,L-1]=-a
    E[L-1,0]=-a
    M=E.conjugate()
    Einv=lin.inv(E)
    
    return x,f,Einv,M,t0,tmax,dt

x,f,Einv,M,t,tmax,dt = init()	# function(x) array (in this line giving the initial condition)
f1= copy.deepcopy(f)	# Copy of the initial condition to compare sums at the end of iteration

# Iterate over time
while t<tmax:
    t+=dt
    f=np.dot(Einv,np.dot(M,f))		# Update main function array

sum0=(sum(np.multiply(f,np.conj(f))))[0,0]	# Integral at the end

sum1=(sum(np.multiply(f,np.conj(f))))[0,0]	# Integral at the beginning

# Plot initial condition and result of iteration
#title="Integral em t=0: %7.4f, integral em tmax: %7.4f" % (sum0, sum1)
title="Integral em t=0: {0}, integral em tmax: {1}".format(sum0, sum1)
plt.title(title)
plt.plot(x,abs(f),'r--',x,abs(f1))
plt.legend()
plt.show()
