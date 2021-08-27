# Apply Crank-Nicholson method on the one dimensional drift equation
import copy
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt

# Function to generate gaussian initial condition
def init():
    L=100   # Size of the array
    D=1.    # Diffusion constant
    v=1.    # Drift velocity
    dt=0.5    # Time increment
    dx=1.   # Spatial increment
    t0=0    # Initial time
    tmax=50.   # Maximum time
    a=1+(D*dt/(dx*dx))
    b=(v*dt/(4*dx))+(D*dt/(2*dx*dx))
    c=(D*dt/(2*dx*dx))-(v*dt/(4*dx))
    d=1-(D*dt/(dx*dx))
    L0=L/3.     
    S = L/10.
    x = np.arange(0,L,1)
    f = np.mat(np.exp((-(x-L0)**2)/(2*S*S))).T
    M = np.zeros(shape=(L,L))
    E = np.zeros(shape=(L,L))
    for j in range(1,L):
        M[j][j]=d
        M[j-1][j]=c
        M[j][j-1]=b
        E[j][j]=a
        E[j-1][j]=-c
        E[j][j-1]=-b
    
    M[0][0]=d
    M[L-1][0]=b
    M[0][L-1]=c
    E[0][0]=a
    E[L-1][0]=-b
    E[0][L-1]=-c
    Einv=lin.inv(E)
    
    return x,f,Einv,M,t0,tmax,dt

x,f,Einv,M,t,tmax,dt = init()	# function(x) array (in this line giving the initial condition)
f1= copy.deepcopy(f)	# Copy of the initial condition to compare sums at the end of iteration

print(M)
print(Einv)

# Iterate over time
while t<tmax:
    t+=dt
    f=np.dot(Einv,np.dot(M,f))		# Update main function array

sum0=sum(f)	# Integral at the end
sum1=sum(f1)	# Integral at the beginning

# Plot initial condition and result of iteration
title="Integral em t=0: %7.4f, integral em tmax: %7.4f" % (sum0, sum1)
plt.title(title)
plt.plot(x,f,'r--',x,f1)
plt.legend()
plt.show()

