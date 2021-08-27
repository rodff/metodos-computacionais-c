# Apply Lax method on the drift equation
import copy
import numpy as np
import matplotlib.pyplot as plt

# Function just to initialize variables
def init():
    L=100;Vx=1.0;dt=1.;dx=1.;t=0;tmax=130.
    k=Vx*dt/(2*dx)
    return L,k,dt,t,tmax

# Function to generate gaussian initial condition
def gauss_init(x,L):
    S = L/10.;
    f0 = np.exp((-(x-L/2.)**2)/(2*S*S))
    return f0

# Initialize model parameters 
L,k,dt,t,tmax = init()

x = np.arange(0,L,1)	# 1D space array
f = gauss_init(x,L)	# function(x) array (in this line giving the initial condition)
f1= gauss_init(x,L)	# Auxiliary array 
f2= gauss_init(x,L)	# Copy of the initial condition to compare sums at the end of iteration

# Iterate over time
while t<tmax:
    t+=dt
    f1[1:L-1]=(1./2)*(f[2:L]+f[0:L-2])-k*(f[2:L]-f[0:L-2])	 # Iteration over the function array
    f1[L-1]=(1./2)*(f[0]+f[L-2])-k*(f[0]-f[L-2])		# Iteration for the boundaries
    f1[0]=(1./2)*(f[1]+f[L-1])-k*(f[1]-f[L-1])		# using cyclic boundary condition
    f=copy.deepcopy(f1)		# Update main function array

sum0=sum(f1)	# Integral at the end
sum1=sum(f2)	# Integral at the beginning

# Plot initial condition and result of iteration
title="Integral em t=0: %7.4f, integral em tmax: %7.4f" % (sum0, sum1)
plt.title(title)
plt.plot(x,f,'r--',x,f2)
plt.show()

# COMMENTS
'''
The iteration is stable when k is equal to 1/2
'''

