#Apply FTCS method on the diffusion equation
import copy
import numpy as np
import matplotlib.pyplot as plt

# Function just to initialize variables
def init():
	L=50;D=1.;dt=0.05;dx=1.;t=0;tmax=100.
	k=D*dt/(dx*dx)    
	return L,k,dt,t,tmax

# Initialize model parameters 
L,k,dt,t,tmax = init()

x = np.arange(0,L,1)	# 1D space array
f = np.zeros(x.shape)	# function(x) array 

# Give initial condition
a = int(L/3)
b = int(2*L/3)
f[a:b]=0.8
f1 = copy.deepcopy(f)	# Auxiliar array
f2 = copy.deepcopy(f)	# Copy of the initial condition to compare sums at the end of iteration

# Iterate over time
while t<tmax:
    t+=dt
    f1[1:L-1]=f[1:L-1]+k*(f[0:L-2]+f[2:L]-2*f[1:L-1])	 # Iteration over the function array
    f1[L-1]=f[L-1]+k*(f[L-2]+f[0]-2*f[L-1])	# Iteration for the boundaries
    f1[0]=f[0]+k*(f[1]+f[L-1]-2*f[0])		# using cyclic boundary condition
    f=copy.deepcopy(f1)		# Update main function array

print("Integration complete.")

sum0=sum(f1)	# Integral at the end
sum1=sum(f2)	# Integral at the beginning

# Plot initial condition and result of iteration
title="Integral em t=0: %7.4f, integral em tmax: %7.4f" % (sum0, sum1)
plt.title(title)
plt.plot(x,f,'r--',x,f2)
plt.show()

'''
Its worth to remeber that from the von Neumann stability analisys
the integration of the diffusion equation through the FTCS method
is stable only for k = D*dt/dx**2 <= 1/2.
'''


