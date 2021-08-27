# Apply Lax method on the drift equation
import copy
import numpy as np
import matplotlib.pyplot as plt

# Function just to initialize variables
def init(Vx=1.0,dt=1.0):
    L=100;dx=1.;t=0;tmax=100.
    k=Vx*dt/(2*dx)
    return L,k,dt,t,tmax

# Function to generate gaussian initial condition
def gauss_init(x,L):
    S = L/10.;
    f0 = np.exp((-(x-L/2.)**2)/(2*S*S))
    return f0

dt_L = np.arange(0.02,2,0.02)	# List of dt's

A_ratio = []		# Mean amplitude ratio of the function
I_ratio = []		# Integral ratio of the function

for dtn in dt_L: 
	print("Integrando para dt={0}".format(dtn))
	# Initialize model parameters 
	L,k,dt,t,tmax = init(dt=dtn)

	x = np.arange(0,L,1)	# 1D space array
	f = gauss_init(x,L)	# function(x) array (in this line giving the initial condition)
	f1= gauss_init(x,L)	# Auxiliary array 
	f2= gauss_init(x,L)	# Copy of the initial condition to compare sums at the end of iteration

	# Iterate over time
	while t<tmax:
	    t+=dt
	    f1[1:L-1]=(1./2)*(f[2:L]+f[0:L-2])-k*(f[2:L]-f[0:L-2])	 # Iteration over the function array
	    f1[L-1]=(1./2)*(f[0]+f[L-2])-k*(f[0]-f[L-2])	# Iteration for the boundaries
	    f1[0]=(1./2)*(f[1]+f[L-1])-k*(f[1]-f[L-1])		# using cyclic boundary condition
	    f=copy.deepcopy(f1)		# Update main function array

	sum1=sum(f1)	# Integral at the end
	sum2=sum(f2)	# Integral at the beginning

	Am=0	# Sum of absolute amplitude ratios
	ni=0	# Number of ignored points
	for i in range(len(f1)):
		if(f1[i]==0 or f2[i]==0):
			ni=ni+1
		else:
			Am=Am+np.abs(f1[i]/f2[i])

	I_ratio.append(sum1/sum2)
	A_ratio.append(Am/(len(f1)-ni))

# Plot stability analisys plots
titulo="Estabilidade de Lax aplicado em eq. da Deriva (tmax={0})".format(tmax)
plt.title(titulo)
plt.plot(dt_L,I_ratio,'r--',label='Integral ratio')
plt.legend()
plt.savefig('integral_ratio_lax_t_{0}.png'.format(tmax))

plt.clf()

titulo="Estabilidade de Lax aplicado em eq. da Deriva (tmax={0})".format(tmax)
plt.title(titulo)
plt.plot(dt_L,A_ratio,'b',label='Amplitude ratio')
plt.legend()
plt.savefig('amplitude_ratio_lax_t_{0}.png'.format(tmax))


# COMMENTS
'''
The integration is not stable for any k = Vx*dt/(2*dx)
'''
