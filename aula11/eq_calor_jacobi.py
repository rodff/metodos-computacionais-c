# Applying FCTS method with Jacobi algorithm to iterate bidimensional heat equation
import copy
import numpy as np
import matplotlib.pyplot as plt

a = 1.0
lbd1 = 0.1
lbd2 = 0.3
L = 30
tmax = 150
t = 0
dx = 1.0
dt = 0.25

# Initial conditions for the entire array
f = np.random.rand(L,L)
f = np.zeros((L,L))

# Boundary conditions
for k,i in enumerate(f[0,:]):
	f[0,k]=np.exp(-lbd1*k)
	f[L-1,k]=np.exp(-lbd2*k)
f[:,0]=1.0
f[:,L-1]=0.0

f0 = copy.deepcopy(f)
f1 = copy.deepcopy(f)

while t<tmax:
	t+=dt
	print('t = {0}'.format(t))
	for i in range(1,f.shape[0]-2):	
		for j in range(1,f.shape[1]-2):
			f1[i,j]=(f[i,j-1]+f[i,j+1]+f[i+1,j]+f[i-1,j])/4.

	f=copy.deepcopy(f1)	


plt.imshow(f)
plt.colorbar()
plt.show()
