import sys
import numpy as np
import matplotlib.pyplot as plt

tmax=300
dt=0.01
raizdt=np.sqrt(dt)
N=int(tmax/dt)	# Number of time iterations
t=np.arange(N)*dt	# Array for t
alfa= 4.0	# Parameter of potential well
beta= 1.5	# Noise parameter

X = []
M = 100	# Number of realizations

for m in range(M):
	print('Realization: {0}'.format(m))
	dW=beta*raizdt*np.random.normal(0,1,N)	# Array with stochastic process differentials for each t
	x=np.zeros(N)	# Array for x(t)
	x[0]=0.		# Initialize particle position

	for i in range(N-1):	# Loop for time evolution
		x[i+1]=x[i]+(alfa-x[i]**2)*x[i]*dt+dW[i]	# x = x + dx = x + x*(alfa-x**2)*dt+dW(t)

	X.append(x)

X = np.array(X)

Xmed = []
Xvar = []

for k in range(N):
	Xmed.append(np.mean(X[:,k]))
	Xvar.append(np.var(X[:,k]))

plt.plot(t,Xmed)
plt.title('{0} Realizations'.format(M))
plt.xlabel('t')
plt.ylabel('<x(t)>')
plt.savefig('media.png')
plt.clf()

plt.plot(t,Xvar)
plt.title('{0} Realizations'.format(M))
plt.xlabel('t')
plt.ylabel(r'$\sigma^2(t)$')
plt.savefig('variancia.png')
plt.clf()


