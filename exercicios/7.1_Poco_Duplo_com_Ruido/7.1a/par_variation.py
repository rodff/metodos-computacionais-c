import sys
import numpy as np
import matplotlib.pyplot as plt

tmax=300
dt=0.01
raizdt=np.sqrt(dt)
N=int(tmax/dt)	# Number of time iterations
t=np.arange(N)*dt	# Array for t
seed = np.random.normal(0,1,N)

# Parameters lists
alfaL= [1.0,4.0,7.0]	# Parameter of potential well
betaL= [0.5,1.5,5.0]	# Noise parameter

for alfa in alfaL:
	beta=1.5
	dW=beta*raizdt*seed	# Array with stochastic process differentials for each t
	x=np.zeros(N)	# Array for x(t)
	x[0]=0.		# Initialize particle position

	for i in range(N-1):	# Loop for time evolution
		x[i+1]=x[i]+(alfa-x[i]**2)*x[i]*dt+dW[i]	# x = x + dx = x + x*(alfa-x**2)*dt+dW(t)

	plt.plot(t,x, label='alfa={0}'.format(alfa))

plt.title('Variation in alfa')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.legend()
plt.savefig('alfa_variation.png')

plt.clf()

for beta in betaL:
	alfa=4.0
	dW=beta*raizdt*seed	# Array with stochastic process differentials for each t
	x=np.zeros(N)	# Array for x(t)
	x[0]=0.		# Initialize particle position

	for i in range(N-1):	# Loop for time evolution
		x[i+1]=x[i]+(alfa-x[i]**2)*x[i]*dt+dW[i]	# x = x + dx = x + x*(alfa-x**2)*dt+dW(t)

	plt.plot(t,x, label='beta={0}'.format(beta))

plt.title('Variation in beta')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.legend()
plt.savefig('beta_variation.png')

