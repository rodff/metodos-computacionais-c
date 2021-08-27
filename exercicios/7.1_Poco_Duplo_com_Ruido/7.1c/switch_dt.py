import sys
import numpy as np
import matplotlib.pyplot as plt

tmax=300

# Time incrments list
dtL= [0.002,0.01,0.04,0.1]

for dt in dtL:
	alfa=4.0	# Parameter of potential well
	beta=1.5	# Noise parameter
	raizdt=dt
	N=int(tmax/dt)	# Number of time iterations
	dW=beta*raizdt*np.random.normal(0,1,N)	# Array with stochastic process differentials for each t
	x=np.zeros(N)	# Array for x(t)
	t=np.arange(N)*dt	# Array for t
	x[0]=0.		# Initialize particle position

	for i in range(N-1):	# Loop for time evolution
		x[i+1]=x[i]+(alfa-x[i]**2)*x[i]*dt+dW[i]	# x = x + dx = x + x*(alfa-x**2)*dt+dW(t)

	plt.plot(t,x, label='dt={0}'.format(dt))

plt.title('Variation of dt')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.legend()
plt.savefig('dt_variation.png')

