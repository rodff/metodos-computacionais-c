import numpy as np
import matplotlib.pyplot as plt

tmax=300	# Maximum time of integration
dt=0.01		# Time increment
alfa=4.		# Parameter of potential well
beta=1.5	# Noise parameter

raizdt=np.sqrt(dt)
N=int(tmax/dt)+1	# Number of time iterations
dW=beta*raizdt*np.random.normal(0,1,N)	# Array with stochastic process differentials for each t
x=np.zeros(N)	# Array for x(t)
t=np.arange(N)*dt	# Array for t
x[0]=0.		# Initialize particle position
M=1
for j in range(M):
    for i in range(N-1):	# Loop for time evolution
        x[i+1]=x[i]+(alfa-x[i]**2)*x[i]*dt+dW[i]	# x = x + dx = x + x*(alfa-x**2)*dt+dW(t)

plt.plot(t,x)
plt.show()
