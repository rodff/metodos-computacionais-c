# Apply Crank-Nicholson method on the one dimensional Schrodinger Equation (with harmonic oscillator potential)
import copy
import numpy as np
import numpy.linalg as lin
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

def V_harm_potential(x,L,V0):
    V = V0*(x-L/2.)**2/(L/2.)**2
    return V

# Function to generate gaussian initial condition
def animate():
    L=100   # Size of the array
    dt=0.5    # Time increment
    dx=1.   # Spatial increment
    t=0    # Initial time
    tmax=3500.   # Maximum time
    k=1.	#
    V0 = 1.
    a=dt/4.j	# convenient constant
    b=1+2*a	# convenient constant
    z = 1.j 	# Short variable for imaginary unit
    L0=L/2.    # Mean of the gaussian initialization
    S = 5.	# StdDev of the gaussian initialization

    x = np.arange(0,L,1)	# Space discretization
    # Wave function initialized as a gaussian with diffusion
    f = np.mat(np.exp(((-(x-L0)**2)/(2*S*S))+k*z*x)/(np.sqrt(S*np.sqrt(np.pi)))).T

    E = np.identity(L)*b
    ax.set_autoscale_on(False)
    line, = ax.plot(x,abs(f))
    Vx=V_harm_potential(x,L,V0)
    ax.plot(x,Vx)

    for i in range(1,L):
        E[i-1,i]=-a
        E[i,i-1]=-a

    E[0,L-1]=-a
    E[L-1,0]=-a

    for i in range(0,L):
        E[i,i]=1+2*a*(1+V_harm_potential(x[i],L,V0))

    M=E.conjugate()
    Einv=lin.inv(E)

    # Iterate over time
    while t<tmax:
        t+=dt
        f=np.dot(Einv,np.dot(M,f))		# Iteration step
        line.set_ydata(abs(f[0:L+1]))
        sumf = (sum(np.multiply(f,np.conj(f))))[0,0]	# Calculate integral to check consistency
        title="Integral: {0}, Tempo: {1}".format(sumf,t)
        ax.set_title(title)
        fig.canvas.draw()

win = fig.canvas.manager.window
fig.canvas.manager.window.after(1000,animate)
plt.axis([0,100,0,1.])
plt.show()

