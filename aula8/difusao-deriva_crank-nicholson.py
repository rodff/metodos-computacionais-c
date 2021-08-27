# Apply Crank-Nicholson method on the one dimensional drift equation
import copy
import numpy as np
import numpy.linalg as lin
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

# Function to generate gaussian initial condition
def animate():
    L=100   # Size of the array
    D=1.    # Diffusion constant
    v=1.    # Drift velocity
    dt=0.5    # Time increment
    dx=1.   # Spatial increment
    t=0    # Initial time
    tmax=1000.   # Maximum time
    a=1+(D*dt/(dx*dx))
    b=(v*dt/(4*dx))+(D*dt/(2*dx*dx))
    c=(D*dt/(2*dx*dx))-(v*dt/(4*dx))
    d=1-(D*dt/(dx*dx))
    L0=L/3.     
    S = L/10.
    x = np.arange(0,L,dx)
    f = np.mat(np.exp((-(x-L0)**2)/(2*S*S))).T
    M = np.zeros(shape=(L,L))
    E = np.zeros(shape=(L,L))

    ax.set_autoscale_on(False)
    line, = ax.plot(x,f)

    for j in range(1,L):
        M[j][j]=d
        M[j-1][j]=b
        M[j][j-1]=c
        E[j][j]=a
        E[j-1][j]=-b
        E[j][j-1]=-c
    
    M[0][0]=d
    M[L-1][0]=b
    M[0][L-1]=-c
    E[0][0]=a
    E[L-1][0]=b
    E[0][L-1]=-c
    Einv=lin.inv(E)

    # Iterate over time
    while t<tmax:
        t+=dt
        f=np.dot(Einv,np.dot(M,f))		# Update main function array
        line.set_ydata(f[0:L+1])
        sumf = sum(f)
        title="Integral: {0}, Tempo: {1}".format(sumf,t)
        ax.set_title(title)
        fig.canvas.draw()

win = fig.canvas.manager.window
fig.canvas.manager.window.after(1000,animate)
plt.axis([0,100.,0,1.])
plt.show()

