# Applying FCTS method with Jacobi algorithm to iterate bidimensional heat equation
import copy
import numpy as np
import numpy.linalg as lin
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)


# Initialize
L = 50
tmax = 5000
a = 2.0
lbd1 = 0.1
lbd2 = 0.3

def animate():
	global a
	global lbd1
	global lbd2
	global L
	global tmax 
	t = 0
	dx = 1.0
	dt = 0.25
	tol = 10**(-4)

	# Initial conditions for the entire array
	f = np.random.rand(L,L)
	#f = np.zeros((L,L))
	#f[10:40,10:40]=1.0

	# Boundary conditions
	for k,i in enumerate(f[0,:]):
        	f[0,k]=np.exp(-lbd1*k)
	        f[0,k]=1.0
		f[L-1,k]=np.exp(-lbd2*k)
		f[L-1,k]=0.0
	f[:,0]=1.0
	f[:,L-1]=0.0

	f0 = copy.deepcopy(f)	
	f1 = copy.deepcopy(f)

	img = ax.imshow(f)
	plt.colorbar(img)

	conv=False

	while (t<tmax and conv==False):
	        t+=dt
        	print('t = {0}'.format(t))
	        for i in range(1,f.shape[0]-2):
        	        for j in range(1,f.shape[1]-2):
                	        f1[i,j]=(f[i,j-1]+f[i,j+1]+f[i+1,j]+f[i-1,j])/4.
		
		delta_f=(f-f1)**2
		conv = np.ndarray.all(delta_f<tol)

	        f=copy.deepcopy(f1)
		
		img.set_array(f)
		fig.canvas.draw()

win = fig.canvas.manager.window
win.after(tmax,animate)
plt.axis([0,L,0,L])
plt.show()
