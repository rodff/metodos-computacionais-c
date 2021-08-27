# Applying FCTS method with Jacobi algorithm to iterate bidimensional heat equation
import copy
import numpy as np
import numpy.linalg as lin
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# Initialize
L = 100
tmax = 1500
a = 1.0
lbd1 = 0.1
lbd2 = 0.3
dx = 1.0
dt = 0.25

fig = plt.figure()
ax = fig.add_subplot(111)

def init(tipo,L):
	if(tipo=='random'):
		f = np.random.rand(L,L)
	elif(tipo=='zero'):
		f = np.zeros((L,L))
	elif(tipo=='half'):
		f = np.zeros((L,L))
		f[0:int(L/2.),:] = 1.0
	elif(tipo=='peak'):
		f = np.zeros((L,L))
		f[int(L/2.),int(L/2.)]=1.0
	elif(tipo=='box'):
		f = np.zeros((L,L))
		ctr = int(L/2.)
		l = int(L/4.)
		f[(ctr-l):(ctr+l),(ctr-l):(ctr+l)]=1.0
	
	return f
		

def time_ev(f,f1,tipo):
    if(tipo=='jacobi'):
        for i in range(1,f.shape[0]-1):
        	        for j in range(1,f.shape[1]-1):
                	        f1[i,j]=(f[i,j-1]+f[i,j+1]+f[i+1,j]+f[i-1,j])/4.
    elif(tipo=='gauss-seidel'):
        for i in range(1,f.shape[0]-1):
        	        for j in range(1,f.shape[1]-1):
                	        f1[i,j]=(f1[i,j-1]+f[i,j+1]+f[i+1,j]+f1[i-1,j])/4.
    elif(tipo=='super-relax'):
        alfa = 0.8
        for i in range(1,f.shape[0]-1):
        	        for j in range(1,f.shape[1]-1):
                	        f1[i,j]=-alfa*f[i,j]+(1+alfa)*(f1[i,j-1]+f[i,j+1]+f[i+1,j]+f1[i-1,j])/4.
    
def animate():
        global a
        global lbd1
	global lbd2
	global L
	global tmax 
	global dx
	global dt
	t = 0
	tol = 10**(-6)

	# Initial condition
	f = init('box',L)

        # Boundary conditions
        for k,i in enumerate(f[0,:]):
            f[0,k]=np.exp(-lbd1*k)
            #f[0,k]=0.0
            f[L-1,k]=np.exp(-lbd2*k)
            #f[L-1,k]=0.0
	f[:,0]=1.0
	f[:,L-1]=0.0

	f0 = copy.deepcopy(f)	
	f1 = copy.deepcopy(f)

	img = ax.imshow(f)
	plt.colorbar(img)

	conv=False

	while (t<tmax):
        	print('t = {0}'.format(t))
	        time_ev(f,f1,'super-relax')
				
		delta_f=(f-f1)**2
		conv = np.ndarray.all(delta_f<tol)

	        f=copy.deepcopy(f1)
		
		img.set_array(f)
		fig.canvas.draw()

		if(t%100==0):
			plt.savefig('snapshot_{0}.png'.format(int(t)))

		t+=dt

		if(conv==True):
			print('Convergiu em t = {0}'.format(t))
			break

win = fig.canvas.manager.window
win.after(tmax,animate)
plt.axis([0,L,0,L])
plt.show()
