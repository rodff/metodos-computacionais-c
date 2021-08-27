import numpy as np
from numpy.linalg import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from pylab import *

fig = plt.figure()
ax = fig.add_subplot(111)

class particle():
	def __init__(self,x,y,vx,vy):
		self.r = np.array([x,y])
		self.v = np.array([vx,vy])

	def evol(self,dt):
		self.r=self.r+self.v*dt 

	def reflex(self,L):
		if self.r[0] <= 0 : self.v[0]=abs(self.v[0])
		if self.r[0] >= L : self.v[0]=-abs(self.v[0])
		if self.r[1] <= 0 : self.v[1]=abs(self.v[1])
		if self.r[1] >= L : self.v[1]=-abs(self.v[1])

def init():
	l=3
	l2=l*l
	L=10
	N=200
	dt=0.05
	tmax=100
	T=0.1
	SQRT=np.sqrt(T)
	x=np.random.uniform(0,L,N)
	y=np.random.uniform(0,L,N)
	vx=np.random.normal(0,SQRT,N)
	vy=np.random.normal(0,SQRT,N)
	T=(np.sum(vx*vx+vy*vy)/(2.*N))
	dt2=dt*dt
	mdt2=dt2/2.
	h=dt/2.
	itmax=int(tmax/dt)
	
	return L,N,dt,tmax,T,x,y,vx,vy,l,l2,dt2,mdt2,h,itmax

def animate():
	L,N,dt,tmax,T,x,y,vx,vy,l,l2,dt2,mdt2,h,itmax = init()
	part = list(particle(x[i],y[i],vx[i],vy[i]) for i in range(N))
	ax.set_autoscale_on(False)
	line, = ax.plot(x,y,'o')
	t=0
	for it in range(itmax):
		x=[]
		y=[]
		vx=[]
		vy=[]
		
		list(map(lambda i:i.evol(dt),part))
		list(map(lambda i:i.reflex(L),part))
		list(map(lambda i:x.append(i.r[0]),part))
		list(map(lambda i:y.append(i.r[1]),part))
		line.set_data(x,y)
		t=t+dt
		ax.set_title(t)
		fig.canvas.draw()
	return

win = fig.canvas.manager.window
fig.canvas.manager.window.after(100,animate)
plt.axis([0,10,0,10])
plt.show()
