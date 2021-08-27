import numpy as np
from numpy.linalg import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from pylab import *

fig = plt.figure()
ax = fig.add_subplot(111)

fname='time_energy.dat'

def init():
	N=64
	L=np.sqrt(N)
	dt=0.005
	tmax=100
	T=1.0
	SQRT=np.sqrt(T)
	x=np.arange(N)%L+0.5    # Nao pode ser totalmente aleatorio
	y=np.arange(N)/L+0.5
	vx = np.random.normal(0,SQRT,N)
	vy = np.random.normal(0,SQRT,N)
	T=(np.sum(vx*vx+vy*vy))/(2.*N)
	itmax=int(tmax/dt)

	return L,N,dt,tmax,T,x,y,vx,vy,itmax

def animate():
	class particle():
		def __init__(self,x,y,vx,vy,ident,N):
			self.r=np.array([x,y])
			self.v=np.array([vx,vy])
			self.ident = ident
			self.Force = np.zeros(2)    # 2 Because it its a 2-dimensional system
			self.Ec = 0.0

		def my_Ec(self):
			self.Ec=np.linalg.norm(self.v)**2/2

		def evol0(self,dt):
			self.Force=self.forces_between_particles(N)
			self.r=self.r+self.v*dt+self.Force*dt**2/2
			self.v = self.v+self.Force*dt/2

		def evol1(self,dt):
			self.evol0(dt)
			self.Force=self.forces_between_particles(N)
			self.v = self.v+self.Force*dt/2

		def reflex(self,L):
			if self.r[0] <= 0 : self.v[0]=abs(self.v[0])          
			if self.r[0] >= L : self.v[0]=-abs(self.v[0])
			if self.r[1] <= 0 : self.v[1]=abs(self.v[1])
			if self.r[1] >= L : self.v[1]=-abs(self.v[1])

		def force(self,dr):
			dr2=np.linalg.norm(dr)**2
			f=12*(1/dr2**7-1/dr2**4)*dr     # Lennard-Jones
                        return f

		def forces_between_particles(self,N):
			self.Force = np.zeros(2)
			for i in range(N):
				if(i != self.ident):
					dr = self.r-part[i].r
					self.Force+=self.force(dr)
			return self.Force

	L,N,dt,tmax,T,x,y,vx,vy,itmax = init()
	part = list(particle(x[i],y[i],vx[i],vy[i],i,N) for i in range(N))
	
	ax.set_autoscale_on(False)
	line, = ax.plot(x,y,'o')
	t=0

        global fname
        arq=open(fname,'w')

	for it in range(itmax):
		x=[]
		y=[]
		vx=[]
		vy=[]
		list(map(lambda i:i.evol1(dt), part))
		list(map(lambda i:i.reflex(L), part))
		t=t+dt
		if(it%5==0):
			list(map(lambda i:x.append(i.r[0]), part))
			list(map(lambda i:y.append(i.r[1]), part))
			list(map(lambda i:i.my_Ec(), part))
			Ec = sum(list(map(lambda i:i.Ec, part)))
			V=0
			for i in range(N):
				for j in range(i+1,N):
					dr2 = np.linalg.norm(part[i].r-part[j].r)**2
					V+=(1/dr2**6-2/dr2**3)
                        
                        arq.write('{0}  {1} {2} {3}\n'.format(t,Ec,V,(Ec+V)))
                        
                        line.set_data(x,y)
			ax.set_title("Ec=%f V=%f Et=%f t=%f"%(Ec,V,Ec+V,t))
			fig.canvas.draw()

        arq.close()

	return

L,N,dt,tmax,T,x,y,vx,vy,itmax = init()
win = fig.canvas.manager.window
fig.canvas.manager.window.after(1,animate)
ax.axis([0,L,0,L])
plt.show() 

t_l,E_l,Ec_l,V_l = [], [], [], []
arq=open('time_energy.dat','r')
data = arq.readlines()
for l in data:
    t_l.append(float(l.split()[0]))
    Ec_l.append(float(l.split()[1]))
    V_l.append(float(l.split()[2]))
    E_l.append(float(l.split()[3]))
arq.close()

fig2 = plt.figure()  
ax2 = fig2.add_subplot(211)
ax3 = fig2.add_subplot(212)
ax2.plot(t_l,E_l)
ax2.axis([0,max(t_l),(0.9*mean(E_l)),(1.1*mean(E_l))])
plt.ylabel('Energy')
plt.xlabel('Time')
ax3.plot(t_l,abs(np.array(Ec_l)), label='|K|')
ax3.plot(t_l,abs(np.array(V_l)), label='|V|')
ax3.legend()
fig2.savefig('time_energy.png')


