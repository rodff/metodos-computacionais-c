import numpy as np
from numpy.linalg import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from pylab import *

fig=plt.figure()
ax=fig.add_subplot(111)

def init():
    N=49
    L=np.sqrt(N)
    dt=0.001
    tmax=100
    T=1.
    SQRT=np.sqrt(T)
    x=np.arange(N)%L+0.5#nao pode ser aleatorio
    y=np.arange(N)/L+0.5
    vx=np.random.normal(0,SQRT,N)
    vy=np.random.normal(0,SQRT,N)
    T=(np.sum(vx*vx+vy*vy))/2./N
    itmax=int(tmax/dt)
    return L,N,dt,tmax,T,x,y,vx,vy,itmax

def animate():
    class particle():
        def __init__(self,x,y,vx,vy,ident,N):
            self.r=np.array([x,y])
            self.v=np.array([vx,vy])
            self.ident=ident
            self.force=np.zeros(2)
            self.Ec=0.0
            #continua aq

        def my_Ec(self,dt):
            self.ec=np.linalg.norm(self.v)**2/2

        def evol(self,dt):
           self.force=self.force_b_p(N) 
           self.r=self.r+self.v*dt+self.force*dt**2/2
           self.v=self.v+self.force*dt/2
           self.force=self.force_b_p(N)
           self.v=self.v+self.force*dt/2
           
        def reflete(self,x,y,vx,vy,L):
           if self.r[0]<=0:self.v[0]=abs(self.v[0])
           if self.r[0]>=L:self.v[0]=-abs(self.v[0])
           if self.r[1]<=0:self.v[1]=abs(self.v[1])
           if self.r[1]>=L:self.v[1]=-abs(self.v[1])
           return vx,vy
       #...?
       
        def force_b_p(self,N):
            self.force=np.zeros(2)
            for i in range(N):
                if i!= self.ident:
                    dr=self.r-part[i].r
                    self.force=self.force+self.f(dr)
            return self.force

L,N,dt,tmax,T,x,y,vx,vy,itmax=init()
part=list(particle(x[i],y[i],vx[i],vy[i],i,N) for i in range(N))

ax.set_autoscale_on(False)
line,=ax.plot(x,y,'o')
t=0

for it in range(itmax):
    x=[];y=[];vx=[];vy=[]
    list(map(lambda i:i.evol(dt),part))
    list(map(lambda i:i.reflete(L),part))
    t=t+dt
    if it%10==0:
        
      list(map(lambda i:x.append(i.r[0]),part))
      list(map(lambda i:y.append(i.r[1]),part))
      list(map(lambda i:i.my_Ec(),part))
      Ec=(sum(list(map(lambda i :i**2,vx)))+sum(list(map(lambda i:i**2, vy))))/2.
      V=0
      for i in range(N):
          for j in range (i+1,N):
              dr2=np.linalg.norm(part[i].r-part[j].r)**2
              V=V+(1/dr2**6-2/dr2**3)
      line.set_data(x,y)        
      ax.set_title("Ec=%f V=%f Et=%f t=%d"%(Ec,V,Ec+V,int(t)))
      fig.canvas.draw()
      #return

L,N,dt,tmax,T,x,y,vx,vy,itmax=init()
win=fig.canvas.manager.window
fig.canvas.manager.window.after(1,animate)
plt.axis([0,L,0,L])
plt.show()

