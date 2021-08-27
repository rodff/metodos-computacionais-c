import numpy as np
from numpy.linalg import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from pylab import *

fig=plt.figure()
ax=fig.add_subplot(111)

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
    T=(np.sum(vx*vx+vy*vy))/2./N
    dt2=dt*dt
    mdt2=dt2/2.
    h=dt/2.
    itmax=int(tmax/dt)
    return L,N,dt,tmax,T,x,y,vx,vy,l,l2,dt2,mdt2,h,itmax

def reflete(x,y,vx,vy,L):
    vx[x<=0.]=abs(vx[x<=0.])
    vx[x>=L]=-abs(vx[x>=L])
    vy[y<=0.]=abs(vy[y<=0.])
    vy[y>=L]=-abs(vy[y>=L])
    return vx,vy
           
def integra(x,y,vx,vy,dt,dt2,mdt2,h,N,l2,L):
    for i in range(N):
        x[i]=x[i]+vx[i]*dt
        y[i]=y[i]+vy[i]*dt
    vx,vy=reflete(x,y,vx,vy,L)
    return x,y,vx,vy

def animate():
    L,N,dt,tmax,T,x,y,vx,vy,l,l2,dt2,mdt2,h,itmax=init()
    ax.set_autoscale_on(False)
    line,=ax.plot(x,y,'o')
    t=0
    for it in range(itmax):
        x,y,vx,vy=integra(x,y,vx,vy,dt,dt2,mdt2,h,N,l2,L)
        t=t+dt
        if it%1==0:
            line.set_data(x,y)
            Ec=(sum(list(map(lambda i :i**2,vx)))+sum(list(map(lambda i:i**2, vy))))/2.
            ax.set_title("t=%d Ec=%f"%(int(t),Ec))
            fig.canvas.draw()
            
win=fig.canvas.manager.window
fig.canvas.manager.window.after(100,animate)
plt.axis([0,10,0,10])
plt.show()
