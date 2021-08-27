# N particulas em uma caixa com paredes refletoras de dim. L*L
# com potencial entre pares tipo Lennard-Jones
import numpy as np
from numpy.linalg import *
import matplotlib
matplotlib.use('TkAgg') # do this before importing pylab
import matplotlib.pyplot as plt
from pylab import *

fig = plt.figure()
ax = fig.add_subplot(111)
        
def init():
    global N,L,lbox,nb,dt,nb2,t,Tau,R0,v0,mu,eta,nse_lim
    N=36
    L=np.int64(ceil(np.sqrt(N)))
    lbox = 1
    nb = np.int64(ceil(L/lbox))
    nb2=nb**2
    tmax=100
    T=1.0
    mu = 1.0    # Mobility 
    v0 = 1.0    # Self-propelling velocity magnitude
    Tau = 1.0   # Relaxation time
    eta = 0.6   # Noise
    R0 = 1.0               
    dt = 0.05*R0/v0
    nse_lim = eta/(2*np.sqrt(dt))    

    x = np.arange(N)%L+0.5
    y = np.int64(np.arange(N)/L)+0.5
    vx = np.full(N,v0)
    vy = np.full(N,v0)
    T=(np.sum(vx*vx+vy*vy))/2./N
    itmax=int(tmax/dt)

    return tmax,T,x,y,vx,vy,itmax


def animate():
    class particle():
        def __init__(self,x,y,vx,vy,ident,N):
            theta0 = np.random.uniform(0,2*np.pi)  
            self.r = np.array([x,y])
            self.v = np.array([vx*np.cos(theta0),vy*np.sin(theta0)])
            self.theta = theta0
            self.n = np.array([np.cos(theta0),np.sin(theta0)])
            self.ident = ident
            self.Force = np.zeros(2)
            self.Ec = 0.0
            self.Mybox = int(self.r[0]/lbox)+nb*int(self.r[1]/lbox)

        def mybox(self): #Each particle calculates the box it is in
            if math.isnan(self.r[0]):
                print("Found Nan in r[0], t=%f, particle_id=%d"%(t,self.ident))
                exit()
            if math.isnan(self.r[1]):
                print("Found Nan in r[1], t=%f, particle_id=%d"%(t,self.ident))
                exit()
                     
            j=int(self.r[0]/lbox)+nb*int(self.r[1]/lbox)
            return int(j)
            
        def my_Ec(self):
            self.Ec=np.linalg.norm(self.v)**2/2
            
        def evol0(self,dt):
            self.r=self.r+self.v*dt
            #self.v = v0*self.n+mu*self.Force*dt/2           

        def evol1(self,dt):
            nse = np.random.uniform(-nse_lim,nse_lim)
            v_ort = np.array([-self.v[1],self.v[0]])
            sarg=np.dot(self.n*self.v,v_ort/np.linalg.norm(v_ort))
            dtheta = np.arcsin(sarg)/Tau+(nse*np.sqrt(dt))
            self.theta = self.theta+nse*dtheta
            self.n[0] = np.cos(self.theta) 
            self.n[1] = np.sin(self.theta)
            self.v = v0*self.n+mu*self.Force*dt/2

        def reflex(self,L):
            h=0.05
            if self.r[0] <= h: self.v[0]=abs(self.v[0])
            if self.r[0] >= L-h : self.v[0]=-abs(self.v[0])
            if self.r[1] <= h : self.v[1]=abs(self.v[1])
            if self.r[1] >= L-h : self.v[1]=-abs(self.v[1])

        def reflex_force(self,L):
            Fwall = 0.5
            if (self.r[0] < R0): 
                self.Force+=-Fwall*np.exp(-2*(np.abs(self.r[0]))/R0)*np.array([1,0])
            if(self.r[0] > L-R0): 
                self.Force+=-Fwall*np.exp(-2*(np.abs(L-self.r[0]))/R0)*np.array([-1,0])
            if(self.r[1] < R0): 
                self.Force+=-Fwall*np.exp(-2*(np.abs(self.r[1]))/R0)*np.array([0,1])
            if(self.r[1] > L-R0): 
                self.Force+=-Fwall*np.exp(-2*(np.abs(L-self.r[1]))/R0)*np.array([0,-1])
            self.v+=self.Force

        def force(self,dr):
            Frep = 30.
            Fadh = 0.75
            Req = 5/6.

            dij=np.linalg.norm(dr)
            eij=dr/dij
            
            if(dij<Req):
                f=Frep*(dij-Req)/Req
            elif(dij>R0):
                f=0.
            else:
                f=Fadh*(dij-Req)/(R0-Req)

            return eij*f            

        def forces_between_particles(self):
            for i in box[self.Mybox].mylist: #forces inside mybox
                if(self.ident!=i):  
                    dr=part[i].r-self.r
                    self.Force+=self.force(dr)
            for i in box[self.Mybox].neighboxlist: #forces in neighboring boxes
                dr=part[i].r-self.r
                self.Force+=self.force(dr)
                part[i].Force-=self.force(dr)

        def zero_force(self):
            self.Force=np.array([0.,0.])


        def changebox(self):
            newbox=self.mybox()
            if(newbox!=self.Mybox): #verify particle box change 
                box[self.Mybox].mylist.remove(self.ident) #this is safe since particles ident is unique
                box[newbox].mylist.append(self.ident)
                self.Mybox=newbox
            return self.Mybox


    class boite:
        def __init__(self,index):
            self.index = index
            self.mylist = []
            self.neighboxlist = []

        def neighbor_list(self):
            self.neighboxlist=[]
            #for reflective walls
            zz1=np.int64(self.index/nb)*nb+(self.index+1)%nb  #right box
            if zz1%nb != 0 : 
                self.neighboxlist.extend(box[zz1].mylist)
            zz2=np.int64((self.index+nb)/nb)*nb+(self.index+nb-1)%nb #lower left box
            if self.index%nb != 0 and zz2 < nb2 :
                self.neighboxlist.extend(box[zz2].mylist)
            zz3=(np.int64((self.index+nb)/nb)*nb+(self.index+nb)%nb)%nb2 #lower center box
            if zz3 < nb2 :
                self.neighboxlist.extend(box[zz3].mylist)
            zz4=np.int64((self.index+nb)/nb)*nb+(self.index+nb+1)%nb #lower right box
            if zz4%nb != 0 and zz4 < nb2 :
                self.neighboxlist.extend(box[zz4].mylist)
            '''
            #for periodic boundary conditions
            zz1=(np.int64(self.index/nb)*nb+(self.index+1)%nb)%nb2
            zz2=(np.int64((self.index+nb)/nb)*nb+(self.index+nb-1)%nb)%nb2
            zz3=(np.int64((self.index+nb)/nb)*nb+(self.index+nb)%nb)%nb2
            zz4=(np.int64((self.index+nb)/nb)*nb+(self.index+nb+1)%nb)%nb2
            self.neighboxlist.extend(box[zz1].mylist)
            self.neighboxlist.extend(box[zz2].mylist)
            self.neighboxlist.extend(box[zz3].mylist)
            self.neighboxlist.extend(box[zz4].mylist)
            ''' 
        
    tmax,T,x,y,vx,vy,itmax = init()
    part = list(particle(x[i],y[i],vx[i],vy[i],i,N) for i in range(N))
    box=list(boite(i) for i in range(nb2))
    
    #construct list of particles in each box

    for i in range(N):
        part[i].Mybox=part[i].mybox()
        box[np.int64(part[i].Mybox)].mylist.append(i)
    # Construct list of particles in neighboring boxes

    list(map(lambda i:i.neighbor_list(), box))

    ax.set_autoscale_on(False)
    line, = ax.plot(x,y,'o')
    t=0
    for it in range(itmax):
        x=[]
        y=[]
        vx=[]
        vy=[]
        list(map(lambda i:i.evol0(dt), part))
        #Reset forces to zero
        list(map(lambda i:i.zero_force(), part))
        list(map(lambda i:i.forces_between_particles(), part))
        list(map(lambda i:i.evol1(dt), part))
        list(map(lambda i:i.reflex(L), part))
        #Find the newboxes 
        list(map(lambda i:i.changebox(), part))
        #Construct the list of particles in neighboring boxes
        list(map(lambda i:i.neighbor_list(), box))
        #Reset forces to zero
        list(map(lambda i:i.zero_force(), part))

        t=t+dt
        if it%1 == 0: 
            list(map(lambda i:x.append(i.r[0]), part))
            list(map(lambda i:y.append(i.r[1]), part))
            list(map(lambda i:i.my_Ec(), part))
#            dr=np.linalg.norm(part[0].r-part[1].r)
            Ec = sum(list(map(lambda i:i.Ec, part)))
            for i in range(N):
                for j in range(i+1,N):
                    dr2=np.linalg.norm(part[i].r-part[j].r)**2
            line.set_data(x,y)
            ax.set_title("Ec=%f t=%d "%(Ec,int(t)))
            fig.canvas.draw()
    return

tmax,T,x,y,vx,vy,itmax = init()
win = fig.canvas.manager.window
fig.canvas.manager.window.after(1, animate)
plt.axis([0,L,0,L])
plt.show()

