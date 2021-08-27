import numpy as np

class particle():
    def __init__(self,x,y,vx,vy):
        self.x=x
        self.y=y
        self.vx=vx
        self.vy=vy


    def evol(self,dt):
        self.x+=self.vx*dt
        self.y+=self.vy*dt

    def contorno(self):
        if self.x<0:
            self.vx=abs(self.vx)
        if self.x>lx:
            self.vx=-1*abs(self.vx)
        if self.y<0:
            self.vy=abs(self.vy)
        if self.y>ly:
            self.vy=-1*abs(self.vy)
       
sigmav=10
N=100
dt=0.1
lx=100
ly=100
t=0
tmax=50
part=[]
for i in range (N):  
    part.append(particle(np.random.random()*lx, np.random.random()*ly, np.random.normal(0,sigmav),np.random.normal(0,sigmav)))
#i=0
while t<tmax:
    list(map(lambda i:i.evol(dt),part))
    list(map(lambda i:i.contorno(),part))
    t=t+dt
    #i=i+1
   # if i%10:
      #ec=list(reduce())      
