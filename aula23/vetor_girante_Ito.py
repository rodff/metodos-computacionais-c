from numpy import *
import matplotlib.pyplot as plt

tmax=30.
dt=0.01
dtg=dt
beta=.5
w0=1.
nu=int(tmax/dt)
num=nu+1
srdt=sqrt(dt)
b2=.5*beta**2
w0dt=w0*dt
xi=zeros(num)
yi=zeros(num)
xe=zeros(num)
ye=zeros(num)
fase=zeros(num)
dW=srdt*random.normal(0,1,nu)
xe[0],xi[0]=1.,1.
fase[0]=0.
W=0.
for j in range(nu):
    W=W+dW[j]
    fase[j+1]=w0dt*j+beta*W
#    print beta*W
#  solucao por Ito
    xij=xi[j]
    yij=yi[j]
    xi[j+1]=xij-(w0*yij+b2*xij)*dt-yij*beta*dW[j]
    yi[j+1]=yij+(w0*xij-b2*yij)*dt+xij*beta*dW[j]
# solucao exata
    xe[j+1]=cos(fase[j+1])
    ye[j+1]=sin(fase[j+1])

t=arange(0,tmax+dt,dt)
#print t
#for i in range(nu):
#    print t[i],xi[i],xe[i]

#print len(xi)
#print len(xe)
#print len(t)
#plt.plot(t,xi)
#plt.show()
plt.plot(t,xi, 'r--', t, yi, 'b--', t,xe,t,ye)


#plt.subplot(311)
#plt.plot(t, xi, 'r--',t,xe)
#plt.title('Comparando x exato e  integrado por Ito')
#plt.subplot(312)
#plt.plot(t, yi, 'r--',t,ye)
#plt.title('Comparando y exato e  integrado por Ito')
#plt.subplot(313)
#plt.plot(t, fase)
#plt.title('Evolucao da fase')
#plt.xlabel('time (s)')
plt.show()
#print t
#print Intito
#print Intp
#print Intg
    
