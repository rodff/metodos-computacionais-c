from numpy import *
import matplotlib.pyplot as plt

tmax=10.
dt=0.05
dtg=dt
dtp=dt/10.
ntg=int(tmax/dtg)
print ntg
sqdtp=sqrt(dtp)
Wp=0
Wg=0
W=zeros(ntg+1)
Intp=zeros(ntg+1)
Intg=zeros(ntg+1)
Intep=0
t=arange(0,tmax+dt,dtg)
for j in range(ntg):
    dWg=0.
    for k in range(10):
        dWp=sqdtp*random.normal(0,1,1)
        dWg=dWg+dWp
        Intep=Intep+Wp*Wp*dWp+Wp*dtp
        Wp=Wp+dWp

    Intp[j+1]=Intep
    Intg[j+1]=Intg[j]+Wg*Wg*dWg+Wg*dtg
    Wg=Wg+dWg
    W[j+1]=Wg

Intito=W*W*W/3.
plt.plot(t, Intito,'r--', t, Intp, 'bs', t, Intg, 'g^')
plt.legend(("Exato","dtp","dtg"),'upper left',shadow=True, fancybox=True)
plt.show()
    
