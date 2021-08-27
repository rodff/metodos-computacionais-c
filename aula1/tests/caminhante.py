import random as p
import matplotlib.pyplot as m

class part():   # Particula
    x = 0.0 
    def mov(zz):
        zz.x=zz.x+(p.random()-0.5)

ob1=part()
ob2=part()

ob = list(part() for i in range(1000))

for i in range(100):
    for j in range(1000):
            ob[j].mov()

a = map(lambda y:y.x, ob)
b = m.hist(a)

for i in range(10000): 
    map(lambda x:x.mov(),ob)

c = map(lambda z:z.x, ob)
d = m.hist(c,bins=33)

m.show()
