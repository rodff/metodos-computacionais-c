#Random walker with spatial restriction
import random as p

class walker():
    def __init__(self,ident,x):
        self.x = x
        self.newx=x
        self.ident=ident

    def mov(self):
        self.newx = self.newx + int(10*(p.random()-0.5))
        for i in wk:
            if i.ident != self.ident:
		if i.x == self.newx :
			break
	else:
		self.x=self.newx

# List of walkers objects
wk = list(walker(i,2*i) for i in range(10))
'''
Above line:
First arguments is the identity of the walker,
Second argument is the initial position of the walker.
'''

print(list(map(lambda i:i.x,wk)))

for i in range(100):
	
