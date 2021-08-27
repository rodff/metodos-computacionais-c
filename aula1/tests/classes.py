class MyClass:
    i = 2
    j = 3
    def FazIsso(self,coisa):
        self.k = self.i + self.j + coisa

MyObject = MyClass()

MyObject.FazIsso(1) 

tst = MyObject.k

print(tst)
