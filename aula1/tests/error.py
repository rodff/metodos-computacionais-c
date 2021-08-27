import math as mt

a = -1

try:
    k = mt.log10(a)
except ValueError:
    print("Argumento deve ser positivo!")

