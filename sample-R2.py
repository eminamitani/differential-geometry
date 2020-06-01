import Ricci
from sympy import *
from sympy.abc import phi

r=Symbol('r')
coordinate=[r,phi]
g=Matrix([[1,0],[0,r*r]])

R2=Ricci.Ricci(coordinate,g)
CF = R2.ChristoffelSymbol()
print(CF)