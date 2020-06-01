import Ricci
from sympy import *
from sympy.abc import phi,theta

coordinate=[theta,phi]
g=Matrix([[1,0],[0,sin(theta)*sin(theta)]])

S2=Ricci.Ricci(coordinate,g)
CF = S2.ChristoffelSymbol()
print(CF)