from typing import List
from sympy import *
from sympy.abc import theta, phi
from sympy import Matrix

init_printing(use_unicode=True)


class Ricci:

## coordinate: list, metric Matrix
    def __init__(self, coordinate:List,  metric:Matrix):
        self.coordinate=coordinate
        self.metric=metric

    def inverseMetric(self):
        return self.metric.inv()

    def ChristoffelSymbol(self):
        dim = len(self.coordinate)
        inv=self.metric.inv()
        total=[]
        #Evaluate Gamma^i_{j,k}
        for i in range(dim):
            container_i=[]
            for j in range(dim):
                container_j=[]
                for k in range(dim):
                    symbol=0
                    for dummy in range(dim):
                        symbol=symbol+(inv[i,dummy]*(diff(self.metric[dummy,j],self.coordinate[k])) \
                                +inv[i,dummy]*(diff(self.metric[dummy,k],self.coordinate[j])) \
                                -inv[i,dummy]*(diff(self.metric[j,k],self.coordinate[dummy])))*1/2
                        #print(symbol)
                    container_j.append(simplify(symbol))
                container_i.append(container_j)
            total.append(container_i)
        return Array(total)

    def RiemannTensor(self):
        CF=self.ChristoffelSymbol()
        #print(CF[1,0,0])
        dim=len(self.coordinate)
        #print(self.coordinate)
        Rie=[]
        #R^i_{jkl}
        for i in range(dim):
            container_i=[]
            for j in range(dim):
                container_j=[]
                for k in range(dim):
                    container_k=[]
                    for l in range(dim):
                        Rijkl=diff(CF[i,l,j],self.coordinate[k])-diff(CF[i,k,j],self.coordinate[l])
                        #print(Rijkl)
                        for dummy in range(dim):
                            Rijkl=Rijkl+CF[dummy,l,j]*CF[i,k,dummy]-CF[dummy,k,j]*CF[i,l,dummy]
                        #print(Rijkl)
                        container_k.append(simplify(Rijkl))
                    #print(container_l)
                    container_j.append(container_k)
                container_i.append(container_j)
            Rie.append(container_i)

        return Array(Rie)

    def RicciTensor(self):
        Rie=self.RiemannTensor()
        dim=len(self.coordinate)
        Ric=[]
        for i in range(dim):
            container_i=[]
            for j in range(dim):
                #Ric_{ij}
                Ricij=0
                for dummy in range(dim):
                    Ricij=Ricij+Rie[dummy,i,dummy,j]
                container_i.append(simplify(Ricij))
            Ric.append(container_i)
        return Array(Ric)

    def RicciScalar(self):
        Ric=self.RicciTensor()
        dim=len(self.coordinate)
        R=0
        for i in range(dim):
            R=R+Ric[i,i]

        return(simplify(R))






if __name__=='__main__':
    #example Schwarzschild metric. In this metric, Ricci Tensor & Ricci Scalar should be exactly 0
    t=Symbol('t')
    r=Symbol('r')
    G=Symbol('G')
    M=Symbol('M')
    f=1-2*G*M/r
    print(simplify(f))
    g=Matrix([[-f,0,0,0],[0,1/f,0,0],[0,0,r*r,0],[0,0,0,r*r*sin(theta)*sin(theta)]])
    coordinate=[t,r,theta,phi]

    Ric = Ricci(coordinate, g)
    print(Ric.inverseMetric())
    CF=Ric.ChristoffelSymbol()
    print(CF)
    print(CF.shape)

    Rie=Ric.RiemannTensor()
    print(Rie)

    RT=Ric.RicciTensor()
    print(RT)

    RS=Ric.RicciScalar()
    print(RS)




