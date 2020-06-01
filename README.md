# About this script

This is simple python script to calculate several tensors appears in differential geometry.

## Ricci.py
`class Ricci` include modules to evaluate Christoffel Symbol, 
Riemann Tensor, Ricci Tensor, and Ricci Scalar for given metric. 
When create the instance of `Ricci`, metric and the index of local coordinate should be given as arguments.
The metric should be NxN `Matrix` of Sympy and index of local coordinate should be a list of `Symbol` of Sympy.

Christoffel Symbol, Ricci Tensor and Ricci Scalar has the form of `Array` in Sympy.
Here I comment on the order of index of `Array` and mathematical form.

A instance of `Ricci` is constructed under some metric g and coordinate. 
```python
Ric = Ricci(coordinate, g)
```
When Christoffel Symbol is defined by
```python
CF=Ric.ChristoffelSymbol()
```
$Gamma^i_{jk}$ can be referred by
```python
CF[i,j,k]
```

When Riemann Tensor is defined by
```python
Rie=Ric.RiemannTensor()
```
$R^i_{jkl}$ can be referred by

```python
Rie[i,j,k,l]
```

When Ricci Tensor is defined by
```python
RT=Ric.RicciTensor()
```
$Ric_{ij}$ can be referred by

```python
Ric[i,j]
```

The `__main__` is an example for Schwarzschild metric.