from math import pow

S = 0.75
p = 2
t_max = 40.0
alfa = 5
tol_array = [1e-2, 1e-5]


def f(v):
    return v


def g(x, v, alfa):
    return alfa*(1-x*x)*v - x


def generate_array_of_arrays(n):
    return [[] for i in range(n)]
