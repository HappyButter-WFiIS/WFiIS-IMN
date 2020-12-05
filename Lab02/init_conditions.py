# init conditions
beta = 0.001
N = 500
gamma = 0.1
t_min = 0
t_max = 100
dt = 0.1
u_0 = 1
tol = 10e-6
mi_max = 21 # in range() stop value is omitted
alfa = beta*N - gamma

# function being used in all methods
def f(u):
    return alfa*u - beta*u*u
