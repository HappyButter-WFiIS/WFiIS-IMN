import numpy as np
from numba import njit
from math import pow, exp 

epsilon = 1.0
delta = 0.1
nx = 150
ny = 100
V1 = 10
V2 = 0
x_max = delta*nx
y_max = delta*ny
sigma_x = 0.1*x_max
sigma_y = 0.1*y_max
tol = 1e-8

@njit
def get_rho1(x,y):
    return 2.0 * exp( (-pow(x-0.35*x_max, 2)/pow(sigma_x,2) ) - (pow(y-0.25*y_max, 2)/pow(sigma_y,2)) )

@njit
def get_rho2(x,y):
    return -2.0 * exp( (-pow(x-0.65*x_max, 2)/pow(sigma_x,2) ) - (pow(y-0.25*y_max, 2)/pow(sigma_y,2)) )

@njit
def get_rho3(x,y):
    return -1.0 * exp( (-pow(x-0.35*x_max, 2)/pow(sigma_x,2) ) - (pow(y-0.75*y_max, 2)/pow(sigma_y,2)) )

@njit
def get_rho4(x,y):
    return 1.0 * exp( (-pow(x-0.65*x_max, 2)/pow(sigma_x,2) ) - (pow(y-0.75*y_max, 2)/pow(sigma_y,2)) )


@njit
def get_rho(x,y):
    return get_rho1(x,y) + get_rho2(x,y) + get_rho3(x,y) + get_rho4(x,y)
