from init_conditions import *
from math import sin, pow, pi
from numba import njit
from plot_generator import plot, plot_map
import numpy as np 


@njit
def set_vertical_V_boarders(V):
    for i in range(ny+1):
        V[0][i] = sin(pi * i * delta / y_max)
        V[nx][i] = sin(pi * i * delta/ y_max)


@njit
def set_horizontal_V_boarders(V):
    for i in range(nx+1):
        V[i][ny] = -1.0 * sin(2.0 * pi * i * delta / x_max)
        V[i][0] = 1.0 * sin(2.0 * pi * i * delta / x_max)


@njit
def set_V_boarders(V):
    set_vertical_V_boarders(V)
    set_horizontal_V_boarders(V)


@njit
def first_step(V,k):
    for i in range(k, nx+1 - k, k):
        for j in range(k, ny+1 - k, k):
            V[i][j] = 0.25 * (V[i + k][j] + V[i - k][j] + V[i][j + k] + V[i][j - k])


@njit
def count_S(V, k):
    summary = 0.0
    for i in range(0, nx+1 -k, k):
        for j in range(0, ny+1 -k, k):
            summary += 0.5 * pow(k*delta, 2) \
                * (pow( (V[i+k][j] - V[i][j]) / (2*k*delta) + (V[i+k][j+k] - V[i][j+k]) / (2*k*delta), 2)\
				+ pow((V[i][j+k] - V[i][j]) / (2*k*delta) + (V[i+k][j+k] - V[i+k][j]) / (2*k*delta), 2))
    return summary


@njit
def upscale_V(V, k):
    
    k_half = int(k/2)

    for i in range(0, nx+1-k, k):
        for j in range(0, ny+1-k, k):
            V[i+k_half][j+k_half] = 0.25 * (V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k])
            if(i!=nx-k):
                V[i+k][j+k_half] = 0.5 * (V[i+k][j] + V[i+k][j+k])
            if(j!=ny-k):
                V[i+k_half][j+k] = 0.5 * (V[i][j+k] + V[i+k][j+k])
            if(j!=0):
                V[i+k_half][j] = 0.5 * (V[i][j] + V[i+k][j])
            if(i!=0):
                V[i][j+k_half] = 0.5 * (V[i][j] + V[i][j+k])
            

@njit
def fill_V_current(V, V_current, k):
    for i in range(0, nx+1 -k, k):
        for j in range(0, ny+1 -k, k):
            V_current[int(i/k)][int(j/k)] = V[i][j]


def prepare_data_and_plot_map(V, k):
    number_of_steps_x = int(nx/k +1)
    number_of_steps_y = int(ny/k +1)

    x_min = y_min = 0.0

    x = np.linspace(x_min, x_max, number_of_steps_x, endpoint=True)
    y = np.linspace(y_min, y_max, number_of_steps_y, endpoint=True)

    V_current = np.zeros((number_of_steps_x, number_of_steps_y))
    fill_V_current(V, V_current, k)

    plot_map(x,y,np.transpose(V_current),k,title='relaksacja', zlabel='V(x,y)')


def relaksacja_wielosiatkowa():
    
    k = 16
    it_sum = 0
    
    V = np.zeros((nx+1, ny+1))
    set_V_boarders(V)

    S_array = []
    k_array = []

    while(k>=1):
        S_prev = 0
        it = 1
        S_current_array = []

        while(True):
            
            first_step(V, k)
            S = count_S(V, k)
            
            if (it > 1 and abs((S-S_prev) / S_prev) < tol):
                break

            S_current_array.append(S)
            S_prev = S
            it+=1

        it_sum += it
        
        prepare_data_and_plot_map(V, k)
        upscale_V(V, k)
        
        S_array.append(S_current_array)
        k_array.append(k)
        
        k = int(k/2)

    plot(k_array, S_array, 'relaksacja', 'it', 'S(it)')
