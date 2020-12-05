from init_conditions import *
from plot_generator import plot, plot_map

omega_L_array = [1.0, 1.4, 1.8, 1.9]

@njit
def first_step(V, rho, omega_L):
    for i in range(1, nx):
        for j in range(1, ny):
            V[i][j] = (1.0-omega_L)*V[i][j] + \
                (omega_L/4.0)*(V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + pow(delta,2)/epsilon*rho[i][j])

@njit
def WB_von_Neumann(V):
    for j in range(1, ny+1):
        V[0][j] = V[1][j]
        V[nx][j] = V[nx-1][j]


@njit
def count_S(V, rho):
    S_sum = 0.0
    for i in range(nx):
        for j in range(ny):
            S_sum += pow(delta, 2) * (0.5 * pow((V[i+1][j] - V[i][j])/delta,2) + 0.5 * pow((V[i][j+1] - V[i][j])/delta,2) - rho[i][j]*V[i][j])
    return S_sum


def generate_array_of_arrays(n):
    return [[] for i in range(n)]


def  local_relaxation():
    
    S = generate_array_of_arrays(len(omega_L_array))
    it_array = generate_array_of_arrays(len(omega_L_array))
    for it_omega, omega_L in enumerate(omega_L_array):
        
        number_of_steps_x = nx + 1 
        number_of_steps_y = ny + 1 
        x_min = y_min = 0.0

        x = np.linspace(x_min, x_max, number_of_steps_x, endpoint=True)
        y = np.linspace(y_min, y_max, number_of_steps_y, endpoint=True)

        rho = np.zeros((nx+1, ny+1))
        V = np.zeros((nx+1, ny+1))

        # filling rho matrix and boundary conditions to V_* matrixes
        for i, x_current in enumerate(x):
            
            V[i][0] = V1
            
            for j, y_current in enumerate(y):
                rho[i][j] =  get_rho(x_current, y_current)

        it = 0
        
        while(True):

            first_step(V, rho, omega_L)
            WB_von_Neumann(V)
                
            S_sum = count_S(V, rho)
            S[it_omega].append(S_sum)
            it_array[it_omega].append(it)
            
            if(it > 0 and abs((S[it_omega][it]-S[it_omega][it-1])/S[it_omega][it-1]) < tol):
                break

            it+=1 

    plot(it_array, S, omega_L_array, title='local', xlabel='it', ylabel='S(it)')