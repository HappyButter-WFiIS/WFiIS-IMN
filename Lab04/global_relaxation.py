from init_conditions import *
from plot_generator import plot, plot_map

omega_G_array = [0.6, 1.0]

@njit
def first_step(V_n, V_s, rho):
    for i in range(1, nx):
        for j in range(1, ny):
            V_n[i][j] = 0.25 * (V_s[i+1][j] + V_s[i-1][j] + V_s[i][j+1] + V_s[i][j-1] + (pow(delta,2)*rho[i][j])/epsilon)


@njit
def second_step(V_n):
    for j in range(1, ny+1):
        V_n[0][j] = V_n[1][j]
        V_n[nx][j] = V_n[nx-1][j]


@njit
def third_step(V_n, V_s, omega_G):
    for i in range(nx+1):
        for j in range(ny+1):
            V_s[i][j] = (1.0 - omega_G)*V_s[i][j] + omega_G*V_n[i][j]


@njit
def count_S(V_n, rho):
    S_sum = 0.0
    for i in range(nx):
        for j in range(ny):
            S_sum += pow(delta, 2) * (0.5 * pow((V_n[i+1][j] - V_n[i][j])/delta,2) + 0.5 * pow((V_n[i][j+1] - V_n[i][j])/delta,2) - rho[i][j]*V_n[i][j])
    return S_sum


@njit
def count_err(V_n, V_err, rho):
    for i in range(1,nx):
        for j in range(1,ny):
            V_err[i][j] = ((V_n[i+1][j] - 2.0*V_n[i][j] + V_n[i-1][j])/pow(delta,2) \
            + (V_n[i][j+1] - 2.0*V_n[i][j] + V_n[i][j-1])/pow(delta,2)) \
            + rho[i][j]/epsilon


def  global_relaxation():
    
    S = [[],[]]
    it_array = [[],[]]
    for it_omega, omega_G in enumerate(omega_G_array):
        
        number_of_steps_x = nx + 1 
        number_of_steps_y = ny + 1 
        x_min = y_min = 0.0

        x = np.linspace(x_min, x_max, number_of_steps_x, endpoint=True)
        y = np.linspace(y_min, y_max, number_of_steps_y, endpoint=True)

        rho = np.zeros((nx+1, ny+1))
        V_s = np.zeros((nx+1, ny+1))
        V_n = np.zeros((nx+1, ny+1))
        V_err = np.zeros((nx+1, ny+1))

        # filling rho matrix and boundary conditions to V_* matrixes
        for i, x_current in enumerate(x):
            
            V_s[i][0] = V1
            V_n[i][0] = V1
            
            for j, y_current in enumerate(y):
                rho[i][j] =  get_rho(x_current, y_current)

        it = 0
        
        while(True):

            first_step(V_n, V_s, rho)
            second_step(V_n)
            third_step(V_n, V_s, omega_G)
                
            S_sum = count_S(V_n, rho)
            S[it_omega].append(S_sum)
            it_array[it_omega].append(it)
            
            if(it > 0 and abs((S[it_omega][it]-S[it_omega][it-1])/S[it_omega][it-1]) < tol):
                break

            it+=1 

        count_err(V_n, V_err, rho)
        plot_map(x,y,np.transpose(V_err), omega_G,title='global', zlabel='Error')
        plot_map(x,y,np.transpose(V_n), omega_G,title='global', zlabel='V(x,y)')
    plot(it_array,S,omega_G_array, title='global',xlabel='it',ylabel='S(it)')