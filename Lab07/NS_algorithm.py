from WB import *
from chart_generator import plot_map, plot_contour_map
omega = 0

@njit
def calculate_current_psi(i, j, psi, zeta):
    return 0.25 * (psi[i+1][j] + psi[i-1][j] \
        + psi[i][j+1] + psi[i][j-1] \
        - delta*delta*zeta[i][j])

@njit
def calculate_current_zeta(i, j, psi, zeta, omega):
    return 0.25 * (zeta[i+1][j] + zeta[i-1][j] \
        + zeta[i][j+1] + zeta[i][j-1]) \
        - (omega*rho/(16.0*mi)) \
        * ( (psi[i][j+1] - psi[i][j-1]) \
        * (zeta[i+1][j] - zeta[i-1][j]) \
        - (psi[i+1][j] - psi[i-1][j]) \
        * (zeta[i][j+1] - zeta[i][j-1]))

@njit
def calculate_main_NS(psi, zeta, omega):
    for i in range(1,nx):
        for j in range(1,ny):
            if( (i <= i_1 and j > j_1) or (i > i_1) ):
                psi[i][j] = calculate_current_psi(i, j, psi, zeta)
                zeta[i][j] = calculate_current_zeta(i, j, psi, zeta, omega)

@njit
def err_control(psi,zeta):
    gamma = 0.0
    j_2 = j_1 + 2
    for i in range(1,nx):
        gamma += (psi[i+1][j_2] + psi[i-1][j_2] \
            + psi[i][j_2+1] + psi[i][j_2-1] \
            - 4.0 * psi[i][j_2] \
            - pow(delta,2)*zeta[i][j_2])
    return gamma

@njit
def calculate_u_v(u, v, psi):
    for i in range(1, nx):
        for j in range(1, ny):
            if (i > i_1 or j > j_1):
                u[i][j] = (psi[i][j+1] - psi[i][j-1]) \
                    / (2.0*delta)
                v[i][j] = -(psi[i+1][j] - psi[i-1][j]) \
                    / (2.0*delta)



def calculate_NS(Qwe):
    x = np.linspace(0.0, (nx+1)*delta, nx+1, endpoint=True)
    y = np.linspace(0.0, (ny+1)*delta, ny+1, endpoint=True)

    psi = np.zeros((nx+1,ny+1))
    zeta = np.zeros((nx+1,ny+1))

    psi[0:i_1,0:j_1] = np.nan
    zeta[0:i_1,0:j_1] = np.nan

    set_WB_psi(psi, Qwe, y)
    # set_WB_zeta(zeta, psi, Qwe, y)

    for it in range(1,IT_MAX+1):
        if(it < 2000):
            omega = 0
        else:
            omega = 1
        
        calculate_main_NS(psi, zeta, omega)
        set_WB_zeta(zeta, psi, Qwe, y)
        # print('it '+ str(it) + ': ' + str(err_control(psi,zeta)))

    u = np.zeros((nx + 1, ny + 1))
    v = np.zeros((nx + 1, ny + 1))

    calculate_u_v(u, v, psi)

    if(Qwe != 4000):
        plot_map(x, y, np.transpose(u), Qwe, 'u(x,y)')
        plot_map(x, y, np.transpose(v), Qwe, 'v(x,y)')

        plot_contour_map(x, y, np.transpose(zeta), Qwe, 'zeta(x,y)')
    
    plot_contour_map(x, y, np.transpose(psi), Qwe, 'psi(x,y)')