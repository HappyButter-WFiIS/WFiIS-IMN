from init_conditions import *


@njit
def get_Qwy(Qwe, y):
    return Qwe * (pow(y[ny],3) - pow(y[j_1],3) \
        - 3.0*y[j_1]*pow(y[ny],2) \
        + 3.0*pow(y[j_1],2) * y[ny]) \
        / pow(y[ny],3)

# --------------- psi ------------------ #

@njit
def set_WB_psi_A(psi, Qwe, y):          # wejscie
    for j in range(j_1, ny+1):
        psi[0][j] = ( Qwe*(pow(y[j],3)/3.0 \
            - pow(y[j],2)*0.5*(y[j_1]+y[ny]) \
            + y[j]*y[j_1]*y[ny]) ) \
            / (2.0*mi)

@njit                           
def set_WB_psi_C(psi, Qwe, y):
    for j in range(0, ny+1):          #wyjscie
        psi[nx][j] = get_Qwy(Qwe, y) * (pow(y[j],3)/3.0 \
            - pow(y[j],2)*0.5*y[ny]) / (2.0*mi) \
            + Qwe * pow(y[j_1],2)*(-y[j_1]+3.0*y[ny]) \
            / (12.0*mi)  

@njit
def set_WB_psi_B(psi):
    for i in range(1,nx):
        psi[i][ny] = psi[0][ny]

@njit
def set_WB_psi_D(psi):
    for i in range(i_1, nx):
        psi[i][0] = psi[0][j_1]

@njit
def set_WB_psi_E(psi):
    for j in range(1, j_1+1):
        psi[i_1][j] = psi[0][j_1]

@njit
def set_WB_psi_F(psi):
    for i in range(1, i_1+1):
        psi[i][j_1] = psi[0][j_1]

@njit
def set_WB_psi(psi, Qwe, y):
    set_WB_psi_A(psi, Qwe, y)
    set_WB_psi_C(psi, Qwe, y)
    set_WB_psi_B(psi)
    set_WB_psi_D(psi)
    set_WB_psi_E(psi)
    set_WB_psi_F(psi)

# --------------- zeta ------------------ #

@njit
def set_WB_zeta_A(zeta, Qwe, y):
    for j in range(j_1,ny+1):
        zeta[0][j] = Qwe * (2.0*y[j]-y[j_1]-y[ny]) \
            / (2.0*mi)

@njit
def set_WB_zeta_C(zeta, Qwe, y):
    for j in range(0,ny+1):
        zeta[nx][j] = get_Qwy(Qwe, y) * (2.0*y[j]-y[ny]) \
            / (2.0*mi)

@njit
def set_WB_zeta_B(zeta, psi):
    for i in range(1, nx):
        zeta[i][ny] = (2.0/pow(delta,2)) \
            * (psi[i][ny-1] - psi[i][ny])

@njit
def set_WB_zeta_D(zeta, psi):
    for i in range(i_1+1, nx):
        zeta[i][0] = (2.0/pow(delta,2)) \
            * (psi[i][1] - psi[i][0])

@njit
def set_WB_zeta_E(zeta, psi):
    for j in range(1, j_1):
        zeta[i_1][j] = (2.0/pow(delta,2)) \
            * (psi[i_1+1][j] - psi[i_1][j])

@njit
def set_WB_zeta_F(zeta, psi):
    for i in range(1, i_1+1):
        zeta[i][j_1] = (2.0/pow(delta,2)) \
            * (psi[i][j_1+1] - psi[i][j_1])


@njit
def set_WB_zeta(zeta, psi, Qwe, y):
    set_WB_zeta_A(zeta, Qwe, y)
    set_WB_zeta_C(zeta, Qwe, y)
    set_WB_zeta_B(zeta, psi)
    set_WB_zeta_D(zeta, psi)
    set_WB_zeta_E(zeta, psi)
    set_WB_zeta_F(zeta, psi)

    zeta[i_1][j_1] = 0.5*(zeta[i_1-1][j_1] + zeta[i_1][j_1-1])




