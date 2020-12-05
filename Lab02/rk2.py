from plot_generator import plot, np
from init_conditions import *


# Butcher's array
a = [[0.25, 0.25 - (np.sqrt(3)/6.0)], \
    [0.25 + (np.sqrt(3)/6.0), 0.25]]
b = [0.5, 0.5]

# help function
def F1(U1,U2,un):
    return U1 - un - dt*( a[0][0]*f(U1) + a[0][1]*f(U2) )


def F2(U1,U2,un):
    return U2 - un - dt*( a[1][0]*f(U1) + a[1][1]*f(U2) )


def calculate_m(U1, U2):
    m00 = 1 - dt*a[0][0]*(alfa-2*beta*U1)
    m01 = -dt*a[0][1]*(alfa-2*beta*U2)
    m10 = -dt*a[1][0]*(alfa-2*beta*U1)
    m11 = 1 - dt*a[1][1]*(alfa-2*beta*U2)

    return [[m00, m01],\
            [m10, m11]]


def calculate_delata_U1(U1,U2,un):
    m = calculate_m(U1, U2)
    return (F2(U1,U2,un) * m[0][1] - F1(U1,U2,un) * m[1][1])\
        / ( m[0][0] * m[1][1] - m[0][1] * m[1][0] )

        
def calculate_delata_U2(U1,U2,un):
    m = calculate_m(U1, U2)
    return (F1(U1,U2,un) * m[1][0] - F1(U1,U2,un) * m[0][0])\
        / ( m[0][0] * m[1][1] - m[0][1] * m[1][0] )



# main rk2 function
def rk2():

    # settings 
    n_max = int(t_max / dt) + 1
    t = np.linspace(t_min, t_max, n_max, endpoint=True)
    u_confirmed = np.ones(n_max)

    # solution
    for i in range(n_max-1):
        
        U1_current = u_confirmed[i]
        U2_current = u_confirmed[i]
 
        for mi in range(mi_max):

            delta_U1 = calculate_delata_U1(U1_current, U2_current, u_confirmed[i])

            delta_U2 = calculate_delata_U2(U1_current, U2_current, u_confirmed[i])
            
            U1_next = U1_current + delta_U1
            U2_next = U2_current + delta_U2

            if ( abs(U1_next - U1_current) < tol or \
                 abs(U2_next - U2_current) < tol):
                break

            U1_current = U1_next
            U2_current = U1_next

        u_confirmed[i+1] = u_confirmed[i] + dt*(b[0]*f(U1_current) + b[1]*f(U2_current))

    plot(t, u_confirmed, title='Rk2', xlabel='t', ylabel='u(t), z(t)')