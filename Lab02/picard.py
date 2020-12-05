from plot_generator import plot, np
from init_conditions import *


def picard():

    # settings 
    n_max = int(t_max / dt) + 1
    t = np.linspace(t_min, t_max, n_max, endpoint=True)
    u_confirmed = np.ones(n_max)

    u_current_mi = u_0
    u_next_mi = 0

    # solution
    for i in range(n_max-1):
        u_current_mi = u_confirmed[i]
 
        for mi in range(mi_max):
            
            u_next_mi = u_confirmed[i] + (dt/2)*(f(u_confirmed[i]) + f(u_current_mi))
 
            if ( abs(u_next_mi - u_current_mi) < tol):
                break
            u_current_mi = u_next_mi

        u_confirmed[i+1] = u_next_mi

    plot(t, u_confirmed, title='Picard', xlabel='t', ylabel='u(t), z(t)')




