from plot_generator import *
from analytical import math

def V(t, omega):
    return 10.0 * math.sin(t * omega)


def rrz2():
    # settings
    labels = ['0.5 * w0', '0.8 * w0', '1.0 * w0', '1.2 * w0']

    ax1, ax2 = generate_plot_area(title1='Q(t) calculated with RRZ2 method', title2='I(t) calculated with RRZ2 method', \
                                        xlabel1='t', xlabel2='t', \
                                        ylabel1='Q(t)', ylabel2='I(t)')


    # the initial conditions
    R = 100.0
    L = 0.1
    C = 0.001

    dt = math.pow(10,-4)
    w0 = 1.0 / math.sqrt(L * C)
    T0 = (2.0 * math.pi) / w0
    t_min = 0.0
    t_max = 4.0 * T0
    wV_array = [.5 * w0, .8 * w0, 1.0 * w0, 1.2 * w0]

    number_of_steps = int(t_max / dt) + 1   # this time number_of_steps is constant 
    t = np.linspace(t_min, t_max, number_of_steps, endpoint=True)   


    for i, omega in enumerate(wV_array):

        Q = np.zeros(number_of_steps)
        I = np.zeros(number_of_steps)

        for j in range(0, number_of_steps):
            
            current_time = j * dt
            current_time_half = (j + 0.5) * dt
            current_time_one = (j + 1.0) * dt
        
            Vt = V(current_time, omega)
            Vt_half = V(current_time_half, omega)
            Vt_one = V(current_time_one, omega)

            k1_Q = I[j-1]
            k1_I = (Vt/L) - ((1.0/(L * C)) * Q[j-1]) - ((R/L) * I[j-1])

            k2_Q = I[j-1] + ((dt/2.0) * k1_I)
            k2_I = (Vt_half / L) - ((1.0 / (L * C)) * (Q[j-1] + (dt / 2.0) * k1_Q)) - ((R / L) * (I[j-1] + (dt / 2.0) * k1_I))

            k3_Q = I[j-1] + ((dt/2.0) * k2_I )
            k3_I = (Vt_half/L) - ( (1.0/(L * C)) * (Q[j-1] + (dt/2.0) * k2_Q) ) - ((R/L) * (I[j-1] + (dt/2.0) * k2_I))

            k4_Q = I[j-1] + (dt * k3_I)
            k4_I = (Vt_one/L) - ( (1.0/( L* C)) * (Q[j-1] + (dt * k3_Q))) - ((R/L) * (I[j-1] + (dt * k3_I)))

            Q[j] = Q[j-1] + (dt/6.0)*(k1_Q + 2.0 * k2_Q + 2.0 * k3_Q + k4_Q)
            I[j] = I[j-1] + (dt/6.0)*(k1_I + 2.0 * k2_I + 2.0 * k3_I + k4_I)
            
        ax1.plot(t, Q, label=labels[i])
        ax2.plot(t, I, label=labels[i])


    ax1.grid(True)
    ax2.grid(True)
    ax1.legend(loc="upper right")
    ax2.legend(loc="lower right")
    plt.savefig("rrz2.png")
    plt.close()
