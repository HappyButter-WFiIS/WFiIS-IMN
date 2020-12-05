from analytical import analytical 
from plot_generator import *

def eulerMethod():
    # settings
    plot_types = [',', 'x', '*']
    labels = ['dt = 0.01', 'dt = 0.1', 'dt = 1.0']
    
    ax1, ax2 = generate_plot_area(title1='Eulers method', title2='Error', \
                                        xlabel1='t', xlabel2='t', \
                                        ylabel1='y(t)', ylabel2='ynum(t) - yanalytical(t)')


    # the initial conditions
    t_min = 0.0
    t_max = 5.0
    lam = -1.0
    dt_array = [0.01, 0.1, 1.0]

    for i, dt in enumerate(dt_array):

        number_of_steps = int(t_max / dt) + 1
        t = np.linspace(t_min, t_max, number_of_steps, endpoint=True)
        
        y = np.ones(number_of_steps)
        y_err = np.zeros(number_of_steps)


        for j in range(1, number_of_steps):
            y[j] = y[j-1] + dt*lam*y[j-1]
            y_err[j] = y[j] - analytical(t[j], lam)
        
        ax1.plot(t,y,plot_types[i], label=labels[i])
        ax2.plot(t,y_err,plot_types[i], label=labels[i])
        
    # analytical solution
    number_of_steps = 10000
    t = np.linspace(t_min, t_max, number_of_steps, endpoint=True)
    y = np.exp(t*lam)

    ax1.plot(t,y)
    ax1.plot(t,y, label='e^(t*lambda)')


    ax1.grid(True)
    ax2.grid(True)
    ax1.legend(loc="upper right")
    ax2.legend(loc="lower right")
    plt.savefig("euler.png")
    plt.close()
