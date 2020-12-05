from init_conditions import *
from plot_generator import plot

def time_control(method_procedure, method_name=''):

    x_array = generate_array_of_arrays(2)
    v_array = generate_array_of_arrays(2)
    t_array = generate_array_of_arrays(2)
    dt_array = generate_array_of_arrays(2)  

    for i, tol in enumerate(tol_array):

        x = 0.01
        v = 0.0
        t = 0.0
        dt = 1.0

        while t < t_max:
            (x_next_two_steps, v_next_two_steps) = method_procedure(x, v, dt, alfa)
            (x_next_two_steps, v_next_two_steps) = method_procedure(x_next_two_steps, v_next_two_steps, dt, alfa)

            (x_next_one_step, v_next_one_step) = method_procedure(x, v, 2*dt, alfa)

            err_x = (x_next_one_step - x_next_two_steps) / (pow(2, p) - 1)
        
            err_v = (v_next_one_step - v_next_two_steps) / (pow(2, p) - 1)

            if (max(abs(err_x), abs(err_v)) < tol):
                t = t + 2*dt
                x = x_next_two_steps
                v = v_next_two_steps

                x_array[i].append(x)
                v_array[i].append(v)
                t_array[i].append(t)
                dt_array[i].append(dt)  
                  
            dt = pow((S*tol) / (max(abs(err_x), abs(err_v))), 1.0/(p+1.0))*dt


    plot(t_array,v_array,title=method_name,xlabel='t',ylabel='v(t)')
    plot(t_array,x_array,title=method_name,xlabel='t',ylabel='x(t)')
    plot(t_array,dt_array,title=method_name,xlabel='t',ylabel='dt(t)')
    plot(x_array,v_array,title=method_name,xlabel='x',ylabel='v(x)')