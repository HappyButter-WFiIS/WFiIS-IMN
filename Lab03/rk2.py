from init_conditions import *

def count_rk2(x_current, v_current, dt, alfa): 
    
    k1x = f(v_current)
    k1v = g(x_current, v_current, alfa)

    k2x = f(v_current + dt*k1v)
    k2v = g(x_current + dt*k1x, v_current + dt*k1v, alfa)

    x_next = x_current + (dt/2.0)*(k1x+k2x)
    v_next = v_current + (dt/2.0)*(k1v + k2v)

    return (x_next, v_next)

