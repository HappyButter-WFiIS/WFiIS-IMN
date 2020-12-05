from init_conditions import *

def count_trapezoid(x_current, v_current, dt, alfa):
    gamma = 1e-10

    x_next = x_current
    v_next = v_current
    
    while True:

        F = x_next - x_current - (dt/2.0)*( f(v_current) + f(v_next) )
        G = v_next - v_current - (dt/2.0)*( g(x_current, v_current, alfa) + g(x_next, v_next, alfa) )

        a00 = 1.0
        a01 = -dt/2.0
        a10 = (-dt/2.0) * (-2.0*alfa*x_next*v_next - 1)
        a11 = 1 - (dt/2.0)*alfa*(1-(x_next*x_next)) 

        dx = ((-F)*a11 - (-G)*a01) / (a00*a11 - a01*a10)
        dv = (a00*(-G) - a10*(-F)) / (a00*a11 - a01*a10) 

        x_next += dx
        v_next += dv

        if(dx < gamma and dv < gamma):
            break

    return (x_next, v_next)

