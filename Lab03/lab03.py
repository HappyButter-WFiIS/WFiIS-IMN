from time_control import time_control
from trapezoid import count_trapezoid
from rk2 import count_rk2

if __name__ == "__main__":
    
    time_control(method_procedure=count_trapezoid,\
                 method_name='trapezoid')
    
    time_control(method_procedure=count_rk2, \
                 method_name='rk2')