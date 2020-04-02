import scipy.integrate as integrate
from numpy import pi, sin, cos, divide
import numpy as np

def my_function(x):
    return np.log(abs( (x-x_k) / (x_p-x) ))

# Solution using Midpoint Rule
def midpoint_rule(x_k, x_p):
    # Function 1
    # Inputs: 
    #   h: int
    # Outputs:
    #   solution to my_function using midpoint rule
    h = 1/3

    x = np.arange(a,b,h)
    midpoint_x = []
    y_a = []
    for i in range(len(x)-1):
        midpoint_x.append((x[i+1]+x[i])/2)
        y_a.append(my_function(midpoint_x[i]))
    return np.sum(y_a)*h

a = 12/6
b = 16/6

x_p = 18/6
x_k = -18/6


y, err = integrate.quad(my_function,a,b)

print(y)

h = 1/3

# a) Midpoint Rule
print('--------------Midpoint Rule--------------')
print('N=1: ',str(midpoint_rule(x_k, x_p)))