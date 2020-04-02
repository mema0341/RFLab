import scipy.integrate as integrate
from numpy import pi, sin, cos, divide
import numpy as np

# Homework 1 Problem 1
def integral_to_evaluate(x):
    return x*cos(10*x)*sin(20*x)

# Solution using Midpoint Rule
def midpoint_rule(h, ):
    x = np.arange(a,b,h)
    midpoint_x = []
    y_a = []
    for i in range(len(x)-1):
        midpoint_x.append((x[i+1]+x[i])/2)
        y_a.append(integral_to_evaluate(midpoint_x[i]))
    return np.sum(y_a)*h

# Solution using Trapezoidal Rule
def trapezoidal_rule(h, ):
    x = np.arange(a,b,h)
    f = []
    y_b = []
    for i in range(len(x)):
        y_b.append(integral_to_evaluate(x[i]))
    return h/2*(y_b[0]+y_b[len(y_b)-1])+np.sum(y_b[1:len(y_b)-1])*h

# Solution using Simpson's Rule
def simpsons_rule(n, a, b):
    N = 2 * n
    h = divide((b-a),N)
    x = np.arange(a,b,h)
    y_c_even = []
    y_c_odd = []
    for i in range(len(x)):
        if i == 0:
            y0 = integral_to_evaluate(x[i])
        elif i == len(x)-1:
            yN = integral_to_evaluate(x[i])
        elif i % 2:
            y_c_even.append(integral_to_evaluate(x[i]))
        else:
            y_c_odd.append(integral_to_evaluate(x[i]))       
    return h/3*(y0 + np.sum(4*y_c_even) + np.sum(2*y_c_odd) + yN)

a = 0
b = 2*pi

y, err = integrate.quad(integral_to_evaluate,a,b)

print(y)

N = [10.0, 50.0, 100.0, 500.0]
h = divide((b-a),N)

# a) Midpoint Rule
solution_a = []
print('--------------Midpoint Rule--------------')
for i in range(len(h)):
    solution_a.append(midpoint_rule(h[i]))
    print('N='+str(N[i])+': ',str(solution_a[i]))
    # print('N='+str(N[i])+' Error: ',str(abs(solution_a[i]-y)))
# For fun used h = 0.00001
# print('N=100,000.00: ', midpoint_rule(0.0001))
print('N=100,000.00 Error: ', abs(y-midpoint_rule(0.0001)))

# Trapezoidal Rule
solution_b = []
print('\n--------------Trapezoidal Rule--------------')
for i in range(len(h)):
    solution_b.append(trapezoidal_rule(h[i]))
    print('N='+str(N[i])+': ',str(solution_b[i]))
    # print('N='+str(N[i])+' Error: ',str(abs(solution_b[i]-y)))

# For fun used h = 0.00001
print('N=100,000.00: ', trapezoidal_rule(0.0001))
# print('N=100,000.00 Error: ', abs(trapezoidal_rule(0.0001)-y))

# Simpson's Rule
solution_c = []
print('\n--------------Simpsons Rule--------------')
for i in range(len(h)):
    solution_c.append(simpsons_rule(N[i], a, b))
    print('N='+str(N[i])+' Error: ',str(abs(solution_c[i]-y)))
    # print('N='+str(N[i])+': ',str(solution_c[i]))
# For fun used h = 0.00001
print('N=100,000.00: ', simpsons_rule(1/0.0001, a, b))
# print('N=100,000.00 Error: ', abs( simpsons_rule(1/0.0001, a, b)-y))
