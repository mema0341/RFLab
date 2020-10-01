import numpy as np 
import sympy
import scipy.integrate as integrate
import scipy.special as special 
#  Important parameters
N = 3
lam = 1
length = 0.1*lam
a = 1e-4*lam 
epsilon_naught = 8.854 *1e-12
val = (4*np.pi*epsilon_naught)
delta = np.round(length/N,3)

#################################################
# Symbols
y_prime = sympy.symbols('y_prime') 

def function(ym, y_prime):
    return 1/sympy.sqrt((ym + y_prime)**2 + a**2)

print('N = %s' % N)
print('Delta = %s' % delta)
old_delta = 0
for y in range(N):
    rm = delta*(y+1)
    # print(rm)
    lower_bound = delta*y
    upper_bound = delta*(y+1)
    ym = (upper_bound+lower_bound)/2
    print(ym)
    # print('y%s: %s to %s; rm = %s, ym = %s' % (str(y+1), str(lower_bound),str(upper_bound),str(rm), str(ym)))
    # print(function(ym, y_prime))
    ans, err = integrate.quad(lambda y_prime: function(ym, y_prime), lower_bound, upper_bound)
    # print(ans,err)
    # print('####################################################')
