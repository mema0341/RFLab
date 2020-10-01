import numpy as np
from homework4helperfunctions import *
from numpy import round, arange, pi, zeros, real, imag
from sympy import symbols, sqrt, ln, sin, cos, exp
import scipy.integrate as integrate
import scipy

# This is to create the given impedance matrix
def draw_impedance_matrix():
    N = 3

    one = 0; two = 1; three = 2; four = 3;

    z12 = -0.493889 - 1577.1j; z21 = z12; z23 = z12; z32 = z12;
    z13 = -0.490242 - 132.354j; z31 = z13
    z11 = -0.495107 + 3426.99j; z22 = z11; z33 = z11;

    impedance_matrix = np.zeros((N,N),dtype='complex')

    impedance_matrix[one,one] = z11
    impedance_matrix[one,two] = z12
    impedance_matrix[one,three] = z13
    impedance_matrix[two, one] = z21
    impedance_matrix[two,two] = z22
    impedance_matrix[two,three] = z23
    impedance_matrix[three,one] = z31
    impedance_matrix[three,two] = z32
    impedance_matrix[three,three] = z33
    
    return impedance_matrix

# Definitions
N = 3
lam = 1
length = 0.1*lam
a = 1e-4*lam 
delta_z = length/N
k0 = 2*pi/lam

n = arange(1,N+1,1); n_idx = 1
m = arange(1,N+1,1); m_idx = 1

z_n = -length/2 + n*delta_z
z_m = -length/2 + m*delta_z

upper_bound = z_m[m_idx]+delta_z
lower_bound = z_m[m_idx]-delta_z

# Grab Impedance Matrix
Z = draw_impedance_matrix()

# Calculation of Vm
V = zeros((N))
# Because we have a delta gap excitation (1V exists at the feed gap only)
V[1], err = integrate.quad(lambda z: -(1/delta_z)*sin(k0*(delta_z-abs(z-z_m[m_idx])))/sin(k0*delta_z), lower_bound, upper_bound)

# Obtain current by solving
I = np.dot(np.linalg.inv(Z), V)

print(I[0])
print(I[1])
print(I[2])

# Zin = V0 / I_(N/2)
Zin = 1/I[1]
print(Zin)
print("Real: "+str(np.real(Zin)))
print("Imaginary "+str(np.imag(Zin)))


