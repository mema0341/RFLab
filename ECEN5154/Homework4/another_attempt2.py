from numpy import round, arange, pi, zeros, real, imag
from sympy import symbols, sqrt, ln, sin, cos, exp
import scipy.integrate as integrate
import scipy

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

V = zeros((N))
V[1], err = integrate.quad(lambda z: -(1/delta_z)*sin(k0*(delta_z-abs(z-z_m[m_idx])))/sin(k0*delta_z), lower_bound, upper_bound)

Zmn = zeros((N))
Zmn[1], err = -30j/sin(k0*delta_z) * integrate.quad(lambda z: real(sin(k0*(delta_z-abs(z-z_m[m_idx])))/sin(k0*delta_z)*(exp(-1j*k0*sqrt(a**2 + (z-z_n[n_idx-1])**2))/sqrt(a**2 + (z-z_n[n_idx-1])**2)-2*cos(k0*delta_z*exp(-1j*k0*sqrt(a**2 + (z-z_n[n_idx])**2))/sqrt(a**2 + (z-z_n[n_idx])**2)+exp(-1j*k0*sqrt(a**2 + (z-z_n[n_idx+1])**2))/sqrt(a**2 + (z-z_n[n_idx+1])**2)))), lower_bound, upper_bound)
print(Zmn[1])