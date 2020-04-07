import numpy as np
from homework4helperfunctions import *

# Reference: http://www.iitg.ac.in/engfac/krs/public_html/lectures/ee340/2014/10_slides.pdf
#   Slide ~92
# delta_z = 0.1*lambda

Z = draw_impedance_matrix()

# Because we have a delta gap excitation (1V exists at the feed gap only)
V = [0,-1,0]

I = np.dot(np.linalg.inv(Z), V)

print(I[0])
print(I[1])
print(I[2])

# Zin = 2*V0 / I_(N/2)
Zin = 2*1/I[1]
print(Zin)
print(np.real(Zin))
print(np.imag(Zin))
print(np.abs(Zin))
