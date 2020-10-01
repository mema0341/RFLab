import numpy as np

N = 3

one = 0; two = 1; three = 2; four = 3;

z12 = 1; z21 = z12; z23 = z12; z32 = z12;
z13 = 2; z31 = z13
z11 = 3; z22 = z11; z33 = z11;

z = np.zeros((N,N))
z[one,one] = z11
z[one,two] = z12
z[one,three] = z13
z[two, one] = z21
z[two,two] = z22
z[two,three] = z23
z[three,one] = z31
z[three,two] = z32
z[three,three] = z33

# Voltage of gap 
V = [0,-1,0]

# I = Z^(-1)*V
I = np.linalg.inv(z,V)


deltaZ_div_lam = 0.1/N
a_div_lam = 1e-4/N
Z0 = 50
Zmn = -1j * Z0 * deltaZ_div_lam * ( (1 + 2j*np.pi*a_div_lam) * (2-3) + np.power(2*np.pi*a_div_lam,2) ) * (np.cos(2*np.pi*a_div_lam) - 1j*np.sin(2*np.pi*a_div_lam) ) / (8*np.pi)

# Compressed Sparse Row
A = np.asmatrix([[3, 0, 0, 4, 5],[7, 0, 4, 0, 2], [4, 0, 7, 0, 0], [0, 0, 8, 0 ,0], [9, 7, 0, 0,0]])

Value = A[np.where(A!=0)]
ColumnPointer = np.column_stack(np.where(A!=0))[:,1]

RowPointer = np.column_stack(np.where(A!=0))[:,1]