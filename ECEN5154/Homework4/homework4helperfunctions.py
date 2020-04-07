import numpy as np

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
