# Problem 2
import scipy.integrate as integrate
from numpy import pi, sin, cos, divide
import numpy as np
import pylab as plt
import random
import matplotlib.cm as cm
from sympy import symbols
from numpy import shape

#############################################################################################
# Variable Declaration
voltage = 10

# This must be changed to 0.05 at some point
step_size = 0.1

# Heights
a = 1
b = a
h1 = 0.4
h2 = 0.2
h3 = 0.4

# define DKs
dk_sub = 2.2
dk_air = 1

w = 0.2
d = 0

#############################################################################################
# Function Declaration

# def parameter_extraction(x):
#     # We gotta find where the solution for the stripline is
#     strip_rows,strip_cols = np.where(x==voltage)  
#     # strip_row = strip_rows[0]
#     stripIndY = round(0.6/step_size+1)
#     stripWidth = round(w/step_size+1)
#     stripIndL = strip_cols[0]
#     stripIndR = strip_cols[len(strip_cols)-1]
    
#     # C = 0.5*(dk_sub+dk_air)*(x[stripInd:stripIndR,stripIndY-1])+
#         x[stripIndR-1,stripIndY]+
#         dk_air*(sum(x[strip_rows,stripIndY+1]))+
#         dk_sub*(sum(x[strip_rows,stripIndY-1]))-
#         (V*((w/step_size)+1)*0.5*(dk_air+dk_sub)

#     Co = (dk_air)*(x[stripIndL-1,stripIndY]+x[stripIndR-1,stripIndY])+
#         dk_air*sum(x[stripIndL:stripIndR, stripIndY+1])+
#         dk_air*sum(x[stripIndL:stripIndR,stripIndY-1])

    # C = C/V
    # Co = Co/V
    # Zc = 1/(3e8)/np.sqrt(Co*C)/8.85*10^-12

def remapping_A(ii,jj,k_vals, node_vals,stencil, km, Am, ax, bv, b_value):
# center, right, up, left, down
# Middle
    i = ii
    j = jj

    ii = ii + 1
    jj = jj + 1

    mx = ((jj-1)*ax+ii)-1
    my = ((jj-1)*ax+ii)-1
    print(mx, my)
    print(k_vals[0])
  
    lx = ((jj-1)*ax+ii)-1
    ly = ((jj-1)*ax+ii-1)-1
    print(lx,ly)
    print(k_vals[3])
    
    rx = ((jj-1)*ax+ii)-1
    ry = ((jj-1)*ax+ii+1)-1
    print(rx,ry)
    print(k_vals[1])
    
    ux = ((jj-1)*ax+ii)-1
    uy = ((jj)*ax+ii)-1
    print(ux,uy)
    print(k_vals[2])
    
    dx = ((jj-1)*ax+ii)-1
    dy = ((jj-2)*ax+ii)-1
    print(dx,dy)
    print(k_vals[4])
    
    b_idx = (((jj-1)*ax+ii)-1)-1
    print(b_idx)

    # # center, right, up, left, down
    # a_idxs = [(mx,my), (rx,ry), (ux,uy), (lx,ly), (dx,dy)]
    # return a_idxs, b_idx

    bx = (jj-1)*ax+ii - 1
    #print(bx)

    # if b_value != 0:
    # #print("**************************FOUND STRIPLINE INTERACTION**************************")
    # center
    if k_vals[0] != -1 and node_vals[0] != 0:
        print("Center: A[",str(mx),str(my),"] =",str(stencil[0]))
        Am[mx, my] = stencil[0]

    # Right
    if k_vals[1] != -1 and node_vals[1] != 0:
        print("Right: A[",str(rx),str(ry),"] =",str(stencil[1]))
        Am[rx, ry] = stencil[1]

    # Up
    if k_vals[2] != -1 and node_vals[2] != 0:
        print("Up: A[",str(ux),str(uy),"] =",str(stencil[2]))

        Am[ux, uy] = stencil[2]

    if k_vals[3] != -1 and node_vals[3] != 0:    
        print("Left: A[",str(lx),str(ly),"] =",str(stencil[3]))
        Am[lx, ly] = stencil[3]

    if k_vals[4] != -1 and node_vals[4] != 0:
        print("Down: A[",str(dx),str(dy),"] =",str(stencil[4]))
        Am[dx, dy] = stencil[4]

    # Place B
    print("B[",str(bx),"] = ",str(b_value))
    bv[bx] = b_value

    return Am, bv

def get_stencil(flag, dk_array, stripline_involved, on_pec_boundary):
    # center, right, up, left, down
    epsT = dk_array[2]
    epsB = dk_array[4]
    
    if stripline_involved:
        # I have decided that you're on the boundary since your stripline is infinitely small
        stencil = stencils[int(np.squeeze(np.where(stencil_types=='air_up_stencil')))]    
        
        # Nevermind
        # # I have decided to temporarily trust Matt
        # stencil = stencils[int(np.squeeze(np.where(stencil_types=='air_up_stencil')))]   

# stencil_types = np.asarray(["air_hom_stencil", "sub_hom_stencil", "air_up_stencil", "sub_up_stencil"])
    elif on_pec_boundary:
        # If you're on the ceiling/floor you're in air
        if flag == "corner":
            stencil = stencils[int(np.squeeze(np.where(stencil_types=='air_hom_stencil')))]    
            flag = flag+"_homo_air"        
        elif flag == "upper_boundary" or flag == "floor_boundary":
            stencil = stencils[int(np.squeeze(np.where(stencil_types=='air_hom_stencil')))]    
            flag = flag+"_homo_air"
        else: # Otherwise you can check for nonhomo and homo
            dk_array = dk_array[np.where(dk_array!=0)]
            udk_array = np.unique(dk_array)
            if len(udk_array) == 1:
                # Case 1: Homogeneous
                homogeneous = True
                if udk_array[0] == dk_air:
                    stencil = stencils[int(np.squeeze(np.where(stencil_types=='air_hom_stencil')))]
                    flag = flag+"_homo_air"
                else:
                    stencil = stencils[int(np.squeeze(np.where(stencil_types=='sub_hom_stencil')))]                
                    flag = flag+"_homo_sub"
            else: # Case 2: Nonhomogeneous - find out if air is top or bottom
                homogeneous = False
                if epsT == dk_air:
                    stencil = stencils[int(np.squeeze(np.where(stencil_types=='air_up_stencil')))]    
                    flag = flag+"_nonhomo_air_up"
                else:
                    stencil = stencils[int(np.squeeze(np.where(stencil_types=='sub_up_stencil')))]    
                    flag = flag+"_nonhomo_sub_up"

    else: # No Stripline or PEC, so check if you're homo or nonhomo
        dk_array = dk_array[np.where(dk_array!=0)]
        udk_array = np.unique(dk_array)
        if len(udk_array) == 1:
            # Case 1: Homogeneous
            homogeneous = True
            if udk_array[0] == dk_air:
                stencil = stencils[int(np.squeeze(np.where(stencil_types=='air_hom_stencil')))]
                flag = flag+"_homo_air"
            else:
                stencil = stencils[int(np.squeeze(np.where(stencil_types=='sub_hom_stencil')))]                
                flag = flag+"_homo_sub"
        else: # Case 2: Nonhomogeneous - find out if air is top or bottom
            homogeneous = False
            if epsT == dk_air:
                stencil = stencils[int(np.squeeze(np.where(stencil_types=='air_up_stencil')))]    
                flag = flag+"_nonhomo_air_up"
            else:
                stencil = stencils[int(np.squeeze(np.where(stencil_types=='sub_up_stencil')))]    
                flag = flag+"_nonhomo_sub_up"

    return flag, np.asarray(stencil)

def reduce_matrix(dks, reduction_vals=[0,voltage]):
    for rr in range(len(reduction_vals)):
        dks = dks[np.where(dks!=reduction_vals[rr])]
    return dks

# center, right, up, left, down
def FD_coefficients_simplified(epsT, epsB):
    phi0,phi1,phi2,phi3,phi4 = symbols('phi0 phi1 phi2 phi3 phi4')
    expr = ( ((epsT+epsB)/2) * (phi1-phi0) + epsT*(phi2-phi0) + ((epsT+epsB)/2) * (phi3-phi0) + epsB*(phi4-phi0) )
    # Down Left Center Right Up
    return np.asarray([expr.coeff(phi0), expr.coeff(phi1), expr.coeff(phi2), expr.coeff(phi4), expr.coeff(phi4)], dtype=float)

def plot_matrix(matrix, title):
    # plt.figure()
    fig, my_ax = plt.subplots()
    pos = my_ax.imshow(matrix, cmap='jet') # pos = ax1.imshow(Zpos, cmap ='jet', interpolation='none')
    fig.colorbar(pos, ax=my_ax)
    plt.title(title)
    plt.show()

def draw_matrix():
    # Define size of array cols/rows
    rows = np.arange(0,a,step_size)
    cols = rows

    # Find where your substrate boundaries are
    idx_h1 = int(np.squeeze(np.where(rows==h1)))
    idx_h2 = int(np.squeeze(np.where(rows==h1+h2)))

    # Find where your stripline is 
    idx_h3 = int(np.squeeze(np.where(rows==h1+h2)))
    tmp = int(np.squeeze(np.where(rows==a/2)))
    # sidx_col = np.arange(tmp-int(w/(2*step_size))-1,tmp+int(w/(2*step_size)))
    sidx_col = np.arange(tmp-int(w/(2*step_size)),tmp+int(w/(2*step_size)))

    # Define Matrix
    matrix = np.zeros((len(rows),len(cols)))

    # Fill in air
    matrix[0:len(rows),:] = dk_air

    # Fill in substrate
    matrix[idx_h1:idx_h2,:] = dk_sub

    # Add stripline
    matrix[idx_h3-1,sidx_col] = voltage

    if d > 0:
        matrix = np.pad(matrix,int(d/step_size))
        idx_h1 = idx_h1+int(d/step_size)
        idx_h2 = idx_h2+int(d/step_size)
        matrix[idx_h1:idx_h2,0] = dk_sub
        matrix[idx_h1:idx_h2,len(matrix)-1] = dk_sub
        
    # Add padding
    padded_matrix = np.pad(matrix,1)

    return matrix, padded_matrix

# Stencils
stencil_types = np.asarray(["air_hom_stencil", "sub_hom_stencil", "air_up_stencil", "sub_up_stencil"])
stencils = [FD_coefficients_simplified(dk_air, dk_air), FD_coefficients_simplified(dk_sub, dk_sub), FD_coefficients_simplified(dk_air, dk_sub), FD_coefficients_simplified(dk_sub, dk_air)]
