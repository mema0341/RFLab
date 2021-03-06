# Problem 2
import scipy.integrate as integrate
from numpy import pi, sin, cos, divide
import numpy as np
import pylab as plt
import random
import matplotlib.cm as cm
from sympy import symbols
from numpy import shape
from Homework2HelperFunctions import *

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

# Stencils
stencil_types = np.asarray(["air_hom_stencil", "sub_hom_stencil", "air_up_stencil", "sub_up_stencil"])
stencils = [FD_coefficients_simplified(dk_air, dk_air), FD_coefficients_simplified(dk_sub, dk_sub), FD_coefficients_simplified(dk_air, dk_sub), FD_coefficients_simplified(dk_sub, dk_air)]

nonpadded_matrix, matrix = draw_matrix()

# plot_matrix(np.flip(matrix),'Does d exist?')

# Find how big your actual matrix is with 0 padding 
pmr, pmc = np.shape(matrix)
unpmr, unpmc = np.shape(nonpadded_matrix)

# fm = np.ndarray((unpmr, unpmc), dtype=object) # Flag Matrix
# fmc = np.ndarray((unpmr,unpmc)) # Flag matrix (color version!)

ax = (unpmr+1)
ay = (unpmc+1)

Am = np.zeros(((ax)*(unpmr+1), (unpmc+1)*(unpmc+1))) # A matrix
np.fill_diagonal(Am, -4)

#print(np.shape(Am))
bv = np.zeros(((ax)*(ay))) # b vector)

# For mapping purposes
km = np.zeros((ax, ay)) # k matrix (for mapping purposes)
for kr in range(ax):
    for kc in range(ay):
        km[kr,kc] = unpmr*(kc)+kr

# Now we want to loop through our padded matrix and figure out all of the flags
for ii in range(ax):
    for jj in range(ay):
        # Init section
        print("##########################################################################################")
        print(str(ii+1),',',str(jj+1))

        b_value = 0
        # Set Flag to 0
        flag = ""
        # This is just for sanity: make i and j equal to where the non-padded solve-domain is
        i = ii-1 # This is to index the inner matrix 
        j = jj-1 # This is to index the inner matrix 

        center_node = matrix[ii,jj] # center node

        # now, are we on the boundary?
        if center_node == 0:
            flag = "Outer_PEC_Boundary"
        else:
            # center, right, up, left, down
            center_k = km[ii,jj] # center node
            right_node = matrix[ii,jj+1] # down node (inverted because matrix is upside down)
            up_node = matrix[ii+1,jj] # left node (inverted because matrix is upside down)
            left_node = matrix[ii,jj-1] # up node (inverted because matrix is upside down)
            down_node = matrix[ii-1,jj] # right node (inverted because matrix is upside down)

            # if right_node == 0:
            #     #print("bitches love titles")
            
            # Weirdly this is allowing negative indexing, so try/except won't work in this case
            if down_node != 0:
                down_k = km[i-1,j] # right node (inverted because matrix is upside down)
            else:
                down_k = -1
                flag = "floor_boundary"
            
            if left_node != 0:
                left_k = km[i,j-1] # up node (inverted because matrix is upside down)
            else:
                left_k = -1
                flag = "left_boundary"

            if up_node != 0: 
                up_k = km[i+1,j] # left node (inverted because matrix is upside down)
            else:
                up_k = -1
                flag = "upper_boundary"
            
            if right_node != 0:
                right_k = km[i,j+1] # down node (inverted because matrix is upside down)
            else:
                right_k = -1
                flag = "right_boundary"
    
            # Unnecessary but could be used for a sanity check
            k_vals = np.asarray([center_k, right_k, up_k, left_k, down_k],dtype=int)

            original_k_vals = np.asarray([center_k, right_k, up_k, left_k, down_k],dtype=int)

            node_vals = np.asarray([center_node, right_node, up_node, left_node, down_node])
            
            # Now we have our flags for the PEC. Finish those up first.
            if len(flag) == 0:
                if len(node_vals[np.where(node_vals==voltage)]) != 0:
                    if center_node == voltage:
                        flag = 'stripline'
                        # This could be put in the function but I'm being lazy
                        k_vals[1:5] = -1
                        stencil = np.asarray([1, -1, -1, -1, -1])
                        b_idx = k_vals[0]
                        b_value = voltage
                    else:
                        flag = 'interacting_stripline_'
                        flag, stencil = get_stencil(flag,node_vals,stripline_involved=True, on_pec_boundary = False)
                        b_idx = k_vals[0]
                        k_vals[np.where(node_vals==voltage)] = -1
                        b_value = -voltage
                else:
                    flag = "non_boundary"
                    flag, stencil = get_stencil(flag,node_vals,stripline_involved=False, on_pec_boundary = False)
            
            # center, right, up, left, down
            else: # Flag: Boundary
                # First we can check for corners
                if len(node_vals[np.where(node_vals == 0)]) > 1:
                    flag = "corner"
                    
                # Step 1: Figure out if it's homogeneous
                flag, stencil = get_stencil(flag,node_vals,stripline_involved=False, on_pec_boundary = True)
        #print(flag)

        # if flag == "ground_boundary":
            # # if flag == "stripline" or flag == "interacting_stripline_":
    # Useful debug
        # if flag != "non_boundary_homo_air" and flag != "corner_homo_air" and flag != "floor_boundary_homo_air" and flag != "Outer_PEC_Boundary" and flag != "floor_boundary_homo_air" and flag != "non_boundary__homo_air" and flag != "left_boundary_homo_air" and flag != "right_boundary_homo_air": 
        #     tmp = matrix[ii,jj]
        #     matrix[ii,jj]=50
        #     plot_matrix(np.flip(matrix,0), flag)
        #     matrix[ii,jj]=tmp

        if flag != "Outer_PEC_Boundary": 
            Am, bv = remapping_A(ii,jj,k_vals, node_vals,stencil,Am, ax,bv,b_value)

#print(Am)
#print(bv)

# plot_matrix(Am,'A')

A = np.linalg.inv(Am)
x = A.dot(bv)
x = x.reshape((ax, ay))

plot_matrix(np.flip(x), "Test")

print(np.flip(x))

if w == 0.2 and step_size == 0.1:
    np.savetxt('bv_02_01.txt',bv)
    np.savetxt('Am_02_01.txt',Am)
elif w == 0.4 and step_size == 0.1:
    np.savetxt('bv_04_01.txt',bv)
    np.savetxt('Am_04_01.txt',Am)
elif w == 0.6 and step_size == 0.1:
    np.savetxt('bv_06_01.txt',bv)
    np.savetxt('Am_06_01.txt',Am)
if w == 0.2 and step_size == 0.05:
    np.savetxt('bv_02_005.txt',bv)
    np.savetxt('Am_02_005.txt',Am)
elif w == 0.4 and step_size == 0.05:
    np.savetxt('bv_04_005.txt',bv)
    np.savetxt('Am_04_005.txt',Am)
else:
    np.savetxt('bv_06_005.txt',bv)
    np.savetxt('Am_06_005.txt',Am)

#print(np.flip(np.round(x,2)))

# strip_rows,strip_cols = np.where(x==voltage)  
# # strip_row = strip_rows[0]
# stripIndY = round(0.6/step_size+1)
# stripWidth = round(w/step_size+1)
# stripIndL = strip_cols[0]
# stripIndR = strip_cols[len(strip_cols)-1]
    
# C = 0.5*(dk_sub+dk_air)*(x[stripIndL:stripIndR,stripIndY-1])+x[stripIndR-1,stripIndY]+dk_air*(sum(x[stripIndL:stripIndR,stripIndY+1]) + dk_sub*(sum(x[stripIndL:stripIndR,stripIndY-1]))) - (voltage*((w/step_size)+1)*0.5*(dk_air+dk_sub))
# Co = (dk_air)*(x[stripIndL-1,stripIndY])+x[stripIndR-1,stripIndY]+dk_air*sum(x[stripIndL:stripIndR, stripIndY+1])+dk_air*sum(x[stripIndL:stripIndR,stripIndY-1])
# C = C/voltage
# Co = Co/voltage
# Zc = 1/(3e8)/np.sqrt(Co*C)