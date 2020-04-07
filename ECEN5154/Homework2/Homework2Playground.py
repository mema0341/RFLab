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

import openpyxl
import pandas as pd

import pandas as pd

import os

save_excel = False
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

# convert arrays to dataframe
df_nonpadded_matrix = pd.DataFrame(np.flip(nonpadded_matrix))
df_matrix = pd.DataFrame(np.flip(matrix))

if save_excel:
    # save arrays to xlsx file
    df_nonpadded_matrix.to_excel(r'Exported Data\nonpadded_matrix.xlsx', index=False)
    df_matrix.to_excel(r'Exported Data\matrix.xlsx', index=False)

# Find how big your actual matrix is with 0 padding 
pmr, pmc = np.shape(matrix)
unpmr, unpmc = np.shape(nonpadded_matrix)

fm = np.ndarray((pmr, pmc), dtype=object) # Flag Matrix
node_array = np.ndarray((pmr, pmc), dtype=object) # Flag Matrix
# fmc = np.ndarray((unpmr,unpmc)) # Flag matrix (color version!)

ax = pmr
ay = pmr

Am = np.zeros((ax*ax, ay*ay)) # A matrix
np.fill_diagonal(Am, -4)

#print(np.shape(Am))
bv = np.zeros(((ax)*(ay))) # b vector)

# For mapping purposes
km = np.zeros((ax, ay)) # k matrix (for mapping purposes)
km_matlab = np.zeros((ax, ay)) # k matrix (for mapping purposes)
km = np.reshape(np.arange(0,ax*ay,1), ((ax,ay)), order='C') # numbering goes down rows
# km = np.reshape(np.arange(0,ax*ay,1), ((ax,ay)), order='F') # F = Fortran (numbering goes down cols)
km_matlab = km + 1 # for debugging purposes

# Now we want to loop through our padded matrix and figure out all of the flags
for ii in range(ax):
    for jj in range(ay):
        # Init section
        print("##########################################################################################")
        print(str(ii+1),',',str(jj+1))

        b_value = 0
        flag = ""
        center_node = matrix[ii,jj] # center node

        # now, are we on the boundary?
        if center_node == 0:
            flag = "Outer_PEC_Boundary"
        else:
            # center, right, up, left, down
            center_k = km[ii,jj] # center node

            right_node = matrix[ii,jj+1] # down node (inverted because matrix is upside down)
            down_node = matrix[ii+1,jj] # left node (inverted because matrix is upside down)
            left_node = matrix[ii,jj-1] # up node (inverted because matrix is upside down)
            up_node = matrix[ii-1,jj] # right node (inverted because matrix is upside down)

            # Weirdly this is allowing negative indexing, so try/except won't work in this case
            if down_node != 0:
                down_k = km[ii+1,jj] # right node (inverted because matrix is upside down)
            else:
                down_k = -1
                flag = "floor_boundary"
            
            if left_node != 0:
                left_k = km[ii,jj-1] # up node (inverted because matrix is upside down)
            else:
                left_k = -1
                flag = "left_boundary"

            if up_node != 0: 
                up_k = km[ii-1,jj] # left node (inverted because matrix is upside down)
            else:
                up_k = -1
                flag = "upper_boundary"
            
            if right_node != 0:
                right_k = km[ii,jj+1] # down node (inverted because matrix is upside down)
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

        # print("Row: %s, Col: %s, K: %s" % (ii+1, jj+1, km_matlab[ii,jj]))
        fm[ii, jj] = flag
        if flag != "Outer_PEC_Boundary": 
            # Am, bv = remapping_A(ii,jj,k_vals, node_vals,stencil,km, Am, ax,bv,b_value)
            # print("Node vals:", str(node_vals))
            # print("K vals: ", str(k_vals))
            # print("Stencil vals: ", str(stencil))
        
            # center, right, up, left, down
            print(flag)
            print("                %s (%s, %s)     " % (node_vals[2], k_vals[0], k_vals[2]))
            print("%s (%s, %s)    %s (%s, %s)    %s (%s, %s)" % (node_vals[3], k_vals[0], k_vals[3], node_vals[0], k_vals[0], k_vals[0], node_vals[1], k_vals[0], k_vals[1]))
            print("                %s (%s, %s)     " % (node_vals[4], k_vals[0], k_vals[4]))            
            
            str_val = ''
            # str_val = flag + '                %s (%s, %s)     ' % (node_vals[2], k_vals[0], k_vals[2])
            # str_val = str_val + '%s (%s, %s)    %s (%s, %s)    %s (%s, %s)' % (node_vals[3], k_vals[0], k_vals[3], node_vals[0], k_vals[0], k_vals[0], node_vals[1], k_vals[0], k_vals[1])
            # str_val = str_val + '                %s (%s, %s)     ' % (node_vals[4], k_vals[0], k_vals[4])
            # node_array[ii,jj] = str_val
            
            # Loop through values and add to A
            # center
            if k_vals[0] != -1 and node_vals[0] != 0:
                print("Center: A[",str(k_vals[0]),str(k_vals[0]),"] =",str(stencil[0]))
                Am[k_vals[0], k_vals[0]] = stencil[0]
                str_val = str_val + "Center: A["+str(k_vals[0])+','+str(k_vals[0])+"] ="+str(stencil[0])

            # Right
            if k_vals[1] != -1 and node_vals[1] != 0:
                print("Right: A[",str(k_vals[0]),str(k_vals[1]),"] =",str(stencil[1]))
                Am[k_vals[0], k_vals[1]] = stencil[1]
                str_val = str_val + "Right: A["+str(k_vals[0])+','+str(k_vals[1])+"] ="+str(stencil[1])

            # Up
            if k_vals[2] != -1 and node_vals[2] != 0:
                print("Up: A[",str(k_vals[0]),str(k_vals[2]),"] =",str(stencil[2]))
                str_val = str_val + "Up: A["+str(k_vals[0])+','+str(k_vals[2])+"] ="+str(stencil[2])
                Am[k_vals[0], k_vals[2]] = stencil[2]

            if k_vals[3] != -1 and node_vals[3] != 0:    
                print("Left: A[",str(k_vals[0]),str(k_vals[3]),"] =",str(stencil[3]))
                Am[k_vals[0], k_vals[3]] = stencil[3]
                str_val = str_val + "Left: A["+str(k_vals[0])+','+str(k_vals[3])+"] ="+str(stencil[3])

            if k_vals[4] != -1 and node_vals[4] != 0:
                print("Down: A["+str(k_vals[0])+str(k_vals[4])+"] ="+str(stencil[4]))
                Am[k_vals[0], k_vals[4]] = stencil[4]
                str_val = str_val + "Down: A["+str(k_vals[0])+','+str(k_vals[4])+"] ="+str(stencil[4])
            
            # Place B
            if b_value != 0:
                print("B["+str(k_vals[0])+"] = "+str(b_value))
                bv[k_vals[0]] = b_value
                str_val = str_val + "B["+str(k_vals[0])+"] = "+str(b_value)
            
            node_array[ii,jj] = str_val

A = np.linalg.inv(Am)
x = A.dot(bv)
# x = Am.dot(bv)
x = x.reshape((ax, ay))

# save_excel = True

df_node_array = pd.DataFrame(node_array)
df_node_array.to_excel(r'Exported Data\df_node_array.xlsx', index=False) # save arrays to xlsx file


if save_excel:
    # convert flags to dataframe and export
    df_fm = pd.DataFrame(fm)
    df_fm.to_excel(r'Exported Data\fm.xlsx', index=False) # save arrays to xlsx file

    # convert A to dataframe and export
    df_Am = pd.DataFrame(Am)
    df_Am.to_excel(r'Exported Data\A.xlsx', index=False) # save arrays to xlsx file

    # convert b to dataframe and export
    df_bv = pd.DataFrame(bv)
    df_bv.to_excel(r'Exported Data\bv.xlsx', index=False) # save arrays to xlsx file

    # convert b to dataframe and export
    df_x = pd.DataFrame(x)
    df_x.to_excel(r'Exported Data\x.xlsx', index=False) # save arrays to xlsx file

plot_matrix(np.flipud(x), 'Test?')

for i in range(ax):
    for j in range(ay):
        # print("A[i, j] = %s ; A[j,i] = %s" % (A[i,j], A[j, i]))
        if Am[i,j] != Am[j,i]:
            print('hey')