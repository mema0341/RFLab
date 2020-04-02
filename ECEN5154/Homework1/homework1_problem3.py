import scipy.integrate as integrate
from numpy import pi, sin, cos, divide
import numpy as np
import matplotlib.pyplot as plt
import random
import matplotlib.cm as cm

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

r1 = 14.2
r2 = 31.8
r3 = 184.1
r4 = 203.2
r5 = 304.8
r6 = 323.8
R = 425.5

# radii = [r2, r3]
radii = [r1, r2, r3, r4, r5, r6, R]

nshift_radii = np.round((radii/np.max(radii))*255+255)

# circle_points = [50, 25, 1, 3, 1, 2, 0.5]
circle_points = [50, 25, 10, 30, 10, 20, 15]
wedge_points = np.arange(1,21)
# Randomize Wedge Points
random.shuffle(wedge_points)
# print(wedge_points)

# Wire Width
w = 1.5 # * mm

pie_angle = 18
num_wedges = 360/pie_angle
wedge = np.arange(0, 256, 1)


sampling_size = 512
theta = np.arange(0, np.pi, np.pi/sampling_size)
# theta = np.append(np.arange(0.0, np.pi, np.pi/511,0)) 

# circ_coords = np.ndarray((len(theta), 2))
circ_coords = np.ndarray((len(theta), 2, len(radii)))
for ii in range(len(radii)):
    x = radii[ii]*np.cos(theta)
    y = radii[ii]*np.sin(theta)
    circ_coords[:,0,ii] = x
    circ_coords[:,1,ii] = y    
    # circ_coords[:,0] = x
    # circ_coords[:,1] = y  
#     plt.plot(x,y) 
# plt.show()

maxim = np.max(circ_coords)
n_circ_coords = circ_coords/maxim
nshift_circ_coords = np.round((n_circ_coords*255)+255,2)

new_x = np.arange(0,512,1)

# Translate to Matrix
dartboard_matrix = np.ndarray((512, 512))
weight_matrix = np.ndarray((512, 512))
for ee in range(len(radii)):
    y1 = []
    y2 = []
    for value in range(512):
        idx = (np.abs(nshift_circ_coords[:,0,ee]-value)).argmin()
        y1.append(nshift_circ_coords[idx,1, ee])
        y2.append(-nshift_circ_coords[idx,1, ee]+512)

        new_value = int(nshift_circ_coords[idx,1, ee])
        if new_value != 255:
            dartboard_matrix[value, new_value]  = circle_points[ee]
            dartboard_matrix[value, int(-nshift_circ_coords[idx,1, ee]+512)] = circle_points[ee]
            if ee == len(radii)-1:
                weight_matrix[value, new_value]  = circle_points[ee]
                weight_matrix[value, int(-nshift_circ_coords[idx,1, ee]+512)] = circle_points[ee]  

##################################################
# weight_matrix = np.ndarray((512, 512))
n = np.arange(int(num_wedges))
x = np.sin((18*n-9)*np.pi/180)*256
y = np.cos((18*n-9)*np.pi/180)*256
# Wedge Coordinates
for ii in range(20):
    try:
        wedge_line_x = np.arange(0,x[ii],x[ii]/256)+256
    except:
        wedge_line_x = np.ones(256)*x[ii]+256
    try:
        wedge_line_y = np.arange(0,y[ii],y[ii]/256)+256
    except:
        wedge_line_y = np.ones(256)*y[ii]+256
    weight_matrix[wedge_line_x.astype(int),wedge_line_y.astype(int)] = wedge_points[ii]
    
    # if ii == 1:
    #     for r in range(512):
    #         weight_matrix[r, prev_wedge_line_y[r].astype(int):wedge_line_y[r].astype(int)] = ii + 1
    prev_wedge_line_x = wedge_line_x
    prev_wedge_line_y = wedge_line_y
    # plt.plot(wedge_line_x,prev_wedge_line_y)
####################################################
# plt.contour(weight_matrix)
# plt.show()

wedge_coords = np.ndarray((len(theta), 2, len(radii)))
n = np.arange(int(num_wedges))

x = np.sin((18*n-9)*np.pi/180)*256
y = np.cos((18*n-9)*np.pi/180)*256

# # Plot Circles
# for ii in range(len(radii)): 
#     plt.plot(nshift_circ_coords[:,0,ii], nshift_circ_coords[:,1,ii])


# Wedge Coordinates
for ii in range(20):
    try:
        wedge_line_x = np.arange(0,x[ii],x[ii]/256)+256
    except:
        wedge_line_x = np.ones(256)*x[ii]+256
    try:
        wedge_line_y = np.arange(0,y[ii],y[ii]/256)+256
    except:
        wedge_line_y = np.ones(256)*y[ii]+256
    weight_matrix[wedge_line_x.astype(int),wedge_line_y.astype(int)] = wedge_points[ii]
    
    # if ii == 1:
    #     for r in range(512):
    #         weight_matrix[r, prev_wedge_line_y[r].astype(int):wedge_line_y[r].astype(int)] = ii + 1
    prev_wedge_line_x = wedge_line_x
    prev_wedge_line_y = wedge_line_y
    # plt.plot(wedge_line_x,prev_wedge_line_y)

plt.contour(weight_matrix)
plt.show()

# Fill in weight matrix
for cc in range(512):
    cl_idx = np.squeeze(np.where(weight_matrix[:,cc]>0))
    # print(cl_idx)
    try:
        l_idx=len(cl_idx)
    except:
        l_idx=0
    for ee in range(1,l_idx):
        # print(str(cl_idx[ee-1]), '-', str(cl_idx[ee]), '=', weight_matrix[cl_idx[ee],cc])
        weight_matrix[cl_idx[ee-1]+1:cl_idx[ee], cc] = weight_matrix[cl_idx[ee], cc]
    # print(weight_matrix[:, cc])
# plt.contour(weight_matrix)
# plt.show()

# Start Filling in Circles
for dd in range(512):
    idx = np.squeeze(np.where(dartboard_matrix[dd,:] > 0))
    # print(idx)
    try:
        l_idx=len(idx)
    except:
        l_idx=0
    for ee in range(1,l_idx):
        # print(str(idx[ee-1]), '-', str(idx[ee]), '=', dartboard_matrix[dd,idx[ee]])
        dartboard_matrix[dd,idx[ee-1]+1:idx[ee]] = dartboard_matrix[dd,idx[ee]]
    # print(dartboard_matrix[dd,:])

dartboard_matrix[:,0:257] = 0 
dartboard_matrix[:,257] = dartboard_matrix[:,258] 
dartboard_matrix[:,256] = dartboard_matrix[:,258] 

dartboard_matrix[:,0:256] = np.flip(dartboard_matrix[:,256:512]) 

fig, my_ax = plt.subplots()
pos = my_ax.imshow(dartboard_matrix) # pos = ax1.imshow(Zpos, cmap ='Blues', interpolation='none')
fig.colorbar(pos, ax=my_ax)
fig.show()

test = np.multiply(weight_matrix, dartboard_matrix)
fig, my_ax = plt.subplots()
pos = my_ax.imshow(weight_matrix) # pos = ax1.imshow(Zpos, cmap ='Blues', interpolation='none')
fig.colorbar(pos, ax=my_ax)
fig.show()

test = np.multiply(weight_matrix, dartboard_matrix)
fig, my_ax = plt.subplots()
pos = my_ax.imshow(test, cmap ='plasma') # pos = ax1.imshow(Zpos, cmap ='Blues', interpolation='none')
fig.colorbar(pos, ax=my_ax)
fig.show()




wm = np.ndarray((512, 512))

wm[np.where(dartboard_matrix==0)] = -50


n = np.arange(int(num_wedges))
x = np.sin((18*n-9)*np.pi/180)*256
y = np.cos((18*n-9)*np.pi/180)*256
# Wedge Coordinates
for ii in range(20):
    try:
        wedge_line_x = np.arange(0,x[ii],x[ii]/256)+256
    except:
        wedge_line_x = np.ones(256)*x[ii]+256
    try:
        wedge_line_y = np.arange(0,y[ii],y[ii]/256)+256
    except:
        wedge_line_y = np.ones(256)*y[ii]+256
    wm[wedge_line_x.astype(int),wedge_line_y.astype(int)] = wedge_points[ii]
    
    # if ii == 1:
    #     for r in range(512):
    #         weight_matrix[r, prev_wedge_line_y[r].astype(int):wedge_line_y[r].astype(int)] = ii + 1
    prev_wedge_line_x = wedge_line_x
    prev_wedge_line_y = wedge_line_y
    # plt.plot(wedge_line_x,prev_wedge_line_y)
####################################################
# plt.contour(weight_matrix)
# plt.show()
 
# Fill in weight matrix
for cc in range(512):
    cl_idx = np.squeeze(np.where(wm[:,cc]>0))
    # print(cl_idx)
    try:
        l_idx=len(cl_idx)
    except:
        l_idx=0
    for ee in range(1,l_idx):
        # print(str(cl_idx[ee-1]), '-', str(cl_idx[ee]), '=', weight_matrix[cl_idx[ee],cc])
        wm[cl_idx[ee-1]+1:cl_idx[ee], cc] = wm[cl_idx[ee], cc]
    # print(weight_matrix[:, cc])
# plt.contour(weight_matrix)
# plt.show()


fig, my_ax = plt.subplots()
pos = my_ax.imshow(wm, cmap ='plasma') # pos = ax1.imshow(Zpos, cmap ='Blues', interpolation='none')
fig.colorbar(pos, ax=my_ax)
fig.show()