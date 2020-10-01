import numpy as np 
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pylab as plt

def plot_matrix(matrix, title):
    fig, my_ax = plt.subplots()
    im = my_ax.imshow(matrix, cmap='jet') # pos = ax1.imshow(Zpos, cmap ='jet', interpolation='none')
    plt.title(title)
    plt.show()

# Definitions
mu0 = 4*np.pi*1e-7
eps0 = 8.85e-12
epsr = 1

# Open EXCEL File
fname = r"FDMSolutions\wb2.5_stepsize0.05.xlsx"
potential = np.asanyarray(pd.read_excel(fname))
x,y = np.shape(potential)

# Plot Potenital
plot_matrix(potential, "Potential Distribution\nW/b = 1.0 Nx = %s, Ny = %s" % (str(x), str(y)))

# Calculate capacitance and characteristic impedance
Nx=40
strip_x, strip_y  = np.where(potential==1)
above_strip = np.sum(potential[strip_x+2, strip_y]) - np.sum(potential[strip_x+1, strip_y])
below_strip = np.sum(potential[strip_x-2, strip_y]) - np.sum(potential[strip_x-1, strip_y])
left_strip = np.sum(potential[[strip_x[0]-1, strip_x[0], strip_x[0]+1], strip_y[0]-2])-np.sum(potential[[strip_x[0]-1, strip_x[0], strip_x[0]+1], strip_y[0]-1])
right_strip = np.sum(potential[[strip_x[0]-1, strip_x[0], strip_x[0]+1], strip_y[len(strip_y)-1]+2])-np.sum(potential[[strip_x[0]-1, strip_x[0], strip_x[0]+1], strip_y[len(strip_y)-1]+1])
total = above_strip+below_strip+left_strip+right_strip
c = abs(total*eps0)/Nx
c0 = abs(total*eps0)
print(c)
c = 3e8
Z0 = 1/(c*c0)
print(Z0)
