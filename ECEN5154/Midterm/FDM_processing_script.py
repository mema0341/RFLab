import numpy as np 
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pylab as plt

# Definitions
mu0 = 4*pi*1e-7
eps0 = 8.85e-12
epsr = 1


def plot_matrix(matrix, title):
    # plt.figure()
    fig, my_ax = plt.subplots()
    im = my_ax.imshow(matrix, cmap='jet') # pos = ax1.imshow(Zpos, cmap ='jet', interpolation='none')
    # divider = make_axes_locatable(my_ax)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # fig.colorbar(im, cax=cax)
    plt.title(title)
    # plt.tight_layout()
    plt.show()

# fname = "FDMSolutions\wb1.0_stepsize0.1.xlsx"
# data = np.asanyarray(pd.read_excel(fname))
# print(np.shape(data))
# plot_matrix(data, "Potential Distribution\nW/b = 1.0 Nx = 20, Ny = 80")

# fname = "FDMSolutions\wb1.5_stepsize0.1.xlsx"
# data = np.asanyarray(pd.read_excel(fname))
# print(np.shape(data))
# plot_matrix(data, "Potential Distribution\nW/b = 1.5 Nx = 20, Ny = 80")

# fname = "FDMSolutions\wb2.0_stepsize0.1.xlsx"
# data = np.asanyarray(pd.read_excel(fname))
# print(np.shape(data))
# plot_matrix(data, "Potential Distribution\nW/b = 2.0 Nx = 20, Ny = 80")

# fname = "FDMSolutions\wb2.5_stepsize0.1.xlsx"
# data = np.asanyarray(pd.read_excel(fname))
# print(np.shape(data))
# plot_matrix(data, "Potential Distribution\nW/b = 2.0 Nx = 20, Ny = 80")

    
fname = r"FDMSolutions\wb1.5_stepsize0.05.xlsx"
data = np.asanyarray(pd.read_excel(fname))
print(np.shape(data))
plot_matrix(data, "Potential Distribution\nW/b = 2.0 Nx = 20, Ny = 80")

fname = r"FDMSolutions\wb2.0_stepsize0.05.xlsx"
data = np.asanyarray(pd.read_excel(fname))
print(np.shape(data))
plot_matrix(data, "Potential Distribution\nW/b = 2.0 Nx = 20, Ny = 80")

fname = r"FDMSolutions\wb2.5_stepsize0.05.xlsx"
data = np.asanyarray(pd.read_excel(fname))
print(np.shape(data))
plot_matrix(data, "Potential Distribution\nW/b = 2.0 Nx = 20, Ny = 80")
