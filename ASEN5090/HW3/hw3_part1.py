import numpy as np
import matplotlib.pyplot as plt
import pyproj
from pptx import Presentation
from pptx.util import Inches
from hw3helpers import send_to_table

# Global constants 
MU = 3.986005e14 # m^3 * s^-2
EARTH_ROT_RATE = 7.2921151467e-5 # rad / sec

ecef_coords = np.asarray([-1288398, -4721697, 4078625])

# Goal: Make a table of ECEF and Lat/Long/Height Positions (wrt WGS84 Ellispoid)
#   Express in meters or degrees with 1m precision 

# Convert ECEF to LLA 
x = -1288398
y = -4721697
z = 4078625

lambda_ = np.arctan(y/x)
print(lambda_)

row = np.sqrt(np.power(x,2)+np.power(y,2))
print(row)

r = np.sqrt(np.power(x,2)+np.power(y,2)+np.power(z,2))

init_guess = np.arcsin(z/r)

c_plus = r_plus