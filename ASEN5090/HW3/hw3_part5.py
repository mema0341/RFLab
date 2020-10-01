import numpy as np
import matplotlib.pyplot as plt
import pyproj
from pptx import Presentation
from pptx.util import Inches
from hw3helpers import *

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

ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
nist_lon, nist_lat, nist_alt = pyproj.transform(ecef, lla, x, y, z, radians=False)

nistECEF = [x,y,z]

# Make sure we get the NIST ECEF coords back
test = pyproj.transform(lla, ecef, nist_lon, nist_lat, nist_alt)
# print(test)

# Get EQUA Things
equa_lat = 0
equa_lon = nist_lon
equa_alt = 0 
equaECEF = np.asarray(pyproj.transform(lla, ecef, equa_lon, equa_lat, equa_alt))
# print(equaECEF)


# Test func
conversion_matrix = ecef2enu(nist_lon, nist_lat)
enu_coords = np.dot([x, y, z], conversion_matrix)
# print(enu_coords)

# Get 89no 
no89_lat = 89 # Wait what is N89 deg
no89_lon = nist_lon
no89_alt = 0
no89ECEF = np.asarray(pyproj.transform(lla, ecef, no89_lon, no89_lat, no89_alt))
# print(no89ECEF)

file_ = open("data\gpspos20200901.txt", "r")

nist_table = []
enu_table = []
no89_table = []
for line in file_.readlines():    
    data = line.split()
    prn_code = float(data[0])
    sat_x = float(data[1])*1000
    sat_y = float(data[2])*1000
    sat_z = float(data[3])*1000

    satECEF = [sat_x, sat_y, sat_z]

    print(sat_x)
    print(sat_y)
    print(sat_z)
    
    # NIST COORDS
    [az_, el_, range_] = compute_azelrange(np.asarray(nistECEF), np.asarray(satECEF))
    nist_table.append([int(prn_code), np.degrees(az_), np.degrees(el_),range_]) # Append values to a list of lists

    # EQUA COORDS
    [az_, el_, range_] = compute_azelrange(np.asarray(equaECEF), np.asarray(satECEF))
    enu_table.append([int(prn_code), np.degrees(az_), np.degrees(el_),range_]) # Append values to a list of lists

    # 89no COORDS
    [az_, el_, range_] = compute_azelrange(np.asarray(no89ECEF), np.asarray(satECEF))
    no89_table.append([int(prn_code), np.degrees(az_), np.degrees(el_),range_]) # Append values to a list of lists 

file_.close()

# NIST: Make the Lists of lists variables into arrays
nist_table = np.asarray(nist_table)
enu_table = np.asarray(enu_table)
no89_table = np.asarray(no89_table)

# Find the proper elevation values 
el_vals1 = nist_table[:,2]>0
el_vals2 = nist_table[:,2]<180
el_idx = el_vals1 & el_vals2
print(nist_table[el_idx,:])
print('\n')

# Reduce the table 
nist_table = nist_table[el_idx,:]

# Make the table print prettily
num_cols, num_rows = np.shape(nist_table)
for cc in range(num_cols):
    print(int(nist_table[cc, 0]), np.round(nist_table[cc, 1],6), np.round(nist_table[cc, 2],6), np.round(nist_table[cc, 3]))

# EQUATOR: Make the Lists of lists variables into arrays
el_vals1 = enu_table[:,2]>0
el_vals2 = enu_table[:,2]<180
el_idx = el_vals1 & el_vals2
print(enu_table[el_idx,:])
print('\n')

# Reduce the table 
enu_table = enu_table[el_idx,:]

# Make the table print prettily
num_cols, num_rows = np.shape(enu_table)
for cc in range(num_cols):
    print(int(enu_table[cc, 0]), np.round(enu_table[cc, 1],6), np.round(enu_table[cc, 2],6), np.round(enu_table[cc, 3]))
print('\n')

# 89 degrees north: Make the Lists of lists variables into arrays
el_vals1 = no89_table[:,2]>0
el_vals2 = no89_table[:,2]<180
el_idx = el_vals1 & el_vals2

# Reduce the table 
no89_table = no89_table[el_idx,:]

# Make the table print prettily
num_cols, num_rows = np.shape(no89_table)
for cc in range(num_cols):
    print(int(no89_table[cc, 0]), np.round(no89_table[cc, 1],6), np.round(no89_table[cc, 2],6), np.round(no89_table[cc, 3]))

