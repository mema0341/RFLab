import numpy as np
import matplotlib.pyplot as plt
import pyproj
from hw3helpers import *

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

print(nist_lat)
print(nist_lon)
print(nist_alt)

# Test func

# Make sure we get the NIST ECEF coords back
test = pyproj.transform(lla, ecef, nist_lon, nist_lat, nist_alt)
# print(test)

# Get EQUA Things
equa_lat = 0
equa_lon = nist_lon
equa_alt = 0 
equa_ecef_coords = np.asarray(pyproj.transform(lla, ecef, equa_lon, equa_lat, equa_alt))
print(equa_ecef_coords)


# Test func
conversion_matrix = ecef2enu(nist_lon, nist_lat)
enu_coords = np.dot([x, y, z], conversion_matrix)
print(enu_coords)

# Get 89no 
no89_lat = 89 # Wait what is N89 deg
no89_lon = nist_lon
no89_alt = 0
no89_ecef_coords = np.asarray(pyproj.transform(lla, ecef, no89_lon, no89_lat, no89_alt))
print(no89_ecef_coords)




