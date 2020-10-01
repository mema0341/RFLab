import numpy as np
import pyproj
from hw3helpers import ecef2enu

ECEF = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
LLA = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

# Convert ECEF to LLA 
x = -1288398
y = -4721697
z = 4078625

nist_lon, nist_lat, nist_alt = pyproj.transform(ECEF, LLA, x, y, z, radians=True)

# conversion_matrix = ecef2enu(nist_lon, nist_lat)
conversion_matrix = ecef2enu(nist_lon, nist_lat)
print(conversion_matrix)

print(np.dot(conversion_matrix,[x,y,z]))

# test = ecef2enu(0,0)
# test = ecef2enu(90,0) # North Pole
# test = ecef2enu(0,90) # North Pole