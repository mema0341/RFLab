import numpy as np
import pyproj
from hw3helpers import ecef2enu

ECEF = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
LLA = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

# Write a function that takes as input ECEF position of a user on the ground and the ECEF position of a GPS
# satellite, and outputs the LOS unit vector pointing from the user to the satellite, expressed in ENU coordinates
# wrt the user location on the Earth. Make up an example to check your function.

# COMPUTE_LOS_ENU takes in the user position and satellite position in ECEF and outputs the LOS unit vector pointing from the user to the satellite in ENU coords
def compute_LOS_ENU(userECEF, satECEF):
    # Part A: Compute the vector pointing from user location to GPS
    vector = satECEF - userECEF
    # Part B: Normalize it to a unit vector (LOS_ECEF)
    ehat_los = vector / np.sqrt(np.power(vector[0], 2) \
         + np.power(vector[1], 2) + np.power(vector[2], 2))
    # Part C: Transform the Unit Vector from ECEF to ENU using the function you wrote in part 2
    # First Grab Lat and lon values 
    lon_, lat_, alt_ = pyproj.transform(ECEF, LLA, userECEF[0], userECEF[1], userECEF[2], radians=True)
    # Now grab conversion matrix
    conversion_matrix = ecef2enu(lon_, lat_)
    # Now take the dot product to get the ENU coords
    enu_coords = np.dot(conversion_matrix,ehat_los)
    # Return ENU Coords
    return enu_coords

# COMPUTE_AZELRANGE function takes ECEF coords of a satellite and a user, converts
# to get the LOS vector in ENU, and returns az, el, and range of the satellite
def compute_azelrange(userECEF, satECEF):
    # First get the LOS in ENU coords
    x_east, x_north, x_up = compute_LOS_ENU(userECEF, satECEF)

    # Compute Azimuth using LOS 
    az_ = np.arctan2(x_east,x_north)

    # Compute Elevation using LOS 
    el_ = np.arcsin(x_up / np.sqrt(np.power(x_east, 2)\
        +np.power(x_north, 2)+np.power(x_up, 2)))

    # Compute Range
    diff = userECEF-satECEF
    range = np.sqrt(np.power(diff[0],2) + np.power(diff[1],2) + np.power(diff[2],2))
    # Returns Az [degrees], El [degrees], and rabg [m]
    return az_, el_, range

# Part 4
# [az_, el_, range_] = compute_azelrange(userECEF, satECEF)

# Convert NIST ECEF coordinates  to LLA 
x = -1288398
y = -4721697
z = 4078625

# Obtain nist lat/long/altitude 
nist_lon, nist_lat, nist_alt = pyproj.transform(ECEF, LLA, x, y, z, radians=False)

# Put the user at NIST location and obtain the ECEF coords of the satellite 
satECEF = pyproj.transform(LLA, ECEF, nist_lon, nist_lat, nist_alt+400e3)
userECEF = np.asarray([-1288398, -4721697, 4078625])

# Compute LOS Unit Vector
los_enu = compute_LOS_ENU(userECEF, satECEF)

print(los_enu)

# # Obtain azimuth, elevation, and range of s
# [az_, el_, range_] = compute_azelrange(userECEF, satECEF)

# print(np.degrees(az_))
# print(np.degrees(el_))
# print(range_)
