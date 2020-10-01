import numpy as np
from pptx import Presentation
from pptx.util import Inches
import pyproj

ECEF = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
LLA = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

R_PLUS = 6378.137 # km
F_P = 1/298.2572235363
E_PLUS = 2*F_P - np.power(F_P,2) # eccentricity

# def lat2ecef(lat, lon, alt):

#     C_plus = R_PLUS / np.sqrt(1-np.power(E_PLUS,2)*np.power(np.sin(lat),2))

#     x = (C_plus+h) * np.cos(lat)*np.cos(lon)
#     y = (C_plus+h) * np.cos(lat)*np.sing(lon)
#     z = (C_plus*(1-E_PLUS)+alt)*np.sin(lat)

#     return x, y, z

# def ecef2lat(x, y, z):
#     lon_ = np.arctan(y/x)
#     row = np.sqrt(np.power(x,2)+np.power(y,2))

#     r = np.sqrt(np.power(x,2)+np.power(y,2)+np.power(z,2))

#     gd_old = np.arcsin(z/r)

#     C_plus = r_plus / np.sqrt(1-E_PLUS*np.power(np.sin(gd_old),2))

#     gd_new = np.arctan( (z + C_plus*E_PLUS*np.sin(gd_old)) / row)


def ecef2lat(x, y, z):
    lat_, lon_, alt_ = pyproj.transform(ecef, lla, x, y, z)

# ECEF3ENU function takes in lat and long as arguments and returns
# the ECEF coordinates using the conversion matrix found on slide 35 of Lectures 6&7
def ecef2enu(lon_, lat_):
    lambda_ = lon_
    phi_ = lat_

    # Define the conversion matrix
    conversion_matrix = \
        [[-np.sin(lambda_), np.cos(lambda_), 0],\
         [-np.sin(phi_)*np.cos(lambda_), -np.sin(phi_)*np.sin(lambda_), np.cos(phi_)],\
         [np.cos(phi_)*np.cos(lambda_), np.cos(phi_)*np.sin(lambda_), np.sin(phi_)]]
    
    # To get the solution, take the dot product of the ECEF coords
    # and the conversion matrix
    return conversion_matrix

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
    # Returns Az [degrees], El [degrees], and X_up [m]
    return az_, el_, range