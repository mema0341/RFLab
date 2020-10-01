import numpy as np
import pyproj
from readyuma import read_GPSyuma

# ephs, ephsdf = read_GPSyuma("data\YUMA245.ALM")

# print(ephs[0])

def lat2ecef(lat, lon, alt):

    R_plus = 6378.137 # km
    f = 1/298.2572235363
    eplus_squared = 2*f - np.power(f,2) # eccentricity

    C_plus = R_plus / np.sqrt(1-np.power(eccentricity,2)*np.power(np.sin(geodetic_lat),2))

    x = (C_plus+h) * np.cos(lat)*np.cos(lon)
    y = (C_plus+h) * np.cos(lat)*np.sing(lon)
    z = (C_plus*(1-eplus_squared)+alt)*np.sin(lat)

    return x, y, z