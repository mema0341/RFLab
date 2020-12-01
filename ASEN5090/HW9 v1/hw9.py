from hw9_helpers import * 

import numpy as np
from matplotlib import pyplot as plt

import pyproj
import pandas as pd
from scipy.io import loadmat

from read_ephemeris import parse_rinex
from rinex_utilities import parse_rinex_nav_file, parse_rinex_obs_file
from readyuma import mean2eccentric, broadcast2pos, read_GPSyuma

from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)

ECEF = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
LLA = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

# The .mat file contains 20 ms of 8-bit real samples. T
data = loadmat("HW9data.mat")['gpsdata']

# 40° 0' 37.54877" N, 105° 14' 37.10560" W, 1613.976 m
roof_lat = 40.010430
roof_lon = -105.243640
roof_alt = 1613.976

roof_ECEF = pyproj.transform(LLA, ECEF, roof_lon, roof_lat, roof_alt)

head, a = parse_rinex_nav_file('brdc2450.20n')
timevec = np.arange(172800, 172800+86401, 31)
weeks = np.ones(len(timevec))*2121
 
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, polar = True)
ax1.set_rlim(bottom = 90, top = 0)
ax1.set_theta_zero_location("N")
ax1.set_theta_direction(-1)
ax1.set_title("GPS Passes at NIST on Sept 1 2020")
# plotting through the day using the broadcast2pos function.
for i in a:
	#print(i)
	# print(np.array([[somezeros],[times]]))
	[h, poses]=broadcast2posM(a,np.array([weeks,timevec]).transpose(),i)
	aer1 = pd.DataFrame(azelrange(roof_ECEF, poses), columns = ['Azimuth', 'Elevation', 'Range'])
	
	ax1.plot(np.deg2rad(aer1['Azimuth']), aer1['Elevation'], label = "PRN {}".format(i))
	
	#aer1 = pd.DataFrame(azelrange(userECEFs[2], poses), columns = ['Azimuth', 'Elevation', 'Range'])

ax1.legend()
plt.show()