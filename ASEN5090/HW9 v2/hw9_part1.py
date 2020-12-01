import numpy as np
import matplotlib.pyplot as plt
import pyproj
from hw9helpers import *
from readyuma import *

from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)

# ECEF and LLA conversions
ECEF = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
LLA = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')


# The .mat file contains 20 ms of 8-bit real samples. T
# data = loadmat("HW9data.mat")['gpsdata']

# 40° 0' 37.54877" N, 105° 14' 37.10560" W, 1613.976 m
roof_lat = 40.010430
roof_lon = -105.243640
roof_alt = 1613.976

roofECEF = pyproj.transform(LLA, ECEF, roof_lon, roof_lat, roof_alt)

# Grab satellite data
ephs, ephsdf = read_GPSyuma(r'data\YUMA310.ALM')
# print(type(ephsdf))
# Create time for entire day
sampling_frequency = 5e6 # 5 MHz
start_time = 20+21/60
time = np.linspace(start_time,start_time+20e-3,11)*60*60+86400*6
print(time)
# new Time
# t = np.arange(0,time_max,1/sampling_frequency)
tow = 24*6*3600+20*3600+21*60
print(tow)
print(roofECEF)

t_input = np.zeros((len(time), 2))

# Create input based off of time
t_input[:,0] = np.ones(np.shape(time))*310
t_input[:,1] = time

# Grab PRN Values
prns = ephsdf['prn'].to_numpy()
print(ephsdf['prn'])
roof_table = []
for pp in range(len(prns)):
    health, x = broadcast2pos(ephsdf, t_input, prns[pp])

    # print(np.shape(x))
    print('PRN: %s' % prns[pp])

    for rr in range(len(time)):
        satECEF = x[rr,:]

        [az_, el_, range_] = compute_azelrange(np.asarray(roofECEF), np.asarray(satECEF))
		
        print("Az: %s" % np.degrees(az_))
        print("Ek: %s" % np.degrees(el_))
        print("Range: %s" % range_)

        # print(az_, el_, range_)
        roof_table.append([int(prns[pp]), np.degrees(az_), np.degrees(el_),range_]) # Append values to a list of lists
        # print(current_table)

# NIST: Make the Lists of lists variables into arrays
roof_table = np.asarray(roof_table)
print(roof_table)
# Find the proper elevation values 
el_vals1 = roof_table[:,2]>0
el_vals2 = roof_table[:,2]<180
el_idx = el_vals1 & el_vals2
print(roof_table[el_idx,:])
print('\n')

# Reduce the table 
roof_table = roof_table[el_idx,:]

# Make the table print prettily
num_cols, num_rows = np.shape(roof_table)
for cc in range(num_cols):
    print(int(roof_table[cc, 0]), np.round(roof_table[cc, 1],6), np.round(roof_table[cc, 2],6), np.round(roof_table[cc, 3]))


fig1 = plt.figure()
ax1 = fig1.add_subplot(111, polar = True)
ax1.set_rlim(bottom = 90, top = 0)
ax1.set_theta_zero_location("N")
ax1.set_theta_direction(-1)
ax1.set_title("GPS Passes at Roof on Nov 7th 2020")

for prn in prns:
    print("PRN: %s" % prn)
	
    idx = np.where(roof_table[:,0]==prn)
    az_vals = np.squeeze(roof_table[idx, 1])
    el_vals = np.squeeze(roof_table[idx, 2])
    label = prn
    try:
        if len(az_vals) > 0:
            ax1.plot(np.deg2rad(az_vals), el_vals, label=label)
            # print("Az Vals: %s" % az_vals)
            print("El Vals: %s" % el_vals)
            print("###########################################")
            ax1.annotate(str(prn), np.deg2rad(az_vals[0]), el_vals[[0]] )
            plt.show()

    except:
        pass
# ax1.legend()

plt.show()
