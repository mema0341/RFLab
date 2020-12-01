import numpy as np
import matplotlib.pyplot as plt
import pyproj
from hw9helpers import *
from readyuma import *

from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)

def broadcast_eph2pos(eph, timevec, prn):
    tow = timevec[0][1]

    muE = 3.986005e14     # WGS-84 value, m^3/s^2
    OmegaE  = 7.2921151467e-5 # WGS-84 value, rad/s 

    # Get satellite data
    M0 = ephsdf['M0'][prn]
    dn = ephsdf['delta_n'][prn]
    ecc = ephsdf['ecc'][prn]
    sqrta = ephsdf['sqrt_a'][prn]
    Loa = ephsdf['Loa'][prn]
    inc = ephsdf['incl'][prn]
    w = ephsdf['perigee'][prn]
    ra_rate = ephsdf['ra_rate'][prn]
    i_rate = ephsdf['i_rate'][prn]
    Cuc = ephsdf['Cuc'][prn]
    Cus = ephsdf['Cus'][prn]
    Crc = ephsdf['Crc'][prn]
    Crs = ephsdf['Crs'][prn]
    Cic = ephsdf['Cic'][prn]
    Cis = ephsdf['Cis'][prn]

    if prn == 30:
        Toe = ephsdf['Toa'][prn]
    else:
        Toe = ephsdf['toa'][prn]

    a = np.power(sqrta,2)
    n0 = np.sqrt(muE/np.power(a,3))
    dt = tow - Toe
    n = n0 + dn

    M   = M0 + n * dt #mean anomaly
    E = mean2eccentric(M,ecc)
    cosE = np.cos(E)
    sinE = np.sin(E)

    nu = np.arctan2(np.sqrt(1. - ecc**2) * sinE,  cosE - ecc) 
    phi = nu + w
    du = Cus*np.sin(2*phi) + Cuc*np.cos(2*phi)
    dr = Crs*np.sin(2*phi) + Crc*np.cos(2*phi)
    di = Cis*np.sin(2*phi) + Cic*np.cos(2*phi)
    u = phi + du
    r = a*(1-ecc*np.cos(E)) + dr
    i = inc + di + i_rate*dt

    xp = r*np.cos(u)
    yp = r*np.sin(u)

    Omega = Loa + (ra_rate-OmegaE)*dt - OmegaE*Toe

    # x,y,z
    x = xp*np.cos(Omega) - yp*np.cos(i)*np.sin(Omega)
    y = xp*np.sin(Omega) + yp*np.cos(i)*np.cos(Omega)
    z = yp*np.sin(i)

    # Output
    pos = [x,y,z]
    health = ephsdf['health'][prn]

    return health, pos

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

start_time = 20+21/60
time = [start_time*60*60+86400*6]
t_input = np.zeros((len(time), 2))

# Create input based off of time
t_input[:,0] = np.ones(np.shape(time))*310
t_input[:,1] = time

roof_table = []
prns = ephsdf['prn'].to_numpy()
for pp in range(len(prns)):
    # print('PRN: %s' % prns[pp])

    health, satECEF = broadcast_eph2pos(ephsdf, t_input, pp)

    for rr in range(len(time)):
        # satECEF = x[rr,:]

        [az_, el_, range_] = compute_azelrange(np.asarray(roofECEF), np.asarray(satECEF))
		
        # print("Az: %s" % np.degrees(az_))
        # print("El: %s" % np.degrees(el_))
        # print("Range: %s" % range_)
        # print("#######################################3")
        # print(az_, el_, range_)
        roof_table.append([int(prns[pp]), np.degrees(az_), np.degrees(el_),range_]) # Append values to a list of lists
        # print(current_table)
# NIST: Make the Lists of lists variables into arrays
roof_table = np.asarray(roof_table)

# Find the proper elevation values 
el_vals1 = roof_table[:,2]>0
el_vals2 = roof_table[:,2]<180
el_idx = el_vals1 & el_vals2


# Reduce the table 
roof_table = roof_table[el_idx,:]

# # Make the table print prettily
# num_cols, num_rows = np.shape(roof_table)
# for cc in range(num_cols):
#     print(int(roof_table[cc, 0]), np.round(roof_table[cc, 1],6), np.round(roof_table[cc, 2],6), np.round(roof_table[cc, 3]))

fig1 = plt.figure()
ax1 = fig1.add_subplot(111, polar = True)
ax1.set_rlim(bottom = 90, top = 0)
ax1.set_theta_zero_location("N")
ax1.set_theta_direction(-1)
ax1.set_title("GPS Visible Satellites at Roof\nNov 7th 2020 at 20:21")
ax1.set_yticks(np.arange(0,90,20))
for prn in roof_table[:,0]:
    print("PRN: %s" % prn)

    idx = np.where(roof_table[:,0]==prn)
    az_val = roof_table[idx, 1][0][0]
    el_val = roof_table[idx, 2][0][0]

    ax1.scatter(np.deg2rad(az_val), el_val, label=prn)
    print("Az Vals: %s" % az_val)
    print("El Vals: %s" % el_val)
    print("###########################################")
    ax1.annotate(str(int(prn)), [np.deg2rad(az_val), el_val] )
    # plt.show()
plt.tight_layout()
# ax1.legend()
plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111, polar = True)
# ax.set_rlim(bottom = 90, top = 0)
# # ax1.set_theta_zero_location("N")
# # ax1.set_theta_direction(-1)
# # ax1.set_title("GPS Passes at Roof on Nov 7th 2020")
# ax.scatter([45,45,45], [0, 45, 90], label='what')
# plt.show()