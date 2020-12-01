from hw6_helpers import *

from read_ephemeris import parse_rinex
from rinex_utilities import parse_rinex_nav_file, parse_rinex_obs_file
import numpy as np
from matplotlib import pyplot as plt
from readyuma import mean2eccentric, broadcast2pos, read_GPSyuma
import pandas as pd
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)

c = 299792458 
##############################################################################################################
# Part 1: Use the given NIST location as the initial guess of the user position. 
# init_guess = get_nistECEF() # Grabs ECEF coords of NIST Location

# Part 7: Now, instead of using the correct location for NIST, use an initial estimate that is way off - for example try a
# location in Washington DC, or the center of the earth. Compute the A matrix and the pre-fit residuals based on this
# assumed position. Compute the least squares correction and iterate your solution 5 times. For each iteration print
# the position correction, and your new estimate. Does it converge to the correct answer?
# washington_lat = 38 #°53′42″ N
# washington_lon = 77 #°02′10″ W
# washington_h = 6
# init_guess = np.asarray(llatoECEF(washington_lat, washington_lon, washington_h))
init_guess = np.asarray([1114105.32957155, -4844619.21156437, 3982733.05752002])

##############################################################################################################

# header, almanac_data = parse_rinex_nav_file('brdc2450.20n') # Grab data 
_, obs_data = parse_rinex_obs_file("nist2450.20o") # Grab even more data 

timeN = 172800


# Part 2: Read the first epoch (set) of observation records. (One time or epoch, all the observed satellites.) For each of the
# observed satellites print and check the following values:
# a. the ionosphere free pseudorange observable PIF (using C1 and P2) in meters
pif = []
# b. the expected range (from the last homework) in meters
exp_r = []
# c. the elevation and azimuth in degrees (based on the adjusted satellite location at time of transmission)
az_vals = []
el_vals = []
# d. the satellite clock correction (from the a0 and a1 terms in the nav message) in meters
clock_corr = []
# e. the relativistic correction (from the orbital elements or r,v of the GPS satellite) in meters
rel_corr = []
# f. the simple tropospheric model (used in HW5) in meters
tropo_vals = []

ecef_vals = []

base = np.datetime64("2020-09-01T00:00:00")

epoch = np.datetime64("2020-09-01T00:00:00")
for my_satt in obs_data:
    for time in obs_data[my_satt].time:
        print(time)

#     base = np.datetime64("2020-09-01T00:00:00") # Weird Time Thing
#     seconds = (obs_data[my_satt].time - base)/1e6 + 172800 # Weird Time Thing
#     weekComp = np.ones(len(seconds))*2121 # Weird Time Thing
#     simplified_time = np.array([weekComp,seconds.astype(float)]).transpose()

#     # Loop Through PRNs
#     _, sattPos, _, _ = broadcast2posM(almanac_data,simplified_time,alm_key) # Grab Satellite Position Using Revised Broadcast2Pos

#     # rxPos = init_guess # this is redundant
#     azelrange_rx2satt = pd.DataFrame(azelrange(init_guess, sattPos), columns = ['Azimuth', 'Elevation', 'Range']) # Az/El/Range between sat and rx

#     # Compute the expected range for each of these times.
#     seconds = [timeN]
#     satExpectedRange = np.zeros(len(seconds))
#     bsv = np.zeros(len(seconds))
#     relsv = np.zeros(len(seconds))
#     azs = np.zeros(len(seconds))
#     els = np.zeros(len(seconds))
#     ecefs = np.zeros((len(seconds), 3))
#     for i in range(len(seconds)):
#         satExpectedRange[i], azs[i], els[i], ecefs[i,:], bsv[i], relsv[i] = ExpectedRange(sattPos[i,:], init_guess,seconds[i], almanac_data, alm_key)

#     C1 = obs_data[my_satt].signals["L1"].pr
#     P2 = obs_data[my_satt].signals["L2"].pr
#     L1 = 1575.42e6 
#     L2 = 1227.60e6

#     zd = 2
#     tropo = zd/np.sin(np.deg2rad(azelrange_rx2satt['Elevation'])) 
#     PRIF, ionocorr = Ionofree(C1, L1, P2, L2)

#     dPR0 = C1 - satExpectedRange
#     dPR1 = C1 - (satExpectedRange - bsv)
#     dPR2 = C1 - (satExpectedRange - bsv - relsv)
#     dPR3 = C1 - (satExpectedRange - bsv - relsv + tropo)
#     dPR4 = PRIF - (satExpectedRange - bsv - relsv + tropo)

#     # Grab the values you need
#     pif.append(PRIF[0])
#     exp_r.append(satExpectedRange[0])
#     clock_corr.append(bsv[0])
#     rel_corr.append(relsv[0])
#     tropo_vals.append(tropo[0])
#     az_vals.append(azs[0])
#     el_vals.append(els[0])

#     ecef_vals.append(ecefs[0])
#     # Get az/el values
#     rev_azelrange_rx2satt = pd.DataFrame(azelrange(init_guess, sattPos), columns = ['Azimuth', 'Elevation', 'Range']) # Az/El/Range between sat and rx

# pif = np.asarray(pif)
# exp_r = np.asarray(exp_r)
# az_vals = np.asarray(az_vals)
# el_vals = np.asarray(el_vals)
# clock_corr = np.asarray(clock_corr)
# rel_corr = np.asarray(rel_corr)
# tropo_vals = np.asarray(tropo_vals)
# ecef_vals = np.asarray(ecef_vals)

# print('Iteration %s...' % str(zz+1))
# g_matrix = np.ones((len(ecef_vals), 4))
# g_matrix[:,0] = -(ecef_vals[:,0] - init_guess[0]) / exp_r
# g_matrix[:,1] = -(ecef_vals[:,1] - init_guess[1]) / exp_r
# g_matrix[:,2] = -(ecef_vals[:,2] - init_guess[2]) / exp_r
# # print(g_matrix)

# # Part 4: Compute and print the prefit residuals (dy) correcting for ionospheric delay, satellite clock, relativity and
# # tropospheric delay.
# dy = pif - (exp_r - clock_corr - rel_corr + tropo_vals)
# # print(dy)

# # Part 5: Plot these residuals as a function of the satellite elevation angle. If it looks like there is a dependence on
# # elevation, maybe the assumed zenith tropospheric delay needs to be adjusted.
# # fig, ax = plt.subplots()
# # plt.scatter(el_vals, dy)
# # ax.set_ylabel('Prefit Residuals [m]')
# # ax.set_xlabel('Elevation [Deg.]')
# # plt.title('Prefit Residuals [m]')
# # plt.grid(True)
# # plt.show()
# # fig.savefig('Pics\HW6_Plot4.png')


# # Part 6: Form the least squares solution for “Dx” the state vector comprising the three position components in ECEF Dx,
# # Dy, Dz, and the receiver clock error b. Report each of the components to a precision of 1cm. If the assumed
# # location is correct, this will be the point solution error
# state_vector, _, _, _ = np.linalg.lstsq(g_matrix, dy, rcond=None)

# print('x error correction: %s m' % np.round(state_vector[0],2))
# print('y error correction: %s m' % np.round(state_vector[1],2))
# print('z error correction: %s m' % np.round(state_vector[2],2))
# print('b error correction: %s m' % np.round(state_vector[3],2))

# init_guess[0] = init_guess[0] + state_vector[0]
# init_guess[1] = init_guess[1] + state_vector[1]
# init_guess[2] = init_guess[2] + state_vector[2]
# seconds[0] = seconds[0] + state_vector[3]/c

# print('New Guess: ')
# print('x Component: %s m' % np.round(init_guess[0],2))
# print('Y Component: %s m' % np.round(init_guess[1],2))
# print('Z Component: %s m' % np.round(init_guess[2],2))    
# print("#######################################################################")

# # Compute the transformation matrix from ECEF to ENU for the assumed location. Convert your position error
# # from ECEF to ENU coordinates and report these values to a precision of 1cm
# lat_, lon_, _ = ECEFtolla(state_vector[0:3])
# r_L = ECEF_ENU_Transform(lat_, lon_)

# tilda_r_L = np.zeros((4,4))
# tilda_r_L[0:3,0:3] = r_L
# tilda_r_L[3,3] = 1
# enu_state_vector = np.dot(tilda_r_L, state_vector)
# print(np.round(enu_state_vector,2))