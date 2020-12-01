from hw6_helpers import *

from read_ephemeris import parse_rinex
from rinex_utilities import parse_rinex_nav_file, parse_rinex_obs_file
import numpy as np
from matplotlib import pyplot as plt
from readyuma import mean2eccentric, broadcast2pos, read_GPSyuma
import pandas as pd
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)

# MY_SATT = "G02"
MY_SATT = "G30"
PRN = int(MY_SATT[1:len(MY_SATT)]) 

##############################################################################################################
# Part 1: Use the given NIST location as the initial guess of the user position. 
nist_loc = get_nistECEF() # Grabs ECEF coords of NIST Location

##############################################################################################################
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

header, almanac_data = parse_rinex_nav_file('brdc2450.20n') # Grab data 
_, obs_data = parse_rinex_obs_file("nist2450.20o") # Grab even more data 

base = np.datetime64("2020-09-01T00:00:00") # Weird Time Thing
seconds = (obs_data[MY_SATT].time - base)/1e6 + 172800 # Weird Time Thing
weekComp = np.ones(len(seconds))*2121 # Weird Time Thing


for kk in range(len(almanac_data)):
    # Loop Through PRNs
    simplified_time = np.array([weekComp,seconds.astype(float)]).transpose()
    _, sattPos, _, _ = broadcast2posM(almanac_data,simplified_time,PRN) # Grab Satellite Position Using Revised Broadcast2Pos

    rxPos = get_nistECEF() # Grabs ECEF coords of NIST Location
    azelrange_rx2satt = pd.DataFrame(azelrange(rxPos, sattPos), columns = ['Azimuth', 'Elevation', 'Range']) # Az/El/Range between sat and rx

    # Compute the expected range for each of these times.
    satExpectedRange = np.zeros(len(seconds))
    bsv = np.zeros(len(seconds))
    relsv = np.zeros(len(seconds))
    for i in range(len(seconds)):
        satExpectedRange[i], bsv[i], relsv[i] = ExpectedRange(sattPos[i,:], rxPos,seconds[i].astype(float), almanac_data, PRN)

    C1 = obs_data[MY_SATT].signals["L1"].pr
    P2 = obs_data[MY_SATT].signals["L2"].pr
    L1 = 1575.42e6 
    L2 = 1227.60e6

    zd = 2
    tropo = zd/np.sin(np.deg2rad(azelrange_rx2satt['Elevation'])) 
    PRIF, ionocorr = Ionofree(C1, L1, P2, L2)

    dPR0 = C1 - satExpectedRange
    dPR1 = C1 - (satExpectedRange - bsv)
    dPR2 = C1 - (satExpectedRange - bsv - relsv)
    dPR3 = C1 - (satExpectedRange - bsv - relsv + tropo)
    dPR4 = PRIF - (satExpectedRange - bsv - relsv + tropo)

    # Grab the values you need
    pif.append(PRIF[0])
    exp_r.append(satExpectedRange[0])

    fig, ax = plt.subplots()
    ax.set_title("DPR (m) PRN%s" % PRN)
    plt.scatter((seconds.astype(float) -172800)/3600 , dPR1, 3, label = "dPR1-CLOCK")
    plt.scatter((seconds.astype(float) -172800)/3600 , dPR2, 3, label = "dPR2-RELATIVITY")
    plt.scatter((seconds.astype(float) -172800)/3600 , dPR3, 3, label = "dPR3-TROPO")
    plt.scatter((seconds.astype(float) -172800)/3600 , dPR4, 3, label = "dPR4-IONO")
    plt.legend(markerscale=2.)
    ax.set_xlabel("Time (hours after Midnight)")
    ax.set_ylabel("Difference in Psuedo Range (m)")
    ax.grid(True)
    plt.show()