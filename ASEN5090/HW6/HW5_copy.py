from read_ephemeris import parse_rinex
from rinex_utilities import parse_rinex_nav_file, parse_rinex_obs_file
import numpy as np
from matplotlib import pyplot as plt
from readyuma import mean2eccentric, broadcast2pos, read_GPSyuma
import pandas as pd
from hw6_helpers import *
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)

##############################################################################################################
# Part 1
# Starting point: Choose a good PRN from the RINEX obs file to work with for this assignment. 
header, almanac_data = parse_rinex_nav_file('brdc2450.20n') # Grab data 
# print(almanac_data)
_, obs_data = parse_rinex_obs_file("nist2450.20o") # Grab even more data 

for alm_key in almanac_data.keys():

    print("Computing for PRN%s..." % str(alm_key))

    if alm_key < 10:
        my_satt = "G0%s" % alm_key
    else:
        my_satt = "G%s" % alm_key
    
    base = np.datetime64("2020-09-01T00:00:00") # Weird Time Thing
    seconds = (obs_data[my_satt].time - base)/1e6 + 172800 # Weird Time Thing
    weekComp = np.ones(len(seconds))*2121 # Weird Time Thing

    _, sattPos, _, _ = broadcast2posM(almanac_data,np.array([weekComp,seconds.astype(float)]).transpose(),alm_key) # Grab Satellite Position Using Revised Broadcast2Pos
    rxPos = get_nistECEF() # Grabs ECEF coords of NIST Location
    azelrange_rx2satt = pd.DataFrame(azelrange(rxPos, sattPos), columns = ['Azimuth', 'Elevation', 'Range']) # Az/El/Range between sat and rx

    # Compute the expected range for each of these times.
    satExpectedRange = np.zeros(len(seconds))
    bsv = np.zeros(len(seconds))
    relsv = np.zeros(len(seconds))
    for i in range(len(seconds)):
        satExpectedRange[i], _, _, bsv[i], relsv[i] = ExpectedRange(sattPos[i,:], rxPos,seconds[i].astype(float), almanac_data, alm_key)

    # Compute and plot the difference between the measured C/A code pseudorange and the
    # expected range (R in meters) versus time (hours) for Sept 1, 2020. What should it look like?
    C1 = obs_data[my_satt].signals["L1"].pr
    P2 = obs_data[my_satt].signals["L2"].pr

    dPR0 = C1 - satExpectedRange

    ##############################################################################################################
    # Part 2: Clock error: Write a function to compute the satellite clock correction (meters) based on the a0 and a1 terms in
    # the broadcast epheremis. Alternatively, you can add the SV clock correction capability to broadcast_eph2pos and provide bsv as an
    # additional output

    dPR1 = C1 - (satExpectedRange - bsv)

    ##############################################################################################################
    # Part 3: Relativity: Write a function or augment broadcast_eph2pos to compute the relativistic correction
    # (meters) based on the orbital elements in the ephemeris data
    dPR1 = C1 - (satExpectedRange - bsv)
    dPR2 = C1 - (satExpectedRange - bsv - relsv)

    ##############################################################################################################
    # Part 4: . Troposphere: Write a function to compute the simple tropospheric delay model (meters) based on the satellite
    # elevation and an assumed zenith delay for the ground station location. For the NIST ground station you may assume the 
    # zenith value is 2 m.
    zd = 2
    tropo = zd/np.sin(np.deg2rad(azelrange_rx2satt['Elevation'])) 

    # tropo = tropomodel(np.array([weekComp,seconds.astype(float)]).transpose(), zd=2)
    # Plot the tropo correction (meters) versus time (hours) for Sept 1, 2020. What should it look like?
    dPR1 = C1 - (satExpectedRange - bsv)
    dPR2 = C1 - (satExpectedRange - bsv - relsv)
    dPR3 = C1 - (satExpectedRange - bsv - relsv + tropo)

    ##############################################################################################################
    # Part 5: Ionosphere: For the ionospheric correction, rather than using a model-based correction (like Klobuchar) we
    # will combine GPS pseudorange measurements into an ionosphere free combination. 
    # Write a function to compute the iono-free combination of the pseudoranges reported in the RINEX obs file and
    # also the iono correction based on these two measurements. To allow your function to also be used for L5
    # measurements you can pass into it the carrier frequencies as well as the measurements.
    # [PRIF, iono] = ionocorr (C1, f1, P2, f2) *high precision people tend to call the iono-free combo P3

    L1 = 1575.42e6 
    L2 = 1227.60e6
    PRIF, ionocorr = Ionofree(C1, L1, P2, L2)

    fig, ax = plt.subplots()
    plt.scatter((seconds.astype(float) -172800)/3600 , ionocorr, 3)
    # ax.set_xlim(6, 12)
    # ax.set_ylim(min(tropo), max(tropo))
    ax.set_title("Iono Correction (m) PRN%s" % alm_key)
    ax.set_xlabel("Time (hours after Midnight)")
    ax.set_ylabel("Iono Correction (m)")
    ax.grid(True)
    # fig.savefig('Pics\PRN%s\PRN%sPLOT8.png' % (PRN, PRN))

    dPR1 = C1 - (satExpectedRange - bsv)
    dPR2 = C1 - (satExpectedRange - bsv - relsv)
    dPR3 = C1 - (satExpectedRange - bsv - relsv + tropo)
    dPR4 = PRIF - (satExpectedRange - bsv - relsv + tropo)
    print(PRIF[0])
    print(satExpectedRange[0])
    print(bsv[0])
    print(relsv[0])
    print(tropo[0])
    # Compute and plot dPR3 = C1 – (R – bsv - relsv + tropo)
    fig, ax = plt.subplots()
    ax.set_title("DPR (m) PRN%s" % alm_key)
    plt.scatter((seconds.astype(float) -172800)/3600 , dPR1, 3, label = "dPR1-CLOCK")
    plt.scatter((seconds.astype(float) -172800)/3600 , dPR2, 3, label = "dPR2-RELATIVITY")
    plt.scatter((seconds.astype(float) -172800)/3600 , dPR3, 3, label = "dPR3-TROPO")
    plt.scatter((seconds.astype(float) -172800)/3600 , dPR4, 3, label = "dPR4-IONO")
    plt.legend(markerscale=2.)
    ax.set_xlabel("Time (hours after Midnight)")
    ax.set_ylabel("Difference in Psuedo Range (m)")
    ax.grid(True)
    # fig.savefig('Pics\PRN%s\PRN%sPLOT9.png' % (PRN, PRN))

    plt.show()