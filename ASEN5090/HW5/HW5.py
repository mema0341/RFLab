from read_ephemeris import parse_rinex
from rinex_utilities import parse_rinex_nav_file, parse_rinex_obs_file
import numpy as np
from matplotlib import pyplot as plt
from readyuma import mean2eccentric, broadcast2pos, read_GPSyuma
import pandas as pd
from hw5_helpers import *
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)

MY_SATT = "G02"
# MY_SATT = "G30"
PRN = int(MY_SATT[1:len(MY_SATT)]) 


##############################################################################################################
# Part 1
# Starting point: Choose a good PRN from the RINEX obs file to work with for this assignment. 
header, almanac_data = parse_rinex_nav_file('brdc2450.20n') # Grab data 
# print(almanac_data)
_, obs_data = parse_rinex_obs_file("nist2450.20o") # Grab even more data 

base = np.datetime64("2020-09-01T00:00:00") # Weird Time Thing
seconds = (obs_data[MY_SATT].time - base)/1e6 + 172800 # Weird Time Thing
weekComp = np.ones(len(seconds))*2121 # Weird Time Thing

_, sattPos, _, _ = broadcast2posM(almanac_data,np.array([weekComp,seconds.astype(float)]).transpose(),PRN) # Grab Satellite Position Using Revised Broadcast2Pos
rxPos = get_nistECEF() # Grabs ECEF coords of NIST Location
azelrange_rx2satt = pd.DataFrame(azelrange(rxPos, sattPos), columns = ['Azimuth', 'Elevation', 'Range']) # Az/El/Range between sat and rx

# Create a time vector just including the times for which you have an observation for the satellite.
# visibility_idx = (azelrange_rx2satt['Elevation'] > 0) & (azelrange_rx2satt['Elevation'] < 180)
# seconds = seconds[visibility_idx]
# azelrange_rx2satt = azelrange_rx2satt[visibility_idx]
# sattPos = sattPos[visibsility_idx]
# azelrange_rx2satt.to_excel('revised'+MY_SATT+'.xlsx') # Export for sanity

# Compute the expected range for each of these times.
satExpectedRange = np.zeros(len(seconds))
bsv = np.zeros(len(seconds))
relsv = np.zeros(len(seconds))
for i in range(len(seconds)):
	satExpectedRange[i], bsv[i], relsv[i] = ExpectedRange(sattPos[i,:], rxPos,seconds[i].astype(float), almanac_data, PRN)

# Compute and plot the difference between the measured C/A code pseudorange and the
# expected range (R in meters) versus time (hours) for Sept 1, 2020. What should it look like?
C1 = obs_data[MY_SATT].signals["L1"].pr
P2 = obs_data[MY_SATT].signals["L2"].pr

dPR0 = C1 - satExpectedRange
# dPR0 = C1[visibility_idx] - satExpectedRange
fig, ax = plt.subplots()
ax.set_title("DPR (m) PRN%s" % PRN)
plt.scatter((seconds.astype(float) -172800)/3600 , dPR0, 3, label = "dPR0")
plt.legend(markerscale=2.)
ax.set_xlabel("Time (hours after Midnight)")
ax.set_ylabel("DPR (m)")
ax.grid(True)
fig.savefig('Pics\PRN%s\PRN%sPLOT1.png' % (PRN, PRN))

##############################################################################################################
# Part 2: Clock error: Write a function to compute the satellite clock correction (meters) based on the a0 and a1 terms in
# the broadcast epheremis. Alternatively, you can add the SV clock correction capability to broadcast_eph2pos and provide bsv as an
# additional output

# Plot the SV clock correction (meters) versus time (hours) for Sept 1, 2020. What should it look like?
fig, ax = plt.subplots()
ax.set_title("Clock Correction (m) PRN%s" % PRN)
plt.scatter((seconds.astype(float) -172800)/3600 , bsv, 3)
plt.legend(markerscale=2.)
ax.set_xlabel("Time (hours after Midnight)")
ax.set_ylabel("Clock Correction (m)")
ax.grid(True)
fig.savefig('Pics\PRN%s\PRN%sPLOT2.png' % (PRN, PRN))

dPR1 = C1 - (satExpectedRange - bsv)
fig, ax = plt.subplots()
ax.set_title("DPR (m) PRN%s" % PRN)
plt.scatter((seconds.astype(float) -172800)/3600 , dPR1, 3, label = "dPR1-CLOCK")
plt.legend(markerscale=2.)
ax.set_xlabel("Time (hours after Midnight)")
ax.set_ylabel("Difference in Psuedo Range (m)")
ax.grid(True)
fig.savefig('Pics\PRN%s\PRN%sPLOT3.png' % (PRN, PRN))

##############################################################################################################
# Part 3: Relativity: Write a function or augment broadcast_eph2pos to compute the relativistic correction
# (meters) based on the orbital elements in the ephemeris data
dPR1 = C1 - (satExpectedRange - bsv)
dPR2 = C1 - (satExpectedRange - bsv - relsv)

fig, ax = plt.subplots()
plt.scatter((seconds.astype(float) -172800)/3600 , relsv, 3)
# ax.set_xlim(6, 12)
# ax.set_ylim(min(relsv), max(relsv))
ax.set_title("Relativistic Correction (m) PRN%s" % PRN)
ax.set_xlabel("Time (hours after Midnight)")
ax.set_ylabel("Relativistic Correction (m)")
ax.grid(True)
fig.savefig('Pics\PRN%s\PRN%sPLOT4.png' % (PRN, PRN))

fig, ax = plt.subplots()
ax.set_title("DPR (m) PRN%s" % PRN)
plt.scatter((seconds.astype(float) -172800)/3600 , dPR1, 3, label = "dPR1-CLOCK")
plt.scatter((seconds.astype(float) -172800)/3600 , dPR2, 3, label = "dPR2-RELATIVITY")
plt.legend(markerscale=2.)
ax.set_xlabel("Time (hours after Midnight)")
ax.set_ylabel("Difference in Psuedo Range (m)")
ax.grid(True)
fig.savefig('Pics\PRN%s\PRN%sPLOT5.png' % (PRN, PRN))

##############################################################################################################
# Part 4: . Troposphere: Write a function to compute the simple tropospheric delay model (meters) based on the satellite
# elevation and an assumed zenith delay for the ground station location. For the NIST ground station you may assume the 
# zenith value is 2 m.
zd = 2
tropo = zd/np.sin(np.deg2rad(azelrange_rx2satt['Elevation'])) 

fig, ax = plt.subplots()
plt.scatter((seconds.astype(float) -172800)/3600 , tropo, 3)
# ax.set_xlim(6, 12)
# ax.set_ylim(min(tropo), max(tropo))
ax.set_title("Tropo Correction (m) PRN%s" % PRN)
ax.set_xlabel("Time (hours after Midnight)")
ax.set_ylabel("Tropo Correction (m)")
ax.grid(True)
fig.savefig('Pics\PRN%s\PRN%sPLOT6.png' % (PRN, PRN))

# tropo = tropomodel(np.array([weekComp,seconds.astype(float)]).transpose(), zd=2)
# Plot the tropo correction (meters) versus time (hours) for Sept 1, 2020. What should it look like?
dPR1 = C1 - (satExpectedRange - bsv)
dPR2 = C1 - (satExpectedRange - bsv - relsv)
dPR3 = C1 - (satExpectedRange - bsv - relsv + tropo)

# Compute and plot dPR3 = C1 – (R – bsv - relsv + tropo)
fig, ax = plt.subplots()
ax.set_title("DPR (m) PRN%s" % PRN)
plt.scatter((seconds.astype(float) -172800)/3600 , dPR1, 3, label = "dPR1-CLOCK")
plt.scatter((seconds.astype(float) -172800)/3600 , dPR2, 3, label = "dPR2-RELATIVITY")
plt.scatter((seconds.astype(float) -172800)/3600 , dPR3, 3, label = "dPR3-TROPO")
plt.legend(markerscale=2.)
ax.set_xlabel("Time (hours after Midnight)")
ax.set_ylabel("Difference in Psuedo Range (m)")
ax.grid(True)
fig.savefig('Pics\PRN%s\PRN%sPLOT7.png' % (PRN, PRN))

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
ax.set_title("Iono Correction (m) PRN%s" % PRN)
ax.set_xlabel("Time (hours after Midnight)")
ax.set_ylabel("Iono Correction (m)")
ax.grid(True)
fig.savefig('Pics\PRN%s\PRN%sPLOT8.png' % (PRN, PRN))

dPR1 = C1 - (satExpectedRange - bsv)
dPR2 = C1 - (satExpectedRange - bsv - relsv)
dPR3 = C1 - (satExpectedRange - bsv - relsv + tropo)
dPR4 = PRIF - (satExpectedRange - bsv - relsv + tropo)

# Compute and plot dPR3 = C1 – (R – bsv - relsv + tropo)
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
fig.savefig('Pics\PRN%s\PRN%sPLOT9.png' % (PRN, PRN))

plt.show()