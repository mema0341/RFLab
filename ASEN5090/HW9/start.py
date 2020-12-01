from hw8_helpers import *

from read_ephemeris import parse_rinex
from rinex_utilities import parse_rinex_nav_file, parse_rinex_obs_file
import numpy as np
from matplotlib import pyplot as plt
from readyuma import mean2eccentric, broadcast2pos, read_GPSyuma
import pandas as pd
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)
                              

fig2 ,(a1,a2,a3) = plt.subplots(3,1)  

head, a = parse_rinex_nav_file('brdc2450.20n')
#it would be nice to understand what a SimpleNamespace is before doing this.
#print(head)

timevec = np.arange(172800, 172800+86401, 31)
weeks = np.ones(len(timevec)) *2121

with open("NIST_location.txt", 'r') as file:
	file.readline()
	pos = file.readline()
strings = pos.split()
pos_1 = np.array([float(part) for part in strings])
vec1_ECEF = pos_1.transpose()

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
	aer1 = pd.DataFrame(azelrange(vec1_ECEF, poses), columns = ['Azimuth', 'Elevation', 'Range'])
	
	ax1.plot(np.deg2rad(aer1['Azimuth']), aer1['Elevation'], label = "PRN {}".format(i))
	
	#aer1 = pd.DataFrame(azelrange(userECEFs[2], poses), columns = ['Azimuth', 'Elevation', 'Range'])

ax1.legend()
plt.show()

# h, PRN1pose = broadcast2posM(a,np.array([weeks,timevec]).transpose(),5)
# aer1 = pd.DataFrame(azelrange(vec1_ECEF, PRN1pose), columns = ['Azimuth', 'Elevation', 'Range'])
# a1.plot((timevec/3600 - 48)[aer1["Elevation"]>0],aer1['Azimuth'][aer1["Elevation"]>0], 'bo', label = "Azimuth", markersize = 3 )
# a2.plot((timevec/3600 - 48)[aer1["Elevation"]>0],aer1['Elevation'][aer1["Elevation"]>0],'bo', label = "Elevation", markersize = 3)
# a3.plot((timevec/3600 - 48)[aer1["Elevation"]>0],aer1['Range'][aer1["Elevation"]>0],'bo', label= "Range", markersize = 3)
# a3.set_xlabel("Time (hours after midnight)")
# fig2.suptitle("Az/El/Range from NIST to PRN 5 based on Broadcast Ephemeris")
# #a1.set_title("Azimuth")
# a1.set_ylabel("Azimuth (deg)")

# #a2.set_title("Elevation")
# a2.set_ylabel("Elevation (deg)")

# a3.set_ylabel("Range (m)")

# #h, PRN1pose = broadcast2posM(a,np.array([weeks,timevec]).transpose(),5)
# # print(PRN1pose)
# # plt.plot(timevec, PRN1pose[:,0],'g-')
# # plt.plot(timevec, PRN1pose[:,1], 'b-')
# # plt.plot(timevec, PRN1pose[:,2], 'r-')

# [diction, dfyuma] =read_GPSyuma('YUMA245.ALM')

# h, posAlm = broadcast2pos(dfyuma,np.array([weeks,timevec]).transpose(),5 )
# plt.plot(timevec,posAlm[:,0]-PRN1pose[:,0],'g--')
# plt.plot(timevec,posAlm[:,1]-PRN1pose[:,1], 'b--')
# plt.plot(timevec,posAlm[:,2]-PRN1pose[:,2], 'r--')

# plt.show()
formatter = DateFormatter("%d %h %H:%M ")
fig = plt.figure()
axpr = fig.add_subplot(111)
head, obs_data = parse_rinex_obs_file("nist2450.20o")
print(obs_data["G05"])
base = np.datetime64("2020-09-01T00:00:00")
seconds = (obs_data["G05"].time - base)/1e6 + 172800
weekComp = np.ones(len(seconds))*2121
#print(seconds.astype(int))

h, PRN1pose = broadcast2posM(a,np.array([weekComp,seconds.astype(float)]).transpose(),5)
aer2 = pd.DataFrame(azelrange(vec1_ECEF, PRN1pose), columns = ['Azimuth', 'Elevation', 'Range'])
Ex1 = np.zeros(len(seconds))
for i in range(len(seconds)):
	Ex1[i] = ExpectedRange(PRN1pose[i,:], vec1_ECEF,seconds[i].astype(float), a)
#print(seconds)

#seconds = 172800 + obs_data["G05"].time.hour *3600 + obs_data["G05"].time.minutes * 60 + obs_data["G05"].time.seconds
#print(seconds)
axpr.set_title("Psuedo Range Measurements for PRN 5")
axpr.scatter((seconds.astype(float) -172800)/3600,obs_data["G05"].signals["L1"].pr, 2, label = "L1 PsuedoRange")
axpr.scatter((seconds.astype(float) -172800)/3600,obs_data["G05"].signals["L2"].pr,2, label = "L2 PsuedoRange")
axpr.set_xlabel("Time after Midnight (hours)")
axpr.set_ylabel("Psuedo Range (m)")
# axpr.scatter((seconds.astype(float) -172800)/3600, aer2["Range"], label = "Expected Geometric Range at reception")
# axpr.scatter((seconds.astype(float) -172800)/3600, Ex1, label = "Expected Range, Corrected")
#axpr.xaxis.set_major_formatter(formatter)
plt.legend()
plt.show()
plt.figure()
plt.title("Difference Between Corrected and Uncorrected Psuedoranges")
plt.scatter((seconds.astype(float) -172800)/3600 , obs_data["G05"].signals["L1"].pr - aer2["Range"], 3, label = "Difference Between Corrected Expectation and L1 PsuedoRange")
plt.scatter((seconds.astype(float) -172800)/3600 , obs_data["G05"].signals["L2"].pr - aer2["Range"], 3, label = "Difference Between Corrected Expectation and L2 PsuedoRange")
plt.legend()
plt.xlabel("Time (hours after Midnight)")
plt.ylabel("Difference in Psuedo Range (m)")
#plt.scatter((seconds.astype(float) -172800)/3600 , Ex1 - aer2["Range"], 3)
#plt.scatter((seconds.astype(float) -172800)/3600 , Ex1 - aer2["Range"])

plt.show()
# #there are 13 of these objects in each thing.
# # a note: a[PRN][entry indx].attribute
# h, PRN1pose = broadcast2posM(a,np.array([[2121,2121],[172800,172800]]).transpose(),5)
# rangeN = ExpectedRange(PRN1pose[0,:], vec1_ECEF, 172800, a)
# print(rangeN)fig2 ,(a1,a2,a3) = plt.subplots(3,1)  

head, a = parse_rinex_nav_file('brdc2450.20n')
#it would be nice to understand what a SimpleNamespace is before doing this.
#print(head)

timevec = np.arange(172800, 172800+86401, 31)
weeks = np.ones(len(timevec)) *2121

with open("NIST_location.txt", 'r') as file:
	file.readline()
	pos = file.readline()
strings = pos.split()
pos_1 = np.array([float(part) for part in strings])
vec1_ECEF = pos_1.transpose()
# 
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111, polar = True)
# ax1.set_rlim(bottom = 90, top = 0)
# ax1.set_theta_zero_location("N")
# ax1.set_theta_direction(-1)
# ax1.set_title("GPS Passes at NIST on Sept 1 2020")
# # plotting through the day using the broadcast2pos function.
# for i in a:
# 	#print(i)
# 	# print(np.array([[somezeros],[times]]))
# 	[h, poses]=broadcast2posM(a,np.array([weeks,timevec]).transpose(),i)
# 	aer1 = pd.DataFrame(azelrange(vec1_ECEF, poses), columns = ['Azimuth', 'Elevation', 'Range'])
	
# 	ax1.plot(np.deg2rad(aer1['Azimuth']), aer1['Elevation'], label = "PRN {}".format(i))
	
# 	#aer1 = pd.DataFrame(azelrange(userECEFs[2], poses), columns = ['Azimuth', 'Elevation', 'Range'])

# ax1.legend()
# plt.show()

# h, PRN1pose = broadcast2posM(a,np.array([weeks,timevec]).transpose(),5)
# aer1 = pd.DataFrame(azelrange(vec1_ECEF, PRN1pose), columns = ['Azimuth', 'Elevation', 'Range'])
# a1.plot((timevec/3600 - 48)[aer1["Elevation"]>0],aer1['Azimuth'][aer1["Elevation"]>0], 'bo', label = "Azimuth", markersize = 3 )
# a2.plot((timevec/3600 - 48)[aer1["Elevation"]>0],aer1['Elevation'][aer1["Elevation"]>0],'bo', label = "Elevation", markersize = 3)
# a3.plot((timevec/3600 - 48)[aer1["Elevation"]>0],aer1['Range'][aer1["Elevation"]>0],'bo', label= "Range", markersize = 3)
# a3.set_xlabel("Time (hours after midnight)")
# fig2.suptitle("Az/El/Range from NIST to PRN 5 based on Broadcast Ephemeris")
# #a1.set_title("Azimuth")
# a1.set_ylabel("Azimuth (deg)")

# #a2.set_title("Elevation")
# a2.set_ylabel("Elevation (deg)")

# a3.set_ylabel("Range (m)")

# #h, PRN1pose = broadcast2posM(a,np.array([weeks,timevec]).transpose(),5)
# # print(PRN1pose)
# # plt.plot(timevec, PRN1pose[:,0],'g-')
# # plt.plot(timevec, PRN1pose[:,1], 'b-')
# # plt.plot(timevec, PRN1pose[:,2], 'r-')

# [diction, dfyuma] =read_GPSyuma('YUMA245.ALM')

# h, posAlm = broadcast2pos(dfyuma,np.array([weeks,timevec]).transpose(),5 )
# plt.plot(timevec,posAlm[:,0]-PRN1pose[:,0],'g--')
# plt.plot(timevec,posAlm[:,1]-PRN1pose[:,1], 'b--')
# plt.plot(timevec,posAlm[:,2]-PRN1pose[:,2], 'r--')

# plt.show()
formatter = DateFormatter("%d %h %H:%M ")
fig = plt.figure()
axpr = fig.add_subplot(111)
head, obs_data = parse_rinex_obs_file("nist2450.20o")
print(obs_data["G05"])
base = np.datetime64("2020-09-01T00:00:00")
seconds = (obs_data["G05"].time - base)/1e6 + 172800
weekComp = np.ones(len(seconds))*2121
#print(seconds.astype(int))

h, PRN1pose = broadcast2posM(a,np.array([weekComp,seconds.astype(float)]).transpose(),5)
aer2 = pd.DataFrame(azelrange(vec1_ECEF, PRN1pose), columns = ['Azimuth', 'Elevation', 'Range'])
Ex1 = np.zeros(len(seconds))
for i in range(len(seconds)):
	Ex1[i] = ExpectedRange(PRN1pose[i,:], vec1_ECEF,seconds[i].astype(float), a)
#print(seconds)

#seconds = 172800 + obs_data["G05"].time.hour *3600 + obs_data["G05"].time.minutes * 60 + obs_data["G05"].time.seconds
#print(seconds)
axpr.set_title("Psuedo Range Measurements for PRN 5")
axpr.scatter((seconds.astype(float) -172800)/3600,obs_data["G05"].signals["L1"].pr, 2, label = "L1 PsuedoRange")
axpr.scatter((seconds.astype(float) -172800)/3600,obs_data["G05"].signals["L2"].pr,2, label = "L2 PsuedoRange")
axpr.set_xlabel("Time after Midnight (hours)")
axpr.set_ylabel("Psuedo Range (m)")
# axpr.scatter((seconds.astype(float) -172800)/3600, aer2["Range"], label = "Expected Geometric Range at reception")
# axpr.scatter((seconds.astype(float) -172800)/3600, Ex1, label = "Expected Range, Corrected")
#axpr.xaxis.set_major_formatter(formatter)
plt.legend()
plt.show()
plt.figure()
plt.title("Difference Between Corrected and Uncorrected Psuedoranges")
plt.scatter((seconds.astype(float) -172800)/3600 , obs_data["G05"].signals["L1"].pr - aer2["Range"], 3, label = "Difference Between Corrected Expectation and L1 PsuedoRange")
plt.scatter((seconds.astype(float) -172800)/3600 , obs_data["G05"].signals["L2"].pr - aer2["Range"], 3, label = "Difference Between Corrected Expectation and L2 PsuedoRange")
plt.legend()
plt.xlabel("Time (hours after Midnight)")
plt.ylabel("Difference in Psuedo Range (m)")
#plt.scatter((seconds.astype(float) -172800)/3600 , Ex1 - aer2["Range"], 3)
#plt.scatter((seconds.astype(float) -172800)/3600 , Ex1 - aer2["Range"])

plt.show()
# #there are 13 of these objects in each thing.
# # a note: a[PRN][entry indx].attribute
# h, PRN1pose = broadcast2posM(a,np.array([[2121,2121],[172800,172800]]).transpose(),5)
# rangeN = ExpectedRange(PRN1pose[0,:], vec1_ECEF, 172800, a)
# print(rangeN)