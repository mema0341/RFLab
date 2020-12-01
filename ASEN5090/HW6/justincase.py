from types import SimpleNamespace
import numpy
from numpy import array, nan
from datetime import datetime, timezone
import re
from hw6_helpers import get_nistECEF, broadcast2posM, azelrange, ECEF_ENU_Transform, ECEFtolla, llatoECEF
from read_ephemeris import parse_rinex
from rinex_utilities import parse_rinex_nav_file, parse_rinex_obs_file, parse_rinex_header
import numpy as np
from matplotlib import pyplot as plt
# from readyuma import mean2eclock_correntric, broadcast2pos, read_GPSyuma
import pandas as pd
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)
import pickle 

def Ionofree (pr1, f1, pr2, f2):
	prif = np.zeros(pr1.size)
	correction = np.zeros(pr1.size)
	for i in range(pr1.size):
		prif[i] = (pr1[i]*f1**2 - pr2[i]*f2**2)/(f1**2-f2**2)
		correction[i] = (f2**2/(f1**2 - f2**2)) *(pr2[i]-pr1[i])
	return prif, correction


def TopoCorrect (els, zd):
	Del = np.zeros(els.size)
	i=0
	if els.size>1:
		for el in els:
			el = el* np.pi/180
			md = 1/(np.sin(el))
			Del[i] = zd *md
			i+=1
	else:
		el = els* np.pi/180
		md = 1/(np.sin(el))
		Del= zd *md
	return Del

def ExpectedRange(station, Tr, a, satNum):
	c = 2.99792458e8
	wE  = 7.2921151467e-5
	h, sat, clock_corr, rel = broadcast2posM(a,np.array([[2121,2121],[Tr, Tr ]]).transpose(),satNum)
	Geo = np.linalg.norm(sat[0,:] - station)
	TimeTrans = Tr - Geo / c
	h, posHere, clock_corr, rel = broadcast2posM(a,np.array([[2121,2121],[TimeTrans, TimeTrans ]]).transpose(),satNum)
	phi = wE * (Tr - TimeTrans)
	rot = np.array([[np.cos(phi), np.sin(phi), 0],[-np.sin(phi), np.cos(phi),0],[0,0,1]])
	changed = np.dot(rot,posHere[0,:])
	rangeNew = np.linalg.norm(changed - station)

	while abs(rangeNew - Geo)> 1e-6:

		Geo = rangeNew
		TimeTrans = Tr - Geo / c
		h, posHere,clock_corr, rel = broadcast2posM(a,np.array([[2121,2121],[TimeTrans, TimeTrans ]]).transpose(),satNum)
		phi = wE * (Tr - TimeTrans)
		rot = np.array([[np.cos(phi), np.sin(phi), 0],[-np.sin(phi), np.cos(phi),0],[0,0,1]])
		changed = rot @ posHere[0,:]
		rangeNew = np.linalg.norm(changed - station)
	aer =azelrange(station, changed.reshape(1,3))
	#print(aer)
	az = aer[0,0]
	el = aer[0,1]
	return rangeNew , clock_corr, rel, az, el, changed

def get_prefit_residuals(sat_ids, dt, current_time):
    prefit_residuals = []

    plot = False 

    g_matrix= np.empty((0,4), float)
    pifs = []
    # b. the expected range (from the last homework) in meters
    exp_rs = []
    # c. the elevation and azimuth in degrees (based on the adjusted satellite location at time of transmission)
    az_vals = []
    el_vals = []
    # d. the satellite clock correction (from the a0 and a1 terms in the nav message) in meters
    clock_corrs = []
    # e. the relativistic correction (from the orbital elements or r,v of the GPS satellite) in meters
    rel_corrs = []
    # f. the simple tropospheric model (used in HW5) in meters
    tropo_vals = []
    ecef_vals = []

    dys = []

    PRNs = []

    for my_satt in sat_ids:
        # print("PRN: {}".format(int(my_satt[1:])))
        dt_idx = np.where(obs_data[my_satt].time == dt)[0][0]

        C1 = obs_data[my_satt].signals["L1"].pr[dt_idx:dt_idx+2]
        P2 = obs_data[my_satt].signals["L2"].pr[dt_idx:dt_idx+2]
        L1 = 1575.42e6 
        L2 = 1227.60e6

        exp_r, clock_corr, rel_corr, azimuth, elevation, ecef_vals = ExpectedRange(init_guess, current_time, almanac_data, int(my_satt[1:])) #grabbing_matrixthe integer PRN from the Satellite GXX string
        tropo =  TopoCorrect(elevation, 2)
        PRIF, ionocorr = Ionofree(C1, L1, P2, L2)
        # print("Expected Range: {} m".format(exp_r))
        # print("Iono-Corrected : {} m".format(PRIF[0]))
        # print("Clock Correction: {} m".format(clock_corr))
        # print("rel_corr Correction: {} m".format(rel_corr))
        # print("Azimuth: {} deg,  Elevation: {} deg".format(azimuth, elevation))
        # print("Tropo Correction: {} m".format(tropo))
        
        # print("#######################################################")

        row = [-(ecef_vals[0]- init_guess[0])/ exp_r, -(ecef_vals[1]- init_guess[1])/ exp_r, - (ecef_vals[2] - init_guess[2])/ exp_r, 1]
        g_matrix= np.append(g_matrix, np.array([row]), axis = 0)
        exp_rs = np.append(exp_rs, exp_r)
        pifs = np.append(pifs, PRIF[0])
        clock_corrs = np.append(clock_corrs, clock_corr)
        rel_corrs = np.append(rel_corrs, rel_corr)
        tropo_vals = np.append(tropo_vals, tropo)
        az_vals = np.append(az_vals, azimuth)
        el_vals = np.append(el_vals, elevation)
        PRNs = np.append(PRNs, int(my_satt[1:]))

        dy = PRIF[0] - (exp_r - clock_corr -rel_corr + tropo )
        dys = np.append(dys, dy)

        if elevation > 10:
            prefit_residuals.append(dy)


    if plot:
        plt.scatter(el_vals, dys)
        plt.grid()
        plt.title('Psuedorange Residuals vs Elevation')
        plt.xlabel('Elevation (deg)')
        plt.ylabel('Psuedorange Residual (m)')
        plt.show()

    return prefit_residuals

    # c = 2.99792458e8
    # base = np.datetime64("2020-09-01T00:00:00")
    # Washington = llatoECEF(38.889, -77.049, 0)

    # #Performs the least-squares estimation, but without crappy numpy inverse methods.
    # delt, _, _, _ = np.linalg.lstsq(g_matrix, dys)

    # #this is the correction I would do if problem 6 asked for it.
    # new_r = init_guess + delt[0:3]
    # new_time = current_time + delt[3]/c

    # new_r = np.array(list(Washington))
    # #new_r = init_guess  
    # print("Intitial Guess: {}".format(new_r))
    # for i in range(5):
    #     g_matrix= np.empty((0,4), float)
    #     dys = []
    #     #go through all the satellites.
    #     for my_satt in sat_ids:
    #         #get data and do calculations if it is visible at time = base
    #         C1 = obs_data[my_satt].signals["L1"].pr[dt_idx:dt_idx+2]
    #         P2 = obs_data[my_satt].signals["L2"].pr[dt_idx:dt_idx+2]
    #         L1 = 1575.42e6 
    #         L2 = 1227.60e6

    #         exp_r, clock_corr, rel_corr, azimuth, elevation, ecef_vals = ExpectedRange(new_r, new_time,almanac_data, int(my_satt[1:])) #grabbing_matrixthe integer PRN from the Satellite GXX string
    #         tropo =  TopoCorrect(elevation, 2)
    #         PRIF, ionocorr = Ionofree(C1, L1, P2, L2)

    #         #a new row for the g_matrixmatrix
    #         row = [-(ecef_vals[0]- init_guess[0])/ exp_r, -(ecef_vals[1]- init_guess[1])/ exp_r, - (ecef_vals[2] - init_guess[2])/ exp_r, 1]
    #         g_matrix= np.append(g_matrix, np.array([row]), axis = 0)
    #         #The dysual, we'll keep trying_matrixto change ER and time 
    #         dy = PRIF[0] - (exp_r - clock_corr -rel_corr + tropo )
    #         dys = np.append(dys, dy)
    #     #do least squares and continue
    #     delt, _, _, _ = np.linalg.lstsq(g_matrix, dys, rcond = None)
    #     new_r = new_r + delt[0:3]
    #     new_time = new_time + delt[3]/c
    #     print("iteration: {}".format(i+1))
    #     print()
    #     print("delta : {}".format(delt))
    #     print("Position Est.: {}".format(new_r))
    #     print("Time Est.: {}".format(new_time))
    #     print(" ")

    # print("Reference Point (NIST): {}".format(init_guess))
    # lati,longi,alt = ECEFtolla(init_guess)
    # err = new_r - init_guess
    # print()
    # C = ECEF_ENU_Transform(lati,longi)
    # errENU = C @ np.transpose([err])
    # print("ECEF-ENU Matrix")
    # print(C)
    # print()

    # print(dt)
    # print("Error(E): {:.2f} m".format(errENU[0,0]))
    # print("Error(N): {:.2f} m".format(errENU[1,0]))
    # print("Error(U): {:.2f} m".format(errENU[2,0]))
    # print()
    
    # c_squiggly = np.zeros((4,4)); c_squiggly[0:3, 0:3] = C; c_squiggly[3, 3] = 1
    # g_squiggly = np.matmul(g_matrix, c_squiggly)
    # h_matrix = np.linalg.inv(np.matmul(g_squiggly.transpose(), g_squiggly))
    # print("EDOP: {:.2f} ".format(h_matrix[0,0]))
    # print("NDOP: {:.2f} ".format(h_matrix[1,1]))
    # print("VDOP: {:.2f} ".format(h_matrix[2,2]))
    # print("TDOP: {:.2f} ".format(h_matrix[3,3]))

prefit_residuals = []

c = 299792458 
current_time = 172800
century = 2000

_, obs_data = parse_rinex_obs_file("nist2450.20o") # Grab even more data 
header, almanac_data = parse_rinex_nav_file('brdc2450.20n') # Grab data 

init_guess = get_nistECEF() # Grabs ECEF coords of NIST Location

filepath = 'nist2450.20o'
with open(filepath, 'r') as f:
    lines = list(f.readlines())
for i, line in enumerate(lines):
    if line.find('END OF HEADER') >= 0:
        break
header_lines = lines[:i + 1]
obs_lines = lines[i + 1:]
header = parse_rinex_header(header_lines)
if not hasattr(header, 'obs_types'):
    raise Exception('RINEX header must contain `# / TYPES OF OBS.` and `header` namespace from `parse_rinex_header` must contain corresponding list `obs_types`')
lines = obs_lines
observations = header.obs_types

data = {}  # <sat_id>: {'time': [<dt...>], <obs_id>: [<values...>]}
lines = iter(lines)
try:
    while True:
        # at each epoch, the two-digit year, month, day, hour, minute, and seconds
        # of the measurement epoch are specified, along with the number and ids of
        # the satellites whose measurements are given
        line = next(lines)
        yy = int(line[:4])
        year = century + yy
        month = int(line[4:7])
        day = int(line[7:10])
        hour = int(line[10:13])
        minute = int(line[13:16])
        seconds = float(line[16:25])
        microseconds = int(1e6 * (seconds % 1))
        seconds = int(seconds)
        dt = numpy.datetime64(datetime(year, month, day, hour, minute, seconds, microseconds))
        flag = int(line[25:28])
        num_sats = int(line[29:32])

        # there is space for (80 - 32) / 3 = 16 satellite ids
        # if there are more than 16, then they continue on the next line
        line = line[32:]
        if num_sats > 16:
            line = (line + next(lines).strip()).replace(' ', '')
        line = line.strip()
        # must replace spaces with zeros: e.g. to convert `'G 1'` to `'G01'`
        sat_ids = [line[3*i:3*(i+1)].replace(' ', '0') for i in range(num_sats)]
        # print(dt)
        # print(sat_ids)

        c1_vals = []
        p2_vals = []
        for sat_id in sat_ids:
            # create new entry if `sat_id` is new
            if sat_id not in data.keys():
                data[sat_id] = {'time': []}
                for obs_id in observations:
                    data[sat_id][obs_id] = []
            # append time first, then append obs values
            data[sat_id]['time'].append(dt)
            # each line of observation values contains up to 5 entries
            # each entry is of width 16, starting at index 0
            num_lines_per_sat = 1 + len(observations) // 5
            line = ''
            while num_lines_per_sat > 0:
                line += next(lines).replace('\n', '')
                num_lines_per_sat -= 1
            for i, obs_id in enumerate(observations):
                try:
                    val = float(line[16 * i:16 * (i + 1)])
                except Exception:
                    val = nan
                data[sat_id][obs_id].append(val)
                # print(obs_id)
                # print(data)
                # print('')
        # print("####################################################################################")
        # print()

        prft_res = get_prefit_residuals(sat_ids, dt, current_time)
        prefit_residuals.append(prft_res)

        current_time = current_time + 30

except StopIteration:
    pass

print(np.shape(prefit_residuals))
filename = "stored_prefit_residuals.pickle"

with open(filename, 'wb') as f:
    pickle.dump(prefit_residuals, f)


