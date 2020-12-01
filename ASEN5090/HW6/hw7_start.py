from types import SimpleNamespace
import numpy
from numpy import array, nan
from datetime import datetime, timezone
import re
from hw6_helpers import get_nistECEF, get_usn8ECEF, broadcast2posM, azelrange, ECEF_ENU_Transform, ECEFtolla, llatoECEF
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

def ExpectedRange(base_station, Tr, a, satNum):
	c = 2.99792458e8
	wE  = 7.2921151467e-5
	h, sat, clock_corr, rel = broadcast2posM(a,np.array([[2121,2121],[Tr, Tr ]]).transpose(),satNum)
	Geo = np.linalg.norm(sat[0,:] - base_station)
	TimeTrans = Tr - Geo / c
	h, posHere, clock_corr, rel = broadcast2posM(a,np.array([[2121,2121],[TimeTrans, TimeTrans ]]).transpose(),satNum)
	phi = wE * (Tr - TimeTrans)
	rot = np.array([[np.cos(phi), np.sin(phi), 0],[-np.sin(phi), np.cos(phi),0],[0,0,1]])
	changed = np.dot(rot,posHere[0,:])
	rangeNew = np.linalg.norm(changed - base_station)

	while abs(rangeNew - Geo)> 1e-6:

		Geo = rangeNew
		TimeTrans = Tr - Geo / c
		h, posHere,clock_corr, rel = broadcast2posM(a,np.array([[2121,2121],[TimeTrans, TimeTrans ]]).transpose(),satNum)
		phi = wE * (Tr - TimeTrans)
		rot = np.array([[np.cos(phi), np.sin(phi), 0],[-np.sin(phi), np.cos(phi),0],[0,0,1]])
		changed = rot @ posHere[0,:]
		rangeNew = np.linalg.norm(changed - base_station)
	aer =azelrange(base_station, changed.reshape(1,3))
	#print(aer)
	az = aer[0,0]
	el = aer[0,1]
	return rangeNew , clock_corr, rel, az, el, changed

def get_prefit_residuals(init_guess, sat_ids, dt, current_time, C):
    prefit_residuals = []
    dys = []
    el_vals = []
    g_matrix= np.empty((0,4), float)

    for my_satt in sat_ids:
        # print("PRN: {}".format(int(my_satt[1:])))
        dt_idx = np.where(obs_data[my_satt].time == dt)[0][0]

        # print(obs_data[my_satt].signals["L1"].pr)
        C1 = obs_data[my_satt].signals["L1"].pr[dt_idx:dt_idx+2]
        P2 = obs_data[my_satt].signals["L2"].pr[dt_idx:dt_idx+2]
        L1 = 1575.42e6 
        L2 = 1227.60e6

        exp_r, clock_corr, rel_corr, azimuth, elevation, ecef_vals = ExpectedRange(init_guess, current_time, almanac_data, int(my_satt[1:])) #grabbing_matrixthe integer PRN from the Satellite GXX string
        # tropo =  TopoCorrect(elevation, 2)
        tropo =  TopoCorrect(elevation, TROPO_VAL)
        PRIF, ionocorr = Ionofree(C1, L1, P2, L2)

        dy = PRIF[0] - (exp_r - clock_corr -rel_corr + tropo )

        if elevation > 10:
            el_vals.append(elevation)

            prefit_residuals.append(dy)
            dys.append(dy)
        
            row = [-(ecef_vals[0]- init_guess[0])/ exp_r, -(ecef_vals[1]- init_guess[1])/ exp_r, - (ecef_vals[2] - init_guess[2])/ exp_r, 1]
            g_matrix= np.append(g_matrix, np.array([row]), axis = 0)

    #Performs the least-squares estimation, but without crappy numpy inverse methods.
    relative_position, _, _, _ = np.linalg.lstsq(g_matrix, dys)
    postfit_residuals = prefit_residuals - np.matmul(g_matrix, relative_position)

    c_squiggly = np.zeros((4,4)); c_squiggly[0:3, 0:3] = C; c_squiggly[3, 3] = 1
    errENU = c_squiggly @ np.transpose([relative_position])  
    g_squiggly = np.matmul(g_matrix, c_squiggly.transpose())
    try:
        h_matrix = np.linalg.inv(np.matmul(g_squiggly.transpose(), g_squiggly))
        dops = [h_matrix[0,0], h_matrix[1,1], h_matrix[2,2], h_matrix[3,3]]
    except:
        dops = [-100,-100,-100,-100]
        print('oh no')
        print(dt)

    if dt ==  np.datetime64("2020-09-01T01:00:00"):
        # print('hour mark')
        print(base_station.upper())
        print(dt)
        print("Error(E): {:.2f} m".format(errENU[0,0]))
        print("Error(N): {:.2f} m".format(errENU[1,0]))
        print("Error(U): {:.2f} m".format(errENU[2,0]))
        print()
        print("EDOP: {:.2f} ".format(h_matrix[0,0]))
        print("NDOP: {:.2f} ".format(h_matrix[1,1]))
        print("VDOP: {:.2f} ".format(h_matrix[2,2]))
        print("TDOP: {:.2f} ".format(h_matrix[3,3]))

    # print(sat_ids)
    # print(prefit_residuals)
    # print(postfit_residuals)
    # print(el_vals)
    # print(errENU)
    # print(dops)
    return prefit_residuals, postfit_residuals, el_vals, errENU, dops

TROPO_VAL = 2 # 2.4
base_station = "nist" # nist
# base_station = "usn8" # nist
directory = "%s_pickles\\" % base_station # NIST_pickles
filepath = "%s2450.20o" % base_station
# filepath = 'usn82450.20o'

prefit_residuals = []
postfit_residuals = []
elevation_values = []
relative_positionENU = []
dop_values = []

c = 299792458 
current_time = 172800
century = 2000
header, almanac_data = parse_rinex_nav_file('brdc2450.20n') # Grab data 
# _, obs_data = parse_rinex_obs_file("nist2450.20o") # Grab even more data 
_, obs_data = parse_rinex_obs_file(filepath) # Grab even more data 

init_guess_NIST = get_nistECEF() # Grabs ECEF coords of NIST Location
init_guess_USN8 = get_usn8ECEF() # Grabs ECEF coords of NIST Location

lati, longi, _ = ECEFtolla(init_guess_NIST)
nist_conversion_matrix = ECEF_ENU_Transform(lati,longi)
lati, longi,alt = ECEFtolla(init_guess_USN8)
usn8_conversion_matrix = ECEF_ENU_Transform(lati,longi)

# print(nist_conversion_matrix)
# print(usn8_conversion_matrix)

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
        if base_station == "usn8":
            prft_res, posft_res, elevs, relposENU, dops = get_prefit_residuals(init_guess_USN8, sat_ids, dt, current_time, usn8_conversion_matrix)
        else:
            prft_res, posft_res, elevs, relposENU, dops = get_prefit_residuals(init_guess_NIST, sat_ids, dt, current_time, nist_conversion_matrix)
        # print(prft_res)
        
        prefit_residuals.append(prft_res)
        postfit_residuals.append(posft_res)
        elevation_values.append(elevs)
        relative_positionENU.append(relposENU)
        dop_values.append(dops)
        current_time = current_time + 30

except StopIteration:
    pass

prft_filename = directory + "stored_prefit_residuals_%s.pickle" % base_station
posft_filename = directory + "stored_postfit_residuals_%s.pickle" % base_station
el_filename = directory + "stored_elevation.pickle_%s" % base_station
relposENU_filename = directory + "stored_relposENU_%s.pickle" % base_station
dops_filename = directory + "stored_dops.pickle_%s" % base_station

print(np.shape(dop_values))
with open(prft_filename, 'wb') as f:
    pickle.dump(prefit_residuals, f)

with open(posft_filename, 'wb') as f:
    pickle.dump(postfit_residuals, f)

with open(el_filename, 'wb') as f:
    pickle.dump(elevation_values, f)
    
with open(relposENU_filename, 'wb') as f:
    pickle.dump(relative_positionENU, f)

with open(dops_filename, 'wb') as f:
    pickle.dump(dop_values, f)