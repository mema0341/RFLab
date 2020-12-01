
import numpy as np 
import pyproj

ECEF = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
LLA = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

# ECEF3ENU function takes in lat and long as arguments and returns
# the ECEF coordinates using the conversion matrix found on slide 35 of Lectures 6&7
def ecef2enu(lon_, lat_):
    lambda_ = lon_
    phi_ = lat_

    # Define the conversion matrix
    conversion_matrix = \
        [[-np.sin(lambda_), np.cos(lambda_), 0],\
         [-np.sin(phi_)*np.cos(lambda_), -np.sin(phi_)*np.sin(lambda_), np.cos(phi_)],\
         [np.cos(phi_)*np.cos(lambda_), np.cos(phi_)*np.sin(lambda_), np.sin(phi_)]]
    
    # To get the solution, take the dot product of the ECEF coords
    # and the conversion matrix
    return conversion_matrix


def azelrange(user, sats):
	#Inputs: ECEF coordinates for user and any number of sats
	#user: 3 * 1 np matrix, a n * 3 np matrix
	n = sats.shape[0]
	out = np.empty((0,3), float)
	#sats[]
	#print('ah')
	for i in range(n):
		loshere = compute_LOS_ENU(user.transpose(), sats[i].transpose())
		az = np.arctan2(loshere[0], loshere[1]) * (180/np.pi)
		el = np.arcsin(loshere[2]) * (180/np.pi)
		r = np.linalg.norm(sats[i] - user)
		out = np.append(out, np.array([[az,el,r]]),axis = 0)
		#print(loshere) 
	return out
def azelrange(user, sats):
	#Inputs: ECEF coordinates for user and any number of sats
	#user: 3 * 1 np matrix, a n * 3 np matrix
	n = sats.shape[0]
	out = np.empty((0,3), float)
	#sats[]
	#print('ah')
	for i in range(n):
		loshere = compute_LOS_ENU(user.transpose(), sats[i].transpose())
		az = np.arctan2(loshere[0], loshere[1]) * (180/np.pi)
		el = np.arcsin(loshere[2]) * (180/np.pi)
		r = np.linalg.norm(sats[i] - user)
		out = np.append(out, np.array([[az,el,r]]),axis = 0)
		#print(loshere) 
	return out
    
def broadcast2posM(ephem_all, t_input, prn):
	"""health, x = broadcast2pos(ephem_all, t_input, prn)

	 Calculates the position from an ephemeris 
	 matrix (see read_GPSyuma function).  The input ephem_all can 
	 be generated by the read_GPSyuma function.


	Modified by P. Axelrad 9/10/2018 to remove extra functionality
	Author: Ben K. Bradley
	Date: 07/19/2009


	INPUT:               Description                                  Units

	ephem_all    - matrix of gps satellite orbit parameters           (nx25)
  
				  col1: prn, PRN number of satellite
				  col2: M0, mean anomaly at reference time, rad
				  col3: delta_n, mean motion difference from computed value, rad/s
				  col4: ecc, eccentricity of orbit
				  col5: sqrt_a, square root of semi-major axis, m^0.5                  col6: Loa, longitude of ascending node of orbit plane at weekly epoch, rad
				  col7: incl, inclination angle at reference time, rad
				  col8: perigee, argument of perigee, rad
				  col9: ra_rate, rate of change of right ascension, rad/s
				 col10: i_rate, rate of change of inclination angle, rad/s
				 col11: Cuc, amplitude of the cosine harmonic correction term to the argument of latitude
				 col12: Cus, amplitude of the sine harmonic correction term to the argument of latitude
				 col13: Crc, amplitude of the cosine harmonic correction term to the orbit radius
				 col14: Crs, amplitude of the sine harmonic correction term to the orbit radius
				 col15: Cic, amplitude of the cosine harmonic correction term to the angle of inclination
				 col16: Cis, amplitude of the cosine harmonic correction term to the angle of inclination
				 col17: Toe, reference time ephemeris (seconds into GPS week)
				 col18: IODE, issue of data (ephemeris) 
				 col19: GPS_week, GPS Week Number (to go with Toe)
				 col20: Toc, time of clock
				 col21: Af0, satellite clock bias (sec)
				 col22: Af1, satellite clock drift (sec/sec)
				 col23: Af2, satellite clock drift rate (sec/sec/sec)
				 col24: Timing Group Delay (TGD), seconds
				 col25: health, satellite health (0=good and usable)

	t_input      - GPS times to calculate values at                 [WN TOW] (nx2)
	prn          - PRN to compute values for (one satellite only)                       



	OUTPUT:       
	
		health       - health of satellite (0=good)                              (nx1)
		x            - position of satellite (ECEF)                  [x y z]   m (nx3)
									 


	 Coupling:

		 mean2eccentric

	 References:
 
	   [1] Interface Control Document: IS-GPS-200D
			 < http://www.navcen.uscg.gov/gps/geninfo/IS-GPS-200D.pdf >
	
	   [2] Zhang, J., et.all. "GPS Satellite Velocity and Acceleration
			 Determination using the Broadcast Ephemeris". The Journal of
			 Navigation. (2006), 59, 293-305.
				< http://journals.cambridge.org/action/displayAbstract;jsess ...
					ionid=C6B8C16A69DD7C910989C661BAB15E07.tomcat1?fromPage=online&aid=425362 >
	
	   [3] skyplot.cpp by the National Geodetic Survey
			  < http://www.ngs.noaa.gov/gps-toolbox/skyplot/skyplot.cpp >


	 Last Updated:
	
	  2015/01/22  B.K. Bradley - the capability to look for updated ephem
								  entries that occur at odd times within each
								  2hr window has been commented out in this 
								  function and added to read_GPSbroadcast.m
								  instead. This moves the computational
								  overhead to the reading which only occurs
								  once.
	2019/08/20   Chelsea Thangavelu - convert matlab function to python 
				 function"""


	# NOTE: Numbered equations in the code (e.g., Eq. 21) correspond to 
	#  equations in the [2] reference.
	
	#==========================================================================
	# Load GPS Accepted WGS-84 Constants 
	#==========================================================================
	muE = 3.986005e14     # WGS-84 value, m^3/s^2
	wE  = 7.2921151467e-5 # WGS-84 value, rad/s 
	c   = 2.99792458e8    # GPS acceptd speed of light, m/s
	
	#==========================================================================
	# Initialize Output Variables for Speed 
	#==========================================================================
	#print(t_input[:,0])
	sz         = t_input.shape[0]
	x          = np.full((sz, 3), np.nan)
	health     = np.full((sz, 1), np.nan)
	
	
	#==========================================================================
	# Pull Out Correct Ephemerides 
	#==========================================================================
	
	# Pull out ephemerides for PRN in question
	sat_ephem = ephem_all[prn]

	
	# No matching PRN found, returning data will be NaNs
	#otherwise continue
	if sat_ephem:
		
		#==========================================================================
		# Start Main Calculation Loop 
		#==========================================================================
		
		# Compute elapsed times of each ephemeris epoch wrt first entry, seconds
		toes = []
		GPS_weeks = []
		for i in range(len(sat_ephem)):
			time = sat_ephem[i].t_oe
			week = sat_ephem[i].week
			toes.append(time)
			GPS_weeks.append(week)
		toe0 = toes[0]
		GPS_week0 = GPS_weeks[0]
		#print(GPS_weeks)
		#print(toes)
		toes = np.array(toes)
		GPS_weeks = np.array(GPS_weeks)

		# GPS_weeks = sat_ephem["GPS_week"].values
		# GPS_week0 = sat_ephem["GPS_week"].values[0]
		# toes = sat_ephem["toa"].values
		# toe0 = sat_ephem["toa"].values[0]
		dt_ephem = (GPS_weeks - GPS_week0) * 604800 + (toes - toe0)
		#print(dt_ephem)
		
		# Compute elapsed times of each input time wrt first ephemeris entry, seconds
		dt_input = (t_input[:,0] - GPS_week0) * 604800 + (t_input[:, 1] - toe0)
		#print(dt_input)
		

		
		for tt in range(sz): # loop through all input times
		
		
			# Pull out most recent ephemeris values
		#     jj = max( find(dt_input(tt) >= dt_ephem) ); # sat_ephem(:,17) = toe (sec into GPS week) of each entry
														# jj = row of specific sat. ephem. data with epoch closest to input time
			
			# Pull out nearest ephemeris values                                                                                        
			min_t_diff = np.abs( dt_input[tt] - dt_ephem ).min()
			jj = np.where(np.abs(dt_input[tt] - dt_ephem) == min_t_diff) #index of minimum
				
			
			if jj[0].size != 1: # no matching ephemeris time found. continue to next input time 
				continue
			jj = int(jj[0])
			#print(jj)
		
			# Pull out common variables from the ephemeris matrix
			#======================================================================
			#toe = sat_ephem(jj,17);           # time of ephemeris
			dt  = dt_input[tt] - dt_ephem[jj] # seconds difference from epoch

			#print(dt)
			a = sat_ephem[jj].sqrt_a ** 2 # semimajor axis
			ecc = sat_ephem[jj].e #eccentricity
			n0 = np.sqrt(muE/a**3) #nominal mean motion (rad/s)
			n   = n0 + sat_ephem[jj].delta_n  #corrected mean motion
			M   = sat_ephem[jj].m_0 + n * dt #mean anomaly
		
#           Compute perigee, true and eccentric anomaly...
#           ======================================================================
		#omega
#            Load argument of perigee to a local variable and add perigee rate, rad
			perigee = sat_ephem[jj].omega
#            % Compute Eccentric Anomaly, rad
			
			E    = mean2eccentric(M,ecc)
			cosE = np.cos(E)
			sinE = np.sin(E)
#        
#            Compute true anomaly, rad
			nu = np.arctan2(np.sqrt(1. - ecc**2) * sinE,  cosE - ecc) 
#        
#            Compute the argument of latitude, rad 
			u = nu + perigee  # true anomaly + argument of perigee

			#correxiones
			du = sat_ephem[jj].c_us * np.sin(2*u) + sat_ephem[jj].c_uc * np.cos(2 * u)
			dr = sat_ephem[jj].c_rs * np.sin(2*u) + sat_ephem[jj].c_rc * np.cos(2 * u)
			di = sat_ephem[jj].c_is * np.sin(2*u) + sat_ephem[jj].c_ic * np.cos(2 * u)
#            Compute radius and inclination
#            ======================================================================
		
			r   = a * (1. - ecc * cosE) + dr                      # corrected radius  
			inc = sat_ephem[jj].i_0 + di + sat_ephem[jj].i_dot *dt  #  inclination 
		
			cosu = np.cos(u)    
			sinu = np.sin(u)  
		
#            Compute satellite position in orbital plane (Eq. 13)
#            ======================================================================
			xo = r * cosu    # satellite x-position in orbital plane
			yo = r * sinu    # satellite y-position in orbital plane
#        
#            Corrected longitude of ascending node for node rate and Earth rotation
#            ======================================================================
			node = sat_ephem[jj].omega_0 + \
			(sat_ephem[jj].omega_dot - wE) * dt - \
			(wE * sat_ephem[jj].t_oe)
#        
#            Calculate GPS Satellite Position in ECEF (m)
#            ======================================================================
			cosi = np.cos(inc)
			sini = np.sin(inc)
			coso = np.cos(node)
			sino = np.sin(node)
		
		
#            Satellite position in ECEF (m)
			x[tt,0] = xo * coso - yo * cosi * sino  #x-position  
		
			x[tt, 1] = xo * sino + yo * cosi * coso;  #y-position 
#        
			x[tt,2] = yo * sini                 #z-position
			
		
#            Keep track of health of each satellite
#            ======================================================================      
			health[tt,0] = sat_ephem[jj].health # satellite health (0.00 is useable)
		
		
		
		
		
#        END of t_input loop =================================================
#        ==========================================================================    

	return health, x

def ExpectedRange(sat,station, Tr, a):
	c = 2.99792458e8
	wE  = 7.2921151467e-5
	Geo = np.linalg.norm(sat - station)
	TimeTrans = Tr - Geo / c
	h, posHere = broadcast2posM(a,np.array([[2121,2121],[TimeTrans, TimeTrans + 3]]).transpose(),5)
	phi = wE * (Tr - TimeTrans)
	rot = np.array([[np.cos(phi), np.sin(phi), 0],[-np.sin(phi), np.cos(phi),0],[0,0,1]])
	changed = np.dot(rot,posHere[0,:])
	rangeNew = np.linalg.norm(changed - station)

	while abs(rangeNew - Geo)> 1e-6:

		Geo = rangeNew
		TimeTrans = Tr - Geo / c
		h, posHere = broadcast2posM(a,np.array([[2121,2121],[TimeTrans, TimeTrans + 3]]).transpose(),5)
		phi = wE * (Tr - TimeTrans)
		rot = np.array([[np.cos(phi), np.sin(phi), 0],[-np.sin(phi), np.cos(phi),0],[0,0,1]])
		changed = np.dot(rot,posHere[0,:])
		rangeNew = np.linalg.norm(changed - station)

	return rangeNew