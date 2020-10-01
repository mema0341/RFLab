# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 16:15:41 2019

@author: chels
functions provided to 
1. parse through almanac files to find ephemeris info
2. read in ephemeris and calculate ECI position (km)

"""
import pandas as pd
import numpy as np
from scipy import optimize
from constants import MU, EARTH_ROT_RATE

def read_GPSyuma(yumafilename):
    """Parses through almanac file and returns ephemeris

    inputs:
        yumafilename (str) --> name of almanac file

    outputs:
        ephs (list[dict]) --> list of dictionaries. separate dictionary for each PRN
        ephsdf (pandas.Dataframe) --> dataframe with columns 
                    "prn"
                    "health"
                    "ecc"
                    "toa" (s)
                    "incl" (rad)
                    "ra_rate" (r/s)
                    "sqrt_a" (m ^ 1/2)
                    "Loa" (rad)
                    "perigee" (rad)
                    "M0" (rad)
                    "Af0" (s)
                    "Af1" (s/s)
                    "GPS_week"
                almanac files do not contain the information below, so these
                column values are set to 0 as their existence is needed
                in the broadcast2pos function (ephemeris --> ECI position)
                    "delta_n"    
                    "i_rate"
                    "Cuc"  
                    "Cus"
                    "Crc"    
                    "Crs"
                    "Cic"    
                    "Cis"
                    "IODE"    
                    "Toc"
                    "Af2"
    example use:
        ephs, ephsdf = read_GPSyuma("240.ALM")
        ephs[0] is a dictionary containing info about one PRN
        ephsdf["prn"] is a series of all the PRNs
        ephs["ecc"] is a series of all the eccentricities
    
    """
    with open(yumafilename) as f: 
        contents = f.readlines()
        ephs = [] #initialize ephemeris list
        for tline in contents:
#            if tline[0] == "*": #start of new line
#                eph = {} #initialize dictionary
            if tline[:2] == "ID":
                prn = int(tline[3:])
            elif tline[:6] == "Health":
                health = int(tline[7:])
            elif tline[:12] == "Eccentricity":
                ecc = float(tline[13:])
            elif tline[:24] == "Time of Applicability(s)":
                toa = float(tline[25:])
            elif tline[:24] == "Orbital Inclination(rad)":
                incl = float(tline[25:])
            elif tline[:24] == "Rate of Right Ascen(r/s)":
                ra_rate = float(tline[25:])
            elif tline[:16] == "SQRT(A)  (m 1/2)":
                sqrt_a = float(tline[17:])
            elif tline[:24] == "Right Ascen at Week(rad)":
                Loa = float(tline[25:])
            elif tline[:24] == "Argument of Perigee(rad)":
                perigee = float(tline[25:])
            elif tline[:14] == "Mean Anom(rad)":
                M0 = float(tline[15:])
            elif tline[:6] == "Af0(s)":
                Af0 = float(tline[7:])
            elif tline[:8] == "Af1(s/s)":
                Af1 = float(tline[9:])
            elif tline[:4] == "week":
                GPS_week = float(tline[5:])
                if GPS_week < 1024:
                    GPS_week += 2048 #Penny added to avoid modulo 1024
            elif tline == '\n':
                # Variables not included in the YUMA almanac (but are in broadcast ephemeris)
                delta_n = 0    
                i_rate = 0
                Cuc = 0    
                Cus = 0
                Crc = 0    
                Crs = 0
                Cic = 0    
                Cis = 0
                IODE = 0    
                Toc = 0 #Time of clock
                Af2 = 0
                eph = {"prn": prn, "M0": M0, "delta_n": delta_n, "ecc": ecc, 
                       "sqrt_a": sqrt_a, "Loa": Loa, "incl": incl,
                       "perigee": perigee, "ra_rate": ra_rate, "i_rate": i_rate,
                       "Cuc": Cuc, "Cus": Cus, "Crc": Crc, "Crs": Crs, "Cic": Cic,
                       "Cis": Cis, "toa": toa, "IODE": IODE, "GPS_week": GPS_week,
                       "Toc": Toc, "Af0": Af0, "Af1": Af1, "Af2": Af2, "0": 0, 
                       "health": health}
                ephs.append(eph)
                continue
        delta_n = 0    
        i_rate = 0
        Cuc = 0    
        Cus = 0
        Crc = 0    
        Crs = 0
        Cic = 0    
        Cis = 0
        IODE = 0    
        Toc = 0 #Time of clock
        Af2 = 0
        eph = {"prn": prn, "M0": M0, "delta_n": delta_n, "ecc": ecc, 
               "sqrt_a": sqrt_a, "Loa": Loa, "incl": incl,
               "perigee": perigee, "ra_rate": ra_rate, "i_rate": i_rate,
               "Cuc": Cuc, "Cus": Cus, "Crc": Crc, "Crs": Crs, "Cic": Cic,
               "Cis": Cis, "Toa": toa, "IODE": IODE, "GPS_week": GPS_week,
               "Toc": Toc, "Af0": Af0, "Af1": Af1, "Af2": Af2, "0": 0, 
               "health": health}
        ephs.append(eph) #append last eph (no '\n' at end of file)
        ephsdf = pd.DataFrame(ephs)
    
    return ephs, ephsdf


def mean2eccentric(M, ecc):
    """E = mean2eccentric(M, e)
    Calculates the eccentric anomaly from the mean anomaly using Newton's 
    method. The tolerance is set to 1e-12.
 
    Author: Ben K. Bradley
    Last Revision Date: 26 October 2010
    
    Edited by Chelsa Thangavelu and Hunter Mellema to convert matlab script to python
    August 28, 2019
    
    inputs:
        M (float) --> mean anomaly (radians)
        ecc (float) --> eccentricity orbit
        
    outputs:
        E (float) --> eccentric anomaly (radians)
    
    References:
    [1] Curtis, H.D. "Orbital Mechanics for Engineering Students".
          Elsevier Ltd. 2005."""

    if ecc > .999 or ecc < 0.:
        M = np.nan
        print("Not circular or elliptical")

    #==========================================================================
    
    # Set tolerance
    twopi   = 2. * np.pi
    tol     = 1.0e-12
    maxiter = 20
    
    
    # Make sure mean anomaly is positive and within 2pi
    M = np.remainder(M, twopi)
    
    M = (M < 0.) * twopi + M
    
    
    #  Make an initial guess for E, Eccentric Anomaly
    sinM = np.sin(M)
    E = M + (ecc * sinM) / (1. - np.sin(M + ecc) + sinM)
    
    
    # Initialize iteration count and error
    count = 1
    err  = 1
    
    
    # Iterate to find E, eccentric anomaly
    while abs(err) > tol and count <= maxiter:
        err  = (E - ecc * np.sin(E) - M) / (1. - ecc * np.cos(E))
    
        E = E - err
        
        count += 1
        
        if (count > maxiter):
            print('Iterations maxed out in mean2eccentric')
    return E

def broadcast2pos(ephem_all, t_input, prn):
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
    sz         = t_input.shape[0]
    x          = np.full((sz, 3), np.nan)
    health     = np.full((sz, 1), np.nan)
    
    
    #==========================================================================
    # Pull Out Correct Ephemerides 
    #==========================================================================
    
    # Pull out ephemerides for PRN in question
    sat_ephem = ephem_all.loc[ephem_all["prn"] == prn]

    
    # No matching PRN found, returning data will be NaNs
    #otherwise continue
    if not sat_ephem.empty:
        
        #==========================================================================
        # Start Main Calculation Loop 
        #==========================================================================
        
        # Compute elapsed times of each ephemeris epoch wrt first entry, seconds
        GPS_weeks = sat_ephem["GPS_week"].values
        GPS_week0 = sat_ephem["GPS_week"].values[0]
        toes = sat_ephem["toa"].values
        toe0 = sat_ephem["toa"].values[0]
        dt_ephem = (GPS_weeks - GPS_week0) * 604800 + (toes - toe0)
        # print(dt_ephem)
        
        # Compute elapsed times of each input time wrt first ephemeris entry, seconds
        dt_input = (t_input[:,0] - GPS_week0) * 604800 + (t_input[:, 1] - toe0)
        

        
        for tt in range(sz): # loop through all input times
        
        
            # Pull out most recent ephemeris values
        #     jj = max( find(dt_input(tt) >= dt_ephem) ); # sat_ephem(:,17) = toe (sec into GPS week) of each entry
                                                        # jj = row of specific sat. ephem. data with epoch closest to input time
            
            # Pull out nearest ephemeris values                                                                                        
            min_t_diff = np.abs( dt_input[tt] - dt_ephem ).min()
            jj = np.where(np.abs(dt_input[tt] - dt_ephem) == min_t_diff) #index of minimum
                
            
            if jj[0].size < 1: # no matching ephemeris time found. continue to next input time 
                continue
            
        
        
            # Pull out common variables from the ephemeris matrix
            #======================================================================
            #toe = sat_ephem(jj,17);           # time of ephemeris
            dt  = dt_input[tt] - dt_ephem[jj] # seconds difference from epoch
            
            a = sat_ephem["sqrt_a"].iloc[jj].values ** 2 # semimajor axis
            ecc = sat_ephem["ecc"].iloc[jj].values #eccentricity
            n0 = np.sqrt(muE/a**3) #nominal mean motion (rad/s)
            n   = n0 + sat_ephem["delta_n"].iloc[jj].values * dt #corrected mean motion
            M   = sat_ephem["M0"].iloc[jj].values + n * dt #mean anomaly
        
#           Compute perigee, true and eccentric anomaly...
#           ======================================================================
        
#            Load argument of perigee to a local variable and add perigee rate, rad
            perigee = sat_ephem["perigee"].iloc[jj].values
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
        
        
#            Compute radius and inclination
#            ======================================================================
        
            r   = a * (1. - ecc * cosE)                        # corrected radius  
            inc = sat_ephem["incl"].iloc[jj].values   #  inclination 
        
            cosu = np.cos(u)    
            sinu = np.sin(u)  
        
#            Compute satellite position in orbital plane (Eq. 13)
#            ======================================================================
            xo = r * cosu    # satellite x-position in orbital plane
            yo = r * sinu    # satellite y-position in orbital plane
#        
#            Corrected longitude of ascending node for node rate and Earth rotation
#            ======================================================================
            node = sat_ephem["Loa"].iloc[jj].values + \
            (sat_ephem["ra_rate"].iloc[jj].values - wE) * dt - \
            (wE * sat_ephem["toa"].iloc[jj].values)
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
            health[tt,0] = sat_ephem["health"].iloc[jj].values # satellite health (0.00 is useable)
        
        
        
        
        
#        END of t_input loop =================================================
#        ==========================================================================    

    return health, x