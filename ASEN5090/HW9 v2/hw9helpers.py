import numpy as np
from pptx import Presentation
from pptx.util import Inches
import pyproj
import matplotlib.pyplot as plt

ECEF = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
LLA = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

R_PLUS = 6378.137 # km
F_P = 1/298.2572235363
E_PLUS = 2*F_P - np.power(F_P,2) # eccentricity

def get_ca(prn, N=1023, shift=0, convert=False):
    # Step 1: Load G1 and G2 with all ones
    shift_reg1 = np.ones(10)
    shift_reg2 = np.ones(10)

    # Step 2: Compute sums
    tapped1 = [2, 9]
    tapped2 = [1, 2, 5, 7, 8, 9]

    ca_code = ""
    ca_code_plot = np.ndarray((N))
    output_plot = np.ndarray((N))

    for i in range(N):
        # Grab Output
        output1 = shift_reg1[len(shift_reg1)-1]
        output2 = shift_reg2[len(shift_reg2)-1]

        # Grab Input
        input1 = np.mod(np.sum(shift_reg1[tapped1]),2)
        input2 = np.mod(np.sum(shift_reg2[tapped2]),2)

        g2i_output = np.mod(np.sum(shift_reg2[prn]),2)

        ca_code += str(int(np.mod(np.sum([output1, g2i_output]),2)))
        ca_code_plot[i] = np.mod(np.sum([output1, g2i_output]),2)

        # Shift
        shift_reg1 = np.roll(shift_reg1,1)
        shift_reg2 = np.roll(shift_reg2,1)

        shift_reg1[0] = input1
        shift_reg2[0] = input2

    # Convert 0 --> 1 and 1 --> -1
    if convert:
        idx0 = np.where(ca_code_plot==0)
        idx1 = np.where(ca_code_plot==1)

        ca_code_plot[idx0] = 1
        ca_code_plot[idx1] = -1
        
    
    return np.roll(ca_code_plot, shift)

def plot16(ca_code, prn_val, N=1023, show_plot=True):
    string_ca_code = ""
    for x in ca_code:
        string_ca_code = string_ca_code + str(int(x))

    fig, (ax1, ax2) = plt.subplots(1,2)

    fig.set_figwidth(fig.get_figwidth()*2)

    # Plot First 16 chips
    ax1.plot(np.arange(1,17,1), ca_code[0:16])
    # ax1.set_title('1023-chip C/A-code PRN%s, Epochs 1024-2046\nFirst 16 Chips: %s' % (str(prn_val), str(hex(int(string_ca_code[0:16],2)))))
    ax1.set_xlabel('Chips')
    ax1.set_ylabel('C/A-code PRN%s' % str(prn_val))
    ax1.set_xlim(1,16)
    ax1.set_ylim(-2, 2)
    ax1.grid(True)

    # Plot Last 16 Chips
    ax2.plot(np.arange(1,17,1), ca_code[int(N/2-16-1):int(N/2-1)])
    # ax2.set_title('1023-chip C/A-code PRN%s, Epochs 1024-2046\nLast 16 Chips: %s' % (str(prn_val), str(hex(int(string_ca_code[int(N/2-16-1):int(N/2-1)],2)))))
    ax2.set_xlabel('Chips')
    ax2.set_ylabel('C/A-code PRN%s' % str(prn_val))
    ax2.set_xlim(1,16)
    ax2.set_ylim(-2, 2)
    ax2.grid(True)

    if show_plot:
        plt.show()

# def lat2ecef(lat, lon, alt):

#     C_plus = R_PLUS / np.sqrt(1-np.power(E_PLUS,2)*np.power(np.sin(lat),2))

#     x = (C_plus+h) * np.cos(lat)*np.cos(lon)
#     y = (C_plus+h) * np.cos(lat)*np.sing(lon)
#     z = (C_plus*(1-E_PLUS)+alt)*np.sin(lat)

#     return x, y, z

# def ecef2lat(x, y, z):
#     lon_ = np.arctan(y/x)
#     row = np.sqrt(np.power(x,2)+np.power(y,2))

#     r = np.sqrt(np.power(x,2)+np.power(y,2)+np.power(z,2))

#     gd_old = np.arcsin(z/r)

#     C_plus = r_plus / np.sqrt(1-E_PLUS*np.power(np.sin(gd_old),2))

#     gd_new = np.arctan( (z + C_plus*E_PLUS*np.sin(gd_old)) / row)


def ecef2lat(x, y, z):
    lat_, lon_, alt_ = pyproj.transform(ecef, lla, x, y, z)

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

def compute_LOS_ENU(userECEF, satECEF):
    # Part A: Compute the vector pointing from user location to GPS
    vector = satECEF - userECEF
    # Part B: Normalize it to a unit vector (LOS_ECEF)
    ehat_los = vector / np.sqrt(np.power(vector[0], 2) \
         + np.power(vector[1], 2) + np.power(vector[2], 2))
    # Part C: Transform the Unit Vector from ECEF to ENU using the function you wrote in part 2
    # First Grab Lat and lon values 
    lon_, lat_, alt_ = pyproj.transform(ECEF, LLA, userECEF[0], userECEF[1], userECEF[2], radians=True)
    # Now grab conversion matrix
    conversion_matrix = ecef2enu(lon_, lat_)
    # Now take the dot product to get the ENU coords
    enu_coords = np.dot(conversion_matrix,ehat_los)
    # Return ENU Coords
    return enu_coords

# COMPUTE_AZELRANGE function takes ECEF coords of a satellite and a user, converts
# to get the LOS vector in ENU, and returns az, el, and range of the satellite
def compute_azelrange(userECEF, satECEF):
    # First get the LOS in ENU coords
    x_east, x_north, x_up = compute_LOS_ENU(userECEF, satECEF)

    # Compute Azimuth using LOS 
    az_ = np.arctan2(x_east,x_north)

    # Compute Elevation using LOS 
    el_ = np.arcsin(x_up / np.sqrt(np.power(x_east, 2)\
        +np.power(x_north, 2)+np.power(x_up, 2)))

    # Compute Range
    diff = userECEF-satECEF
    range = np.sqrt(np.power(diff[0],2) + np.power(diff[1],2) + np.power(diff[2],2))
    # Returns Az [degrees], El [degrees], and X_up [m]
    return az_, el_, range

