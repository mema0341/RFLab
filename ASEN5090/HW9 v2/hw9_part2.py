import numpy as np
from scipy.io import loadmat
from hw9helpers import *

# Load the sample data set – this is the received signal sR.
sr_data = np.squeeze(loadmat("HW9data.mat")['gpsdata']) # 20 ms of 8-bit real samples.

# Create a time vector at intervals of ∆ts corresponding to 1ms worth of samples:
sampling_frequency = 5e6 # 5 MHz
t_i = np.arange(0,1e-3,1/sampling_frequency)
num_samples = len(t_i)

# Create a vector of the corresponding PRN31 C/A code values (+1/-1) for each element of time vector
prn31 = [2, 7]
x = get_ca(prn31, convert=True)

# Sample C/A Code 1024
f0 = 1.023e6
y = []
Tc = 1/f0
idx = np.mod(np.floor(t_i/Tc), 1023)
y = []
for tt in range(len(t_i)):
    y.append(x[int(idx[tt])])
y = np.asarray(y)

# Create a vector of the corresponding IF carrier phase for the given time vector and Doppler.
if_freq = 1.25e6 # IF = 1.25 MHz 

# Shift Tests
doppler_freq = [0,1000] 
tau = [9, 2943]
for dt in range(len(doppler_freq)):
    theta_i = 2*np.pi*(if_freq + doppler_freq[dt])*t_i
    sr_shift = np.roll(sr_data,-tau[dt])[0:num_samples]
    summation = 0
    for ii in range(len(t_i)):
        sfdt = sr_shift[ii]*y[ii]*np.exp(-1j*theta_i[ii])
        summation = summation + sfdt
    print("Tau: %s samples" % tau[dt])
    print("Doppler Frequency: %s Hz" % doppler_freq[dt])
    print("Complex Correlator: %s" % summation)
    print('')
