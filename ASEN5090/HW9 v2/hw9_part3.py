import numpy as np
from scipy.io import loadmat
from hw9helpers import *

NUM_SAMPLES = 5000

# Load the sample data set – this is the received signal sR.
sr_data = np.squeeze(loadmat("HW9data.mat")['gpsdata']) # 20 ms of 8-bit real samples.

# Create a time vector at intervals of ∆ts corresponding to 1ms worth of samples:
sampling_frequency = 5e6 # 5 MHz
t_i = np.arange(0,1e-3,1/sampling_frequency)

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

# Set up the delay axis of the grid in steps corresponding to the sample interval.
delay_axis = np.arange(0,1e-3,1/sampling_frequency)

taus = np.arange(0,5000,1)

# Set up the Doppler axis of the grid in steps of 1kHz / (2 X integration time in msec).
# Tint = 1e-3 # integration is 1 ms
# doppler_step = 1/(2*Tint) # 500 Hz 
# doppler_axis = np.arange(0,NUM_SAMPLES+1,1)*2*doppler_step
doppler_axis = np.linspace(-6000,6000,25)


# Compute and display a 3D mesh plot of the magnitude of S(F,t) for PRN31 using 1ms of data.
sfts = np.zeros((len(doppler_axis),len(delay_axis)))
for dd in range(len(doppler_axis)):
    for tt in range(len(delay_axis)):
        theta_i = 2*np.pi*(if_freq + doppler_axis[dd])*t_i
        sr_shift = np.roll(sr_data,-taus[tt])[0:NUM_SAMPLES]
        summation = 0
        for ii in range(len(t_i)):
            sfdt = sr_shift[ii]*y[ii]*np.exp(-1j*theta_i[ii])
            summation = summation + sfdt
        sfts[dd,tt] = abs(summation)
sfts = sfts/np.max(sfts)

fig, ax = plt.subplots() 
im = ax.contour(doppler_axis, delay_axis, sfts)
plt.show()