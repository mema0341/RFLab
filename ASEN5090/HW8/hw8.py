from plot_scope_spectrum_analyzer import plot_scope_spectrum_analyzer
import numpy as np 
from matplotlib import pyplot as plt
from scipy import signal
from hw8_helpers import get_ca
import numpy.matlib
import pandas as pd

time_max = 10e-3

# Part 1
sampling_frequency = 50e6 # 50 MHz
t = np.arange(0,time_max,1/sampling_frequency)
nyquist_frequency = sampling_frequency / 2

# # Get frequency resolution
# res = sampling_frequency/len(t)
# print(res)

# Part 2
f = 200000
y = np.sin(2*np.pi*f*t)

# plot_scope_spectrum_analyzer(t, y, zoomt=4e-5, zoomf=[0,5e5])

# fourier_transform = np.fft(y)
# plt.plot(fourier_transform)
# plt.show() 

# Part 3
f = 10000
y = signal.square(2*np.pi*f*t)
# plot_scope_spectrum_analyzer(t, y, zoomt=4e-4, zoomf=[0,5e5])

# Part 4
f0hip = 1e-3

N = 1023
prn_val = 5
prn5 = [0, 8]
ca_code5, g1 = get_ca(prn5, N, convert=True)

f0 = 1.023e6
y = []
Tc = 1/f0
tidx = 1

ind = np.floor(t/Tc)
idx = np.mod(ind, N)

# y = []
# for tt in range(len(t)):
#     y.append(g1[int(idx[tt])])
# y = np.asarray(y)
# plot_scope_spectrum_analyzer(t, y, zoomt=4e-4, zoomf=[0,10000], zoomp=[-100, 0])

# print(t[485:495])
# print(y[485:495])

print(np.max(idx))
y = []
for tt in range(len(t)):
    y.append(ca_code5[int(idx[tt])])
y = np.asarray(y)
# plot_scope_spectrum_analyzer(t, y, zoomt=4e-4, zoomf=[0,5e6], zoomp=[-100, 0])

# print(t[485:495])
# print(y[485:495])

# Part 6 Use your maximal length G1 code to modulate a carrier signal (cosine wave) with frequency 5 times the
# chipping rate above.
car_sig = np.cos(2*np.pi*5*f0*t)
print(np.shape(car_sig))
mod_car_sig = car_sig*y
# plot_scope_spectrum_analyzer(t, mod_car_sig, zoomt=0.00001, zoomf=[f0/2,9*f0/2], zoomp=[-100, 0])

# Part 7 Implement a BOC(1,1) code using the G1 for the PRN code and the same carrier, code chipping rate, and
# time vectors used above.
square_wave = signal.square(2*np.pi*f0*t)
boc_code = square_wave*y*car_sig
plot_scope_spectrum_analyzer(t, boc_code, zoomt=0.00001, zoomp=[-100, 0], zoomf=[f0/2,f0*2])

# # Add white noise with standard deviation of 1V to the modulated signal you created in 6a. This represents
# # the signal arriving at the receiving antenna. 
# noise = np.random.normal(0,1,len(t))
# noisy_mod_car_sig = mod_car_sig + noise
# # plot_scope_spectrum_analyzer(t, noisy_mod_car_sig, zoomt=0.00001, zoomf=[1*f0,9*f0], zoomp=[-100, 0])

# # Multiply the noisy received signal by a perfectly aligned code replica.
# new_sig = noisy_mod_car_sig*y
# plot_scope_spectrum_analyzer(t, new_sig, zoomt=0.00001, zoomf=[1*f0,9*f0], zoomp=[-100, 0])
