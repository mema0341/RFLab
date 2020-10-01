# A ground station observes an aircraft flying by at constant speed (v0), overhead at
# constant altitude (h0), traveling from horizontal positions –x0 to +x0. Write a wellorganized program/code to simulate this motion,
# compute and plot the ideal zenith angle, range, and range rate observations as a function of the of the horizontal position of the aircraft (x).
# Plot the measurements for the following values: v0=50m/s, x0=250m, h0=100 m.
# Assume that the aircraft height is known perfectly and the measurements are to be
# used to estimate the horizontal position of the aircraft, x.
# Compare the value/utility of the different measurement types as the aircraft goes
# by. Do the observations get better if the aircraft is flying at a lower altitude? Explain.

import numpy as np
import matplotlib.pyplot as plt

v0 = 50 # m/s
h0 = 100
x0 = 250

x = np.arange(-x0, x0+10, 10)
t = x/v0
r = np.sqrt(np.power(x,2)+np.power(h0,2))

r_rate = []
for tt in range(1,len(x)):
    r_rate.append( (r[tt] - r[tt-1]) / (t[tt] - t[tt-1]) )

gamma = np.rad2deg(np.arccos(h0/r))

fig, ax = plt.subplots()
ax.plot(x, r)
ax.set_title('Range [meters]')
ax.set_xlabel('X-Location [meters]')
ax.set_ylabel('Range [meters]')
ax.grid(True)
plt.xlim(min(x),max(x))
plt.ylim(min(r),max(r))
plt.grid(True)

fig, ax = plt.subplots()
ax.plot(x[1:len(x)], r_rate, color='r')
ax.set_title('Range [meters]')
ax.set_xlabel('X-Location [meters]')
ax.set_ylabel('Range Rate [meters/second]')
plt.xlim(min(x),max(x))
plt.ylim(min(r_rate),max(r_rate))
plt.grid(True)

fig, ax = plt.subplots()
ax.plot(x, gamma, color='g')
ax.set_title('Zenith Angle [Degrees]\nWhat')
ax.set_xlabel('X-Location [meters]')
ax.set_ylabel('Zenith Angle [Degrees]')
plt.xlim(min(x),max(x))
plt.ylim(min(gamma),max(gamma))
plt.grid(True)

plt.show()