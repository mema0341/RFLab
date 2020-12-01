import numpy as np
import matplotlib.pyplot as plt
from helper_functions import get_ca, get_correlation, plot_correlation #, plot16, get_correlation, plot_partf

# Define a code and the recieved a code
a_code =  [ 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, -1, -1]
a_recieved =  [-1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1]

# Obtain the auto-correlation for a-code and cross correlation with the recieved a-code
auto_correlation = get_correlation(a_code, a_code)
cross_correlation = get_correlation(a_code, a_recieved)

# Plot Auto Correlation and Cross Correlation
fig, ax = plt.subplots()
ax.plot(cross_correlation)
ax.plot(auto_correlation)
ax.set_title('Correlation')
ax.set_xlabel('LAG')
ax.set_ylabel('CORRELATION')
ax.legend(['Cross-Correlation', 'Auto-Correlation'])
ax.set_ylim(-5, 15)
ax.grid(True)
# plt.show()

shift = len(a_code) - 3
a_shift =  np.roll(a_recieved, shift)
cross_correlation = get_correlation(a_code, a_shift)

# Plot Auto Correlation and Revised Cross Correlation
fig, ax = plt.subplots()
ax.plot(cross_correlation)
ax.plot(auto_correlation,'--')
ax.set_title('Correlation, A-Recieved Shifted %s Bits' % str(len(a_code)-3))
ax.set_xlabel('LAG')
ax.set_ylabel('CORRELATION')
ax.legend(['Cross-Correlation', 'Auto-Correlation'])
ax.set_ylim(-5, 15)
ax.grid(True)
# plt.show()

# If the chipping rate of the A-code is 1 MHz, and the delay is found to be 5 chips (this is not the right
# answer to part c), what is the delay in seconds? In meters?
chip_rate = 3e6 # 1/s
# chips = 3 # num_chips
chips = len(a_code) # num_chips
c = 3e8 # m/s

delay_sec = chips/chip_rate
print(delay_sec)

delay_m = delay_sec * c 
print(delay_m)