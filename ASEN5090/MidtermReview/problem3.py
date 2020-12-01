import numpy as np
import matplotlib.pyplot as plt
from helper_functions import get_ca, get_correlation, plot_correlation #, plot16, get_correlation, plot_partf

N = 15
STAGE = 4

shift_reg1 = np.ones(STAGE)

# Step 2: Compute sums
tapped1 = [3, 0]

a_code = []
for i in range(N):
    # Grab Output
    output1 = shift_reg1[len(shift_reg1)-1]
    input1 = np.mod(np.sum(shift_reg1[tapped1]),2)

    # Shift
    shift_reg1 = np.roll(shift_reg1,1)
    shift_reg1[0] = input1
    
    a_code.append(output1)
print(a_code)

ca_code = np.asarray(a_code)
convert = True
if convert:
    idx0 = np.where(ca_code==0)
    idx1 = np.where(ca_code==1)

    ca_code[idx0] = 1
    ca_code[idx1] = -1

print(ca_code)

auto_correlation = get_correlation(ca_code, ca_code, N)

b_code = [-0.8, 0.4, 0.6, 0.9, -0.6, 0.3, -0.2, -1.7, -0.8, 1.2, -0.8, 0.9, -0.9, -1.2, -0.6]

cross_correlation = get_correlation(ca_code, b_code, N)

fig, ax = plt.subplots()
ax.plot(auto_correlation)
ax.plot(cross_correlation)
ax.set_title('Correlation')
ax.set_xlabel('Chip Delay')
ax.set_ylabel('Correlation')
ax.legend(['Autocorrelation', 'Cross-Correlation'])
ax.set_ylim(-5, 15)
ax.grid(True)
# plt.show()

# If the chipping rate of the A-code is 1 MHz, and the delay is found to be 5 chips (this is not the right
# answer to part c), what is the delay in seconds? In meters?
chip_rate = 1e6 # 1/s
chips = 5 # num_chips
c = 3e8 # m/s

delay_sec = chips/chip_rate
print(delay_sec)

delay_m = delay_sec * c 
print(delay_m)