import numpy as np

L1 = 1.57542
L2 = 1.2276

# scale_factor = np.sqrt(np.power(L1,2)/np.power(L2,2))
scale_factor = np.power(L1,2)/np.power(L2,2)

# L1delay = 1.5 # m
L1delay = 0.2 # m
L2delay = L1delay*scale_factor

print(L2delay)