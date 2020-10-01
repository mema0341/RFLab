import numpy as np
import matplotlib.pyplot as plt
from hw2_helpers import get_ca, plot16, get_correlation, plot_correlation, plot_partf

N = 1023

# 2-1 Part A
prn_val = 19
prn19 = [2, 5]
ca_code19 = get_ca(prn19, N, convert=True)
# plot16(ca_code19, N, prn_val, show_plot=False)

# 2-1 Part B
N = 1024*2
ca_code = get_ca(prn19, N)
new_ca_code = ca_code[int(N/2-1):N-1]
# plot16(new_ca_code, N, prn_val)

# # 2-1 Part C
N = 1024
prn_val = 25
prn25 = [4, 6]
ca_code25 = get_ca(prn25, N)
# plot16(ca_code25, N, prn_val)

# # 2-1 Part D
prn_val = 5
prn5 = [0, 8]
ca_code5 = get_ca(prn5, N)
# plot16(ca_code5, N, prn_val)

# # 2-2 Part A
prn19_autocorrelation = get_correlation(ca_code19, ca_code19, N)
# plot_correlation(prn19_autocorrelation, 'Auto-Correlation PRN 19')

# # 2-2 P/art B
ca_code19_shifted = get_ca(prn19, N, shift=200, convert=True) # C/A code
cross_correlation_19shifted = get_correlation(ca_code19, ca_code19_shifted, N)
# plot_correlation(cross_correlation_19shifted, 'Cross-Correlation PRN 19 and PRN 19 Shifted')

# # 2-2 Part C
cross_correlation19_25 = get_correlation(ca_code19, ca_code25, N)
# plot_correlation(cross_correlation19_25, 'Cross-Correlation PRN %s and PRN %s' % (str(19), str(25)))

# 2-2 Part D
cross_correlation19_5 = get_correlation(ca_code19, ca_code5, N)
# plot_correlation(cross_correlation19_5, 'Cross-Correlation PRN %s and PRN %s' % (str(19), str(5)))

# 2-2 Part E
x1 = get_ca(prn19, N, shift=350, convert=True) # C/A Code using PRN19 shifted 350 chips
x2 = get_ca(prn25, N, shift=905, convert=True) # C/A Code usingPRN 25 delayed by 905 chips
x3 = get_ca(prn5, N, shift=75, convert=True) # C/A Code using PRN 5 Delayed by 75 chips
sum_x = x1 + x2 + x3
cross_correlation19_sumx = get_correlation(ca_code19, sum_x, N)
# plot_correlation(cross_correlation19_sumx, 'Cross-Correlation PRN 19 and Sum of x1, x2, and x3')

# 2-2 Part F
noise = 4*np.random.randn(N)
# plot_partf(x1, x2, x3, noise, N)

# 2-2 Part G
sum_x_noise = x1 + x2 + x3 + noise
cross_correlation19_sumxnoise_19 = get_correlation(ca_code19, sum_x_noise, N)
plot_correlation(cross_correlation19_sumxnoise_19, 'Cross-Correlation PRN 19 and Sum of x1, x2, x3, and Noise')
