import numpy as np

c = 299792458

total_distance = 1000
distance_a = 0
distance_b = total_distance - distance_a

time_a = distance_a / c * 1e6
time_b = distance_b / c * 1e6

# print(time_a)
# print(time_b)

clock_difference = 1.334

new_time_a = (time_a + clock_difference) /  1e6
new_time_b = (time_b + clock_difference) /  1e6

new_distance_a = new_time_a * c
new_distance_b = new_time_b * c

print(np.int(np.round(new_distance_a)), np.int(np.round(new_distance_b)))
