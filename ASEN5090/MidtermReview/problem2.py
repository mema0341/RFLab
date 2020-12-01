
import numpy as np

weekdays = ["Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"]
id_ = 99
health = 0
eccentricty = 0
toa = 270000
orbital_inclination = 0.96 # rad
rate_of_ascension = 0
sqrt_a = 5153.65
right_ascen_at_week = 0.1e001
af0 = 0.3e-3
af1 = 0.3e-11
week = 2

mod_num = 24*60*60

day_of_week = int(np.mod(toa, mod_num)/3600)
time_of_week = np.remainder(toa,mod_num)
print(weekdays[day_of_week])
