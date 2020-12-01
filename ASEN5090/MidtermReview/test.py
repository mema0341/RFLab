from numpy import sqrt, arctan, sin, cos, degrees, deg2rad, tan, power

# # At North Pole
R_Earth = 6380
R_Gps = 26560
Inc_Gps = 55

# z = deg2rad(90 - Inc_Gps)

# z_prime = degrees(arctan(R_Gps*sin(z)/(R_Gps*cos(z)-R_Earth)))
# el_prime = 90 - z_prime
# print(el_prime)

mu = 398600.4
v_t = sqrt(mu/R_Gps)*1000 # m/s
# print(v_t)
c = 3e8 # m/s
f_t = 1575.42e6 # transmit frequency
f_d = v_t * f_t / c
# print(f_t-f_d)
# print(f_t+f_d)

# print(R_Gps-R_Earth)

x = arctan(deg2rad(90)) * R_Gps
# print(x)
dist = sqrt(power((R_Gps-R_Earth),2)+power(x, 2))
# print(dist)
# print(dist*1000/c)

clock_error = 300*1e-9
# clock_error = 1.334*1e-6
print(clock_error*c)
