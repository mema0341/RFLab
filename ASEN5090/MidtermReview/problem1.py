import numpy as np
import pyproj
import helper_functions

ECEF = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
LLA = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

# For this problem, except where specifically indicated, assume the GPS orbits are circular with radius
# RGPS = 26,560,000 m and inclination= 55 deg, and the Earth is spherical with radius RE = 6,378,137 m. The
# table below givs locations of some satellites in the sky (not necessarily GPS satellites). 

# a) An observer (A) is located at latitude N 30 deg, longitude E 90 deg, height above the ellipsoid h=0.
# Assuming a spherical Earth with radius given above, compute the WGS-84 ECEF coordinates of A in
# meters. Express this result (x0, y0, z0) to 1 m precision.'

lat = 30
lon = 90
alt = 0 # above ellipsoid

# Testnt(pyproj.transform(ECEF, LLA, ECEFcoords[0], ECEFcoords[1], ECEFcoords[2]))

# b) Recompute the ECEF coordinates of observer A, accounting for the actual WGS-84 ellipsoid; that is,
# assume that N 30 deg is the geodetic latitude, with longitude and ellipsoid heights as given above
# and use the a and e values of the WGS-84 ellipsoid. Express this result (x1, y1, z1) to 1 m precision
# and also give the differences in the coordinates from part (a), i.e. (x1-x0, y1-y0, z1-z0). Does your
# result make sense?
a = 6378137


# c) Compute the magnitude of the velocity for each of the satellites in the table, assuming they are all in
# circular orbits. Which satellite would have the largest possible range rate (and largest magnitude
# Doppler shift) when tracked by observers on the Earth? 
userECEF = pyproj.transform(LLA, ECEF, lon, lat, alt)

satECEFs = [np.asarray([894992,  5075750, 4324763]), np.asarray([-3994192, 22652188, 13280000]), np.asarray([21189069, -36700543, 0]), np.asarray([0, 2314857, 26458931]), np.asarray([5953261 -22217873 -13280000])]

az_vals = []
el_vals = []
range_vals = []
for satECEF in satECEFs:
    # Sat 1
    [az_, el_, range_] = helper_functions.compute_azelrange(userECEF, satECEF)
    print("Az: %2.8f, El: %2.8f, Range: %2.8f" % (np.degrees(az_), np.degrees(el_), range_))
    az_vals.append(az_)
    el_vals.append(el_)
    range_vals.append(range_)
print(max(range_vals))
