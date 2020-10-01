N = 3

l_dipole = 0.1
bottom_dipole = -l_dipole/2
top_dipole = l_dipole/2

delta_z = l_dipole/N

a = 1e-4

segments = np.arange(bottom_dipole, top_dipole+delta_z/2, delta_z)


# N = 1, m = 1
# Center of the segment m = 1 
point_m = segments[nn+1]+segments[nn])/2

lower_bound = segments[nn+1]
upper_bound = segments[nn+1]



