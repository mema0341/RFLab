from numpy import sqrt, pi, zeros, exp, mod, log, shape, real, imag

N = 3;N_ext = 2*N-1 # This is five

k0 = 2*pi; eps_0 = 120*pi # 377 Impedance of free space
l_dipole = 0.1;r_dipole = 1e-4

half_dipole = l_dipole/2
idx_center = (N+1)/2 # center of diple

dz = 2*half_dipole/N
zm = half_dipole-dz/2
a=-0.5*dz; b=0.5*dz

simpsons_reps = 79
hd = (b-a)/(simpsons_reps+1)
lmbda = 1

# Impedance Matrices
zmnw = zeros((N_ext,N_ext),dtype='complex')
zmn = zeros((N_ext),dtype='complex')
I = zeros((N),dtype='complex') # Current Vector

# Choice of Magnetic Frill or Delta Gap
magnetic_frill = True
delta_gap = False
# magnetic_frill = False
# delta_gap = True

# Calculate Zmn Using Approximation of Integral since can't use scipy.quad (Reference: https://www.youtube.com/watch?v=7EqRRuh-5Lk)
for nn in range(N):
    # Part 1 (first term)
    zn = half_dipole - (nn+0.5)*dz; zal = zn - zm + a; R_a = sqrt(r_dipole**2 + zal**2)
    PockIntEq_a = exp(-1j*k0*R_a)/(4*pi*R_a**5)*((1+1j*k0*R_a)*(2*R_a**2-3*r_dipole**2) + (k0*r_dipole*R_a**2)**2)

    # Part 2 (second term)
    zbl = zn - zm + b; R_b = sqrt(r_dipole**2 + zbl**2)
    PockIntEq_b = exp(-1j*k0*R_b)/(4*pi*R_b**5)*((1+1j*k0*R_b)*(2*R_b**2-3*r_dipole**2) + (k0*r_dipole*R_b**2)**2)
    
    # Combination
    PockIntEq_t = PockIntEq_a + PockIntEq_b

    # Add to Impedance Matrix
    zmn[nn] = PockIntEq_t

    # Intermediate terms
    for rr in range(simpsons_reps):
        xk = a + (rr+1)*hd; zxl = zn-zm+xk
        
        R_x = sqrt(r_dipole**2 + zxl**2)
        PockIntEq_c = exp(-1j*k0*R_x)/(4*pi*R_x**5)*((1+1j*k0*R_x)*(2*R_x**2-3*r_dipole**2) + (k0*r_dipole*R_x**2)**2)
        # print(PockIntEq_c)
        # Check Out Simpson's Rule: https://www.youtube.com/watch?v=7EqRRuh-5Lk
        if mod(rr, 2):
            PockIntEq_t = PockIntEq_t + 2*PockIntEq_c
            # print(PockIntEq_t)
        else:
            PockIntEq_t = PockIntEq_t + 4*PockIntEq_c
            # print(PockIntEq_t)
    
    # Combination Step
    PockIntEq_t = PockIntEq_t*hd/3
    zmn[nn] = PockIntEq_t
    if nn != 0:
        zmn[N+nn-1] = PockIntEq_t

# Next we need to calculate the Current
rb = 2.3*r_dipole # need for magnetic frill
tlab = 2*log(2.3) # need for magnetic frill

# This is grabbing the current 
for i in range(N):
    zi = half_dipole - (i+0.5)*dz
    r1 = 2*pi*sqrt(zi**2 + r_dipole**2)
    r2 = 2*pi*sqrt(zi**2 + rb**2)

    # If magnetic frill, perform calculation
    if magnetic_frill:
        I[i] = -1j*(k0**2)/(eps_0*tlab)* (exp(-1j*r1)/r1-exp(-1j*r2)/r2)
    # If it's the delta gap, find the center value and calculate (The rest are zeros)
    if delta_gap:
        if i+1 == idx_center:
            I[i] = -1j*k0/(eps_0*dz)
        else:
            I[i] = 0

# Final Step, Solve System Matrix A*x=b (T-Matrix Toeplitz Variation)
m=N;a=zmn;b=I;a1=a;a2=a[m:len(a)];t=b;c1=zeros((N_ext),dtype='complex');c2=c1[m-2:len(a)-1];r1=a1[0];t[0]=b[0]/r1
if m != 0:
    for n in range(1,N):
        n1 = n-1;n2 = n-2;r5 = a2[n1];r6 = a1[n]
        if n != 1:
            c1[n1] = r2 
            for i1 in range(n2+1):
                i2=n-i1-1;r5=r5+a2[i1]*c1[i2];r6=r6+a1[i1+1]*c2[i1]
        r2=-r5/r1;r3=-r6/r1;r1=r1+r5*r3
        if n != 1:
            r6 = c2[0]
            c2[n1] = 0
            for i1 in range(1,n1):
                r5 = c2[i1];c2[i1] = c1[i1]*r3+r6;c1[i1] = c1[i1]+r6*r2
        c2[0] = r3;r5 = 0
        for i1 in range(n1+1):
            i2 = n-(i1+1);r5 = r5+a2[i1]*t[i2]
        r6 = (b[n]-r5)/r1
        for i1 in range(n1+1):
            t[i1] = t[i1] + c2[i1]*r6
        t[n] = r6
I = t

# Calculate Input Impedance
zin = 1/I[1]
print("###### Input Impedance ######\nReal(Zin): %s\nImag(Zin): %s" % (real(zin), imag(zin)))
