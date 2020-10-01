from numpy import pi, sinh, arange, mod, zeros, ones, exp, shape, linalg, sqrt, median, average
import matplotlib.pyplot as plt

# Definitions
mu0 = 4*pi*(10**-7) 
eps0 = 8.854*(10**-12)
epsr = 1


w = 1.5
b = 1

# Best Val seems to be 45
N = 30
M = N

hx = w/M

m = arange(1,M+1,1)
n = arange(1,N+1,2)
xm = -w/2 + hx*(m-0.5)
xp = -w/2 + hx*(m-0.5)

Ipm = zeros((M,M))
for pp in range(M):
    for mm in range(M):
        
        Ipm[pp,mm] = (1/(pi*eps0*epsr))*(2*b/pi)*sum((1/n**2)*exp(-n*pi*abs(xp[pp]-xm[mm])/b)*sinh(((n*pi)/(2*b))*hx))

V = ones((M))
Ipm = linalg.inv(Ipm)
alpha = Ipm.dot(V)

# Plot Charge Distribution
fig, ax = plt.subplots()
ax.plot(arange(-w/2,w/2,hx)/b, alpha*1e11,'*')
ax.set_title("Charge Distribution")
ax.grid(True)
ax.set_ylim(2,14)
ax.set_xlim(-0.8,0.8)
ax.set_xlabel("x/b")
ax.set_ylabel("Charge Density")
plt.show()


Z0 = sqrt(mu0*eps0*epsr)/sum(alpha)

print("Z0 = %s" % str(Z0))
print("Median C0 = %s" % str(average(median(alpha))))
print("Avg. C0 = %s" % str(average(alpha)))