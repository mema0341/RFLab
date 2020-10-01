from numpy import pi, sinh, arange, mod, zeros, ones, exp, shape, linalg, sqrt, median, average
import matplotlib.pyplot as plt

w = 2.5
b = 1

# Definitions
mu0 = 4*pi*1e-7
eps0 = 8.85e-12
epsr = 1

N = 30 # Terms in expansion
M = 9 # M = Segments

hx = w/M # hx is length of segments 

m = arange(1,M+1,1)
xm = -w/2+hx*(m-0.5)
xp = -w/2+hx*(m-0.5)

Ipm = zeros((M,M))

n = arange(1,N,2)

pp = 0
mm = 0
for pp in range(M):
    for mm in range(M):
        Ipm[pp, mm] = (1/(pi*eps0*epsr))*((2*b)/pi)*sum((1/(n**2))*exp((-n*pi*abs(xp[pp]-xm[mm])/b)*sinh(((n*pi)/(2*b)*hx))))

V = ones((M))
Ipm = linalg.inv(Ipm)
alpha = Ipm.dot(V)

# # Plot Charge Distribution
# fig, ax = plt.subplots()
# ax.plot(arange(-w/2,w/2,hx)/b, alpha*1e12,'*')
# ax.set_title("Charge Distribution")
# ax.grid(True)
# ax.set_ylim(0,14)
# ax.set_xlim(-w/(2*b)-hx/2,w/(2*b)+hx/2)
# ax.set_xlabel("x/b")
# ax.set_ylabel("Charge Density")
# plt.show()

print(alpha)
print(sum(alpha))
Z0 = sqrt(mu0*eps0*epsr)/sum(alpha)
print(Z0)
