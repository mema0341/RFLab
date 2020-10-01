from numpy import pi, sinh, arange, mod, zeros, ones, exp, shape, linalg, sqrt, median, average
import matplotlib.pyplot as plt

w = 1.5
b = 1

# Definitions
mu0 = 4*pi*1e-7
eps0 = 8.85e-12
epsr = 2

N = 30 # Terms in expansion
M = 30 # M = Segments

hx = w/M # hx is length of segments 

m = arange(1,M+1,1)
xm = -w/2+hx*(m-0.5)
xp = -w/2+hx*(m-0.5)

Ipm = zeros((M,M))

n = arange(1,N*2,2)

pp = 0
mm = 0
for pp in range(M):
    for mm in range(M):
        first_term = 1/(pi*eps0*epsr)
        second_term = (2*b)/pi

        third_term = 1/(n**2)

        tmp0 = abs(xp[pp]-xm[mm])
        fourth_term = exp((-n*pi*tmp0)/b)

        tmp1 = (n*pi)/(2*b)
        fifth_term = sinh(tmp1*hx)

        Ipm[pp, mm] = first_term*second_term*sum(third_term*fourth_term *fifth_term)

V = ones((M))
Ipm = linalg.inv(Ipm)
alpha = Ipm.dot(V)

# # Plot Charge Distribution
# fig, ax = plt.subplots()
# ax.plot(arange(-w/2,w/2,hx)/b, alpha*1e12,'*')
# ax.set_title("Charge Distribution")
# ax.grid(True)
# ax.set_ylim(2,14)
# ax.set_xlim(-0.8,0.8)
# ax.set_xlabel("x/b")
# ax.set_ylabel("Charge Density")
# plt.show()

print(sum(alpha))
Z0 = sqrt(mu0*eps0*epsr)/sum(alpha)
print(Z0)
