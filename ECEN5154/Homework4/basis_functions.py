import numpy as np
# Domain Size 0 <= x < = 1

N = 3

hx = 1/(N+1)

m = np.arange(1,5,1)

xm = m/(N+1)

print(xm)

for x in range(N):
    print('Xm = %s, Lower Bound = %s, Upper Bound: %s' % (str(xm[x]), str(xm[x]-hx/2), str(xm[x]+hx/2)))
