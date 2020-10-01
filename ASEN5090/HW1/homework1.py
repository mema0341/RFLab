import os
import numpy as np

file = open("HW1data.dat", "r")

column0 = []
column1 = []
for f in file:
    line = f.split()
    for col in range(len(line)):
        column0.append(np.float(line[0]))
        column1.append(np.float(line[1]))

print(np.mean(column0))
print(np.std(column0))

print(np.mean(column1))
print(np.std(column1))