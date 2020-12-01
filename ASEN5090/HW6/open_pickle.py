import pickle
import numpy as np 
from matplotlib import pyplot as plt

# directory = "NIST_pickles" # USN8_pickles
base_station = "NIST" # NIST
directory = "%s_pickles\\" % base_station # NIST_pickles
prft_filename = directory + "stored_prefit_residuals_%s.pickle" % base_station.lower()
posft_filename = directory + "stored_postfit_residuals_%s.pickle" % base_station.lower()
el_filename = directory + "stored_elevation.pickle_%s" % base_station.lower()
relposENU_filename = directory + "stored_relposENU_%s.pickle" % base_station.lower()
dops_filename = directory + "stored_dops.pickle_%s" % base_station.lower()

nist_relposENU_filename = "NIST_pickles\stored_relposENU_nist.pickle"

with open(prft_filename, 'rb') as f:
    prefit_residuals = pickle.load(f)

with open(posft_filename, 'rb') as f:
    postfit_residuals = pickle.load(f)

with open(el_filename, 'rb') as f:
    elevation_values = pickle.load(f)

with open(relposENU_filename, 'rb') as f:
    relative_positionENU = pickle.load(f)

with open(nist_relposENU_filename, 'rb') as f:
    nist_relative_positionENU = pickle.load(f)

with open(dops_filename, 'rb') as f:
    dops = pickle.load(f)

flat_prefit_residuals = []
for sublist in prefit_residuals:
    for item in sublist:
        flat_prefit_residuals.append(item)

flat_postfit_residuals = []
for sublist in postfit_residuals:
    for item in sublist:
        flat_postfit_residuals.append(item)
        
flat_elevation = []
for sublist in elevation_values:
    for item in sublist:
        flat_elevation.append(item)

num_epochs = len(prefit_residuals)

num_sats = []
epoch = []
for ee in range(num_epochs):
    epoch.append(ee*30/3600)
    # print(prefit_residuals[ee])
    num_sats.append(len(prefit_residuals[ee]))

    # if num_sats[ee] < 4:
    #     print(prefit_residuals[ee])

dops = np.asarray(dops)
relative_positionENU = np.reshape(relative_positionENU, (num_epochs, 4))
nist_relative_positionENU = np.reshape(nist_relative_positionENU, (num_epochs, 4))

fig, ax = plt.subplots()
plt.scatter(epoch, nist_relative_positionENU[:,3], s=1.5)
plt.scatter(epoch, relative_positionENU[:,3], s=1.5)
ax.grid(True)
plt.legend(['NIST', 'USN8'])
ax.set_title('NIST and USN8 Clock Solutions')
ax.set_ylabel('Clock Corrections [m]')
ax.set_xlabel('Time [hrs]')
plt.show()
plt.savefig('PicsHW7\clock_solutions')

fig, ax = plt.subplots()
plt.scatter(relative_positionENU[:,0], relative_positionENU[:,1], s=1.5)
ax.grid(True)
ax.set_title('dNorth vs. dEast [m] %s' % base_station)
ax.set_ylabel('North Relative Position [m]')
ax.set_xlabel('East Relative Position [m]')
plt.savefig('PicsHW7\%s\dNorthvsdEast %s' % (base_station, base_station))
ax.set_xlim(-10,10)
ax.set_ylim(-10,10)
plt.show()

fig, ax = plt.subplots()
plt.scatter(dops[:,0], nist_relative_positionENU[:,0], s=1.5)
ax.grid(True)
plt.legend(['NIST', 'USN8'])
ax.set_title('NIST and USN8 Clock Solutions')
ax.set_ylabel('Clock Corrections [m]')
ax.set_xlabel('Time [hrs]')
# plt.show()
plt.savefig('PicsHW7\clock_solutions')

fig, (ax1, ax2, ax3) = plt.subplots(3,1)
ax1.scatter(epoch, dops[:,0], marker='o', color='C0', alpha=1.0, s=1.5)
ax2.scatter(epoch, dops[:,1], marker='o', color='C1', alpha=1.0, s=1.5)
ax3.scatter(epoch, dops[:,2], marker='o', color='C2', alpha=1.0, s=1.5)

ax1.set_ylabel('EDOP [m]')
ax2.set_ylabel('NDOP [m]')
ax3.set_ylabel('VDOP [m]')

ax1.grid(True)
ax2.grid(True)
ax3.grid(True)

ax1.set_ylim(0,10)
ax2.set_ylim(0,10)
ax3.set_ylim(0,10)

ax1.set_title('DOP Values %s [m]' % base_station)
ax3.set_xlabel('Time [hrs]')

# plt.show()
plt.savefig("PicsHW7\%s\DOPs %s"% (base_station, base_station))

# print(relative_positionENU)
fig, (ax1, ax2, ax3) = plt.subplots(3,1)
ax1.scatter(epoch, relative_positionENU[:,0], marker='o', color='C0', alpha=1.0, s=1.5)
ax2.scatter(epoch, relative_positionENU[:,1], marker='o', color='C1', alpha=1.0, s=1.5)
ax3.scatter(epoch, relative_positionENU[:,2], marker='o', color='C2', alpha=1.0, s=1.5)

ax1.set_ylabel('EAST [m]')
ax2.set_ylabel('North [m]')
ax3.set_ylabel('Up [m]')

ax1.set_title('Solution Errors %s [m]' % base_station)
ax3.set_xlabel('Time [hrs]')

ax1.grid(True)
ax2.grid(True)
ax3.grid(True)

ax1.set_ylim(-10,10)
ax2.set_ylim(-10,10)
ax3.set_ylim(-10,10)
plt.savefig("PicsHW7\%s\Solution Errors %s"% (base_station, base_station))
# plt.show()

fig, ax = plt.subplots() 
plt.scatter(epoch, num_sats, color='C0', s=3)
plt.grid()
plt.title('Number of Sats Used in Point Solution %s' % base_station)
plt.xlabel('Time (hrs)')
plt.ylabel('Num Sats')
# plt.show()
plt.savefig("PicsHW7\%s\HW7_NumSats_%s" % (base_station, base_station))

fig, ax = plt.subplots() 
for ee in range(num_epochs):
    plt.scatter(elevation_values[ee], prefit_residuals[ee], color='C0', s=1.5)
plt.grid()
plt.title('Prefit Residuals vs Elevation %s' % base_station)
plt.xlabel('Elevation (deg)')
plt.ylabel('Prefit ResidualsPrefit Residuals (m)')
# plt.show()
plt.savefig("PicsHW7\%s\HW7_PrefitvsElevation_%s" % (base_station, base_station))

fig, ax = plt.subplots() 
for ee in range(num_epochs):
    plt.scatter(elevation_values[ee], postfit_residuals[ee], color='C0', s=1.5)
plt.grid()
plt.title('Postfit Residuals vs Elevation %s' % base_station)
plt.xlabel('Elevation (deg)')
plt.ylabel('Postfit ResidualsPrefit Residuals (m)')
plt.savefig("PicsHW7\%s\HW7_PostfitvsElevation_%s" % (base_station, base_station))

fig, ax = plt.subplots() 
for ee in range(num_epochs):
    plt.scatter(np.ones(np.shape(prefit_residuals[ee]))*ee*30/3600, prefit_residuals[ee], color='C0', s=1.5)
plt.grid()
plt.title('Prefit Residuals vs Time %s' % base_station)
plt.xlabel('Time (hrs)')
plt.ylabel('Prefit Residuals (m)')
plt.savefig("PicsHW7\%s\HW7_PrefitfitvsTime_%s" % (base_station, base_station))

fig, ax = plt.subplots() 
for ee in range(num_epochs):
    plt.scatter(np.ones(np.shape(postfit_residuals[ee]))*ee*30/3600, postfit_residuals[ee], color='C0', s=1.5)
plt.grid()
plt.title('Postfit Residuals vs Time %s' % base_station)
plt.xlabel('Time (hrs)')
plt.ylabel('Postfit Residuals (m)')
plt.savefig("PicsHW7\%s\HW7_PostfitvsTime_%s" % (base_station, base_station))

# print(np.shape(relative_positionENU))
eastENU = relative_positionENU[:,0]
northENU = relative_positionENU[:,1]
upENU = relative_positionENU[:,2]

# Comput RMS Values
N = len(flat_prefit_residuals)
rms_prefit = np.sqrt(np.sum(np.power(flat_prefit_residuals,2)/N))
rms_postfit = np.sqrt(np.sum(np.power(flat_postfit_residuals,2)/N))

rms_east = np.sqrt(np.sum(np.power(eastENU,2)/N))
rms_north = np.sqrt(np.sum(np.power(northENU,2)/N))
rms_up = np.sqrt(np.sum(np.power(upENU,2)/N))

print(base_station)
print("East Coord STD: %s" % np.round(np.std(eastENU),2))
print("East Coord RMS: %s m" % np.round(rms_east,2))
print()
print("North Coord STD: %s" % np.round(np.std(northENU),2))
print("North Coord RMS: %s m" % np.round(rms_north,2))
print()
print("Up Coord STD: %s" % np.round(np.std(upENU),2))
print("Up Coord RMS: %s m" % np.round(rms_up))
print()
print("Pre-fit Residuals STD: %s" % np.round(np.std(np.squeeze(flat_prefit_residuals)),2))
print("Pre-fit Residuals RMS: %s m" % np.round(rms_prefit,2))
print()
print("Post-fit Residuals STD: %s" % np.round(np.std(np.squeeze(flat_postfit_residuals)),2))
print("Post-fit Residuals RMS: %s m" % np.round(rms_postfit,2))

