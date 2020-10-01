from readyuma import *

# Grab satellite data
ephs, ephsdf = read_GPSyuma(r'data\YUMA245.ALM')

# Create time for entire day
time = np.arange(0,24*60+1,1)*60+172800
t_input = np.zeros((len(time), 2))

# Create input based off of time
t_input[:,0] = np.ones(np.shape(np.arange(0,24*60+1,1)*60+172800))*73
t_input[:,1] = time

# Grab PRN Values
prns = ephsdf['prn'].to_numpy()

for pp in range(len(prns)):
    health, x = broadcast2pos(ephsdf, t_input, prns[pp])
    print(health)
    print(x)
