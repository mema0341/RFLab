import numpy as np
import matplotlib.pyplot as plt

# Obtain range using aircraft speed and height. 
def get_range(v0, h0, step):
    x = np.arange(-x0, x0+10, 10)
    t = x/v0
    r = np.sqrt(np.power(x,2)+np.power(h0,2))
    return x, t, r

# Obtain range rate range and time to get time derivative of range [m/s] 
def get_range_rate(x, t, r):
    r_rate = []
    for tt in range(1,len(x)):
        r_rate.append( (r[tt] - r[tt-1]) / (t[tt] - t[tt-1]) )
    return r_rate

# Obtain range rate using aircraft height and range [degrees]
def get_zenith(h0, r):
    z_angle = np.rad2deg(np.arccos(h0/r))
    return z_angle

# Main Function
if __name__ == "__main__":

    # Define variables
    v0 = 50 # speed [m/s]
    h0 = 100 # aircraft height [m]
    x0 = 250 # horizontal position of aircraft [m]
    step = 0.1 # Step size for graphing purposes

    # Obtain -x0 to x0, time, and range
    x, t, r = get_range(v0, h0, step)

    # Obtain range rate [m/s]
    r_rate = get_range_rate(x, t, r)

    # Obtain zenith angle [degrees]
    z_angle = get_zenith(h0, r)

    # Plot Range
    fig, ax = plt.subplots()
    ax.plot(x, r)
    ax.set_title('Range [meters]\nv0 = %dm/s, x0 = %dm, h0 = %dm' % (int(v0), int(x0), int(h0)))
    ax.set_xlabel('X-Location [meters]')
    ax.set_ylabel('Range [meters]')
    ax.grid(True)
    plt.xlim(min(x),max(x))
    plt.ylim(min(r),max(r))
    plt.grid(True)

    # Plot Range Rate
    fig, ax = plt.subplots()
    ax.plot(x[1:len(x)], r_rate, color='r')
    ax.set_title('Range Rate [meters/second]\nv0 = %dm/s, x0 = %dm, h0 = %dm' % (int(v0), int(x0), int(h0)))
    ax.set_xlabel('X-Location [meters]')
    ax.set_ylabel('Range Rate [meters/second]')
    plt.xlim(min(x),max(x))
    plt.ylim(min(r_rate),max(r_rate))
    plt.grid(True)

    # Plot Zenith Angle
    fig, ax = plt.subplots()
    ax.plot(x, z_angle, color='g')
    ax.set_title('Zenith Angle [Degrees]\nv0 = %dm/s, x0 = %dm, h0 = %dm' % (int(v0), int(x0), int(h0)))
    ax.set_xlabel('X-Location [meters]')
    ax.set_ylabel('Zenith Angle [Degrees]')
    plt.xlim(min(x),max(x))
    plt.ylim(min(z_angle),max(z_angle))
    plt.grid(True)

    # Comparison to Lower Altitude
    h02 = 20 # aircraft 2 height [m]
    x2, t2, r2 = get_range(v0, h02, step) # Obtain -x0 to x0, time, and range
    r_rate2 = get_range_rate(x2, t2, r2) # Obtain range rate [m/s]
    z_angle2 = get_zenith(h0, r) # Obtain zenith angle [degrees
    # Plot Range Rate Again
    fig, ax = plt.subplots()
    ax.plot(x[1:len(x)], r_rate)
    ax.plot(x[1:len(x)], r_rate2)
    ax.set_title('Range Rate Comparison [meters/second]')
    ax.set_xlabel('X-Location [meters]')
    ax.set_ylabel('Range Rate [meters/second]')
    plt.legend(["h0 = %d [m]" % h0, "h02 = %d [m]" % h02])
    plt.grid(True)

    # Plot Range Comparison
    fig, ax = plt.subplots()
    ax.plot(x, r)
    ax.plot(x, r2)
    ax.set_title('Range Comparison [meters]')
    ax.set_xlabel('X-Location [meters]')
    ax.set_ylabel('Range [meters]')
    ax.grid(True)
    plt.legend(["h0 = %d [m]" % h0, "h02 = %d [m]" % h02])
    plt.grid(True)

    plt.show()