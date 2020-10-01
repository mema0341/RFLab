import numpy as np

def min2meters(min):
    return min*1852.0

def meters2mins(meters):
    return meters/1852.0

def min2degree(min):
    return min/60

def degree2min(degree):
    return degree*60

def degree2meter(degree):
    return min2meters(degree2min(degree))

def meter2degree(meter):
    return min2degree(meters2mins(meter))

def meter2sec(meter):
    return min2sec(meters2mins(meter))

def sec2meter(sec):
    return min2meters(sec2min(sec))

def sec2min(sec):
    return sec/60

def min2sec(min):
    return min*60

mymeter = sec2meter(60)
print(mymeter,"meters")
print(meter2sec(round(mymeter,2)), "seconds")
# meter2degree(
# def meters2degrees(meters):
#     return meters*0.0000089992803

# def degrees2meters(degrees):
#     return degrees/0.00000899928,2

# def seconds2meters(seconds):
#     return seconds*1582/60

# num_meters = 111120
# # print(degrees2meters(1),"degrees")

# print(seconds2meters(60))

