from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import sin, cos, sqrt, pi

def to_hammer_xy(lat, lon):
    lat_rad = np.abs(np.multiply(lat, pi/180))
    lon_rad = np.multiply(lon, pi/180)
    denominator = np.sqrt(1+np.cos(lat_rad)*np.cos(np.divide(lon_rad, 2)))
    x = np.divide(2*np.sqrt(2)*np.cos(lat_rad)*np.sin(np.divide(lon_rad, 2)), denominator)
    y = np.divide(np.sqrt(2)*np.sin(lat_rad), denominator)
    return np.round(np.multiply(x, 180/pi)), np.round(np.multiply(y, 180/pi))

def to_hammer_z(values):
    vals_shape = values.shape
    new_values = np.empty((90,360))
    new_values[:] = np.nan
    for i in range(81):
        lat_rad = abs(i-90)*pi/180
        for j in range(361):
            lon_rad = (j-180)*pi/180

            # project from lat/lon to x/y in Hammer
            denominator = sqrt(1+cos(lat_rad)*cos(lon_rad/2))
            x = (2*sqrt(2)*cos(lat_rad)*sin(lon_rad/2))/denominator
            y = (sqrt(2)*sin(lat_rad))/denominator

            # convert back from rads to degrees
            x = int(round(x*180/pi)) + 180
            y = int(abs(round(y*180/pi) - 90))

            # set value in new array in appropriate place
            new_values[y,x] = values[i,j]
    return new_values


lat = np.abs(np.subtract(zen, 90))
vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
        39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52]
f = interpolate.interp2d(azi, zen, vals, kind='cubic')
azinew = np.arange(-180, 181, 1)
zennew = np.arange(0, 81, 1)
znew = f(azinew, zennew)
x_hammer, y_hammer = to_hammer_xy(lat, azi)
y_hammer = np.subtract(np.abs(np.subtract(y_hammer, 90)), 10)
z_hammer = to_hammer_z(znew)

# fill in nan values within hemisphere if any exist
for row in range(z_hammer.shape[0]):
    ind = np.where(~np.isnan(z_hammer[row]))[0]
    if ind.size != 0:
        first, last = ind[0], ind[-1]
        mask = np.isnan(z_hammer[row,:])
        mask[:first] = 0
        mask[last:] = 0
        z_hammer[row, mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), z_hammer[row, ~mask])

fig = plt.figure(1, figsize=(10, 6), dpi=100)
ax = fig.add_subplot(111)
ax.set_facecolor('black')
cmap = cm.jet
cmap.set_bad('black',1)
ax.imshow(z_hammer, extent=(-180, 180, 80, 0), interpolation='nearest', cmap=cmap)
ax.scatter(x=x_hammer, y=y_hammer, c='r', s=10)
plt.show()

