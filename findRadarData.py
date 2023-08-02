import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from pyhdf.SD import SD, SDC
from pydoc import help
from pyhdf.VS import *
from pyhdf.HDF import *
from read_CloudSat import reader
from geopy.distance import geodesic
from utilities import utl


boundaries = [1553, 2546, 163, 989]
GEO_FILE_NAME = 'MYD03.A2016157.0355.061.2018057232219.hdf'
GEO_FILE_NAME_EXTRA = 'MYD03.A2016157.0400.061.2018057232330.hdf'
MOD06_FILE = 'MYD06_L2.A2016157.0355.061.2018059101934.hdf'
MOD06_FILE_EXTRA = 'MYD06_L2.A2016157.0400.061.2018059102206.hdf'

[ts_row, ts_col] = [2049, 519]

# Read Radar Data
fname = '2016157030402_53752_CS_2B-CLDCLASS-LIDAR_GRANULE_P1_R05_E06_F00.hdf'
f = reader(fname)
lon, lat, elv = f.read_geo()
[radar_nums] = lon.shape

# Read MODIS geolocation data
latitude, longitude = utl.getGeolocationInfo(GEO_FILE_NAME)

latitude_extra, longitude_extra = utl.getGeolocationInfo(GEO_FILE_NAME_EXTRA)

total_latitude = np.vstack((latitude, latitude_extra))

total_longtitude = np.vstack((longitude, longitude_extra))

[total_rows, total_cols] = total_latitude.shape

geo_data = pd.read_csv("geolocation.csv", header=None)
ts_around_geo = geo_data.to_numpy()
outside_data = pd.read_csv("outsideGeo.csv", header=None)
outside_geo = outside_data.to_numpy()

matching_list = []

[tsLat, tsLon] = [total_latitude[ts_row][ts_col], total_longtitude[ts_row][ts_col]]

# lat[i] >= ts_lat_min and lat[i] <= ts_lat_max) and (lon[i] >= ts_lon_min and lon[i] <= ts_lon_max
for i in range(radar_nums):
    if lat[i] >= 0 and lon[i] >= 0:
        dist = geodesic((tsLat, tsLon), (lat[i], lon[i])).km
        if 0 <= dist <= 500:
            matching_list.append([lat[i], lon[i], i])
            # print("Distance: ", dist)
            # print("Latitude: ", lat[i])
            # print("Longitude: ", lon[i])

# matching_data = pd.DataFrame(matching_list)
# matching_data.to_csv("matchedData.csv", header=None, index=False)
#print(len(matching_list))
#print(matching_list[0][0], matching_list[0][1])
matched_all_geo = pd.read_csv("matched_all_geo.csv", header=None).to_numpy()

total_cot = np.vstack((utl.readModis(MOD06_FILE, 'Cloud_Optical_Thickness'), utl.readModis(MOD06_FILE_EXTRA, 'Cloud_Optical_Thickness')))
total_ctp = np.vstack((utl.readModis(MOD06_FILE, 'cloud_top_pressure_1km'), utl.readModis(MOD06_FILE_EXTRA, 'cloud_top_pressure_1km')))
ts_around_cot = total_cot
# Update the original cot data
ts_around_ctype = utl.getAroundCloudTypeInAll(ts_around_geo, outside_geo, total_ctp, total_cot)


[TSlat, TSlon] = [total_latitude[ts_row][ts_col], total_longtitude[ts_row][ts_col]]
[lat_m, lon_m] = [np.nanmean(total_latitude), np.nanmean(total_longtitude)]
print(lat_m, lon_m)

'''
# Visualization for around 500km only
m = Basemap(projection='laea', resolution='c', lat_0=lat_m, lon_0=lon_m, llcrnrlat=13, llcrnrlon=113, urcrnrlat=60, urcrnrlon=180, width=3000000, height=2500000)
m.drawcoastlines(linewidth=0.5)
m.drawparallels(np.arange(14., 53., 4.0), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(113., 154., 10.0), labels=[0, 0, 0, 1])
'''

m = Basemap(projection='laea', resolution='c', lat_0=32, lon_0=136, llcrnrlat=20, llcrnrlon=118, urcrnrlat=50, urcrnrlon=156, width=3000000, height=2500000)
m.drawcoastlines(linewidth=0.5)
m.drawparallels(np.arange(16., 49., 3.0), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(120., 159., 4.0), labels=[0, 0, 0, 1])
m.pcolormesh(total_longtitude, total_latitude, ts_around_ctype, latlon=True, cmap=plt.cm.RdBu_r)
# lo, la = np.meshgrid(total_longtitude, total_latitude)
# fig = m.contourf(lo, la, ts_around_ctype, cmap=plt.cm.RdBu_r)
cb = m.colorbar()
cb.set_label("Cloud Types", fontsize=8)
cb.ax.set_yticklabels(['Cu', 'Sc', 'St', 'Ac', 'As', 'Ns', 'Ci', 'Cs', 'Dc'], fontsize=10)

plt.title("Different Cloud Types with 500km around the Tokyo Station")

xpt, ypt = m(TSlon, TSlat)
m.plot(xpt, ypt, 'c*', markersize=5)
#m.plot(matching_list[0][0], matching_list[0][1], 'c*', markersize=5)

for i in range(len(matching_list)):
    radar_xpt, radar_ypt = m(matching_list[i][1], matching_list[i][0])
    m.plot(radar_xpt, radar_ypt, color='r', marker='.', markersize=0.1)

fig = plt.gcf()
pngfile = "test.png"
fig.savefig(pngfile)
