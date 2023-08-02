import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC
from geopy.distance import geodesic
from read_CloudSat import reader

GEO_FILE_NAME = 'MYD03.A2016157.0355.061.2018057232219.hdf'
GEO_FILE_NAME_EXTRA = 'MYD03.A2016157.0400.061.2018057232330.hdf'

# Read dataset.
hdf_geo = SD(GEO_FILE_NAME, SDC.READ)
hdf_geo_extra = SD(GEO_FILE_NAME_EXTRA, SDC.READ)

# Read geolocation dataset from MOD03 product.
lat = hdf_geo.select('Latitude')
latitude = lat[:, :]
lon = hdf_geo.select('Longitude')
longitude = lon[:, :]

print(latitude.shape)
[rows, cols] = latitude.shape

# Read another geolocation from different time zone
lat_extra = hdf_geo_extra.select('Latitude')
latitude_extra = lat_extra[:, :]
lon_extra = hdf_geo_extra.select('Longitude')
longitude_extra = lon_extra[:, :]

print(latitude_extra.shape)
[rows_next, cols_next] = latitude_extra.shape

latLng = np.zeros((rows, cols, 2))
latLng2 = np.zeros((rows_next, cols_next, 2))
# print(latLng.dtype)
for i in range(rows):
    for j in range(cols):
        latLng[i][j][0] = latitude[i][j]
        latLng[i][j][1] = longitude[i][j]

for i in range(rows_next):
    for j in range(cols_next):
        latLng2[i][j][0] = latitude_extra[i][j]
        latLng2[i][j][1] = longitude_extra[i][j]

update_latLng = np.around(latLng, decimals=2)
update_latLng2 = np.around(latLng2, decimals=2)

total_latLng_origin = np.vstack((latLng, latLng2))
total_latLng = np.vstack((update_latLng, update_latLng2))

print("更新后总经纬度数组的长度：", total_latLng.shape)
findLocation = np.where(total_latLng == [35.68, 139.7])
findLat = findLocation
# findLat = np.where(update_lat == 35.68)

a = []
b = []

length = findLocation[0].size
# print(findLocation)
# print(findLat)
# print(findLng)
# print(latLng[findLocation])
for i in range(length):
    x = findLocation[0][i]
    y = findLocation[1][i]
    if total_latLng[x][y][0] == 35.68 and total_latLng[x][y][1] == 139.77:
        a.append(i)
    elif (total_latLng[x][y][0] - 35.68) <= 0.02 and (total_latLng[x][y][1] - 139.77) <= 0.02:
        b.append(i)
    # print(update_latLng[x][y][0], " ", update_latLng[x][y][1], " No.", i)

appropriate = np.array(a)
less_proper = np.array(b)
print("最适点的坐标集的尺寸: ", appropriate.size)
print("次适点的坐标集的尺寸: ", less_proper.size)

if appropriate.size != 0:
    for i in range(appropriate.size):
        print("经纬度信息（离东京站最近的点）：")
        print("x: ", findLocation[0][appropriate[i]], " y: ", findLocation[1][appropriate[i]])
        Tokyo_row = findLocation[0][appropriate[i]]
        Tokyo_col = findLocation[1][appropriate[i]]
else:
    for i in range(less_proper.size):
        print("经纬度信息（离东京站次近的点）：")
        print("x: ", findLocation[0][less_proper[i]], " y: ", findLocation[1][less_proper[i]])

# print("经纬度信息（离东京站最近的点）：")
# print("x: ", findLocation[0][2164], " y: ", findLocation[1][2164])
print("Latitude: ", total_latLng[Tokyo_row][Tokyo_col][0], " Longitude: ", total_latLng[Tokyo_row][Tokyo_col][1])
print("------")

[total_rows, total_cols, nums] = total_latLng.shape
[tsLat, tsLon] = [total_latLng_origin[Tokyo_row][Tokyo_col][0], total_latLng_origin[Tokyo_row][Tokyo_col][1]]
#print(total_rows, total_cols)
ts_aroud_list = []
outside_list = []

for i in range(total_rows):
    for j in range(total_cols):
        lat_temp = total_latLng_origin[i][j][0]
        lon_temp = total_latLng_origin[i][j][1]
        dist = geodesic((tsLat, tsLon), (lat_temp, lon_temp))
        if 0 <= dist <= 500:
            ts_aroud_list.append([lat_temp, lon_temp, i, j])
        else:
            outside_list.append([i, j])


ts_around_data = pd.DataFrame(ts_aroud_list)
ts_around_data.to_csv("geolocation.csv", header=None, index=False)

outside_data = pd.DataFrame(outside_list)
outside_data.to_csv("outsideGeo.csv", header=None, index=False)

# lat_data = pd.DataFrame(total_latLng_origin[:, :, 0])
# lat_data.to_csv("latitude.csv", header=None, index=False)
# lon_data = pd.DataFrame(total_latLng_origin[:, :, 1])
# lon_data.to_csv("longitude.csv", header=None, index=False)

'''
for i in range(total_rows):
    for j in range(total_cols):
        dist = geodesic((total_latLng_origin[Tokyo_row][Tokyo_col][0], total_latLng_origin[Tokyo_row][Tokyo_col][1]),
                    (total_latLng_origin[i][j][0], total_latLng_origin[i][j][1])).km
        if dist <= 500:
            ts_aroud_list.append([total_latLng_origin[i][j][0], total_latLng_origin[i][j][1], i, j])

#print(type(ts_aroud_list[0]))
with open('output.txt', 'w') as file:
    for item in ts_aroud_list:
        file.write(str(item[0]) + ', ' + str(item[1]) + ', ' + str(item[2]) + ', ' + str(item[3]) + '\n')

fname = '2016157030402_53752_CS_2B-CLDCLASS-LIDAR_GRANULE_P1_R05_E06_F00.hdf'
f = reader(fname)
lon, lat, elv = f.read_geo()
[radar_nums] = lon.shape
#print(around_rows)

matching_list = []

for i in range(total_rows):
    dist = geodesic((total_latLng[Tokyo_row][Tokyo_col][0], total_latLng[Tokyo_row][Tokyo_col][1]),
                    (total_latLng[i][Tokyo_col][0], total_latLng[i][Tokyo_col][1])).km
        
    if 500 <= dist <= 501:
            print("Distance: ", dist)
            print("Latitude: ", total_latLng[i][Tokyo_col][0], " Longitude: ", total_latLng[i][Tokyo_col][1])
            print("x: ", i, " y: ", Tokyo_col)
            print("------")
            ts_aroud_list.append(i)    
        

for i in range(total_cols):
    dist = geodesic((total_latLng[Tokyo_row][Tokyo_col][0], total_latLng[Tokyo_row][Tokyo_col][1]),
                    (total_latLng[Tokyo_row][i][0], total_latLng[Tokyo_row][i][1])).km
        
    if 500 <= dist <= 501.5:
            print("Distance: ", dist)
            print("Latitude: ", total_latLng[Tokyo_row][i][0], " Longitude: ", total_latLng[Tokyo_row][i][1])
            print("x: ", Tokyo_row, " y: ", i)
            print("------")
            ts_aroud_list.append(i)

print(ts_aroud_list)

for i in range(rows):
    dist = geodesic((latitude[Tokyo_row][Tokyo_col], longitude[Tokyo_row][Tokyo_col]),
                    (latitude[i][Tokyo_col], longitude[i][Tokyo_col])).km
    if 500 <= dist <= 501:
        print("Distance: ", dist)
        print("The boundary point's index: ", i, " and ", Tokyo_col)
        print("The boundary point's latitude: ", latitude[i][Tokyo_col], " longitude: ", longitude[i][Tokyo_col])

print("------")
for i in range(rows_next):
    dist = geodesic((latitude_extra[i][Tokyo_col], longitude_extra[i][Tokyo_col]), (35.682762, 139.76604)).km
    if 500 <= dist <= 501:
        print("Find the boundary in the extra geolocation")
        print("Distance: ", dist)
        print("The boundary point's index: ", i, " and ", Tokyo_col)
        print("The boundary point's latitude: ", latitude_extra[i][Tokyo_col], " longitude: ",
              longitude_extra[i][Tokyo_col])

print("------")
for i in range(cols):
    dist = geodesic((latitude[Tokyo_row][Tokyo_col], longitude[Tokyo_row][Tokyo_col]),
                    (latitude[Tokyo_row][i], longitude[Tokyo_row][i])).km
    if 500 <= dist <= 502:
        print("Distance: ", dist)
        print("The boundary point's index: ", Tokyo_row, " and ", i)
        print("The boundary point's latitude: ", latitude[Tokyo_row][i], " longitude: ", longitude[Tokyo_row][i])

print("------")
for i in range(cols_next):
    dist = geodesic((latitude_extra[Tokyo_row][i], longitude_extra[Tokyo_row][i]), (35.682762, 139.76604)).km
    if 500 <= dist <= 501:
        print("Find the boundary in the extra geolocation")
        print("Distance: ", dist)
        print("The boundary point's index: ", Tokyo_row, " and ", i)
        print("The boundary point's latitude: ", latitude_extra[Tokyo_row][i], " longitude: ",
              longitude_extra[Tokyo_row][i])

'''
