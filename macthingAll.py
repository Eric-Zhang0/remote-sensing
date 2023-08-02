import numpy as np
import pandas as pd
from read_CloudSat import reader
from utilities import utl


GEO_FILE_NAME = 'MYD03.A2016157.0355.061.2018057232219.hdf'
GEO_FILE_NAME_EXTRA = 'MYD03.A2016157.0400.061.2018057232330.hdf'
MOD06_FILE = 'MYD06_L2.A2016157.0355.061.2018059101934.hdf'
MOD06_FILE_EXTRA = 'MYD06_L2.A2016157.0400.061.2018059102206.hdf'
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

matching_all_list = []

total_latLng = np.zeros((total_rows, total_cols, 2))
radar_latLng = np.zeros((radar_nums, 2))

for i in range(total_rows):
    for j in range(total_cols):
        total_latLng[i][j][0] = total_latitude[i][j]
        total_latLng[i][j][1] = total_longtitude[i][j]

for i in range(radar_nums):
    radar_latLng[i][0] = lat[i]
    radar_latLng[i][1] = lon[i]

'''
print(total_latLng.shape)
print(radar_latLng.shape)
for latLon in total_latLng:
    print("纬度：", latLon[0])
    print("经度：", latLon[1])

for latLon in total_latLng:
    for radar_latLon in radar_latLng:
        print(latLon.shape)
        [lat_temp, lon_temp] = np.around(latLon, decimals=2)
        [lat_r_temp, lon_r_temp] = np.around(radar_latLon, decimals=2)
        if [lat_temp, lon_temp] == [lat_r_temp, lon_r_temp]:
            matching_all_list.append(radar_latLon)
            print("匹配成功的经纬度的位置：", np.where(radar_latLng == radar_latLon))

count = 0
for lat in np.nditer(total_latitude):
    print(count)
    count += 1
    if (count == 2):
        break

for i in range(radar_nums):
    radarLat = update_radar_geo[i][0]
    radarLon = update_radar_geo[i][1]
    find_macthing = np.where(update_modis_geo == [radarLat, radarLon])
    find_macthing = np.asarray(find_macthing).reshape(-1)
    if find_macthing.size != 0:
        matching_all_list.append([radar_latLng[i][0], radar_latLng[i][1], find_macthing[0], find_macthing[1]])

matching_all_data = pd.DataFrame(matching_all_list)
matching_all_data.to_csv("matching_all.csv", header=None, index=False)
print("finished")
'''

update_modis_geo = np.around(total_latLng, decimals=2)
update_radar_geo = np.around(radar_latLng, decimals=2)

find = np.where(update_modis_geo == [-23.98, 153.48])
print(find)
print(update_modis_geo[1][1][1])
# -23.98,153.48
