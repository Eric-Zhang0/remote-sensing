import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.basemap import Basemap
from utilities import utl
from read_CloudSat import reader
from geopy.distance import geodesic


GEO_FILE_NAME = 'MYD03.A2016157.0355.061.2018057232219.hdf'
GEO_FILE_NAME_EXTRA = 'MYD03.A2016157.0400.061.2018057232330.hdf'
MOD06_FILE = 'MYD06_L2.A2016157.0355.061.2018059101934.hdf'
MOD06_FILE_EXTRA = 'MYD06_L2.A2016157.0400.061.2018059102206.hdf'
radarfname = '2016157030402_53752_CS_2B-CLDCLASS-LIDAR_GRANULE_P1_R05_E06_F00.hdf'
[ts_row, ts_col] = [2049, 519]

# Read MODIS geolocation data
latitude, longitude = utl.getGeolocationInfo(GEO_FILE_NAME)

latitude_extra, longitude_extra = utl.getGeolocationInfo(GEO_FILE_NAME_EXTRA)

total_latitude = np.vstack((latitude, latitude_extra))

total_longtitude = np.vstack((longitude, longitude_extra))

ts_around_geo = pd.read_csv("geolocation.csv", header=None).to_numpy()
outside_geo = pd.read_csv("outsideGeo.csv", header=None).to_numpy()
matched_all_geo = pd.read_csv("matched_all_geo.csv", header=None).to_numpy()
matched_all_ctypes = pd.read_csv("matched_all_ctypes.csv", header=None).to_numpy()
matching_list = pd.read_csv("matchedData.csv", header=None).to_numpy()
total_ctypes = pd.read_csv("total_cloud_types.csv", header=None).to_numpy()

cloud_bottom_height = pd.read_csv("retrieved_cloudBottomHeight.csv", header=None).to_numpy()
cloud_bottom_height *= 1000

total_cth = np.vstack((utl.readModis(MOD06_FILE, 'cloud_top_height_1km'), utl.readModis(MOD06_FILE_EXTRA, 'cloud_top_height_1km')))
ts_around_cth = utl.getAroundSDS(total_cth, ts_around_geo)

f = reader(radarfname)
lon, lat, elv = f.read_geo()
[radar_nums] = lon.shape
radar_geo = utl.getRadarGeo(lat, lon)

# Read SDS data
cloudLayerBase = f.read_sds('CloudLayerBase', process=True)
cloudLayerType = f.read_sds('CloudLayerType', process=True)

comparison = np.hstack((cloud_bottom_height, ts_around_cth))

'''
for i in range(len(ts_around_geo)):
    cbh = comparison[i][0]
    cth = comparison[i][1]
    j = int(ts_around_geo[i][2])
    k = int(ts_around_geo[i][3])

    if cbh == 0:
        if cth >= 6000:
            ctype = total_ctypes[j][k]
            if np.isnan(ctype):
                comparison[i][0] = np.nan
                comparison[i][1] = np.nan
            else:
                ctype = ctype.astype(int)
                matched_cloud = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 5)

                up = 0.0
                down = 0.0
                for l in range(len(matched_cloud)):
                    cloud_index = matched_cloud[l][1].astype(int)
                    max_layer = utl.getMaxBaseLayer(cloudLayerBase[cloud_index])
                    one_layer_base = cloudLayerBase[cloud_index]
                    ref_cbh = 0.0
                    if cloudLayerBase[cloud_index][max_layer] < cth:
                        ref_cbh = cloudLayerBase[cloud_index][max_layer]
                    else:
                        one_layer_base[max_layer] = np.ma.masked
                        max_layer2 = utl.getMaxBaseLayer(one_layer_base)
                        if np.isnan(max_layer2):
                            next
                        else:
                            ref_cbh = cloudLayerBase[cloud_index][max_layer2]
                    dist = geodesic((ts_around_geo[i][0], ts_around_geo[i][1]), (lat[cloud_index], lon[cloud_index])).km
                    weighted_dist = utl.weighthedDistance(dist)
                    up += weighted_dist * ref_cbh
                    down += weighted_dist


                comparison[i][0] = (up/down) * 1000
                print("云类型：", ctype)
                print("云顶高度：", cth)
                print("反演的云底高度：", up/down)
        else:
            comparison[i][0] = np.nan
            comparison[i][1] = np.nan
    elif np.isnan(cth):
        comparison[i][0] = np.nan
        comparison[i][1] = np.nan
    elif cbh >= cth:
        comparison[i][0] = np.nan
        comparison[i][1] = np.nan
'''


#pd.DataFrame(comparison[:, 0]).to_csv("update_around_cbh.csv", header=None, index=False)
#pd.DataFrame(comparison[:, 1]).to_csv("update_around_cth.csv", header=None, index=False)
update_around_cth = pd.read_csv("update_around_cth.csv", header=None).to_numpy()
update_around_cbh = pd.read_csv("update_around_cbh.csv", header=None).to_numpy()

ts_all_cbh = utl.getAroundSDSInAll(update_around_cbh, ts_around_geo, outside_geo, total_cth)
ts_all_cth = utl.getAroundSDSInAll(update_around_cth, ts_around_geo, outside_geo, total_cth)
[tsLat, tsLon] = [total_latitude[ts_row][ts_col], total_longtitude[ts_row][ts_col]]
ts_all_thickness = utl.getAroundSDSInAll(utl.getAroundPhysicalThickness(update_around_cbh, update_around_cth), ts_around_geo, outside_geo, total_cth)

pd.DataFrame(utl.getAroundPhysicalThickness(update_around_cbh, update_around_cth)).to_csv("update_around_cpt.csv", header=None, index=False)

m = Basemap(projection='laea', resolution='c', lat_0=32, lon_0=136, llcrnrlat=20, llcrnrlon=118, urcrnrlat=50, urcrnrlon=156, width=3000000, height=2500000)
m.drawcoastlines(linewidth=0.5)
m.drawparallels(np.arange(16., 49., 3.0), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(120., 159., 4.0), labels=[0, 0, 0, 1])
m.pcolormesh(total_longtitude, total_latitude, ts_all_thickness, latlon=True, cmap=plt.cm.RdBu_r)

tick_locator = ticker.MaxNLocator(nbins=9)
cb = m.colorbar()
cb.locator = tick_locator
cb.set_ticks([0, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000])
#cb.set_label("Cloud Bottom Height(m)", fontsize=8)
#cb.set_label("Cloud Top Height(m)", fontsize=8)
cb.set_label("Cloud Physical Thickness", fontsize=8)
# cb.ax.set_yticklabels(['Cu', 'Sc', 'St', 'Ac', 'As', 'Ns', 'Ci', 'Cs', 'Dc'], fontsize=10)

#plt.title("Different Cloud Bottom Height with 500km around the Tokyo Station")
#plt.title("Different Cloud Top Height with 500km around the Tokyo Station")
plt.title("Different Cloud Physical Thickness with 500km around the Tokyo Station")

xpt, ypt = m(tsLon, tsLat)
m.plot(xpt, ypt, 'c*', markersize=5)
#m.plot(matching_list[0][0], matching_list[0][1], 'c*', markersize=5)

'''
for i in range(len(matching_list)):
    radar_xpt, radar_ypt = m(matching_list[i][1], matching_list[i][0])
    m.plot(radar_xpt, radar_ypt, color='r', marker='.', markersize=0.1)
'''

fig = plt.gcf()
pngfile = "ts_around_cpt.png"
fig.savefig(pngfile)

# comparison_data = pd.DataFrame(comparison)
# comparison_data.to_csv("cbh_and_cth.csv", header=None, index=False)

