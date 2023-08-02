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
from utilities import utl
from geopy.distance import geodesic


# Satellite data files names
GEO_FILE_NAME = 'MYD03.A2016157.0355.061.2018057232219.hdf'
GEO_FILE_NAME_EXTRA = 'MYD03.A2016157.0400.061.2018057232330.hdf'
MOD06_FILE = 'MYD06_L2.A2016157.0355.061.2018059101934.hdf'
MOD06_FILE_EXTRA = 'MYD06_L2.A2016157.0400.061.2018059102206.hdf'
MATCHING_FILE = "matchedData.csv"
AROUND_FILE = "geolocation.csv"
radarfname = '2016157030402_53752_CS_2B-CLDCLASS-LIDAR_GRANULE_P1_R05_E06_F00.hdf'

[ts_row, ts_col] = [2049, 519]

# Read several data by matching and around the sensor position in 500km
matched_geo = pd.read_csv("matchedData.csv", header=None).to_numpy()
matched_len = len(matched_geo)
matched_all_geo = pd.read_csv("matched_all_geo.csv", header=None).to_numpy()
matched_all_len = len(matched_all_geo)

outside_geo = pd.read_csv("outsideGeo.csv", header=None).to_numpy()
ts_around_geo = pd.read_csv("geolocation.csv", header=None).to_numpy()

# Read original satellite data
# Read geolocation data from MYD03 and CloudSat
latitude, longitude = utl.getGeolocationInfo(GEO_FILE_NAME)

latitude_extra, longitude_extra = utl.getGeolocationInfo(GEO_FILE_NAME_EXTRA)

total_latitude = np.vstack((latitude, latitude_extra))
total_longtitude = np.vstack((longitude, longitude_extra))
[total_rows, total_cols] = total_latitude.shape
total_latLng = utl.getTotalLatLng(total_latitude, total_longtitude)

f = reader(radarfname)
lon, lat, elv = f.read_geo()
[radar_nums] = lon.shape
radar_geo = utl.getRadarGeo(lat, lon)

# Read SDS data
cloudLayerBase = f.read_sds('CloudLayerBase', process=True)
cloudLayerType = f.read_sds('CloudLayerType', process=True)

[rows, cols] = cloudLayerBase.shape
#cloudLayerBase_data = pd.DataFrame(cloudLayerBase)
#cloudLayerBase_data.to_csv("cloudLayerBase.csv", header=None, index=False)

cloud_bottom_height_list = []

total_cot = np.vstack((utl.readModis(MOD06_FILE, 'Cloud_Optical_Thickness'), utl.readModis(MOD06_FILE_EXTRA, 'Cloud_Optical_Thickness')))
total_ctp = np.vstack((utl.readModis(MOD06_FILE, 'cloud_top_pressure_1km'), utl.readModis(MOD06_FILE_EXTRA, 'cloud_top_pressure_1km')))
total_cth = np.vstack((utl.readModis(MOD06_FILE, 'cloud_top_height_1km'), utl.readModis(MOD06_FILE_EXTRA, 'cloud_top_height_1km')))

ts_length = len(ts_around_geo)

ts_around_cot = utl.getAroundSDS(total_cot, ts_around_geo)
ts_around_ctp = utl.getAroundSDS(total_ctp, ts_around_geo)

[tsLat, tsLon] = [total_latitude[ts_row][ts_col], total_longtitude[ts_row][ts_col]]

# matched_index_list = utl.getMatchedIndex(matched_geo, total_latLng)
# matched_index = np.array(matched_index_list)

# ts_around_ctypes = utl.getAroundCloudType(ts_around_geo, outside_geo, total_ctp, total_cot)
matched_cbh = utl.getMatchedCLoudBaseLayers(matched_all_geo, cloudLayerBase)
matched_cbh_data = pd.DataFrame(matched_cbh)
matched_cbh_data.to_csv("matched_CloudBaseHeight.csv", header=None, index=False)

# total_ctypes = utl.getAllCLoudTypes(total_ctp, total_cot)
total_ctypes = pd.read_csv("total_cloud_types.csv", header=None).to_numpy()
matched_all_ctypes = utl.getMatchedCloudTypes(matched_all_geo, total_ctypes, cloudLayerType)

# 纬度： 38.43033981323242 经度： 141.0393829345703
# pd.DataFrame(matched_all_ctypes).to_csv("matched_all_ctypes.csv", header=None, index=False)
# ts_around_ctypes = utl.getAroundCloudType(ts_around_geo, total_ctp, total_cot)
# ts_around_cbh = utl.cloudBottomHeightRetrieval(matched_all_geo, matched_all_ctypes, ts_around_geo, ts_around_ctypes, radar_geo, cloudLayerBase)
# cbh_data = pd.DataFrame(ts_around_cbh)
# cbh_data.to_csv("retrieved_cloudBottomHeight.csv", header=None, index=False)
# matched_all_geo = utl.getMatchedIndex(radar_geo, total_latLng)
# matched_all_geo_data = pd.DataFrame(matched_all_geo)
# matched_all_geo_data.to_csv("matched_all_geo.csv", header=None, index=False)
'''
for cbh in np.nditer(cloudLayerBase):
    if cbh is np.ma.masked:
        cbh = np.nan
'''
