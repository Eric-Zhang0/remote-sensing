import numpy as np
import pandas as pd
from utilities import utl

GEO_FILE_NAME = 'MYD03.A2016157.0355.061.2018057232219.hdf'
GEO_FILE_NAME_EXTRA = 'MYD03.A2016157.0400.061.2018057232330.hdf'
MOD06_FILE = 'MYD06_L2.A2016157.0355.061.2018059101934.hdf'
MOD06_FILE_EXTRA = 'MYD06_L2.A2016157.0400.061.2018059102206.hdf'
radarfname = '2016157030402_53752_CS_2B-CLDCLASS-LIDAR_GRANULE_P1_R05_E06_F00.hdf'

ts_around_geo = pd.read_csv("geolocation.csv", header=None).to_numpy()

cloud_bottom_height = pd.read_csv("retrieved_cloudBottomHeight4.csv", header=None).to_numpy()

# total_cth = np.vstack((utl.readModis(MOD06_FILE, 'cloud_top_height_1km'), utl.readModis(MOD06_FILE_EXTRA, 'cloud_top_height_1km')))
# total_cot = np.vstack((utl.readModis(MOD06_FILE, 'Cloud_Optical_Thickness'), utl.readModis(MOD06_FILE_EXTRA, 'Cloud_Optical_Thickness')))
# total_ctp = np.vstack((utl.readModis(MOD06_FILE, 'cloud_top_pressure_1km'), utl.readModis(MOD06_FILE_EXTRA, 'cloud_top_pressure_1km')))
total_lwp = np.vstack((utl.readModis(MOD06_FILE, 'Cloud_Water_Path'), utl.readModis(MOD06_FILE_EXTRA, 'Cloud_Water_Path')))

ts_around_ctypes = pd.read_csv("ts_around_ctype.csv", header=None).to_numpy()
ts_around_cth = pd.read_csv("ts_around_cth.csv", header=None).to_numpy()
ts_around_thickness = utl.getAroundPhysicalThickness(ts_around_cth, cloud_bottom_height)

# utl.classifyCloudEtages(ts_around_geo, ts_around_cth, cloud_bottom_height, ts_around_thickness, ts_around_ctypes)

# Debug for the output csv file
# testIn = pd.read_csv("low_etage.csv", header=None).to_numpy()
# print(type(testIn[0][0]))
etage_files = np.array(["low_etage.csv", "middle_etage.csv", "high_etage.csv", "cu_genus_etage.csv"])
etage_list = []
etage_tf = []

try:
    for i in range(len(etage_files)):
        fname = etage_files[i]
        etage = pd.read_csv(fname, header=None).to_numpy()
        etage_list.append(etage)
        etage_tf.append(True)
    
    '''
    low = pd.read_csv("low_etage.csv", header=None).to_numpy()
    middle = pd.read_csv("middle_etage.csv", header=None).to_numpy()
    high = pd.read_csv("high_etage.csv", header=None).to_numpy()
    cu_genus = pd.read_csv("cu_genus_etage.csv", header=None).to_numpy()
    '''
except pd.errors.EmptyDataError:
    etage_list.append(0)
    etage_tf.append(False)

for i in range(4):
    if etage_tf[i] == True:
        etage = etage_list[i]
        cth = etage[:, 0:1]
        cbh = etage[:, 1:2]
        cpt = etage[:, 2:3]
        
        if i == 0:
            print("低云族：")
        elif i == 1:
            print("中云族：")
        elif i == 2:
            print("高云族：")
        else:
            print("直展云族：")
        
        print("云底高度最低：", np.amin(cbh), "最高：", np.amax(cbh))
        print("云顶高度最低：", np.amin(cth), "最高：", np.amax(cth))
        print("云物理厚度最低：", np.amin(cpt), "最高：", np.amax(cpt))
        print("------")
    

validate_height = utl.getValidateHeight(ts_around_geo, ts_around_cth, cloud_bottom_height, ts_around_thickness, ts_around_ctypes)
validate_cth = validate_height[:, 0:1]
validate_cbh = validate_height[:, 1:2]
validate_cpt = validate_height[:, 2:3]

cirrus_list = []
if etage_tf[2] == True:
    high_etage = etage_list[2]
    for i in range(len(high_etage)):
        ctype = high_etage[i][3]
        cbh = high_etage[i][1]
        cth = high_etage[i][0]
        cpt = high_etage[i][2]
        j = int(high_etage[i][4])
        k = int(high_etage[i][5])
        lwp = total_lwp[j][k]
        lwc = lwp / cpt

        if ctype == 8 and lwc >= 0.2:
            print("卷云云底高度：",cbh)
            print("卷云云顶高度：",cth)
            print("卷云的液体水含量：",lwc)
            cirrus_list.append([cbh, cth, cpt, lwc, j, k])
            high_etage[i][3] = 9
        
    pd.DataFrame(high_etage).to_csv("high_etage.csv", header=None, index=False)

cirrus_arr = np.array(cirrus_list)
pd.DataFrame(cirrus_arr).to_csv("cirrus.csv", header=None, index=False)
