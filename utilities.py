import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC
from pyhdf.VS import *
from pyhdf.HDF import *
from geopy.distance import geodesic


class utl:
    def readModis(filename, parameter):
        hdf = SD(filename, SDC.READ)

        raw_data = hdf.select(parameter)
        data = raw_data.get().astype(np.double)
        attrs = raw_data.attributes(full=1)
        aoa = attrs["add_offset"]
        add_offset = aoa[0]
        fva = attrs["_FillValue"]
        _FillValue = fva[0]
        sfa = attrs["scale_factor"]
        scale_factor = sfa[0]
        vra = attrs["valid_range"]
        valid_min = vra[0][0]
        valid_max = vra[0][1]

        invalid = np.logical_or(data > valid_max, data < valid_min)
        invalid = np.logical_or(invalid, data == _FillValue)
        data[invalid] = np.nan
        data = (data - add_offset) * scale_factor
        data = np.ma.masked_array(data, np.isnan(data))

        return data


    def get5kmGeoInfo(geofile):
        hdf_geo = SD(geofile, SDC.READ)
        # Read geolocation dataset from MOD03 product.
        lat = hdf_geo.select('Latitude')
        lat = lat.get().astype(np.double)
        lon = hdf_geo.select('Longitude')
        lon = lon.get().astype(np.double)

        return lat, lon


    def getGeolocationInfo(geofile):
        hdf_geo = SD(geofile, SDC.READ)
        # Read geolocation dataset from MOD03 product.
        lat = hdf_geo.select('Latitude')
        latitude = lat.get().astype(np.double)
        lon = hdf_geo.select('Longitude')
        longitude = lon.get().astype(np.double)

        return latitude, longitude


    def getModisUnit(filename, parameter):
        hdf = SD(filename, SDC.READ)

        raw_data = hdf.select(parameter)
        data = raw_data.get().astype(np.double)
        attrs = raw_data.attributes(full=1)
        ua=attrs["units"]
        units = ua[0]

        return units
    

    def getAroundSDS(total_sds, around_geo):
        length = len(around_geo)
        around_sds = np.zeros((length, 1))
        
        for i in range(length):
            j = int(around_geo[i][2])
            k = int(around_geo[i][3])
            around_sds[i] = total_sds[j][k]
        
        return around_sds
    

    def getAroundSDSInAll(around_sds, around_geo, outside_geo, total_sds):
        [rows, cols] = total_sds.shape
        around_sds_all = np.zeros((rows, cols))

        for i in range(len(around_geo)):
            j = int(around_geo[i][2])
            k = int(around_geo[i][3])
            around_sds_all[j][k] = around_sds[i]
        
        for i in range(len(outside_geo)):
            j = outside_geo[i][0]
            k = outside_geo[i][1]
            around_sds_all[j][k] = np.nan
        
        return around_sds_all

    

    def getMatchedIndex(matched_latLng, total_latLng):
        matched_len = len(matched_latLng)

        update_latLng = np.around(total_latLng, decimals=2)
        update_matched = np.around(matched_latLng, decimals=2)
        latLng_one = np.around(total_latLng, decimals=1)
        matched_one = np.around(matched_latLng, decimals=1)
        matched_index_list = []
        count = 0
        
        for i in range(matched_len):
            matched_lat = matched_latLng[i][0]
            matched_lon = matched_latLng[i][1]
            find = np.where(update_latLng == [update_matched[i][0], update_matched[i][1]])
            find_len = find[0].size
            find_list = []
            comparison_list = []
            for j in range(find_len):
                x = find[0][j]
                y = find[1][j]
                if update_latLng[x][y][0] == update_matched[i][0] and update_latLng[x][y][1] == update_matched[i][1]:
                    # print("匹配到origin的编号, row: ", x, " col: ", y)
                    # print("原本的编号：", i)
                    find_list.append([x, y])
    
            if len(find_list) > 1:
                for m in range(len(find_list)):
                    x = find_list[m][0]
                    y = find_list[m][1]
                    comparison = abs(total_latLng[x][y][0] - matched_lat) + abs(total_latLng[x][y][1] - matched_lon)
                    comparison_list.append(comparison)

                min_index = comparison_list.index(min(comparison_list))
                row_update = find_list[min_index][0]
                col_update = find_list[min_index][1]

                #print("匹配到origin的编号, row: ", row_update, " col: ", col_update)
                #print("原本的编号：", i)
                matched_index_list.append([matched_lat, matched_lon, i, row_update, col_update])
            elif len(find_list) == 1:
                matched_row = find_list[0][0]
                matched_col = find_list[0][1]
                #print("匹配到origin的编号, row: ", matched_row, " col: ", matched_col)
                #print("原本的编号：", i)
                matched_index_list.append([matched_lat, matched_lon, i, matched_row, matched_col])
            else:
                find_near = np.where(latLng_one == [matched_one[i][0], matched_one[i][1]])
                find_near_len = find_near[0].size
                find_near_list = []
                near_comparison = []
                for k in range(find_near_len):
                    row_near = find_near[0][k]
                    col_near = find_near[1][k]
                    if latLng_one[row_near][col_near][0] == matched_one[i][0] and latLng_one[row_near][col_near][1] == matched_one[i][1]:
                        find_near_list.append([row_near, col_near])

                find_near_list_len = len(find_near_list)
                if find_near_list_len > 1:
                    for n in range(find_near_list_len):
                        x = find_near_list[n][0]
                        y = find_near_list[n][1]
                        c = abs(total_latLng[x][y][0] - matched_lat) + abs(total_latLng[x][y][1] - matched_lon)
                        near_comparison.append(c)
                    near_index = near_comparison.index(min(near_comparison))
                    matched_index_list.append([matched_lat, matched_lon, i, find_near_list[near_index][0], find_near_list[near_index][1]])
                    #print("匹配到origin的编号, row: ", find_near_list[near_index][0], "col: ", find_near_list[near_index][1])
                    #print("原本的编号：", i)
                elif find_near_list_len == 1:
                    matched_index_list.append([matched_lat, matched_lon, i, find_near_list[0][0], find_near_list[0][1]])
                    #print("匹配到origin的编号, row: ", find_near_list[0][0], "col: ", find_near_list[0][1])
                    #print("原本的编号：", i)
                else:
                    count += 1
        
        print("无法匹配的点的数量：", count)
        return matched_index_list
    

    def getTotalLatLng(total_latitude, total_longitude):
        [rows, cols] = total_latitude.shape
        total_latLng = np.zeros((rows, cols, 2))

        for i in range(rows):
            for j in range(cols):
                total_latLng[i][j][0] = total_latitude[i][j]
                total_latLng[i][j][1] = total_longitude[i][j]
        
        return total_latLng
    

    def getAroundInTotal(outside, total):
        outside_len = len(outside)

        for i in range(outside_len):
            j = outside[i][0]
            k = outside[i][1]
            total[j][k] = np.nan
        
        return total
    
    
    def getAroundCloudTypeInAll(around, outside, around_ctp, around_cot):
        [rows, cols] = around_cot.shape
        around_cloud_type = utl.getAroundInTotal(outside, np.zeros((rows, cols)))
        around_len = len(around)
        [countCu, countSc, countSt, countAc, countAs, countNs, countCi, countCs, countDc] = [0, 0, 0, 0, 0, 0, 0, 0, 0]

        for i in range(around_len):
            j = int(around[i][2])
            k = int(around[i][3])
            top_pressure = around_ctp[j][k]
            optical_thickness = around_cot[j][k]

            if np.isnan(optical_thickness) or np.isnan(top_pressure):
                around_cloud_type[j][k] = np.nan

            if 0 < optical_thickness <= 3.6:
                if 50 <= top_pressure < 440:
                    around_cloud_type[j][k] = 6
                    countCi += 1
                elif 440 <= top_pressure < 680:
                    around_cloud_type[j][k] = 3
                    countAc += 1
                elif 680 <= top_pressure <= 1000:
                    around_cloud_type[j][k] = 0
                    countCu += 1
            elif 3.6 < optical_thickness <= 23:
                if 50 <= top_pressure < 440:
                    around_cloud_type[j][k] = 7
                    countCs += 1
                elif 440 <= top_pressure < 680:
                    around_cloud_type[j][k] = 4
                    countAs += 1
                elif 680 <= top_pressure <= 1000:
                    around_cloud_type[j][k] = 1
                    countSc += 1
            elif 23 < optical_thickness <= 379:
                if 50 <= top_pressure < 440:
                    around_cloud_type[j][k] = 8
                    countDc += 1
                elif 440 <= top_pressure < 680:
                    around_cloud_type[j][k] = 5
                    countNs += 1
                elif 680 <= top_pressure <= 1000:
                    around_cloud_type[j][k] = 2
                    countSt += 1
        
        print("500km范围内的云类分布")
        print("积云：", countCu)
        print("层积云：", countSc)
        print("层云", countSt)
        print("高积云：", countAc)
        print("高层云：", countAs)
        print("雨层云：", countNs)
        print("卷云：", countCi)
        print("卷层云：", countCs)
        print("深对流：", countDc)
        print()

        return around_cloud_type
    

    def getAroundCloudType(around_geo, total_ctp, total_cot):
        around_len = len(around_geo)
        around_cloud_type = np.zeros(around_len)

        for i in range(around_len):
            j = int(around_geo[i][2])
            k = int(around_geo[i][3])
            top_pressure = total_ctp[j][k]
            optical_thickness = total_cot[j][k]

            if np.isnan(optical_thickness) or np.isnan(top_pressure):
                around_cloud_type[i] = np.nan

            if 0 < optical_thickness <= 3.6:
                if 50 <= top_pressure < 440:
                    around_cloud_type[i] = 6
                elif 440 <= top_pressure < 680:
                    around_cloud_type[i] = 3
                elif 680 <= top_pressure <= 1000:
                    around_cloud_type[i] = 0
            elif 3.6 < optical_thickness <= 23:
                if 50 <= top_pressure < 440:
                    around_cloud_type[i] = 7
                elif 440 <= top_pressure < 680:
                    around_cloud_type[i] = 4
                elif 680 <= top_pressure <= 1000:
                    around_cloud_type[i] = 1
            elif 23 < optical_thickness <= 379:
                if 50 <= top_pressure < 440:
                    around_cloud_type[i] = 8
                elif 440 <= top_pressure < 680:
                    around_cloud_type[i] = 5
                elif 680 <= top_pressure <= 1000:
                    around_cloud_type[i] = 2
        
        return around_cloud_type
    

    def getMatchedCloudTypes(matched_geo, total_ctyps, cloudLayerType):
        matched_len = len(matched_geo)
        matched_ctypes = np.zeros(matched_len)
        [countCu, countSc, countSt, countAc, countAs, countNs, countCi, countCs, countDc] = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        non_count = 0
        ctype = 0
        
        for i in range(matched_len):
            j = int(matched_geo[i][3])
            k = int(matched_geo[i][4])
            l = int(matched_geo[i][2])

            if np.isnan(total_ctyps[j][k]):
                if np.ma.is_masked(cloudLayerType[l].all()):
                    matched_ctypes[i] = np.nan
                    non_count += 1
                else:
                    exist_layer = []
                    for m in range(cloudLayerType[0].size):
                        if np.ma.is_masked(cloudLayerType[l][m]) == False:
                            exist_layer.append(m)
                    max_layer = exist_layer.index(max(exist_layer))
                    ctype = cloudLayerType[l][max_layer].astype(int)
            else:
                matched_ctypes[i] = total_ctyps[j][k]
                ctype = total_ctyps[j][k].astype(int)
            
            if ctype == 2:
                countCu += 1
            elif ctype == 1:
                countSc += 1
            elif ctype == 0:
                countSt += 1
            elif ctype == 5:
                countAc += 1
            elif ctype == 4:
                countAs += 1
            elif ctype == 3:
                countNs += 1
            elif ctype == 8:
                countCi += 1
            elif ctype == 7:
                countCs += 1
            elif ctype == 6:
                countDc += 1

        print("匹配点的云类分布")
        print("积云：", countCu)
        print("层积云：", countSc)
        print("层云", countSt)
        print("高积云：", countAc)
        print("高层云：", countAs)
        print("雨层云：", countNs)
        print("卷云：", countCi)
        print("卷层云：", countCs)
        print("深对流：", countDc)
        print("无云类型：", non_count)

        return matched_ctypes


    def getMatchedCLoudBaseLayers(matched_geo, cloudLayerBase):
        matched_len = len(matched_geo)
        cbh_layers = cloudLayerBase[0].size
        matched_cbh = np.zeros((matched_len, cbh_layers))

        for i in range(matched_len):
            k = int(matched_geo[i][2])
            for j in range(cbh_layers):
                matched_cbh[i][j] = cloudLayerBase[k][j]
        
        return matched_cbh
    

    def getAllCLoudTypes(total_ctp, total_cot):
        [rows, cols] = total_cot.shape
        total_ctyps = np.zeros((rows, cols))
        [countCu, countSc, countSt, countAc, countAs, countNs, countCi, countCs, countDc] = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        non_count = 0

        for i in range(rows):
            for j in range(cols):
                top_pressure = total_ctp[i][j]
                optical_thickness = total_cot[i][j]

                if np.isnan(optical_thickness) or np.isnan(top_pressure):
                    total_ctyps[i][j] = np.nan
                    non_count += 1

                if 0 < optical_thickness <= 3.6:
                    if 50 <= top_pressure < 440:
                        total_ctyps[i][j] = 8
                        countCi += 1
                    elif 440 <= top_pressure < 680:
                        total_ctyps[i][j] = 5
                        countAc += 1
                    elif 680 <= top_pressure <= 1000:
                        total_ctyps[i][j] = 2
                        countCu += 1
                elif 3.6 < optical_thickness <= 23:
                    if 50 <= top_pressure < 440:
                        total_ctyps[i][j] = 7
                        countCs += 1
                    elif 440 <= top_pressure < 680:
                        total_ctyps[i][j] = 4
                        countAs += 1
                    elif 680 <= top_pressure <= 1000:
                        total_ctyps[i][j] = 1
                        countSc += 1
                elif 23 < optical_thickness <= 379:
                    if 50 <= top_pressure < 440:
                        total_ctyps[i][j] = 6
                        countDc += 1
                    elif 440 <= top_pressure < 680:
                        total_ctyps[i][j] = 3
                        countNs += 1
                    elif 680 <= top_pressure <= 1000:
                        total_ctyps[i][j] = 0
                        countSt += 1
        
        print("500km范围内的云类分布")
        print("积云：", countCu)
        print("层积云：", countSc)
        print("层云", countSt)
        print("高积云：", countAc)
        print("高层云：", countAs)
        print("雨层云：", countNs)
        print("卷云：", countCi)
        print("卷层云：", countCs)
        print("深对流：", countDc)
        print("无法匹配的：", non_count)

        return total_ctyps
    

    def getRadarGeo(radar_lat, radar_lon):
        radar_nums = len(radar_lat)
        radar_geo = np.zeros((radar_nums, 2))

        for i in range(radar_nums):
            radar_geo[i][0] = radar_lat[i]
            radar_geo[i][1] = radar_lon[i]
        
        return radar_geo
    

    def weighthedDistance(dist):
        w = (7.8899e-9) * pow(dist, 3) - (9.6364e-6) * pow(dist, 2) + 0.0047 * dist + 1.1649
        w = 1 / pow(w, 2)
        return w
    

    def cloudBottomHeightRetrieval(matched_all_geo, matched_all_ctypes, around_geo, around_ctypes, around_cth, radar_geo, cloudLayerBase):
        around_len = len(around_geo)
        around_cbh = np.zeros(around_len)
        matched_st = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 0)
        matched_sc = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 1)
        matched_cu = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 2)
        matched_ns = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 3)
        matched_as = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 4)
        matched_ac = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 5)
        matched_dc = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 6)
        matched_cs = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 7)
        matched_ci = utl.getOneMatchedCLoudType(matched_all_ctypes, matched_all_geo, 8)


        for i in range(around_len):
            around_ctype = around_ctypes[i]
            around_lat, around_lon = around_geo[i][0], around_geo[i][1]
            cth = around_cth[i]

            if np.isnan(around_ctype):
                around_cbh[i] = np.nan
            else:
                around_ctype = around_ctype.astype(int)

                if around_ctype == 0:
                    matched_cloud = matched_st
                elif around_ctype == 1:
                    matched_cloud = matched_sc
                elif around_ctype == 2:
                    matched_cloud = matched_cu
                elif around_ctype == 3:
                    matched_cloud = matched_ns
                elif around_ctype == 4:
                    matched_cloud = matched_as
                elif around_ctype == 5:
                    matched_cloud = matched_ac
                elif around_ctype == 6:
                    matched_cloud = matched_dc
                elif around_ctype == 7:
                    matched_cloud = matched_cs
                else:
                    matched_cloud = matched_ci

                up = 0.0
                down = 0.0

                for j in range(len(matched_cloud)):
                    cloud_index = matched_cloud[j][1].astype(int)
                    max_layer = utl.getMaxBaseLayer(cloudLayerBase[cloud_index])
                    one_layer_base = cloudLayerBase[cloud_index]
                    ref_cbh = 0.0

                    if  np.isnan(max_layer):
                        next
                    else:
                        if around_ctype == 8:
                            ref_cbh = cloudLayerBase[cloud_index][max_layer]
                        else:
                            if cloudLayerBase[cloud_index][max_layer] < cth:
                                ref_cbh = cloudLayerBase[cloud_index][max_layer]
                            else:
                                one_layer_base[max_layer] = np.ma.masked
                                max_layer2 = utl.getMaxBaseLayer(one_layer_base)
                                if np.isnan(max_layer2):
                                    ref_cbh = cloudLayerBase[cloud_index][max_layer]
                                else:
                                    ref_cbh = cloudLayerBase[cloud_index][max_layer2]

                    dist = geodesic((around_lat, around_lon), (radar_geo[cloud_index][0], radar_geo[cloud_index][1])).km
                    weighted_dist = utl.weighthedDistance(dist)
                    up += weighted_dist * ref_cbh
                    down += weighted_dist
                
                if down == 0:
                    around_cbh[i] = np.nan
                else:
                    around_cbh[i] = (up / down) * 1000
        
        return around_cbh
    

    def getOneMatchedCLoudType(matched_ctypes, matched_geo, ctype):
        matched_x = []
        
        for i in range(len(matched_geo)):
            if np.isnan(matched_ctypes[i]):
                next
            else:
                if matched_ctypes[i].astype(int) == ctype:
                    matched_x.append([matched_ctypes[i], matched_geo[i][2]])

        matched_x = np.array(matched_x)
        return matched_x 


    def getMaxBaseLayer(one_base_layers):
        layers = len(one_base_layers)
        exsit_layers = []
        max_index = 0
        
        if np.ma.is_masked(one_base_layers.all()):
            max_index = np.nan
        else:
            for i in range(layers):
                if np.ma.is_masked(one_base_layers[i]):
                    next
                else:
                    exsit_layers.append(i)
            
            max_index = exsit_layers.index(max(exsit_layers))
        
        return exsit_layers[max_index]
    

    def getAroundPhysicalThickness(around_cth, around_cbh):
        around_thickness = np.zeros((len(around_cth), 1))

        for i in range(len(around_cth)):
            if np.isnan(around_cbh[i]) or np.isnan(around_cth[i]):
                around_thickness[i] = np.nan
            elif around_cbh[i] == 0.0 or around_cth[i] == 0.0:
                around_thickness[i] = np.nan
            else:
                around_thickness[i] = around_cth[i] - around_cbh[i]
                if around_thickness[i] < 0:
                    around_thickness[i] = np.nan
                else:
                    if around_thickness[i] >= 9000:
                        print("Cloud Top Height: ", around_cth[i])
                        print("Cloud Bottom Height: ", around_cbh[i])
                        print("-------")

        return around_thickness
    

    def classifyCloudEtages(around_geo, around_cth, around_cbh, around_cpt, around_ctypes):
        low_list = []
        middle_list = []
        high_list = []
        cu_genus_list = []

        for i in range(len(around_geo)):
            j = around_geo[i][2]
            k = around_geo[i][3]
            cth = around_cth[i][0]
            cbh = around_cbh[i][0]
            thickness = around_cpt[i][0]
            ctype = around_ctypes[i][0]

            if np.isnan(thickness):
                next
            else:
                if cbh < 2500:
                    if ctype != 2:
                        low_list.append([cth, cbh, thickness, ctype, j, k])
                    else:
                        cu_genus_list.append([cth, cbh, thickness, ctype, j, k])
                elif 2500 <= cbh < 6000:
                    middle_list.append([cth, cbh, thickness, ctype, j, k])
                elif cbh >= 6000:
                    high_list.append([cth, cbh, thickness, ctype, j, k])
        
        pd.DataFrame(low_list).to_csv("low_etage.csv", header=None, index=False)
        pd.DataFrame(middle_list).to_csv("middle_etage.csv", header=None, index=False)
        pd.DataFrame(high_list).to_csv("high_etage.csv", header=None, index=False)
        pd.DataFrame(cu_genus_list).to_csv("cu_genus_etage.csv", header=None, index=False)
    
    def getValidateInfo(around_geo, around_cth, around_cbh, around_cpt, around_ctypes, total_lwp, total_reff):
        validate_list = []
        for i in range(len(around_geo)):
            j = int(around_geo[i][2])
            k = int(around_geo[i][3])
            cth = around_cth[i][0]
            cbh = around_cbh[i][0]
            thickness = around_cpt[i][0]
            ctype = around_ctypes[i][0]
            lwp = total_lwp[j][k]
            lwc = lwp / thickness
            reff = total_reff[j][k]

            if np.isnan(thickness):
                next
            else:
                if ctype == 8 and lwc >= 0.2:
                    ctype = 9
                    validate_list.append([cth, cbh, thickness, ctype, lwc, reff, j, k])
                else:
                    validate_list.append([cth, cbh, thickness, ctype, lwc, reff, j, k])

        validate_arr = np.array(validate_list)

        return validate_arr
    

    def getHeightIndex(h, z_line):
        down = 0
        up = 0
        for i in range(len(z_line)):
            if i <= len(z_line) - 2:
                if z_line[i] <= h <= z_line[i+1]:
                    down = i
                    up = i + 1
                    break
        
        if abs(h - z_line[down]) < abs(h - z_line[up]):
            index = down
        else:
            index = up

        return index
    

    def getReffIndex(reff, re_line):
        down = 0
        up = 0
        for i in range(len(re_line)):
            if i <= len(re_line) - 2:
                if re_line[i] <= reff <= re_line[i+1]:
                    down = i
                    up = i + 1
                    break
        
        if abs(reff - re_line[down]) < abs(reff - re_line[up]):
            index = down
        else:
            index = up
        
        return index
    

    def getGridPointIndex(validate_index, x, y):
        matched_index = np.where(validate_index == np.array([x, y]))
        existed = False
        index = 0
        length = matched_index[0].size

        for i in range(length):
            index = matched_index[0][i]
            grid_x = validate_index[index][0]
            grid_y = validate_index[index][1]

            if grid_x == x and grid_y == y:
                existed = True
                break
        
        if existed == True:
            return index
        else:
            return np.nan
    

    def updateGridInfo(lwp, cbh, cth, reff, interpolate_index, validate_index, total_lwp, validate_cbh, validate_cth, validate_reff):
        if np.isnan(interpolate_index):
            next
        else:
            interpolate_x = validate_index[interpolate_index][0]
            interpolate_y = validate_index[interpolate_index][1]
            interpolate_lwp = total_lwp[interpolate_x][interpolate_y]
            interpolate_cbh = validate_cbh[interpolate_index][0]
            interpolate_cth = validate_cth[interpolate_index][0]
            interpolate_reff = validate_reff[interpolate_index][0]

            reff = interpolate_reff * (interpolate_lwp / (interpolate_lwp + lwp)) + reff * (lwp / (interpolate_lwp + lwp))

            lwp += interpolate_lwp

            if interpolate_cbh < cbh:
                cbh = interpolate_cbh
            
            if interpolate_cth > cth:
                cth = interpolate_cth
        
        return lwp, cbh, cth, reff


    def updateUsedFlag(validate_flag, interpolate_index):
        if np.isnan(interpolate_index) == False:
            validate_flag[interpolate_index][0] = 1.0
