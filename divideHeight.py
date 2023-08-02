import numpy as np
import pandas as pd
from utilities import utl

MOD06_FILE = 'MYD06_L2.A2016157.0355.061.2018059101934.hdf'
MOD06_FILE_EXTRA = 'MYD06_L2.A2016157.0400.061.2018059102206.hdf'
'''
ts_around_geo = pd.read_csv("geolocation.csv", header=None).to_numpy()
x_index = ts_around_geo[:, 2:3]
y_index = ts_around_geo[:, 3:4]

x_index_min = np.amin(x_index)
x_index_max = np.amax(x_index)

y_index_min = np.amin(y_index)
y_index_max = np.amax(y_index)

print(x_index_min, x_index_max)
print(y_index_min, y_index_max)
'''
test_line = np.linspace(5.0, 60.0, num=110, endpoint=True)
# print(test_line)
round_line = np.around(test_line, decimals=4)
print(round_line)
# total_reff = np.vstack((utl.readModis(MOD06_FILE, 'Cloud_Effective_Radius'), utl.readModis(MOD06_FILE_EXTRA, 'Cloud_Effective_Radius')))
# total_surt = np.vstack((utl.readModis(MOD06_FILE, 'surface_temperature_1km'), utl.readModis(MOD06_FILE_EXTRA, 'surface_temperature_1km')))
total_lwp = np.vstack((utl.readModis(MOD06_FILE, 'Cloud_Water_Path'), utl.readModis(MOD06_FILE_EXTRA, 'Cloud_Water_Path')))

# pd.DataFrame(total_surt).to_csv("surface_temperature_1km.csv", header=None, index=False)
ts_around_geo = pd.read_csv("geolocation.csv", header=None).to_numpy()
# ts_around_surft = utl.getAroundSDS(total_surt, ts_around_geo)
# print("地表温度最低：", np.nanmin(ts_around_surft))
# print("地表温度最高：", np.nanmax(ts_around_surft))

validate_info = pd.read_csv("validate_info.csv", header=None).to_numpy()
validate_cth = validate_info[:, 0:1]
validate_cbh = validate_info[:, 1:2]
print("最低的云底高度：", np.amin(validate_cbh))
print("最高的云顶高度：", np.amax(validate_cth))
validate_cpt = validate_info[:, 2:3]
validate_ctype = validate_info[:, 3:4]
validate_lwc = validate_info[:, 4:5]
validate_reff = validate_info[:, 5:6]

validate_x = validate_info[:, 6:7].astype(int)
validate_y = validate_info[:, 7:8].astype(int)
x_min = np.amin(validate_x)
y_min = np.amin(validate_y)
validate_index = np.hstack((validate_x, validate_y)).astype(int)
validate_flag = np.zeros((len(validate_info), 1))

divide_z = np.linspace(0.0, 16.55, num=10, endpoint=True)
divide_z = np.around(divide_z, decimals=4)
# print(divide_z)

divide_t = np.linspace(260.0, 298.0, num=10, endpoint=True)
divide_t = np.around(divide_t, decimals=4)
divide_t = divide_t[::-1]
# print(divide_t)

volume_grids = []
volume_grids.append(["3   propgen particle file"])
volume_grids.append([438, 413, 10])
volume_grids.append([2.0000, 2.0000])
lines = []
z_list = []
t_list = []
z_list.append('')
t_list.append('')

for z in divide_z:
    z_list.append(z)
volume_grids.append(z_list)

for t in divide_t:
    t_list.append(t)
volume_grids.append(t_list)

for i in range(len(validate_info)):
    x = validate_index[i][0]
    y = validate_index[i][1]
    cbh = validate_cbh[i][0]
    cth = validate_cth[i][0]
    cpt = validate_cpt[i][0]
    lwp = total_lwp[x][y]
    reff = validate_reff[i][0]
    used_flag = validate_flag[i][0]

    if used_flag == True:
        next
    else:
        if y % 2 == 1:
            for j in range(3):
                if j == 0:
                    interpolate_index = utl.getGridPointIndex(validate_index, x, y-1)
                    lwp, cbh, cth, reff = utl.updateGridInfo(lwp, cbh, cth, reff, interpolate_index, validate_index, total_lwp, validate_cbh, validate_cth, validate_reff)
                    if np.isnan(interpolate_index) == False:
                        validate_flag[interpolate_index][0] = 1.0
                elif j == 1:
                    interpolate_index = utl.getGridPointIndex(validate_index, x+1, y)
                    lwp, cbh, cth, reff = utl.updateGridInfo(lwp, cbh, cth, reff, interpolate_index, validate_index, total_lwp, validate_cbh, validate_cth, validate_reff)
                    if np.isnan(interpolate_index) == False:
                        validate_flag[interpolate_index][0] = 1.0
                else:
                    interpolate_index = utl.getGridPointIndex(validate_index, x+1, y-1)
                    lwp, cbh, cth, reff = utl.updateGridInfo(lwp, cbh, cth, reff, interpolate_index, validate_index, total_lwp, validate_cbh, validate_cth, validate_reff)
                    if np.isnan(interpolate_index) == False:
                        validate_flag[interpolate_index][0] = 1.0
        else:
            for j in range(3):
                if j == 0:
                    interpolate_index = utl.getGridPointIndex(validate_index, x, y+1)
                    lwp, cbh, cth, reff = utl.updateGridInfo(lwp, cbh, cth, reff, interpolate_index, validate_index, total_lwp, validate_cbh, validate_cth, validate_reff)
                    if np.isnan(interpolate_index) == False:
                        validate_flag[interpolate_index][0] = 1.0
                elif j == 1:
                    interpolate_index = utl.getGridPointIndex(validate_index, x+1, y)
                    lwp, cbh, cth, reff = utl.updateGridInfo(lwp, cbh, cth, reff, interpolate_index, validate_index, total_lwp, validate_cbh, validate_cth, validate_reff)
                    if np.isnan(interpolate_index) == False:
                        validate_flag[interpolate_index][0] = 1.0
                else:
                    interpolate_index = utl.getGridPointIndex(validate_index, x+1, y+1)
                    lwp, cbh, cth, reff = utl.updateGridInfo(lwp, cbh, cth, reff, interpolate_index, validate_index, total_lwp, validate_cbh, validate_cth, validate_reff)
                    if np.isnan(interpolate_index) == False:
                        validate_flag[interpolate_index][0] = 1.0
        
        validate_flag[i][0] = 1.0
        down_index = utl.getHeightIndex(cbh / 1000, divide_z)
        up_index = utl.getHeightIndex(cth / 1000, divide_z)
        z_range = np.linspace(down_index, up_index, num=up_index-down_index+1, endpoint=True)
        re_index = utl.getReffIndex(reff, round_line)
    

        cpt = cth - cbh
        lwc = lwp / (cpt * 4)
        lwc = round(lwc, 5)
        re = round_line[re_index]

        for j in range(len(z_range)):
            volume_grids.append(['', int(np.ceil((x-x_min+1)/2)), int(np.ceil((y-y_min+1)/2)), int(z_range[j]), 1, 1, lwc, re])

    if i == 10000:
        break

    

'''
    if i != 0:
        y_current = np.ceil((y-y_min+1)/2)
        y_previous = np.ceil((validate_index[i-1][1]-y_min)/2)
        if y_current == y_previous:
            next
        else:
            for j in range(len(z_range)):
                volume_grids.append(['', int(np.ceil((x-x_min+1)/2)), int(np.ceil((y-y_min+1)/2)), int(z_range[j]), 1, 1, lwc, re])
    else:
        for j in range(len(z_range)):
            volume_grids.append(['', int(np.ceil((x-x_min+1)/2)), int(np.ceil((y-y_min+1)/2)), int(z_range[j]), 1, 1, lwc, re])
'''

    
        

        
# 7.411636149

with open('big_les_cloud_2.part', 'w') as file:
    for item in volume_grids:
        line = ' '.join([str(i) for i in item])
        lines.append(line)
    result = '\n'.join(lines)
    file.write(result)

file.close()

# pd.DataFrame(validate_flag).to_csv("validate_flag.csv", header=None, index=False)
