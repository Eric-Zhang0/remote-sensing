import netCDF4 as nc
import numpy as np
import pandas as pd
from utilities import utl

MOD06_FILE = 'MYD06_L2.A2016157.0355.061.2018059101934.hdf'
MOD06_FILE_EXTRA = 'MYD06_L2.A2016157.0400.061.2018059102206.hdf'
MYD05_FILE = "MYD05_L2.A2016157.0355.061.2018059095657.hdf"
MYD05_FILE_EXTRA = "MYD05_L2.A2016157.0400.061.2018059095928.hdf"

test_line = np.linspace(5.0, 60.0, num=110, endpoint=True)
# print(test_line)
round_line = np.around(test_line, decimals=4)
# print(round_line)
# total_reff = np.vstack((utl.readModis(MOD06_FILE, 'Cloud_Effective_Radius'), utl.readModis(MOD06_FILE_EXTRA, 'Cloud_Effective_Radius')))
# total_surt = np.vstack((utl.readModis(MOD06_FILE, 'surface_temperature_1km'), utl.readModis(MOD06_FILE_EXTRA, 'surface_temperature_1km')))
total_lwp = np.vstack((utl.readModis(MOD06_FILE, 'Cloud_Water_Path'), utl.readModis(MOD06_FILE_EXTRA, 'Cloud_Water_Path')))
# unit: g/m^2
total_water_vapor = np.vstack((utl.readModis(MYD05_FILE, 'Water_Vapor_Near_Infrared'), utl.readModis(MYD05_FILE_EXTRA, 'Water_Vapor_Near_Infrared')))
# unit: cm

# pd.DataFrame(total_surt).to_csv("surface_temperature_1km.csv", header=None, index=False)
ts_around_geo = pd.read_csv("geolocation.csv", header=None).to_numpy()
# ts_around_surft = utl.getAroundSDS(total_surt, ts_around_geo)
# print("地表温度最低：", np.nanmin(ts_around_surft))
# print("地表温度最高：", np.nanmax(ts_around_surft))

validate_info = pd.read_csv("validate_info.csv", header=None).to_numpy()
validate_cth = validate_info[:, 0:1]
validate_cbh = validate_info[:, 1:2]
# print("最低的云底高度：", np.amin(validate_cbh))
# print("最高的云顶高度：", np.amax(validate_cth))
validate_cpt = validate_info[:, 2:3]
validate_ctype = validate_info[:, 3:4]
validate_lwc = validate_info[:, 4:5]
validate_reff = validate_info[:, 5:6]

validate_x = validate_info[:, 6:7].astype(int)
validate_y = validate_info[:, 7:8].astype(int)
x_min = np.amin(validate_x)
x_max = np.amax(validate_x)
y_min = np.amin(validate_y)
y_max = np.amax(validate_y)

validate_index = np.hstack((validate_x, validate_y)).astype(int)
validate_flag = np.zeros((len(validate_info), 1))

x_num = x_max - x_min + 1
y_num = y_max - y_min + 1
num_cell_grid = 0
if x_num > y_num:
    num_cell_grid = x_num
else:
    num_cell_grid = y_num

layer_num = 43
z_start = 0.02
z_interval = 0.04
# create height levels line
divide_z = np.arange(z_start, z_start + z_interval * layer_num, z_interval)
# print(divide_z)

# create temperature levels line
divide_t = np.linspace(215.7000, 294.2000, num=layer_num, endpoint=True)
divide_t = np.around(divide_t, decimals=4)
divide_t = divide_t[::-1]

# create pressure levels line
divide_p = np.linspace(10100.4199, 101300.0000, num=layer_num, endpoint=True)
divide_p = np.around(divide_p, decimals=4)
divide_p = divide_p[::-1]

# Create a netCDF file
nc_file = nc.Dataset("clouds_km.nc", "w")

# Create dimensions
n_l = 875  # Number of grid cells in the W-E direction
n_m = 825  # Number of grid cells in the S-N direction
n_n = len(divide_z)  # Number of vertical levels

nc_file.createDimension("W_E_direction", n_l)
nc_file.createDimension("S_N_direction", n_m)
nc_file.createDimension("vertical_levels", n_n)
nc_file.createDimension("RVT", 1)
nc_file.createDimension("RCT", 1)
nc_file.createDimension("PABST", 1)
nc_file.createDimension("THT", 1)

# Create variables
W_E_direction = nc_file.createVariable("W_E_direction", "f4", ("W_E_direction",))
S_N_direction = nc_file.createVariable("S_N_direction", "f4", ("S_N_direction",))
vertical_levels = nc_file.createVariable("vertical_levels", "f4", ("vertical_levels",))
RVT = nc_file.createVariable("RVT", "f4", ("W_E_direction", "S_N_direction", "vertical_levels", "RVT",))
RCT = nc_file.createVariable("RCT", "f4", ("W_E_direction", "S_N_direction", "vertical_levels", "RCT",))
PABST = nc_file.createVariable("PABST", "f4", ("W_E_direction", "S_N_direction", "vertical_levels", "PABST",))
THT = nc_file.createVariable("THT", "f4", ("W_E_direction", "S_N_direction", "vertical_levels", "THT",))

# Fill in variable data
intervel = 1000
start = 500
W_E_direction[:] = np.arange(start, start + intervel * n_l, intervel)
S_N_direction[:] = np.arange(start, start + intervel * n_m, intervel)
vertical_levels[:] = divide_z*1000
# print("test for np.arange's length: ", len(np.arange(start, start + intervel * n_l, intervel)))

RVT_data = np.ones((n_l, n_m, n_n, 1))*0.01
RCT_data = np.ones((n_l, n_m, n_n, 1))*0.0001
PABST_data = np.zeros((n_l, n_m, n_n, 1))
THT_data = np.zeros((n_l, n_m, n_n, 1))

# The constant for calculating air density next
specific_gas_constant = 287.05

for i in range(n_l):
    for j in range(n_m):
        for k in range(n_n):
            PABST_data[i][j][k][0] = divide_p[k]
            THT_data[i][j][k][0] = divide_t[k]


for m in range(len(validate_info)):
    x = validate_index[m][0]
    y = validate_index[m][1]
    cbh = validate_cbh[m][0]
    cth = validate_cth[m][0]
    cpt = validate_cpt[m][0]
    lwp = total_lwp[x][y]
    water_vapor = total_water_vapor[x][y]
    used_flag = validate_flag[m][0]

    lwc = validate_lwc[m][0]
    lwc = round(lwc, 4)
    down_index = utl.getHeightIndex(cbh / 1000, divide_z)
    up_index = utl.getHeightIndex(cth / 1000, divide_z)
    z_range = np.linspace(down_index, up_index+1)
    air_density = 0.0

    for index in range(len(z_range)):

        if x-x_min <= (n_l - 1) and y-y_min <= (n_m - 1):
            pressure = PABST_data[x-x_min][y-y_min][int(z_range[index])][0]
            temperature = THT_data[x-x_min][y-y_min][int(z_range[index])][0]
            air_density = pressure / (specific_gas_constant * temperature)

            RCT_data[x-x_min][y-y_min][int(z_range[index])][0] = lwc / 1000

            if np.isnan(water_vapor) == True:
                RVT_data[x-x_min][y-y_min][int(z_range[index])][0] = 0.01
            else:
                RVT_data[x-x_min][y-y_min][int(z_range[index])][0] = 0.02   
        '''
        else:
            x_module = (x-x_min) % 200
            y_module = (y-y_min) % 100

            pressure = PABST_data[x_module][y_module][int(z_range[index])][0]
            temperature = THT_data[x_module][y_module][int(z_range[index])][0]
            air_density = pressure / (specific_gas_constant * temperature)

            RCT_data[x_module][y_module][int(z_range[index])][0] += lwc * air_density / 1000

            if np.isnan(water_vapor) == True:
                RVT_data[x_module][y_module][int(z_range[index])][0] += lwc * air_density / 1000
            else:
                RVT_data[x_module][y_module][int(z_range[index])][0] += (water_vapor / 100) * 1000 / air_density
        '''

# Generate dummy data for VLEV
RVT[:, :, :] = RVT_data

RCT[:, :, :] = RCT_data

PABST[:, :, :] = PABST_data

THT[:, :, :] = THT_data

# Set variable attributes
W_E_direction.units = "m"
S_N_direction.units = "m"
vertical_levels.units = "m"
RVT.units = "kg/kg"
RCT.units = "kg/kg"
PABST.units = "pa"
THT.unitis = "K"

# Close the netCDF file
nc_file.close()

# pd.DataFrame(RVT_data.flatten()).to_csv("water_vapor_rvt.csv", index=False, header=None)
# pd.DataFrame(RCT_data.flatten()).to_csv("liquid_water_rct.csv", index=False, header=None)