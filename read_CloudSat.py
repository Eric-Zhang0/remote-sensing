from pyhdf import HDF, VS, SD

import numpy as np
import numpy.ma as ma
import pandas as pd

class reader:
    '''用于读取CloudSat R05 2B产品的类.'''

    def __init__(self, fname):
        '''打开HDF,Vdata table,和scientific dataset.'''

        self.hdf = HDF.HDF(fname, HDF.HC.READ)
        self.vs = self.hdf.vstart()
        self.sd = SD.SD(fname, SD.SDC.READ)

    def attach_vdata(self, varname):
        '''从Vdata table中读取一个变量的所有记录.'''

        vdata = self.vs.attach(varname)
        data = vdata[:]
        vdata.detach()

        return data

    def scale_and_mask(self, data, varname):
        '''依据变量的factor进行放缩,并根据valid_range来mask掉填充值.'''

        factor = self.attach_vdata(f'{varname}.factor')[0][0]
        valid_min, valid_max = \
            self.attach_vdata(f'{varname}.valid_range')[0][0]

        invalid = np.logical_or(data <= valid_min, data >= valid_max)
        data = ma.array(data, mask=invalid)
        data = data / factor

        return data

    def read_geo(self, process=True):
        '''读取经纬度和地形高度.参数process指定是否放缩并mask地形高度.'''

        lon = np.array(self.attach_vdata('Longitude'))[:, 0]
        lat = np.array(self.attach_vdata('Latitude'))[:, 0]
        elv = np.array(self.attach_vdata('DEM_elevation'))[:, 0]

        if process:
            elv = self.scale_and_mask(elv, 'DEM_elevation')

        return lon, lat, elv

    def read_time(self, datetime=True):
        '''读取每个数据点对应的时间.

        datetime=True: 返回所有数据点的日期时间组成的DatetimeIndex.
        datetime=False: 返回所有数据点相对于第一个点所经过的秒数组成的numpy数组.
        '''

        seconds = np.array(self.attach_vdata('Profile_time'))[:, 0]

        if datetime:
            TAI = self.attach_vdata('TAI_start')[0][0]
            start = pd.to_datetime('1993-01-01') + pd.Timedelta(seconds=TAI)
            offsets = pd.to_timedelta(seconds, unit='s')
            time = pd.date_range(start=start, end=start, periods=offsets.size)
            time = time + offsets
            return time
        else:
            return seconds

    def read_sds(self, varname, process=True):
        '''读取scientific dataset.参数process指定是否放缩并mask.'''

        data = self.sd.select(varname)[:]
        if process:
            data = self.scale_and_mask(data, varname)

        return data

    def close(self):
        '''关闭HDF文件.'''

        self.vs.end()
        self.sd.end()
        self.hdf.close()


# 测试数据.
if __name__ == '__main__':
    fname = 'D:/remote-sensing/2016157030402_53752_CS_2B-CLDCLASS-LIDAR_GRANULE_P1_R05_E06_F00.hdf'
    f = reader(fname)
    f.close()
