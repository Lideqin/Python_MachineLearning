# 下面的代码放到colab里面可以运行
#  功能：step1：批量的读取WRF数据（一次预报的数据，总共84个文件，从wrfout UPP后的grib数据）
#        step2：读取站点的数据
#        step3：将WRF数据插值到站点上
#        step4：输出到csv格式的文件

!pip install pygrib
!pip install matplotlib cartopy


# 如果在本地机器上安装了pygrib，则直接把下面的代码存储到本地机器上运行

import pygrib
import numpy as np
import os
import csv
from glob import glob
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.interpolate import RegularGridInterpolator
from google.colab import drive

#设置google云端数据路径
drive.mount('/content/drive')
wrf_grib_file_path = '/content/drive/My Drive/WRF_Data/'
# 创建 GRIB 文件名列表
grib_files = [os.path.join(wrf_grib_file_path, f'BCSY_DB3KM_2024061000_{i:02d}.grib2') for i in range(85)]

station_file_path = '/content/drive/My Drive/WRF_Data/LN_ID-.csv'
output_file_path = '/content/drive/My Drive/WRF_Data/interpolated_BCSY_DB3KM_2024061000.csv'

def read_station_data(station_file_path):
    stations = []
    with open(station_file_path, 'r', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            stations.append({
                'ID': row['ID'],
                'Lat': float(row['Lat']),
                'Lon': float(row['Lon']),
                'Alti': float(row['Alti'])
            })
    return stations

def interpolate_data(data, lats, lons, point):
    interp_func = RegularGridInterpolator((lats[:, 0], lons[0, :]), data)
    return interp_func(point)

def process_grib_file(grib_file, stations, writer):
    try:
        wrf_grbs = pygrib.open(grib_file)
        print(f"Processing file: {grib_file}")
    except OSError as e:
        print(f"Error opening file {grib_file}: {e}")
        return
    #        0  1   2     3
    # 提取时间信息 BCSY_DB3KM_2024061000_11.grib2
    #        012345678901234567890123
    file_name_parts = os.path.basename(grib_file).split('_')
    date_str = file_name_parts[2]  # '2024010100'
    year = date_str[:4]
    month = date_str[4:6]
    day = date_str[6:8]
    hour = date_str[8:10]
    fcst_str = file_name_parts[3]
    Fcst_Hour = fcst_str[0:2]
    print(Fcst_Hour)
    # 初始化变量
    variable_data = {
        'U10': None, 'V10': None, 'T2_celsius': None, 'RH2': None, 'TD2': None, 'Vis2': None,
        'CC': None, 'TCICW': None, 'Precip': None, 'U_1000': None, 'V_1000': None, 'W_1000': None,
        'T_1000': None, 'U_925': None, 'V_925': None, 'W_925': None, 'T_925': None, 'T_850': None,
        'T_500': None, 'FricV': None, 'CAPE': None, 'CIN': None
    }
    lats = lons = None

    for grb in wrf_grbs:
        if grb.name == '10 metre U wind component':
            variable_data['U10'], lats, lons = grb.data()
        elif grb.name == '10 metre V wind component':
            variable_data['V10'], lats, lons = grb.data()
        elif grb.name == '2 metre temperature':
            T2_kelvin, lats, lons = grb.data()
            variable_data['T2_celsius'] = T2_kelvin - 273.15
        elif grb.name == '2 metre relative humidity':
            variable_data['RH2'], lats, lons = grb.data()
        elif grb.name == 'Dew point temperature':
            variable_data['TD2'], lats, lons = grb.data()
        elif grb.name == 'Visibility':
            variable_data['Vis2'], lats, lons = grb.data()
        elif grb.name == 'Total Cloud Cover':
            variable_data['CC'], lats, lons = grb.data()
        elif grb.name == 'Total column-integrated cloud water':
            variable_data['TCICW'], lats, lons = grb.data()
        elif grb.name == 'Total Precipitation':
            variable_data['Precip'], lats, lons = grb.data()

        elif grb.name == 'U component of wind' and grb.level == 1000:
            variable_data['U_1000'], lats, lons = grb.data()
        elif grb.name == 'V component of wind' and grb.level == 1000:
            variable_data['V_1000'], lats, lons = grb.data()
        elif grb.name == 'Vertical velocity' and grb.level == 1000:
            variable_data['W_1000'], lats, lons = grb.data()
        elif grb.name == 'Temperature' and grb.level == 1000:
            variable_data['T_1000'], lats, lons = grb.data()
        elif grb.name == 'U component of wind' and grb.level == 925:
            variable_data['U_925'], lats, lons = grb.data()
        elif grb.name == 'V component of wind' and grb.level == 925:
            variable_data['V_925'], lats, lons = grb.data()
        elif grb.name == 'Vertical velocity' and grb.level == 925:
            variable_data['W_925'], lats, lons = grb.data()
        elif grb.name == 'Temperature' and grb.level == 925:
            variable_data['T_925'], lats, lons = grb.data()
        elif grb.name == 'Temperature' and grb.level == 850:
            variable_data['T_850'], lats, lons = grb.data()
        elif grb.name == 'Temperature' and grb.level == 500:
            variable_data['T_500'], lats, lons = grb.data()
        elif grb.name == 'Frictional velocity':
            variable_data['FricV'], lats, lons = grb.data()
        elif grb.name == 'Convective available potential energy':
            variable_data['CAPE'], lats, lons = grb.data()
        elif grb.name == 'Convective inhibition':
            variable_data['CIN'], lats, lons = grb.data()

    wrf_grbs.close()

    for station in stations:
        lat = station['Lat']
        lon = station['Lon']
        point = (lat, lon)
        interpolated_data = {
            'ID': station['ID'],
            'Lat': lat,
            'Lon': lon,
            'Alti': station['Alti'],
            'Year': year,
            'Month': month,
            'Day': day,
            'Hour': hour,
            'Fcst_Hour': Fcst_Hour
        }

        for var_name in variable_data:
            if variable_data[var_name] is not None:
                interpolated_data[var_name] = interpolate_data(variable_data[var_name], lats, lons, point)

        writer.writerow(interpolated_data)

def main(grib_files, station_file_path, output_file_path):
    stations = read_station_data(station_file_path)

    with open(output_file_path, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = [
            'ID', 'Lat', 'Lon', 'Alti',
            'Year', 'Month', 'Day', 'Hour',
            'Fcst_Hour',
            'U10', 'V10', 'T2_celsius', 'RH2',
            'TD2', 'Vis2', 'CC', 'TCICW',
            'Precip',
            'U_1000', 'V_1000', 'W_1000',
            'T_1000',
            'U_925', 'V_925', 'W_925',
            'T_925',
            'T_850', 'T_500',
            'FricV',
            'CAPE', 'CIN',
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for grib_file in grib_files:
            process_grib_file(grib_file, stations, writer)

if __name__ == '__main__':
    main(grib_files, station_file_path, output_file_path) 
    
    
    
