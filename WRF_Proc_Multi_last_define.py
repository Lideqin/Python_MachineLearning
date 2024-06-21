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
    
    
    

# 下面的代码在colab里面不用运行，只是用作画图的测试

# -------------------------
# 检查是否提取到数据
if T2_celsius is None:
    print("没有找到2米温度数据")
    exit()

# 创建风速和风向数组
wind_speed = (U10**2 + V10**2)**0.5
wind_direction = 180 + (180 / 3.14) * np.arctan2(V10, U10)

# 子采样风场数据以减少箭头数量
subsample_factor = 20
lons_subsampled = lons[::subsample_factor, ::subsample_factor]
lats_subsampled = lats[::subsample_factor, ::subsample_factor]
U10_subsampled = U10[::subsample_factor, ::subsample_factor]
V10_subsampled = V10[::subsample_factor, ::subsample_factor]

# 绘制2米温度图
plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, edgecolor='black')
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAKES, edgecolor='black')
ax.add_feature(cfeature.RIVERS)

# 添加经纬度网格
gl = ax.gridlines(draw_labels=True, linestyle='--')
gl.top_labels = False
gl.right_labels = False
gl.xlocator = plt.FixedLocator(range(-180, 181, 5))  # 经度间隔
gl.ylocator = plt.FixedLocator(range(-90, 91, 5))   # 纬度间隔

# 绘制温度数据
temp_contour = ax.contourf(lons, lats, T2_celsius, cmap='coolwarm', transform=ccrs.PlateCarree())
plt.colorbar(temp_contour, ax=ax, orientation='horizontal', label='Temperature (°C)')

# 绘制风场箭头
#ax.quiver(lons, lats, U10, V10, transform=ccrs.PlateCarree(), color='black')
ax.quiver(lons_subsampled, lats_subsampled, U10_subsampled, V10_subsampled,
          transform=ccrs.PlateCarree(), color='black', scale=500)

# 添加箭头长度比例尺
scale_lon = lons_subsampled.max() - 10
scale_lat = lats_subsampled.min() + 2
ax.quiver(scale_lon, scale_lat, 5, 0, transform=ccrs.PlateCarree(), scale=500, label='5 m/s')
#ax.quiver(scale_lon, scale_lat, 5, 0, transform=ccrs.PlateCarree(), scale=500)
#plt.text(scale_lon + 0.5, scale_lat, '5 m/s', transform=ccrs.PlateCarree(), fontsize=10, color='black')


# 显示图例
plt.legend(loc='lower right', handlelength=2.5, handletextpad=1)

# 添加目标点
#ax.plot(target_lon, target_lat, 'ro', markersize=10, transform=ccrs.PlateCarree())
#plt.text(target_lon + 1, target_lat, f'{T2_at_target:.2f} °C', transform=ccrs.PlateCarree(), color='red')

# 添加标题
plt.title('2 Meter Temperature (°C)')

# 显示图像
plt.show()