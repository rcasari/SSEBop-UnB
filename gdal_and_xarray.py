# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 10:38:20 2019

@author: nrel
"""


from osgeo import gdal, osr,gdal_array
import xarray as xr
import numpy as np
import sys

def GetTiffInfobyName(in_filename, xa_dataarray):
#    Open netCDF file
#    src_ds = gdal.Open(in_filename)
#    if src_ds is None:
#        print('Open failed')
#        sys.exit()
#        
    src_ds_sd = gdal.Open(in_filename)
    # begin to read info of the named variable (i.e., subdataset)
    NDV = src_ds_sd.GetRasterBand(1).GetNoDataValue()
    xsize = src_ds_sd.RasterXSize
    ysize = src_ds_sd.RasterYSize
    GeoT = src_ds_sd.GetGeoTransform()
    Projection = osr.SpatialReference()
#        Projection = src_ds_sd.SpatialReference()
    Projection.ImportFromWkt(src_ds_sd.GetProjectionRef())
    # Close the subdataaet and the whole dataset
    src_ds_sd = None
#    src_ds = None
    data = np.ma.masked_array(xa_dataarray, mask=xa_dataarray==NDV, fill_value = NDV)

    return NDV, xsize, ysize, GeoT, Projection, data



def create_geotiff(suffix, Array, NDV, xsize, ysize, GeoT, Projection):
    '''
    Creates new Geotiff from array
    adapted from https://www.linkedin.com/pulse/convert-netcdf4-file-geotiff-using-python-chonghua-yin
    '''
    DataType = gdal_array.NumericTypeCodeToGDALTypeCode(Array.dtype)
    
    if type(DataType)!=np.int:
        if DataType.startswith('gdal.GDT_') == False:
            DataType = eval('gdal.GDT_'+DataType)
            
    NewFileName = suffix + '.tif'
    
    if len(Array.shape) > 2:
        zsize = Array.shape[0]
        # create a driver
        driver = gdal.GetDriverByName('GTiff')
        # Set nans to the original No Data Value
        Array[np.isnan(Array)] = NDV
        # Set up the dataset with zsize bands
        DataSet = driver.Create(NewFileName, xsize, ysize,zsize, DataType)
        DataSet.SetGeoTransform(GeoT)
        DataSet.SetProjection(Projection.ExportToWkt())
        # Write each slice of the array along the zsize
        for i in range(0,zsize):
            DataSet.GetRasterBand(i+1).WriteArray(Array[i])
            DataSet.GetRasterBand(i+1).SetNoDataValue(NDV)

    else:
        zsize = 1
        # create a driver
        driver = gdal.GetDriverByName('GTiff')
        # Set nans to the original No Data Value
        Array[np.isnan(Array)] = NDV
        # Set up the dataset with zsize bands
        DataSet = driver.Create(NewFileName, xsize, ysize,zsize, DataType)
        DataSet.SetGeoTransform(GeoT)
        DataSet.SetProjection(Projection.ExportToWkt())
        # Write each slice of the array along the zsize

        DataSet.GetRasterBand(1).WriteArray(Array)
        DataSet.GetRasterBand(1).SetNoDataValue(NDV)
            

    DataSet.FlushCache()
    return NewFileName
    


def GetnetCDFInfobyName(in_filename, var_name):
    #Open netCDF file
    src_ds = gdal.Open(in_filename)
    if src_ds is None:
        print('Open failed')
        sys.exit()
        
    if src_ds.GetSubDatasets()>1:
        # If exists more than one var in the NetCDF
        subdataset = 'NETCDF:"'+ in_filename + '":' + var_name 
        src_ds_sd = gdal.Open(subdataset)
        # begin to read info of the named variable (i.e., subdataset)
        NDV = src_ds_sd.GetRasterBand(1).GetNoDataValue()
        xsize = src_ds_sd.RasterXSize
        ysize = src_ds_sd.RasterYSize
        GeoT = src_ds_sd.GetGeoTransform()
        Projection = osr.SpatialReference()
#        Projection = src_ds_sd.SpatialReference()
        Projection.ImportFromWkt(src_ds_sd.GetProjectionRef())
        # Close the subdataaet and the whole dataset
        src_ds_sd = None
        src_ds = None
        # read data using xarray
        xr_ensemble = xr.open_dataset(in_filename)
        data = xr_ensemble[var_name]
        data = np.ma.masked_array(data, mask=data==NDV, fill_value = NDV)
        return NDV, xsize, ysize, GeoT, Projection, data















