# -*- coding: utf-8 -*-
"""
Created on Wed Jan 2 10:38:20 2019
Author: nrel
"""
from osgeo import gdal, osr, gdal_array
import xarray as xr
import numpy as np
import sys

def GetTiffInfobyName(in_filename, xa_dataarray):
    """
    Função para obter informações de um arquivo TIFF e mascarar os dados de um DataArray.

    Args:
        in_filename (str): Caminho para o arquivo TIFF.
        xa_dataarray (xarray.DataArray): DataArray com os dados a serem mascarados.

    Returns:
        tuple: NDV, xsize, ysize, GeoT, Projection, data
    """
    # Abrir o arquivo TIFF
    src_ds_sd = gdal.Open(in_filename)
    if src_ds_sd is None:
        print('Falha ao abrir o arquivo TIFF')
        sys.exit()

    # Obter informações do TIFF
    NDV = src_ds_sd.GetRasterBand(1).GetNoDataValue()
    xsize = src_ds_sd.RasterXSize
    ysize = src_ds_sd.RasterYSize
    GeoT = src_ds_sd.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(src_ds_sd.GetProjectionRef())

    # Fechar o dataset
    src_ds_sd = None

    # Mascarar dados
    data = np.ma.masked_array(xa_dataarray, mask=xa_dataarray == NDV, fill_value=NDV)
    return NDV, xsize, ysize, GeoT, Projection, data

def create_geotiff(suffix, Array, NDV, xsize, ysize, GeoT, Projection):
    """
    Cria um novo arquivo GeoTIFF a partir de um array.

    Args:
        suffix (str): Sufixo para o nome do arquivo.
        Array (np.array): Array com os dados.
        NDV (float): Valor de NoData.
        xsize (int): Largura do raster.
        ysize (int): Altura do raster.
        GeoT (tuple): Geotransformação.
        Projection (osr.SpatialReference): Projeção.

    Returns:
        str: Nome do novo arquivo GeoTIFF.
    """
    driver = gdal.GetDriverByName('GTiff')
    DataType = gdal_array.NumericTypeCodeToGDALTypeCode(Array.dtype)

    if isinstance(DataType, str) and not DataType.startswith('gdal.GDT_'):
        DataType = eval('gdal.GDT_' + DataType)

    NewFileName = suffix + '.tif'
    Array[np.isnan(Array)] = NDV

    if len(Array.shape) > 2:
        zsize = Array.shape[0]
    else:
        zsize = 1

    DataSet = driver.Create(NewFileName, xsize, ysize, zsize, DataType)
    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection(Projection.ExportToWkt())

    for i in range(zsize):
        band = DataSet.GetRasterBand(i + 1)
        band.WriteArray(Array[i] if zsize > 1 else Array)
        band.SetNoDataValue(NDV)

    DataSet.FlushCache()
    return NewFileName

def GetnetCDFInfobyName(in_filename, var_name):
    """
    Função para obter informações de um arquivo NetCDF e mascarar os dados de uma variável.

    Args:
        in_filename (str): Caminho para o arquivo NetCDF.
        var_name (str): Nome da variável de interesse.

    Returns:
        tuple: NDV, xsize, ysize, GeoT, Projection, data
    """
    # Abrir o arquivo NetCDF
    src_ds = gdal.Open(in_filename)
    if src_ds is None:
        print('Falha ao abrir o arquivo NetCDF')
        sys.exit()

    if len(src_ds.GetSubDatasets()) > 1:
        subdataset = f'NETCDF:"{in_filename}":{var_name}'
        src_ds_sd = gdal.Open(subdataset)
    else:
        src_ds_sd = src_ds

    # Obter informações da variável
    NDV = src_ds_sd.GetRasterBand(1).GetNoDataValue()
    xsize = src_ds_sd.RasterXSize
    ysize = src_ds_sd.RasterYSize
    GeoT = src_ds_sd.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(src_ds_sd.GetProjectionRef())

    # Fechar os datasets
    src_ds_sd = None
    src_ds = None

    # Ler dados usando xarray
    xr_ensemble = xr.open_dataset(in_filename)
    data = xr_ensemble[var_name].values
    data = np.ma.masked_array(data, mask=data == NDV, fill_value=NDV)
    return NDV, xsize, ysize, GeoT, Projection, data
