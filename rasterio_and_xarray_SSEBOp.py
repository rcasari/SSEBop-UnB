# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 15:27:06 2018
Author: Leandro Salles
"""
import rasterio
import xarray as xr
import numpy as np
import pandas as pd
import os

def rasterio_to_xarray(fname):
    """
    Converte um arquivo compatível com rasterio para um objeto xarray.DataArray.

    Args:
        fname (str): Nome do arquivo rasterio-compatível a ser lido.

    Returns:
        xarray.DataArray: Objeto xarray.DataArray contendo os dados do arquivo, 
        juntamente com os metadados geográficos relevantes.
    """
    with rasterio.open(fname) as src:
        data = src.read(1)
        data = np.where(data == src.nodata, np.nan, data)
        
        # Dimensões e coordenadas
        nx, ny = src.width, src.height
        x0, y0 = src.bounds.left, src.bounds.top
        dx, dy = src.res[0], -src.res[1]
        
        lat = np.arange(start=y0, stop=(y0 + ny * dy), step=dy)
        lon = np.arange(start=x0, stop=(x0 + nx * dx), step=dx)
        
        coords = {'latitude': lat[:ny], 'longitude': lon[:nx]}
        dims = ('latitude', 'longitude')
        attrs = {}
        
        try:
            aff = src.transform
            attrs['affine'] = aff.to_gdal()
        except AttributeError:
            pass
        
        try:
            c = src.crs
            attrs['crs'] = c.to_string()
            kwargs = src.meta.copy()
            attrs.update(kwargs)
        except AttributeError:
            pass
    
    return xr.DataArray(data, dims=dims, coords=coords, attrs=attrs)

def xarray_to_rasterio(xa, output_filename):
    """
    Converte um objeto xarray.DataArray para um arquivo raster usando rasterio.

    Args:
        xa (xarray.DataArray): Objeto xarray.DataArray a ser convertido.
        output_filename (str): Nome do arquivo GeoTIFF de saída.

    Notes:
        Esta função suporta apenas DataArrays 2D ou 3D e saída em GeoTIFF.
        O DataArray de entrada deve conter atributos (armazenados em xa.attrs)
        especificando metadados geográficos.
    """
    xa = xa.load()
    
    if len(xa.shape) == 2:
        count, height, width = 1, xa.shape[0], xa.shape[1]
        band_indices = [1]
    else:
        count, height, width = xa.shape[0], xa.shape[1], xa.shape[2]
        band_indices = np.arange(count) + 1
    
    processed_attrs = {}
    
    try:
        val = xa.attrs['affine']
        processed_attrs['transform'] = rasterio.Affine.from_gdal(*val)
    except KeyError:
        pass
    
    try:
        val = xa.attrs['crs']
        processed_attrs['crs'] = rasterio.crs.CRS.from_string(val)
    except KeyError:
        pass
    
    with rasterio.open(output_filename, 'w', driver='GTiff', height=height, width=width,
                       dtype=str(xa.dtype), count=count, **processed_attrs) as dst:
        dst.write(xa.values, band_indices)

def xarray_to_rasterio_by_band(xa, output_basename, dim='time', date_format='%Y-%m-%d'):
    """
    Converte um objeto xarray.DataArray para múltiplos arquivos raster, um para cada banda.

    Args:
        xa (xarray.DataArray): Objeto xarray.DataArray a ser convertido.
        output_basename (str): Base do nome do arquivo de saída.
        dim (str): Nome da dimensão a ser usada para separar as bandas.
        date_format (str): Formato de data para os nomes dos arquivos.

    Notes:
        Esta função cria um arquivo GeoTIFF separado para cada banda do DataArray.
    """
    for i in range(len(xa[dim])):
        args = {dim: i}
        data = xa.isel(**args)
        index_value = data[dim].values
        
        if isinstance(index_value, np.datetime64):
            formatted_index = pd.to_datetime(index_value).strftime(date_format)
        else:
            formatted_index = str(index_value)
        
        filename = f"{output_basename}{formatted_index}.tif"
        xarray_to_rasterio(data, filename)
        print(f'Exported {formatted_index}')

