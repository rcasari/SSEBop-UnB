# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 18:45:33 2018

@author: LeandroSalles
"""


#
#import gdal, ogr, osr
#import sys
import matplotlib.pyplot as plt
#import rasterio
import numpy as np
#import math
#import fiona
#import rasterio.mask
import glob
import os
from datetime import datetime as dt
import pandas as pd

#
#path_met_FAL = 'C:/FAL/dados_meteorologicos_FAL'
#
#arquivo_met_FAL = 'unb_2019_11_SSEBop.xlsx'
#planilha_met_FAL = 'outubro'

def estacao_FAL(path_met_FAL,arquivo_met_FAL,planilha_met_FAL):
    
    
    os.chdir(path_met_FAL)
    DB = pd.read_excel(arquivo_met_FAL, sheet_name= planilha_met_FAL, index_col='Data') 
    DB['ano'] = DB.index.year
    DB['mes'] = DB.index.month
    DB['dia'] = DB.index.day

    
    
    
    
    
    
    
    
    
    ## Estacao FAL :
    
    # Graus decimais
    #lat= 15º56' S
    #lon= 47º56´W
    #Graus min seg
    lat=-1*(15+(56/60))
    lon=-1*(47+(56/60))
    Lat=(np.pi/180)*lat #lat latitude em radiano
    Altitude= 1080
    
    ## Constantes:
    cp=1013#*10**-3
    Patm=101.3*((293-(0.0065*Altitude))/293)**5.26
    e_ratio_molec=0.622
    #SBC= #SBC Stefan-Boltzmann constant
    Gsc=0.0820
    SBC=4.903*10**-9
    
    
    
    #dados_str=pd.DataFrame(dtype=np.str)
    dados=pd.DataFrame(dtype=np.float64)
    dados['ano'] = DB['ano']
    dados['mes'] = DB['mes']
    dados['dia'] = DB['dia']
    dados['P'] = DB["Precipitacao (mm)"]#.groupby(pd.Grouper(freq="d")).mean()

    ############ VERIFICAR 
#    dados['Pacum1'] = dados.groupby(['dia'])['P'].cumsum()  
#    dados['Pacum_soma_ate_fim'] = dados['P'].cumsum()#.groupby(pd.Grouper(freq="d")).mean()
    dados['Tmin']=DB["T min (oC)"]#.groupby(pd.Grouper(freq="d")).mean()
    dados['Tmax']=DB["T max (oC)"]#.groupby(pd.Grouper(freq="d")).mean()
#    dados['Insol']=DB["Insolação (hrs)"]#.groupby(pd.Grouper(freq="d")).mean()
    dados['URMED']=DB["UR med (%)"]#.groupby(pd.Grouper(freq="d")).mean()
    dados['URMAX']=DB["UR max (%)"]
    dados['URMIN']=DB["UR min (%)"]
    dados['URmed']=(DB["UR max (%)"]+DB["UR min (%)"])/2
    #dados['UR']=DB["URMED_L"]
    dados['UR']=DB["UR med (%)"]
#    dados['Vento']=DB["Vento med (m/s)"]#.groupby(pd.Grouper(freq="d")).mean()
    dados['Vento']=DB['Vento media BSB (m/s)']
    
    #dados['P'].to_csv('P3.txt', header=True, index=True, sep=',', mode='a')
    #dados['Tmin'].to_csv('Tmin.txt', header=True, index=True, sep=',', mode='a')
    #dados['Insol'].to_csv('Insolacao.txt', header=True, index=True, sep=',', mode='a')
    
    dados['Tmin_K']=dados['Tmin']+273.16
    dados['Tmax_K']=dados['Tmax']+273.16
    
    #Tmed: media temperatura entre min e máximo - Conforme a conveção definida pela FAO 56
    # FAO: dados['Tmed']=(dados['Tmax']+dados['Tmin'])/2
    dados['Tmed']=(dados['Tmax']+dados['Tmin'])/2
    dados['Tmed_p']=DB["T med (oC)"]
    #Jdia: Dia N
    dados['DOY'] = dados.index.dayofyear
    #ea(kPa) - Actual vapour pressure - Tmed para Slope of saturation pressure curve
    dados['ea_Tmed']=0.6108*np.exp((17.27*dados['Tmed'])/(dados['Tmed']+237.3))
    #ea_min(kPa) - Actual vapour pressure minimum SSEBop
    dados['ea_Tmin']=0.6108*np.exp((17.27*dados['Tmin'])/(dados['Tmin']+237.3))
    #ea_max(kPa) - Actual vapour pressure maximum
    dados['ea_Tmax']=0.6108*np.exp((17.27*dados['Tmax'])/(dados['Tmax']+237.3))
    #es (kPa) - Mean saturation Vapor Pressure
    
    # FAO: dados['es']=(dados['ea_Tmax']+dados['ea_Tmin'])/2
    # vFAL
    dados['es']=(dados['ea_Tmax']+dados['ea_Tmin'])/2
    #dados['es']=dados['ea_Tmed']
    #ea(kPa) - Actual vapour pressure - Tmed para Slope of saturation pressure curve
    dados['ea_URmed']=dados['es']*(dados['UR']/100)
    dados['ea']=(dados['ea_Tmin']*(dados['URMAX']/100) + dados['ea_Tmax']*(dados['URMIN']/100))/2

    #ea_defict: Vapour pressure deficit
    dados['ea_defict']=dados['es']-dados['ea']
    # Radiação
    # dr(rad) - inverse relative distance Earth-Sun
    dados['dr']=1+0.033*np.cos((2*np.pi/365)*dados['DOY'])
    # δ: Solar declination
    dados['Sol_decl']=0.409*np.sin(((2*np.pi/365)*dados['DOY'])-1.39)
    #  Ws (rad) - sunset hour angle
    # FAO: dados['Ws']=np.arccos(-np.tan(Lat)*np.tan(dados['Sol_decl']))
#    dados['X']=(1-(np.tan(Lat))**2*np.tan(dados['Sol_decl'])**2)
    #for i in range(0,len(dados['X'])):
    #    
    #    
    #    if dados['X'][i] > 0:
    #       continue      
    ##    dados.loc[dados['X'][i]] = dados['X'][i] 
    #    else:
    #        dados.loc[['X'][i]]=0.00001
    #        
    #
    #dados['Ws']=(np.pi/2)-np.arctan((-np.tan(Lat)*np.tan(dados['Sol_decl']))/dados['X']**0.5)
    
    dados['Ws']=np.arccos(-np.tan(Lat)*np.tan(dados['Sol_decl']))
    # Ra (MJm-2d-1) - Extraterrestrial daily radiation
    dados['Ra']=((24*60)/np.pi)*Gsc*dados['dr']*((dados['Ws']*np.sin(Lat)*np.sin(dados['Sol_decl']))+(np.cos(Lat)*np.cos(dados['Sol_decl'])*np.sin(dados['Ws'])))
    ## N (horas) - Day light hours
    #dados['N']=(24/np.pi)*dados['Ws']


    ## Rs (MJ m-2) - Solar radiation
  
    kRs=0.16
    dados['Rs']=kRs*(dados['Tmax']-dados['Tmin'])**0.5*dados['Ra'] # Atualizacao do Rs usando a equacao (n50) da FAO-56 proposto no paper de Senay,2018 
    #Rso (MJm-2d-1) - Clear Sky solar radiation
    dados['Rso']=(0.75+(2*10**-5*Altitude))*dados['Ra']

    #Rns (MJm-2d-1) - Net shortwave radiation
    # FAO: dados['Rns']=(1-0.23)*dados['Rs']
    dados['Rns']=(1-0.23)*dados['Rs']
    dados['Rnso']=(1-0.23)*dados['Rso']
    # dados['γ']=((cp*Patm)/(e_ratio_molec*dados['λ']))*10**-3
#    dados['Rnl']=SBC*((dados['Tmax_K']**4+dados['Tmin_K']**4)/2)*(0.34-(0.14*dados['ea']**0.5))*((1.35*(dados['Rs']/dados['Rso']))-0.35)
    
    # SSEBop 2018
    dados['Rnl']=SBC*((dados['Tmax_K']**4+dados['Tmin_K']**4)/2)*(0.34-(0.14*dados['ea_Tmin']**0.5))*((1.35*(dados['Rs']/dados['Rso']))-0.35)
    # SSEBop Clear sky condition
    dados['Rnlo']=SBC*((dados['Tmax_K']**4+dados['Tmin_K']**4)/2)*(0.34-(0.14*dados['ea_Tmin']**0.5))*((1.35*(dados['Rso']/dados['Rso']))-0.35)

    # Rn (MJ/m²/dia) - Net Radiation
    dados['Rn']=dados['Rns']-dados['Rnl']   # SSEBop 2018
    dados['Rno']=dados['Rnso']-dados['Rnlo'] # SSEBop Clear sky condition
    
    #Rn (W/m²)
    dados['Rn_Wm-2']=dados['Rn']/0.0864
    #λ Latent Heat of Vaporization
    dados['λ']=2.501-(2.361*10**-3*dados['Tmed'])
    #γ: Constante Psicometrica
    # FAO: dados['γ']=((cp*Patm)/(e_ratio_molec*dados['λ']))*10**-3
    dados['γ']=Patm*0.665*10**-3
    # Slope vapour pressure curve
    dados['delta']=(4098*dados['ea_Tmed'])/(dados['Tmed']+237.3)**2
    #Soil heat Flux
    dados['G']=0
    
    
    # ETo Penman-Montheith Equation
    # Rs_EToPM,Rns_EToPM,Rnl_EToPM and Rn_EToPM are for ETo Penmann-Montheith equation while other Rs are for SSEBop 2013 e 2018
    #dados['Rs']=(0.25+(0.5*(dados['Insol']/dados['N'])))*dados['Ra']
#    dados['Rs_EToPM']=(0.25+(0.5*(dados['Insol']/dados['N'])))*dados['Ra']
        
    dados['Rs_EToPM']=DB["Rad. Global (MJ/m2 d)"] # dado da planilha da FAL: checar essa variavel
    dados['Rns_EToPM']=(1-0.23)*dados['Rs_EToPM']     # ETo Penmann-Montheith equation
    dados['Rnl_EToPM']=SBC*((dados['Tmax_K']**4+dados['Tmin_K']**4)/2)*(0.34-(0.14*dados['ea']**0.5))*((1.35*(dados['Rs_EToPM']/dados['Rso']))-0.35)
    dados['Rn_EToPM']=dados['Rns_EToPM']-dados['Rnl_EToPM']     # Rn for ETo Penmann-Montheith
    dados['ETo']=((0.408*dados['delta']*(dados['Rn_EToPM']-dados['G']))+(dados['γ']*(900/(dados['Tmed']+273))*dados['Vento']*(dados['es']-dados['ea'])))/(dados['delta']+dados['γ']*(1+0.34*dados['Vento']))
    #Coeficientes para alfafa 1600 e 0.38 (http://edis.ifas.ufl.edu/pdffiles/AE/AE45900.pdf e outros)
    dados['ETo_alfafa']=((0.408*dados['delta']*(dados['Rn_EToPM']-dados['G']))+(dados['γ']*(1600/(dados['Tmed']+273))*dados['Vento']*(dados['es']-dados['ea'])))/(dados['delta']+dados['γ']*(1+0.38*dados['Vento']))


#    dados['ETo']=((0.408*dados['delta']*(dados['Rn']-dados['G']))+(dados['γ']*(900/(dados['Tmed']+273))*dados['Vento']*(dados['es']-dados['ea_Tmin'])))/(dados['delta']+dados['γ']*(1+0.34*dados['Vento']))
#    dados['ETo']=((0.408*dados['delta']*(dados['Rn']-dados['G']))+(dados['γ']*(900/(dados['Tmed']+273))*dados['Vento']*(dados['es']-dados['ea'])))/(dados['delta']+dados['γ']*(1+0.34*dados['Vento']))
    
#    #dT - SSEBop
    dados['Tkv']=1.01*(dados['Tmed']+273.16)
    dados['pa']=(1000*Patm)/(dados['Tkv']*287)
#    dados['dT_Rso']=((dados['Rso']/0.0864)*110)/(1.013*10**-3*dados['pa'])
#    dados['dT_Rs']=((dados['Rs']/0.0864)*110)/(1.013*10**-3*dados['pa'])
#    dados['dT']=((dados['Rs']/0.0864)*110)/(1.013*10**-3*dados['pa'])


    cp=1013
    dados['dT_v2018']=((dados['Rn']/0.0864)*110)/(cp*dados['pa']) # SSEBop 2018
    dados['dT_v2013']=((dados['Rno']/0.0864)*110)/(cp*dados['pa']) # SSEBop Clear sky condition
    dados['dT']=((dados['Rno']/0.0864)*110)/(cp*dados['pa']) # SSEBop 2018

    
#    dados['Rso_Wm-2']=dados['Rso']/0.0864

#    dados['dT_Rso_Wm-2']=(dados['Rso_Wm-2']*110)/(1013*dados['pa'])
    return dados
