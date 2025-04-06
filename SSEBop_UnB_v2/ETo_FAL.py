# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 18:45:33 2018
Author: Leandro Salles
"""
import numpy as np
import pandas as pd
import os

def estacao_FAL(path_met_FAL, arquivo_met_FAL, planilha_met_FAL):
    """
    Função que calcula variáveis meteorológicas e estimativa da evapotranspiração (ETo) 
    com base nos dados da estação meteorológica da FAL.

    Args:
        path_met_FAL (str): Caminho do diretório com os dados.
        arquivo_met_FAL (str): Nome do arquivo Excel com os dados.
        planilha_met_FAL (str): Nome da planilha dentro do Excel.

    Returns:
        pd.DataFrame: DataFrame com os dados processados.
    """
    
    # === LEITURA DOS DADOS ============================
    full_path = os.path.join(path_met_FAL, arquivo_met_FAL)
    DB = pd.read_excel(full_path, sheet_name=planilha_met_FAL, index_col='Data')
    
    DB['ano'] = DB.index.year
    DB['mes'] = DB.index.month
    DB['dia'] = DB.index.day

    # === CONSTANTES GEOGRÁFICAS E ATMOSFÉRICAS ========
    # Estação FAL (latitude/longitude decimais)
    lat = -1 * (15 + 59/60)
    lon = -1 * (47 + 62/60)
    Lat_rad = np.radians(lat)
    Altitude = 1030  # metros

    # Constantes físicas
    cp = 1013  # J/(kg.K)
    Gsc = 0.0820  # MJ m⁻² min⁻¹ constante solar
    SBC = 4.903e-9  # MJ K⁻⁴ m⁻² dia⁻¹ (Stefan-Boltzmann)
    Patm = 101.3 * ((293 - (0.0065 * Altitude)) / 293) ** 5.26  # kPa
    kRs = 0.16  # Coef. empírico para radiação
    cp_J = cp  # redefinido para dT

    # === CÁLCULO DAS VARIÁVEIS BASE ====================
    dados = pd.DataFrame(dtype=np.float64)
    dados[['ano', 'mes', 'dia']] = DB[['ano', 'mes', 'dia']]
    dados['P'] = DB["Precipitacao (mm)"]
    dados['Tmin'] = DB["T min (oC)"]
    dados['Tmax'] = DB["T max (oC)"]
    dados['UR'] = DB["UR med (%)"]
    dados['URMAX'] = DB["UR max (%)"]
    dados['URMIN'] = DB["UR min (%)"]
    dados['Vento'] = DB["Vento med (m/s)"]
    dados['Tmed'] = (dados['Tmax'] + dados['Tmin']) / 2
    dados['Tmed_p'] = DB.get("T med (oC)", dados['Tmed'])  # alternativa
    dados['Tmin_K'] = dados['Tmin'] + 273.16
    dados['Tmax_K'] = dados['Tmax'] + 273.16
    dados['DOY'] = dados.index.dayofyear
    dados['URmed'] = (dados['URMAX'] + dados['URMIN']) / 2

    # === PRESSÕES DE VAPOR =============================
    dados['ea_Tmed'] = 0.6108 * np.exp((17.27 * dados['Tmed']) / (dados['Tmed'] + 237.3))
    dados['ea_Tmin'] = 0.6108 * np.exp((17.27 * dados['Tmin']) / (dados['Tmin'] + 237.3))
    dados['ea_Tmax'] = 0.6108 * np.exp((17.27 * dados['Tmax']) / (dados['Tmax'] + 237.3))
    dados['es'] = (dados['ea_Tmin'] + dados['ea_Tmax']) / 2
    dados['ea_URmed'] = dados['es'] * (dados['UR'] / 100)
    dados['ea'] = (
        (dados['ea_Tmin'] * (dados['URMAX'] / 100)) +
        (dados['ea_Tmax'] * (dados['URMIN'] / 100))
    ) / 2
    dados['ea_defict'] = dados['es'] - dados['ea']

    # === RADIÂNCIA SOLAR ===============================
    dados['dr'] = 1 + 0.033 * np.cos((2 * np.pi / 365) * dados['DOY'])
    dados['Sol_decl'] = 0.409 * np.sin(((2 * np.pi / 365) * dados['DOY']) - 1.39)
    dados['Ws'] = np.arccos(-np.tan(Lat_rad) * np.tan(dados['Sol_decl']))
    
    dados['Ra'] = ((24 * 60) / np.pi) * Gsc * dados['dr'] * (
        dados['Ws'] * np.sin(Lat_rad) * np.sin(dados['Sol_decl']) +
        np.cos(Lat_rad) * np.cos(dados['Sol_decl']) * np.sin(dados['Ws'])
    )
    
    dados['Rs'] = kRs * np.sqrt(dados['Tmax'] - dados['Tmin']) * dados['Ra']
    dados['Rso'] = (0.75 + (2e-5 * Altitude)) * dados['Ra']
    dados['Rns'] = 0.77 * dados['Rs']  # (1 - albedo)
    dados['Rnso'] = 0.77 * dados['Rso']

    # === RADIAÇÃO DE ONDA LONGA ========================
    dados['Rnl'] = SBC * ((dados['Tmax_K']**4 + dados['Tmin_K']**4) / 2) * \
        (0.34 - 0.14 * np.sqrt(dados['ea_Tmin'])) * \
        (1.35 * (dados['Rs'] / dados['Rso']) - 0.35)
        
    dados['Rnlo'] = SBC * ((dados['Tmax_K']**4 + dados['Tmin_K']**4) / 2) * \
        (0.34 - 0.14 * np.sqrt(dados['ea_Tmin'])) * \
        (1.35 * (dados['Rso'] / dados['Rso']) - 0.35)

    dados['Rn'] = dados['Rns'] - dados['Rnl']
    dados['Rno'] = dados['Rnso'] - dados['Rnlo']
    dados['Rn_Wm-2'] = dados['Rn'] / 0.0864

    # === ETO PENMAN-MONTEITH ===========================
    dados['λ'] = 2.501 - 0.002361 * dados['Tmed']
    dados['γ'] = Patm * 0.665 * 1e-3
    dados['delta'] = 4098 * dados['ea_Tmed'] / (dados['Tmed'] + 237.3) ** 2
    dados['G'] = 0  # fluxo de calor no solo

    # Usando radiação global da planilha, se existir
    dados['Rs_EToPM'] = DB["Rad. Global (MJ/m2 d)"]
    dados['Rns_EToPM'] = 0.77 * dados['Rs_EToPM']

    dados['Rnl_EToPM'] = SBC * ((dados['Tmax_K']**4 + dados['Tmin_K']**4) / 2) * \
        (0.34 - 0.14 * np.sqrt(dados['ea'])) * \
        (1.35 * (dados['Rs_EToPM'] / dados['Rso']) - 0.35)

    dados['Rn_EToPM'] = dados['Rns_EToPM'] - dados['Rnl_EToPM']

    dados['ETo'] = (
        (0.408 * dados['delta'] * (dados['Rn_EToPM'] - dados['G'])) +
        (dados['γ'] * (900 / (dados['Tmed'] + 273)) * dados['Vento'] * (dados['es'] - dados['ea']))
    ) / (
        dados['delta'] + dados['γ'] * (1 + 0.34 * dados['Vento'])
    )

    # Para cultura de alfafa
    dados['ETo_alfafa'] = (
        (0.408 * dados['delta'] * (dados['Rn_EToPM'] - dados['G'])) +
        (dados['γ'] * (1600 / (dados['Tmed'] + 273)) * dados['Vento'] * (dados['es'] - dados['ea']))
    ) / (
        dados['delta'] + dados['γ'] * (1 + 0.38 * dados['Vento'])
    )

    # === ΔT (SSEBop) ===================================
    dados['Tkv'] = 1.01 * (dados['Tmed'] + 273.16)
    dados['pa'] = (1000 * Patm) / (dados['Tkv'] * 287)

    dados['dT'] = ((dados['Rno'] / 0.0864) * 110) / (cp_J * dados['pa'])
    dados['dT_v2013'] = dados['dT']
    dados['dT_v2018'] = ((dados['Rn'] / 0.0864) * 110) / (cp_J * dados['pa'])

    return dados
