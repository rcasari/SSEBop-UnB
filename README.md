# SSEBop-UnB
Python Script created for SSEBop calculation using satellite images and  UAV Images

📝 Papers:

https://doi.org/10.1590/0102-77863910007

Authorship:

Leandro Salles

Modified by Raphael Casari



📌 SSEBop Evapotranspiration Processor

Script Python para cálculo de evapotranspiração real (ETr) utilizando o modelo SSEBop (Operational Simplified Surface Energy Balance).

🔍 O que este código faz?

Este script automatiza o processamento de:

Dados meteorológicos (temperatura, ETo)
Imagens de satélite (NDVI e LST - Land Surface Temperature)
Cálculo da evapotranspiração real usando o algoritmo SSEBop
Geração de mapas de saída em formato GeoTIFF

📊 Fluxo de Processamento

Entrada de dados meteorológicos (planilha Excel)
Carregamento de imagens NDVI e LST (GeoTIFF)
Pré-processamento (alinhamento espacial, máscaras)
Cálculo do fator de correção (c) baseado em estatísticas
Aplicação do modelo SSEBop para cálculo de ETf e ETr
Geração de mapas de saída (ETf e ETr em GeoTIFF)

⚙️ Tecnologias Utilizadas

Python 3.x
Bibliotecas principais:
NumPy
GDAL
xarray
rasterio
matplotlib

📂 Estrutura de Arquivos

/Processamento-evapo

├── /Meteorologia

│ └── Dadosmeteorologicos2019_atualizado.xlsx

├── /Images

│ ├── /NDVI

│ └── /Thermal

└── SSEBop_processor.py

📋 Pré-requisitos

Dados meteorológicos no formato especificado
Imagens NDVI e LST (Thermal) no mesmo sistema de coordenadas
Bibliotecas Python listadas em requirements.txt

🚀 Como Usar

Configure os paths no início do script
Execute o arquivo principal SSEBop_processor.py
Os resultados serão salvos no mesmo diretório das imagens LST

📄 Saídas Geradas

ETa_YYYY-MM-DD-Thermal.tif (Evapotranspiração real)
ETf_YYYY-MM-DD-Thermal.tif (Fração de evapotranspiração)

📚 Referência Científica

Baseado no modelo SSEBop descrito em:

Senay, G.B., et al. (2013). "Operational Evapotranspiration Mapping Using Remote Sensing and Weather Datasets: A New Parameterization for the SSEB Approach". Journal of the American Water Resources Association.

📝 Notas

Desenvolvido originalmente por Leandro Salles
Modificado e aprimorado por Raphael Casari
Data da última modificação: Agosto 2023
Sugestões adicionais:
Badges (opcional - adicione no topo do README):

requirements.txt (adicione um arquivo com as dependências):

numpy

gdal

xarray

rasterio

matplotlib

openpyxl
