# -*- coding: utf-8 -*-
"""
config.py — Configuração centralizada do pipeline PM2.5
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

Histórico residencial:
  Cada residência aceita date_start / date_end (formato "AAAA-MM-DD",
  preferido) ou year_start / year_end (anos inteiros, compat retroativa).
  Veja a documentação em cumulative.py para detalhes.

  PAC_02 demonstra a feature de fracionamento: muda de Campinas para SP
  em 2017-07-01, e o ano 2017 é dividido proporcionalmente nas duas
  residências (~0.496 em Campinas + ~0.504 em SP, soma = 1.0).
"""
import os

# ============================================================
# CAMINHOS
# ============================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

CETESB_ANNUAL_CSV = os.path.join(BASE_DIR, "dados_pm25_cetesb",
                                  "RESUMO_MEDIAS_ANUAIS_PM25.csv")
CETESB_STATIONS_CSV = os.path.join(BASE_DIR, "dados_pm25_cetesb",
                                    "cetesb_estacoes_coordenadas.csv")
CAMS_DIR = os.path.join(BASE_DIR, "CAMS")
MERRA2_DIR = os.path.join(BASE_DIR, "MERRA2")
BDQUEIMADAS_DIR = os.path.join(BASE_DIR, "BDQueimadas")
MESORREGIOES_SHP = os.path.join(BASE_DIR, "SP_Mesorregioes_2022",
                                 "SP_Mesorregioes_2022.shp")
FROTA_CSV = os.path.join(BASE_DIR, "frota_ativa_sp.csv")

OUTPUT_DIR = os.path.join(BASE_DIR, "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# PARÂMETROS DE ANÁLISE
# ============================================================
BUFFER_KM = [5, 10, 25]           # raios para análise de sensibilidade
CAMS_CONVERSION = 1e9             # kg/m³ → µg/m³
SP_LAT_RANGE = (-26.0, -19.0)
SP_LON_RANGE = (-54.0, -44.0)

# Fórmula NASA MERRA-2 para PM2.5 (fatores de conversão)
MERRA2_FORMULA = {
    "DUSMASS25": 1.0,     # dust PM2.5
    "SSSMASS25": 1.0,     # sea salt PM2.5
    "BCSMASS":   1.0,     # black carbon
    "OCSMASS":   1.4,     # organic carbon → organic matter
    "SO4SMASS":  1.375,   # sulfate → ammonium sulfate
}

# ============================================================
# PACIENTES-MODELO
# ============================================================
# Cada residência pode usar:
#   - date_start / date_end ("AAAA-MM-DD"): ponderação diária precisa
#   - year_start / year_end (int): equivale a 01/jan e 31/dez do ano
TEST_PATIENTS = [
    {
        "id": "PAC_01",
        "name": "Paciente SP – Av. Paulista",
        "residences": [
            {"cep": "01310-100", "lat": -23.5636, "lon": -46.6544,
             "city": "São Paulo",
             "date_start": "2015-01-01", "date_end": "2024-12-31"},
        ],
    },
    {
        "id": "PAC_02",
        "name": "Paciente Campinas → SP (multi-endereço, mudança em meio de ano)",
        "residences": [
            # Campinas: 2012-01-01 a 2017-06-30 → 2017 conta com fração 181/365
            {"cep": "13083-970", "lat": -22.8168, "lon": -47.0688,
             "city": "Campinas",
             "date_start": "2012-01-01", "date_end": "2017-06-30"},
            # SP: 2017-07-01 a 2024-12-31 → 2017 conta com fração 184/365
            {"cep": "01310-100", "lat": -23.5636, "lon": -46.6544,
             "city": "São Paulo",
             "date_start": "2017-07-01", "date_end": "2024-12-31"},
        ],
    },
    {
        "id": "PAC_03",
        "name": "Paciente S.J. Rio Preto",
        "residences": [
            {"cep": "15015-100", "lat": -20.8113, "lon": -49.3758,
             "city": "S.J. Rio Preto",
             "date_start": "2014-01-01", "date_end": "2024-12-31"},
        ],
    },
    {
        "id": "PAC_04",
        "name": "Paciente Itaí – Rural/Queimadas",
        "residences": [
            # CEP 18730-003: Avenida Santo Antônio, Centro, Itaí/SP
            {"cep": "18730-003", "lat": -23.4183, "lon": -49.0917,
             "city": "Itaí",
             "date_start": "2010-01-01", "date_end": "2024-12-31"},
        ],
    },
    {
        "id": "PAC_05",
        "name": "Paciente Santos – Costeiro",
        "residences": [
            {"cep": "11013-100", "lat": -23.9618, "lon": -46.3322,
             "city": "Santos",
             "date_start": "2016-01-01", "date_end": "2024-12-31"},
        ],
    },
]
