# -*- coding: utf-8 -*-
"""
data_loaders.py — Carregamento e padronização de todas as fontes de dados
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP
"""
import os
import glob
import zipfile
import tempfile

import pandas as pd
import geopandas as gpd
import numpy as np
import xarray as xr
from shapely.geometry import Point

import config


def load_cetesb_annual() -> gpd.GeoDataFrame:
    """
    Carrega médias anuais de PM2.5 das estações CETESB.
    Returns: GeoDataFrame com colunas:
        estacao_nome, estacao_codigo, ano, media_anual_pm25_ugm3, geometry (Point)
    """
    print("📡 Carregando CETESB médias anuais...")
    df = pd.read_csv(config.CETESB_ANNUAL_CSV)
    # Remover linhas sem coordenadas
    df = df.dropna(subset=["lat", "lon", "media_anual_pm25_ugm3"])
    geometry = [Point(lon, lat) for lon, lat in zip(df["lon"], df["lat"])]
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")
    print(f"   ✅ {len(gdf)} registros, {gdf['estacao_nome'].nunique()} estações, "
          f"anos {gdf['ano'].min()}-{gdf['ano'].max()}")
    return gdf


def load_cetesb_stations() -> gpd.GeoDataFrame:
    """Carrega coordenadas de todas as estações CETESB."""
    df = pd.read_csv(config.CETESB_STATIONS_CSV)
    df = df.dropna(subset=["lat", "lon"])
    geometry = [Point(lon, lat) for lon, lat in zip(df["lon"], df["lat"])]
    return gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")


def _find_cams_file(year: int) -> str:
    """Encontra o arquivo CAMS NetCDF para um dado ano."""
    # Os arquivos CAMS são nomeados data_sfcYY.nc (08=2008, 09=2009, etc.)
    # data_sfc10.nc é 2010 (maior), outros são anos individuais
    yy = year - 2000 if year >= 2000 else year - 1900
    pattern = os.path.join(config.CAMS_DIR, f"data_sfc{yy:02d}.nc")
    if os.path.exists(pattern):
        return pattern
    # Tentar variantes
    candidates = glob.glob(os.path.join(config.CAMS_DIR, "data_sfc*.nc"))
    for c in candidates:
        basename = os.path.basename(c)
        num = int(basename.replace("data_sfc", "").replace(".nc", ""))
        if num == yy:
            return c
    return None


def load_cams_pm25(year: int) -> xr.DataArray:
    """
    Carrega dados CAMS EAC4 PM2.5 para um ano.
    Converte de kg/m³ para µg/m³ e retorna média anual.

    Returns: xarray DataArray com dimensões (latitude, longitude) em µg/m³
    """
    filepath = _find_cams_file(year)
    if filepath is None:
        return None

    ds = xr.open_dataset(filepath)
    # A variável PM2.5 no CAMS é 'pm2p5' (units: kg/m³)
    if "pm2p5" not in ds:
        print(f"   ⚠ Variável 'pm2p5' não encontrada no CAMS para {year}")
        ds.close()
        return None

    pm25 = ds["pm2p5"] * config.CAMS_CONVERSION  # → µg/m³
    # Média temporal (anual) — CAMS usa 'valid_time' como dim temporal
    time_dim = None
    for d in pm25.dims:
        if "time" in d.lower() or "valid" in d.lower():
            time_dim = d
            break

    if time_dim:
        pm25_annual = pm25.mean(dim=time_dim)
    else:
        pm25_annual = pm25

    pm25_annual.attrs["units"] = "µg/m³"
    pm25_annual.attrs["source"] = "CAMS EAC4"
    pm25_annual.attrs["year"] = year
    return pm25_annual


def _find_merra2_files(year: int) -> list:
    """Encontra todos os 12 arquivos MERRA-2 mensais de um ano."""
    pattern = os.path.join(config.MERRA2_DIR,
                           f"MERRA2_*.tavgM_2d_aer_Nx.{year}*.nc4.nc4")
    files = sorted(glob.glob(pattern))
    # Também procurar variante _401 (meses com reprocessamento)
    return files


def load_merra2_pm25(year: int) -> xr.DataArray:
    """
    Carrega dados MERRA-2 e calcula PM2.5 usando fórmula NASA:
    PM2.5 = DUSMASS25 + SSSMASS25 + BCSMASS + 1.4×OCSMASS + 1.375×SO4SMASS

    Returns: xarray DataArray com dimensões (lat, lon) em µg/m³
    """
    files = _find_merra2_files(year)
    if not files:
        return None

    monthly_pm25 = []
    for f in files:
        ds = xr.open_dataset(f)
        pm25 = None
        for var_name, factor in config.MERRA2_FORMULA.items():
            if var_name in ds:
                component = ds[var_name].squeeze() * factor
                if pm25 is None:
                    pm25 = component
                else:
                    pm25 = pm25 + component
        if pm25 is not None:
            # MERRA-2 já está em µg/m³ (surface mass concentration)
            # Mas precisamos confirmar — as variáveis SMASS são kg/m³? Não!
            # MERRA-2 tavgM_2d_aer_Nx: unidades são kg m-2 (column mass)...
            # Na verdade, *SMASS são "surface mass concentration" em kg/m³
            # Converter para µg/m³:
            pm25_ugm3 = pm25 * 1e9
            monthly_pm25.append(pm25_ugm3)
        ds.close()

    if not monthly_pm25:
        return None

    # Média anual dos 12 (ou menos) meses
    annual = sum(monthly_pm25) / len(monthly_pm25)
    annual.attrs["units"] = "µg/m³"
    annual.attrs["source"] = "MERRA-2"
    annual.attrs["year"] = year
    annual.attrs["n_months"] = len(monthly_pm25)
    return annual


def load_fire_foci(year: int) -> gpd.GeoDataFrame:
    """
    Carrega focos de calor do BDQueimadas para um ano.
    Descompacta o ZIP e lê o CSV contido.

    Returns: GeoDataFrame com localização dos focos de calor
    """
    zip_path = os.path.join(config.BDQUEIMADAS_DIR, f"focos_br_sp_ref_{year}.zip")
    if not os.path.exists(zip_path):
        return None

    tmpdir = tempfile.mkdtemp()
    try:
        with zipfile.ZipFile(zip_path, 'r') as z:
            z.extractall(tmpdir)

        # Encontrar arquivo CSV extraído
        csv_files = glob.glob(os.path.join(tmpdir, "*.csv"))
        if not csv_files:
            # Tentar subdiretórios
            csv_files = glob.glob(os.path.join(tmpdir, "**", "*.csv"), recursive=True)

        if not csv_files:
            return None

        # Ler o CSV de focos
        df = pd.read_csv(csv_files[0], encoding="latin-1")

        # Os campos de coordenadas podem ser 'lat'/'lon' ou 'latitude'/'longitude'
        lat_col = None
        lon_col = None
        for c in df.columns:
            cl = c.lower().strip()
            if cl in ("lat", "latitude"):
                lat_col = c
            elif cl in ("lon", "longitude"):
                lon_col = c

        if lat_col is None or lon_col is None:
            print(f"   ⚠ Colunas de coordenadas não encontradas em BDQueimadas {year}")
            print(f"     Colunas disponíveis: {list(df.columns)}")
            return None

        df = df.dropna(subset=[lat_col, lon_col])
        geometry = [Point(lon, lat) for lon, lat in zip(df[lon_col], df[lat_col])]
        gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")
        return gdf

    except Exception as e:
        print(f"   ⚠ Erro ao carregar BDQueimadas {year}: {e}")
        return None


def load_mesorregioes() -> gpd.GeoDataFrame:
    """Carrega shapefile das Mesorregiões de SP."""
    print("🗺️  Carregando Mesorregiões SP...")
    gdf = gpd.read_file(config.MESORREGIOES_SHP)
    print(f"   ✅ {len(gdf)} mesorregiões carregadas")
    return gdf


# ============================================================
# Cache global para evitar recarregamentos
# ============================================================
class DataCache:
    """Cache simples para evitar recarregar dados pesados."""

    def __init__(self):
        self.cetesb = None
        self.cams = {}      # year -> DataArray
        self.merra2 = {}    # year -> DataArray
        self.fires = {}     # year -> GeoDataFrame
        self.mesorregioes = None

    def get_cetesb(self):
        if self.cetesb is None:
            self.cetesb = load_cetesb_annual()
        return self.cetesb

    def get_cams(self, year):
        if year not in self.cams:
            self.cams[year] = load_cams_pm25(year)
        return self.cams[year]

    def get_merra2(self, year):
        if year not in self.merra2:
            self.merra2[year] = load_merra2_pm25(year)
        return self.merra2[year]

    def get_fires(self, year):
        if year not in self.fires:
            self.fires[year] = load_fire_foci(year)
        return self.fires[year]

    def get_mesorregioes(self):
        if self.mesorregioes is None:
            self.mesorregioes = load_mesorregioes()
        return self.mesorregioes
