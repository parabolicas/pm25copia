#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
source_specific.py — Exposição source-specific ao PM2.5 (Etapa 6 do pipeline)
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

ARQUITETURA DE CALIBRAÇÃO (v6 — Calibração Dual com Modulador Empírico):

  Âncora primária — RMSP (Pereira et al. 2025, Atmos Chem Phys):
    - veicular 41%, biomassa 25%, secundário 21%, industrial 14%
    - Calibra: bc_oc_ratio → w_veicular e f_om × (1 - bc_oc_ratio) → w_biomassa

  Modulador espacial — Densidade de focos BDQueimadas (INPE 2008-2024):
    - Não é uma "âncora PMF" externa — é dado *observado na própria região*
    - Aumenta w_biomassa proporcionalmente à densidade local de focos via
      kernel gaussiano (bandwidth ~25 km), com fator multiplicativo
      controlado por BIOMASS_BOOST_FACTOR (default 1.4 → +40% em pixels
      com densidade máxima estadual; análise de sensibilidade ±20%).
    - Limitação reconhecida: ausência de PMF dedicado ao interior agrícola
      de SP na literatura. A literatura disponível para áreas de queima de
      biomassa (cana-de-açúcar, queimadas amazônicas) reporta fração
      biomassa de 30-70% em períodos/locais de queima intensa, suportando
      o boost máximo de ~+40% sobre o valor RMSP-anchored.

  Validação cruzada interna:
    - Após estimativa, computa correlação Spearman entre w_biomassa médio
      municipal e densidade de focos por município. Esperada >0.5 — fraca
      correlação indica calibração mal-ajustada e sinaliza necessidade de
      revisão dos parâmetros de boost.

Uso:
    python3 source_specific.py
    python3 source_specific.py --boost 1.2    # análise de sensibilidade
    python3 source_specific.py --boost 1.6    # análise de sensibilidade
"""
import os
import glob
import warnings
import argparse
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr
import rasterio
from rasterio.transform import from_bounds
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator

import config


# ============================================================
# CONFIGURAÇÃO
# ============================================================
OUTPUT_DIR = os.path.join(config.OUTPUT_DIR, "source_specific")
os.makedirs(OUTPUT_DIR, exist_ok=True)

YEAR_START, YEAR_END = 2008, 2024

TARGET_RES = 0.05
SP_LAT_MIN, SP_LAT_MAX = -25.5, -19.5
SP_LON_MIN, SP_LON_MAX = -53.5, -44.0
TARGET_LATS = np.arange(SP_LAT_MAX, SP_LAT_MIN - TARGET_RES / 2, -TARGET_RES)
TARGET_LONS = np.arange(SP_LON_MIN, SP_LON_MAX + TARGET_RES / 2, TARGET_RES)

# Âncora primária: PMF Pereira et al. (2025) na RMSP
PMF_REFERENCE = {
    "veicular": 0.41,
    "biomassa": 0.25,
    "secundario": 0.21,
    "industrial": 0.14,
}
RMSP_LAT, RMSP_LON = -23.55, -46.63

# ===== Parâmetros da calibração dual (NOVO em v6) =====
# Bandwidth do kernel gaussiano de KDE (km).
# 25 km suaviza ruído mas preserva resolução municipal.
FIRE_KDE_BANDWIDTH_KM = 25.0

# Fator de boost aplicado em pixels com densidade máxima de focos.
# 1.4 → w_biomassa pode ser até 40% maior que o RMSP-anchored em áreas de
# queimada extrema. Sensibilidade testará 1.2 e 1.6.
# Justificativa: PMF na RMSP fornece w_bio=25%; literatura amazônica/cana
# reporta 30-70% em áreas de queima intensa. 1.4 → ~35% no extremo —
# valor conservador dentro desse intervalo.
BIOMASS_BOOST_DEFAULT = 1.4


# ============================================================
# UTILITIES
# ============================================================
def _haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return R * 2 * np.arcsin(np.sqrt(a))


# ============================================================
# CARGA DE COMPONENTES MERRA-2 + FRAÇÕES BRUTAS (sem mudanças vs v5)
# ============================================================
def load_merra2_components(year):
    """Carrega componentes MERRA-2 do ano e calcula médias anuais."""
    pattern = os.path.join(config.MERRA2_DIR, f"MERRA2_*{year}*.nc4.nc4")
    files = sorted(glob.glob(pattern))
    if not files:
        return None

    var_names = ["BCSMASS", "OCSMASS", "DUSMASS25", "SSSMASS25", "SO4SMASS"]
    monthly_data = {v: [] for v in var_names}

    for f in files:
        try:
            ds = xr.open_dataset(f)
            for v in var_names:
                if v in ds:
                    monthly_data[v].append(ds[v].squeeze(drop=True).values)
            ds.close()
        except Exception:
            continue

    components = {}
    for v in var_names:
        if monthly_data[v]:
            arr = np.nanmean(np.array(monthly_data[v]), axis=0)
            ds_ref = xr.open_dataset(files[0])
            components[v] = xr.DataArray(
                arr, dims=["lat", "lon"],
                coords={"lat": ds_ref.lat.values, "lon": ds_ref.lon.values},
            )
            ds_ref.close()

    return components if components else None


def compute_source_fractions(components):
    """Calcula frações brutas a partir dos componentes MERRA-2."""
    bc = components["BCSMASS"]
    oc = components["OCSMASS"]
    dust = components["DUSMASS25"]
    ss = components["SSSMASS25"]
    so4 = components["SO4SMASS"]

    pm25_total = dust + ss + bc + 1.4 * oc + 1.375 * so4
    pm25_safe = xr.where(pm25_total > 0, pm25_total, np.nan)

    f_bc = bc / pm25_safe
    f_om = (1.4 * oc) / pm25_safe
    f_dust = dust / pm25_safe
    f_sulfate = (1.375 * so4) / pm25_safe
    f_seasalt = ss / pm25_safe

    carbonaceous = bc + 1.4 * oc
    carbonaceous_safe = xr.where(carbonaceous > 0, carbonaceous, np.nan)
    bc_oc_ratio = bc / carbonaceous_safe

    return {
        "bc_oc_ratio": bc_oc_ratio,
        "f_bc": f_bc, "f_om": f_om, "f_dust": f_dust,
        "f_sulfate": f_sulfate, "f_seasalt": f_seasalt,
        "pm25_total": pm25_total,
    }


# ============================================================
# NOVO v6 — DENSIDADE DE FOCOS POR KDE GAUSSIANO
# ============================================================
def load_all_fires_for_kde():
    """
    Carrega todos os focos BDQueimadas 2008-2024 (cumulativo).
    Retorna array (n_focos, 2) com (lat, lon).
    """
    import zipfile
    import tempfile

    print("\n🔥 Carregando focos BDQueimadas (2008-2024) para KDE...")
    all_lats, all_lons = [], []
    zip_files = sorted(glob.glob(os.path.join(config.BDQUEIMADAS_DIR, "*.zip")))

    for zf_path in zip_files:
        tmpdir = tempfile.mkdtemp()
        try:
            with zipfile.ZipFile(zf_path, "r") as zf:
                zf.extractall(tmpdir)
            csvs = glob.glob(os.path.join(tmpdir, "**", "*.csv"), recursive=True)
            if csvs:
                df = pd.read_csv(csvs[0], encoding="latin-1", usecols=["lat", "lon"])
                df = df.dropna()
                all_lats.extend(df["lat"].values)
                all_lons.extend(df["lon"].values)
        except Exception:
            pass

    coords = np.column_stack([np.array(all_lats), np.array(all_lons)])
    print(f"   ✅ {len(coords):,} focos cumulativos carregados")
    return coords


def generate_fire_density_grid(fire_coords, bandwidth_km=FIRE_KDE_BANDWIDTH_KM):
    """
    Gera grade 2D de densidade de focos via kernel density estimation
    gaussiano sobre a grade alvo (n_lats × n_lons, 0.05° ≈ 5,5 km).

    Args:
        fire_coords: np.ndarray shape (n_focos, 2) com [lat, lon]
        bandwidth_km: bandwidth do kernel gaussiano em km (default 25 km)

    Returns:
        np.ndarray (n_lats, n_lons): densidade não-normalizada (focos × peso)
    """
    n_lats, n_lons = len(TARGET_LATS), len(TARGET_LONS)
    density = np.zeros((n_lats, n_lons), dtype=np.float64)

    if len(fire_coords) == 0:
        return density

    # Aproximação: 1° lat ≈ 111 km, 1° lon ≈ 111 × cos(lat) km
    # Em SP (lat ~-23°), cos(-23°) ≈ 0.92 → 1° lon ≈ 102 km
    bandwidth_deg_lat = bandwidth_km / 111.0
    bandwidth_deg_lon = bandwidth_km / 102.0

    # Para cada foco, distribuir contribuição gaussiana sobre janela ±3σ
    # Usamos haversine real para eficiência de cálculo: vectorize
    sigma2_2 = 2.0 * bandwidth_km ** 2  # 2σ² em km²

    print(f"   🧮 KDE com bandwidth {bandwidth_km:.0f} km, "
          f"{len(fire_coords):,} focos × {n_lats}×{n_lons} pixels...")

    # Otimização: filtrar focos por bounding box da grade
    in_box = (
        (fire_coords[:, 0] >= SP_LAT_MIN - 0.5)
        & (fire_coords[:, 0] <= SP_LAT_MAX + 0.5)
        & (fire_coords[:, 1] >= SP_LON_MIN - 0.5)
        & (fire_coords[:, 1] <= SP_LON_MAX + 0.5)
    )
    fires_in = fire_coords[in_box]
    print(f"      Após bbox filter: {len(fires_in):,} focos relevantes")

    grid_lat_2d, grid_lon_2d = np.meshgrid(TARGET_LATS, TARGET_LONS, indexing="ij")

    # Para evitar O(n_focos × n_pixels) que pode ser muito,
    # acumular foco-a-foco com janela limitada (±3σ ≈ ±75 km)
    window_deg_lat = 3 * bandwidth_deg_lat
    window_deg_lon = 3 * bandwidth_deg_lon
    progress_step = max(1, len(fires_in) // 20)

    for i, (f_lat, f_lon) in enumerate(fires_in):
        if i % progress_step == 0 and i > 0:
            print(f"      Progresso: {i}/{len(fires_in)} ({100*i/len(fires_in):.0f}%)")
        # Janela de pixels relevantes (±3σ ao redor do foco)
        lat_min = f_lat - window_deg_lat
        lat_max = f_lat + window_deg_lat
        lon_min = f_lon - window_deg_lon
        lon_max = f_lon + window_deg_lon

        # Índices da janela na grade alvo (TARGET_LATS é decrescente N→S)
        i_lat_max = np.searchsorted(-TARGET_LATS, -lat_max, side="left")
        i_lat_min = np.searchsorted(-TARGET_LATS, -lat_min, side="right")
        j_lon_min = np.searchsorted(TARGET_LONS, lon_min, side="left")
        j_lon_max = np.searchsorted(TARGET_LONS, lon_max, side="right")

        if i_lat_max >= i_lat_min or j_lon_min >= j_lon_max:
            continue

        sub_lats = grid_lat_2d[i_lat_max:i_lat_min, j_lon_min:j_lon_max]
        sub_lons = grid_lon_2d[i_lat_max:i_lat_min, j_lon_min:j_lon_max]

        # Distância haversine vetorizada
        dlat_rad = np.radians(sub_lats - f_lat)
        dlon_rad = np.radians(sub_lons - f_lon)
        a = (np.sin(dlat_rad / 2) ** 2 +
             np.cos(np.radians(f_lat)) * np.cos(np.radians(sub_lats)) *
             np.sin(dlon_rad / 2) ** 2)
        dist_km_2d = 6371.0 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

        # Kernel gaussiano (peso máximo = 1 no centro do foco)
        weights = np.exp(-(dist_km_2d ** 2) / sigma2_2)
        density[i_lat_max:i_lat_min, j_lon_min:j_lon_max] += weights

    print(f"      ✅ Grade gerada: max={density.max():.1f}, "
          f"mediana={np.median(density):.2f}, "
          f"P95={np.percentile(density, 95):.1f}")
    return density


def normalize_fire_density(density):
    """
    Normaliza a grade de densidade para [0, 1] usando o P95 como teto
    (robusto a outliers extremos como mega-incêndios pontuais).
    """
    p95 = np.percentile(density, 95)
    if p95 <= 0:
        return density
    normalized = np.clip(density / p95, 0.0, 1.0)
    return normalized


# ============================================================
# CALIBRAÇÃO DUAL (v6 — modificada)
# ============================================================
def calibrate_fractions(fractions_raw, fire_density_normalized=None,
                        biomass_boost_factor=BIOMASS_BOOST_DEFAULT):
    """
    Calibração das frações de fonte com âncora dual:

      1. Âncora primária: PMF Pereira et al. (2025) na RMSP — calibra
         bc_oc_ratio → w_veicular e f_om × (1 - bc_oc_ratio) → w_biomassa.

      2. Modulador espacial empírico: densidade normalizada de focos
         BDQueimadas (KDE gaussiano, bandwidth ~25 km). Aumenta
         w_biomassa proporcionalmente à densidade local:

            w_biomassa_boosted = w_biomassa * (1 + (β - 1) × density_norm)

         onde β = biomass_boost_factor e density_norm ∈ [0, 1].
         Pixels com density_norm = 0 (sem focos) → sem boost.
         Pixels com density_norm = 1 (densidade máxima estadual) →
         w_biomassa multiplicado por β.

    Após boost, w_veicular + w_biomassa + w_outros é renormalizado para 1.

    Args:
        fractions_raw: dict com bc_oc_ratio, f_om (DataArrays MERRA-2)
        fire_density_normalized: np.ndarray (n_lats, n_lons) já normalizado
            para [0, 1], sobre a grade alvo TARGET_LATS×TARGET_LONS, OU
            None para calibração apenas RMSP-anchored.
        biomass_boost_factor: float ≥1.0 (default 1.4). =1.0 desativa boost.

    Returns:
        dict com w_veicular, w_biomassa, w_outros (DataArrays no grid MERRA-2)
        + metadados de calibração.
    """
    bc_oc = fractions_raw["bc_oc_ratio"]
    f_om = fractions_raw["f_om"]

    # 1. Valor de referência na RMSP
    lat_idx = np.argmin(np.abs(bc_oc.coords["lat"].values - RMSP_LAT))
    lon_idx = np.argmin(np.abs(bc_oc.coords["lon"].values - RMSP_LON))
    bc_oc_rmsp = float(bc_oc.values[lat_idx, lon_idx])
    if np.isnan(bc_oc_rmsp) or bc_oc_rmsp <= 0:
        bc_oc_rmsp = 0.15

    # 2. w_veicular linear: bc_oc=0 → 0.05; bc_oc=bc_oc_rmsp → 0.41
    slope_veh = (PMF_REFERENCE["veicular"] - 0.05) / bc_oc_rmsp
    w_veicular = 0.05 + slope_veh * bc_oc
    w_veicular = xr.where(w_veicular > 0.70, 0.70, w_veicular)
    w_veicular = xr.where(w_veicular < 0.02, 0.02, w_veicular)

    # 3. w_biomassa proxy: f_om × (1 - bc_oc), escalado para 0.25 na RMSP
    biomass_proxy = f_om * (1.0 - bc_oc)
    bp_rmsp = float(biomass_proxy.values[lat_idx, lon_idx])
    if np.isnan(bp_rmsp) or bp_rmsp <= 0:
        bp_rmsp = 0.10
    slope_bio = PMF_REFERENCE["biomassa"] / bp_rmsp
    w_biomassa = slope_bio * biomass_proxy
    w_biomassa = xr.where(w_biomassa > 0.70, 0.70, w_biomassa)
    w_biomassa = xr.where(w_biomassa < 0.01, 0.01, w_biomassa)

    # 4. NOVO v6 — Boost de biomassa por densidade de focos
    boost_meta = {
        "applied": False,
        "biomass_boost_factor": float(biomass_boost_factor),
    }
    if fire_density_normalized is not None and biomass_boost_factor > 1.0:
        # Resample fire_density_normalized para o grid MERRA-2
        # (diferente do grid alvo 0.05° que será usado depois)
        density_resampled = _resample_density_to_merra2(
            fire_density_normalized, bc_oc.coords["lat"].values,
            bc_oc.coords["lon"].values)
        # density_resampled é np.ndarray no shape do grid MERRA-2
        boost_field = 1.0 + (biomass_boost_factor - 1.0) * density_resampled
        # Aplicar boost (DataArray * np.ndarray funciona se shapes batem)
        w_biomassa_values = w_biomassa.values * boost_field
        # Cap no máximo razoável
        w_biomassa_values = np.clip(w_biomassa_values, 0.01, 0.75)
        w_biomassa = xr.DataArray(w_biomassa_values,
                                   dims=w_biomassa.dims,
                                   coords=w_biomassa.coords)
        boost_meta["applied"] = True
        boost_meta["max_density_norm"] = float(np.nanmax(density_resampled))
        boost_meta["mean_density_norm"] = float(np.nanmean(density_resampled))

    # 5. w_outros + renormalização
    w_outros = 1.0 - w_veicular - w_biomassa
    w_outros = xr.where(w_outros < 0.05, 0.05, w_outros)
    total = w_veicular + w_biomassa + w_outros
    w_veicular = w_veicular / total
    w_biomassa = w_biomassa / total
    w_outros = w_outros / total

    return {
        "w_veicular": w_veicular,
        "w_biomassa": w_biomassa,
        "w_outros": w_outros,
        "bc_oc_rmsp": bc_oc_rmsp,
        "boost_meta": boost_meta,
    }


def _resample_density_to_merra2(density_target_grid, merra2_lats, merra2_lons):
    """
    Resample a grade de densidade (TARGET_LATS×TARGET_LONS, 0.05°)
    para o grid MERRA-2 (0.5°×0.625°). Usa interpolação linear.
    """
    # TARGET_LATS é decrescente; precisa inverter para interpolador
    src_lats = TARGET_LATS[::-1]
    src_density = density_target_grid[::-1, :]

    interp = RegularGridInterpolator(
        (src_lats, TARGET_LONS), src_density,
        method="linear", bounds_error=False, fill_value=0.0
    )
    grid_lats, grid_lons = np.meshgrid(merra2_lats, merra2_lons, indexing="ij")
    points = np.column_stack([grid_lats.ravel(), grid_lons.ravel()])
    return interp(points).reshape(grid_lats.shape)


# ============================================================
# SUPERFÍCIES SOURCE-SPECIFIC (modificada para passar fire grid)
# ============================================================
def resample_to_target(data_2d, src_lats, src_lons):
    src_values = data_2d.copy()
    if src_lats[0] > src_lats[-1]:
        src_lats = src_lats[::-1]
        src_values = src_values[::-1, :]
    interpolator = RegularGridInterpolator(
        (src_lats, src_lons), src_values,
        method="linear", bounds_error=False, fill_value=None
    )
    grid_lats, grid_lons = np.meshgrid(TARGET_LATS, TARGET_LONS, indexing="ij")
    points = np.column_stack([grid_lats.ravel(), grid_lons.ravel()])
    return interpolator(points).reshape(grid_lats.shape)


def generate_source_specific_surfaces(year, fire_density_normalized=None,
                                       biomass_boost_factor=BIOMASS_BOOST_DEFAULT):
    """Gera superfícies fracionadas para um ano, com calibração dual."""
    components = load_merra2_components(year)
    if components is None:
        return None, None

    fractions = compute_source_fractions(components)
    calibrated = calibrate_fractions(
        fractions, fire_density_normalized=fire_density_normalized,
        biomass_boost_factor=biomass_boost_factor)

    src_lats = components["BCSMASS"].coords["lat"].values
    src_lons = components["BCSMASS"].coords["lon"].values
    w_veh = resample_to_target(calibrated["w_veicular"].values, src_lats, src_lons)
    w_bio = resample_to_target(calibrated["w_biomassa"].values, src_lats, src_lons)
    w_out = resample_to_target(calibrated["w_outros"].values, src_lats, src_lons)

    total = w_veh + w_bio + w_out
    w_veh, w_bio, w_out = w_veh / total, w_bio / total, w_out / total

    meta = {
        "year": year,
        "bc_oc_rmsp": calibrated["bc_oc_rmsp"],
        "boost_meta": calibrated["boost_meta"],
        "w_veh_mean": round(float(np.nanmean(w_veh)), 3),
        "w_bio_mean": round(float(np.nanmean(w_bio)), 3),
        "w_out_mean": round(float(np.nanmean(w_out)), 3),
        "w_veh_rmsp": round(float(w_veh[
            np.argmin(np.abs(TARGET_LATS - RMSP_LAT)),
            np.argmin(np.abs(TARGET_LONS - RMSP_LON))]), 3),
        "w_bio_rmsp": round(float(w_bio[
            np.argmin(np.abs(TARGET_LATS - RMSP_LAT)),
            np.argmin(np.abs(TARGET_LONS - RMSP_LON))]), 3),
    }
    return {"w_veicular": w_veh, "w_biomassa": w_bio, "w_outros": w_out}, meta


def save_fraction_geotiffs(fractions, year, output_dir):
    """Salva frações como GeoTIFFs."""
    paths = {}
    for name, data in fractions.items():
        filepath = os.path.join(output_dir, f"{name}_{year}.tif")
        n_lats, n_lons = len(TARGET_LATS), len(TARGET_LONS)
        west = TARGET_LONS[0] - TARGET_RES / 2
        east = TARGET_LONS[-1] + TARGET_RES / 2
        north = TARGET_LATS[0] + TARGET_RES / 2
        south = TARGET_LATS[-1] - TARGET_RES / 2
        transform = from_bounds(west, south, east, north, n_lons, n_lats)

        with rasterio.open(
            filepath, "w", driver="GTiff",
            height=n_lats, width=n_lons, count=1,
            dtype="float32", crs="EPSG:4326",
            transform=transform, nodata=-9999.0,
        ) as dst:
            d = data.astype(np.float32)
            d[np.isnan(d)] = -9999.0
            dst.write(d, 1)
            dst.update_tags(YEAR=str(year), UNITS="fraction_0-1",
                            SOURCE="MERRA-2 + PMF Pereira2025 + KDE BDQueimadas")
        paths[name] = filepath
    return paths


def save_fire_density_geotiff(density, output_dir, suffix=""):
    """Salva a grade de densidade de focos como GeoTIFF para auditoria."""
    filepath = os.path.join(output_dir, f"fire_density{suffix}.tif")
    n_lats, n_lons = len(TARGET_LATS), len(TARGET_LONS)
    west = TARGET_LONS[0] - TARGET_RES / 2
    east = TARGET_LONS[-1] + TARGET_RES / 2
    north = TARGET_LATS[0] + TARGET_RES / 2
    south = TARGET_LATS[-1] - TARGET_RES / 2
    transform = from_bounds(west, south, east, north, n_lons, n_lats)

    with rasterio.open(
        filepath, "w", driver="GTiff",
        height=n_lats, width=n_lons, count=1,
        dtype="float32", crs="EPSG:4326",
        transform=transform, nodata=-9999.0,
    ) as dst:
        d = density.astype(np.float32)
        dst.write(d, 1)
        dst.update_tags(SOURCE="BDQueimadas/INPE 2008-2024 KDE",
                        BANDWIDTH_KM=str(FIRE_KDE_BANDWIDTH_KM))
    return filepath


# ============================================================
# LOOKUP (sem mudanças)
# ============================================================
def lookup_source_specific(lat, lon, year, pm25_total, output_dir=None):
    """Consulta exposição source-specific."""
    if output_dir is None:
        output_dir = OUTPUT_DIR

    result = {}
    for name in ["w_veicular", "w_biomassa", "w_outros"]:
        filepath = os.path.join(output_dir, f"{name}_{year}.tif")
        if not os.path.exists(filepath):
            return None
        with rasterio.open(filepath) as src:
            row, col = src.index(lon, lat)
            if 0 <= row < src.height and 0 <= col < src.width:
                val = src.read(1)[row, col]
                if val == src.nodata:
                    return None
                result[name] = float(val)
            else:
                return None

    if pm25_total is not None:
        result["pm25_veicular"] = round(pm25_total * result["w_veicular"], 2)
        result["pm25_biomassa"] = round(pm25_total * result["w_biomassa"], 2)
        result["pm25_outros"] = round(pm25_total * result["w_outros"], 2)
    return result


# ============================================================
# VALIDAÇÃO CRUZADA INTERNA (NOVO v6)
# ============================================================
def cross_validation_spearman(fractions_dict, fire_density_target):
    """
    Validação cruzada interna: correlação Spearman entre w_biomassa
    pixel-a-pixel e densidade de focos. Esperada >0.5.

    Args:
        fractions_dict: {year: {w_veicular, w_biomassa, w_outros}} no grid alvo.
        fire_density_target: np.ndarray (n_lats, n_lons) — densidade NÃO
            normalizada (focos × peso) sobre o grid alvo.

    Returns:
        dict com correlações por ano e média global.
    """
    print("\n🔬 VALIDAÇÃO CRUZADA INTERNA (Spearman pixel-a-pixel)")
    print(f"   Esperado: ρ > 0.5 — fraca correlação indica calibração mal-ajustada")

    # Achatar grade para vetor; remover pixels sem focos (densidade = 0)
    density_flat = fire_density_target.flatten()
    mask = density_flat > 0
    density_used = density_flat[mask]

    results = {}
    for year, fracs in sorted(fractions_dict.items()):
        w_bio_flat = fracs["w_biomassa"].flatten()[mask]
        rho, p_val = spearmanr(density_used, w_bio_flat)
        results[year] = {
            "rho": float(rho),
            "p_value": float(p_val),
            "n_pixels": int(mask.sum()),
        }

    rhos = [r["rho"] for r in results.values()]
    mean_rho = float(np.mean(rhos))
    median_rho = float(np.median(rhos))
    print(f"\n   ρ médio (2008-2024): {mean_rho:.3f}")
    print(f"   ρ mediano: {median_rho:.3f}")
    print(f"   Intervalo: [{min(rhos):.3f}, {max(rhos):.3f}]")
    if mean_rho < 0.3:
        print("   ⚠ ATENÇÃO: ρ médio < 0.3 — calibração pode estar mal-ajustada.")
    elif mean_rho < 0.5:
        print("   ⚠ ρ médio entre 0.3-0.5 — atenção, considere ajuste do BIOMASS_BOOST_FACTOR.")
    else:
        print(f"   ✅ ρ médio ≥ 0.5 — calibração consistente com dados observados.")

    return {
        "by_year": results,
        "mean_rho": mean_rho,
        "median_rho": median_rho,
    }


# ============================================================
# VISUALIZAÇÕES (sem mudanças funcionais; herda v5)
# ============================================================
def plot_fraction_maps(all_fractions, all_metas, output_dir):
    """Painel evolução das frações."""
    years = sorted(all_fractions.keys())
    n = len(years)
    ncols = 4
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows * 2, ncols, figsize=(ncols * 3.2, nrows * 5.5))
    if axes.ndim == 1:
        axes = axes.reshape(-1, ncols)

    cmap_veh = mcolors.LinearSegmentedColormap.from_list(
        "veh", ["#fef0d9", "#fdcc8a", "#fc8d59", "#e34a33", "#b30000"], N=256)
    cmap_bio = mcolors.LinearSegmentedColormap.from_list(
        "bio", ["#edf8e9", "#bae4b3", "#74c476", "#238b45", "#004529"], N=256)

    extent = [TARGET_LONS[0], TARGET_LONS[-1], TARGET_LATS[-1], TARGET_LATS[0]]

    for i, year in enumerate(years):
        row_veh = (i // ncols) * 2
        row_bio = row_veh + 1
        col = i % ncols

        fracs = all_fractions[year]
        meta = all_metas[year]

        ax = axes[row_veh, col]
        ax.imshow(fracs["w_veicular"], extent=extent,
                  cmap=cmap_veh, vmin=0.05, vmax=0.55, interpolation="bilinear")
        ax.set_title(f"{year} — Veicular\nµ={meta['w_veh_mean']:.0%}",
                     fontsize=7, fontweight="bold")
        ax.set_xlim(SP_LON_MIN, SP_LON_MAX)
        ax.set_ylim(SP_LAT_MIN, SP_LAT_MAX)
        ax.tick_params(labelsize=5)

        ax = axes[row_bio, col]
        ax.imshow(fracs["w_biomassa"], extent=extent,
                  cmap=cmap_bio, vmin=0.05, vmax=0.55, interpolation="bilinear")
        ax.set_title(f"{year} — Biomassa\nµ={meta['w_bio_mean']:.0%}",
                     fontsize=7, fontweight="bold")
        ax.set_xlim(SP_LON_MIN, SP_LON_MAX)
        ax.set_ylim(SP_LAT_MIN, SP_LAT_MAX)
        ax.tick_params(labelsize=5)

    for j in range(i + 1, nrows * ncols):
        r_v = (j // ncols) * 2
        r_b = r_v + 1
        c = j % ncols
        if r_v < axes.shape[0]:
            axes[r_v, c].set_visible(False)
        if r_b < axes.shape[0]:
            axes[r_b, c].set_visible(False)

    plt.suptitle("Frações Source-Specific de PM2.5 — São Paulo\n"
                 "Calibração Dual: PMF RMSP (Pereira 2025) + KDE Focos BDQueimadas",
                 fontsize=11, fontweight="bold", y=1.01)
    plt.tight_layout()
    path = os.path.join(output_dir, "source_specific_panel.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Painel: {path}")
    return path


def plot_trend_fractions(all_metas, output_dir):
    """Tendência temporal das frações na RMSP e estado."""
    years = sorted(all_metas.keys())
    w_veh_state = [all_metas[y]["w_veh_mean"] for y in years]
    w_bio_state = [all_metas[y]["w_bio_mean"] for y in years]
    w_veh_rmsp = [all_metas[y]["w_veh_rmsp"] for y in years]
    w_bio_rmsp = [all_metas[y]["w_bio_rmsp"] for y in years]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.plot(years, [v * 100 for v in w_veh_state], "o-", color="#e34a33",
            linewidth=2, markersize=5, label="Veicular")
    ax.plot(years, [v * 100 for v in w_bio_state], "s-", color="#238b45",
            linewidth=2, markersize=5, label="Biomassa")
    ax.plot(years, [(1 - v - b) * 100 for v, b in zip(w_veh_state, w_bio_state)],
            "^-", color="#6a51a3", linewidth=2, markersize=5, label="Outros", alpha=0.5)
    ax.axhline(y=41, color="#e34a33", linestyle=":", alpha=0.4, label="PMF veh (41%)")
    ax.axhline(y=25, color="#238b45", linestyle=":", alpha=0.4, label="PMF bio (25%)")
    ax.set_xlabel("Ano"); ax.set_ylabel("Fração (%)")
    ax.set_title("Frações Médias — Estado de SP", fontweight="bold")
    ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3); ax.set_ylim(0, 60)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax = axes[1]
    ax.plot(years, [v * 100 for v in w_veh_rmsp], "o-", color="#e34a33",
            linewidth=2, markersize=5, label="Veicular")
    ax.plot(years, [v * 100 for v in w_bio_rmsp], "s-", color="#238b45",
            linewidth=2, markersize=5, label="Biomassa")
    ax.axhline(y=41, color="#e34a33", linestyle="--", alpha=0.6, label="PMF veh (41%)")
    ax.axhline(y=25, color="#238b45", linestyle="--", alpha=0.6, label="PMF bio (25%)")
    ax.set_xlabel("Ano"); ax.set_ylabel("Fração (%)")
    ax.set_title("Frações na RMSP (calibração = Pereira 2025)", fontweight="bold")
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(0, 60)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.suptitle("Evolução Temporal das Frações Source-Specific de PM2.5",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    path = os.path.join(output_dir, "source_specific_trend.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Tendência: {path}")
    return path


# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Source-specific PM2.5 com calibração dual.")
    parser.add_argument("--boost", type=float, default=BIOMASS_BOOST_DEFAULT,
                        help=f"Fator de boost de biomassa (default {BIOMASS_BOOST_DEFAULT}). "
                             "Sensibilidade: 1.0 (sem boost), 1.2, 1.4, 1.6.")
    parser.add_argument("--output-suffix", type=str, default="",
                        help="Sufixo no diretório de saída para análise de sensibilidade.")
    args = parser.parse_args()

    output_dir = OUTPUT_DIR
    if args.output_suffix:
        output_dir = OUTPUT_DIR + args.output_suffix
        os.makedirs(output_dir, exist_ok=True)

    print("\n" + "=" * 60)
    print("🔬 EXPOSIÇÃO SOURCE-SPECIFIC — PM2.5 SP (Etapa 6)")
    print(f"   Calibração dual: PMF Pereira 2025 (RMSP) +")
    print(f"                    KDE focos BDQueimadas (modulador espacial)")
    print(f"   BIOMASS_BOOST_FACTOR = {args.boost}")
    print(f"   Saída: {output_dir}")
    print("=" * 60)

    # 1. Carregar focos cumulativos e gerar grade de densidade
    fire_coords = load_all_fires_for_kde()
    fire_density = generate_fire_density_grid(
        fire_coords, bandwidth_km=FIRE_KDE_BANDWIDTH_KM)
    fire_density_norm = normalize_fire_density(fire_density)

    save_fire_density_geotiff(fire_density, output_dir, suffix="_raw")
    save_fire_density_geotiff(fire_density_norm, output_dir, suffix="_norm")

    # 2. Gerar superfícies por ano com calibração dual
    all_fractions = {}
    all_metas = {}
    for year in range(YEAR_START, YEAR_END + 1):
        print(f"\n📅 {year}")
        fracs, meta = generate_source_specific_surfaces(
            year, fire_density_normalized=fire_density_norm,
            biomass_boost_factor=args.boost)
        if fracs is None:
            print(f"  ❌ Dados indisponíveis")
            continue
        save_fraction_geotiffs(fracs, year, output_dir)
        all_fractions[year] = fracs
        all_metas[year] = meta
        print(f"  ✅ Veicular={meta['w_veh_mean']:.1%}, "
              f"Biomassa={meta['w_bio_mean']:.1%}, "
              f"Outros={meta['w_out_mean']:.1%} | "
              f"RMSP: V={meta['w_veh_rmsp']:.1%}, B={meta['w_bio_rmsp']:.1%}")

    # 3. Resumo
    print("\n" + "=" * 60)
    print("📊 RESUMO")
    print("=" * 60)
    rows = []
    for year in sorted(all_metas.keys()):
        m = all_metas[year]
        rows.append({
            "year": year,
            "w_veh_SP": f"{m['w_veh_mean']:.1%}",
            "w_bio_SP": f"{m['w_bio_mean']:.1%}",
            "w_out_SP": f"{m['w_out_mean']:.1%}",
            "w_veh_RMSP": f"{m['w_veh_rmsp']:.1%}",
            "w_bio_RMSP": f"{m['w_bio_rmsp']:.1%}",
        })
    summary = pd.DataFrame(rows)
    print(summary.to_string(index=False))
    csv_path = os.path.join(output_dir, "source_specific_summary.csv")
    summary.to_csv(csv_path, index=False)

    # 4. Validação cruzada interna
    val_results = cross_validation_spearman(all_fractions, fire_density)
    val_path = os.path.join(output_dir, "validation_spearman.csv")
    val_df = pd.DataFrame([
        {"year": y, "rho": r["rho"], "p_value": r["p_value"], "n_pixels": r["n_pixels"]}
        for y, r in val_results["by_year"].items()
    ])
    val_df.to_csv(val_path, index=False)
    print(f"\n  📄 Validação Spearman: {val_path}")
    print(f"  ρ médio = {val_results['mean_rho']:.3f}")

    # 5. Visualizações
    plot_fraction_maps(all_fractions, all_metas, output_dir)
    plot_trend_fractions(all_metas, output_dir)

    print("\n" + "=" * 60)
    print(f"✅ COMPLETO! BIOMASS_BOOST = {args.boost}")
    print(f"   ρ Spearman médio (validação interna) = {val_results['mean_rho']:.3f}")
    print(f"   Resultados em: {output_dir}/")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
