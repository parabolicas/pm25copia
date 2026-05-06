#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
source_specific.py — Exposição source-specific ao PM2.5 (Etapa 5)
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

Metodologia (conforme Protocolo v2.0):
  1. Calcular fração veicular via MERRA-2:
     f_veh = BCSMASS / (BCSMASS + 1.4 × OCSMASS)
  2. Calibrar com percentuais PMF de referência (Pereira et al. 2025):
     RMSP: veicular 41%, biomassa 25%, secundário 21%, industrial 14%
  3. Fracionar PM2.5 total em componentes:
     PM2.5_veicular  = PM2.5_total × w_veicular(lat, lon, ano)
     PM2.5_biomassa  = PM2.5_total × w_biomassa(lat, lon, ano)
     PM2.5_outros    = PM2.5_total × (1 - w_veicular - w_biomassa)

Uso:
    python3 source_specific.py
"""
import os
import glob
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr
import rasterio
from rasterio.transform import from_bounds
from scipy.interpolate import RegularGridInterpolator
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

# Grid alvo (igual ao pm25_surface.py)
TARGET_RES = 0.05
SP_LAT_MIN, SP_LAT_MAX = -25.5, -19.5
SP_LON_MIN, SP_LON_MAX = -53.5, -44.0
TARGET_LATS = np.arange(SP_LAT_MAX, SP_LAT_MIN - TARGET_RES / 2, -TARGET_RES)
TARGET_LONS = np.arange(SP_LON_MIN, SP_LON_MAX + TARGET_RES / 2, TARGET_RES)

# Percentuais PMF de referência — Pereira et al. (2025), RMSP
# Usados como âncora de calibração para a RMSP
PMF_REFERENCE = {
    "veicular": 0.41,       # 41% — tráfego (HDV + LDV)
    "biomassa": 0.25,       # 25% — queima de biomassa
    "secundario": 0.21,     # 21% — aerossol secundário (sulfato, nitrato)
    "industrial": 0.14,     # 14% — poeira + industrial
}

# Coordenadas de referência RMSP (para calibração)
RMSP_LAT, RMSP_LON = -23.55, -46.63


# ============================================================
# CÁLCULO DE FRAÇÕES POR COMPONENTE MERRA-2
# ============================================================
def load_merra2_components(year):
    """
    Carrega componentes individuais do MERRA-2 para um ano e
    calcula médias anuais de cada variável (kg/m³).

    Returns:
        dict {variável: DataArray 2D (lat, lon)}
    """
    merra_dir = config.MERRA2_DIR
    pattern = os.path.join(merra_dir, f"MERRA2_*{year}*.nc4.nc4")
    files = sorted(glob.glob(pattern))

    if not files:
        return None

    components = {}
    var_names = ["BCSMASS", "OCSMASS", "DUSMASS25", "SSSMASS25", "SO4SMASS"]

    # Acumular médias mensais
    monthly_data = {v: [] for v in var_names}

    for f in files:
        try:
            ds = xr.open_dataset(f)
            for v in var_names:
                if v in ds:
                    data = ds[v].squeeze(drop=True)
                    monthly_data[v].append(data.values)
            ds.close()
        except Exception as e:
            continue

    # Média anual
    for v in var_names:
        if monthly_data[v]:
            arr = np.nanmean(np.array(monthly_data[v]), axis=0)
            # Usar coordenadas do último dataset lido
            ds_ref = xr.open_dataset(files[0])
            components[v] = xr.DataArray(
                arr,
                dims=["lat", "lon"],
                coords={"lat": ds_ref.lat.values, "lon": ds_ref.lon.values},
            )
            ds_ref.close()

    return components if components else None


def compute_source_fractions(components):
    """
    Calcula frações de fonte a partir dos componentes MERRA-2.

    Proxy veicular: BCSMASS / (BCSMASS + 1.4 × OCSMASS)
      - BC (black carbon) é primariamente emitido por combustão diesel
      - OC (organic carbon) vem tanto de veículos quanto de biomassa
      - A razão BC/(BC+OC) é alta em áreas urbanas (diesel) e baixa em áreas rurais (biomassa)

    Proxy biomassa: 1.4 × OCSMASS / PM2.5_total (fração orgânica corrigida)
      - OC é o traçador principal de queima de biomassa
      - OC alto + BC baixo → dominância de biomassa

    Returns:
        dict com DataArrays 2D: f_veh_raw, f_bio_raw, f_dust, f_sulfate, f_seasalt
    """
    bc = components["BCSMASS"]
    oc = components["OCSMASS"]
    dust = components["DUSMASS25"]
    ss = components["SSSMASS25"]
    so4 = components["SO4SMASS"]

    # PM2.5 total (fórmula NASA)
    pm25_total = dust + ss + bc + 1.4 * oc + 1.375 * so4

    # Evitar divisão por zero
    pm25_safe = xr.where(pm25_total > 0, pm25_total, np.nan)

    # Frações brutas por componente
    f_bc = bc / pm25_safe                   # Black carbon
    f_om = (1.4 * oc) / pm25_safe           # Organic matter
    f_dust = dust / pm25_safe               # Dust
    f_sulfate = (1.375 * so4) / pm25_safe   # Sulfate
    f_seasalt = ss / pm25_safe              # Sea salt

    # Proxy veicular: razão BC/carbonáceo total
    carbonaceous = bc + 1.4 * oc
    carbonaceous_safe = xr.where(carbonaceous > 0, carbonaceous, np.nan)
    bc_oc_ratio = bc / carbonaceous_safe  # 0→puro OC(biomassa), 1→puro BC(diesel)

    return {
        "bc_oc_ratio": bc_oc_ratio,     # Proxy veicular bruto (0-1)
        "f_bc": f_bc,                    # Fração black carbon
        "f_om": f_om,                    # Fração organic matter
        "f_dust": f_dust,                # Fração poeira
        "f_sulfate": f_sulfate,          # Fração sulfato
        "f_seasalt": f_seasalt,          # Fração sal marinho
        "pm25_total": pm25_total,        # PM2.5 total (kg/m³)
    }


def calibrate_fractions(fractions_raw, source_class_df=None):
    """
    Calibra as frações brutas MERRA-2 usando os percentuais PMF
    de Pereira et al. (2025) como âncora na RMSP.

    Método:
    1. Extrair bc_oc_ratio na RMSP → valor de referência
    2. Mapear bc_oc_ratio → w_veicular usando relação linear calibrada:
       - Na RMSP: bc_oc_ratio_rmsp → w_veh = 0.41 (PMF)
       - Em áreas rurais com queima: bc_oc_ratio → 0 → w_veh ≈ 0.05 (mínimo veicular)
    3. w_biomassa = w_om_total - w_veicular × f_om_in_veh - w_industrial × f_om_in_ind
       Simplificação: w_biomassa ≈ f_om × (1 - bc_oc_ratio) × calibration_factor

    Returns:
        dict com DataArrays: w_veicular, w_biomassa, w_outros
    """
    bc_oc = fractions_raw["bc_oc_ratio"]
    f_om = fractions_raw["f_om"]
    f_bc = fractions_raw["f_bc"]

    # 1. Valor de referência na RMSP
    lat_idx = np.argmin(np.abs(bc_oc.coords["lat"].values - RMSP_LAT))
    lon_idx = np.argmin(np.abs(bc_oc.coords["lon"].values - RMSP_LON))
    bc_oc_rmsp = float(bc_oc.values[lat_idx, lon_idx])

    if np.isnan(bc_oc_rmsp) or bc_oc_rmsp <= 0:
        bc_oc_rmsp = 0.15  # Fallback

    # 2. Calibração linear: bc_oc_ratio → w_veicular
    # Na RMSP: bc_oc_rmsp → 0.41 (PMF)
    # Em área limpa: bc_oc→0 → 0.05 (background veicular)
    # Slope: (0.41 - 0.05) / bc_oc_rmsp
    slope_veh = (PMF_REFERENCE["veicular"] - 0.05) / bc_oc_rmsp
    w_veicular = 0.05 + slope_veh * bc_oc
    w_veicular = xr.where(w_veicular > 0.70, 0.70, w_veicular)  # Cap máximo
    w_veicular = xr.where(w_veicular < 0.02, 0.02, w_veicular)  # Floor mínimo

    # 3. Fração biomassa: proporcional ao OC excedente (não explicado pelo BC/veicular)
    # Na RMSP: w_biomassa ≈ 0.25 (PMF)
    # Em áreas rurais com queimadas: w_biomassa pode chegar a 0.60
    # Proxy: f_om × (1 - bc_oc_ratio) → alta quando OC domina e BC é baixo
    biomass_proxy = f_om * (1.0 - bc_oc)
    bp_rmsp = float(biomass_proxy.values[lat_idx, lon_idx])

    if np.isnan(bp_rmsp) or bp_rmsp <= 0:
        bp_rmsp = 0.10

    slope_bio = PMF_REFERENCE["biomassa"] / bp_rmsp
    w_biomassa = slope_bio * biomass_proxy
    w_biomassa = xr.where(w_biomassa > 0.70, 0.70, w_biomassa)
    w_biomassa = xr.where(w_biomassa < 0.01, 0.01, w_biomassa)

    # 4. Outros (secundário + industrial + poeira + sal)
    w_outros = 1.0 - w_veicular - w_biomassa
    w_outros = xr.where(w_outros < 0.05, 0.05, w_outros)

    # Renormalizar para somar 1.0
    total = w_veicular + w_biomassa + w_outros
    w_veicular = w_veicular / total
    w_biomassa = w_biomassa / total
    w_outros = w_outros / total

    return {
        "w_veicular": w_veicular,
        "w_biomassa": w_biomassa,
        "w_outros": w_outros,
        "bc_oc_rmsp": bc_oc_rmsp,
    }


# ============================================================
# SUPERFÍCIES SOURCE-SPECIFIC
# ============================================================
def resample_to_target(data_2d, src_lats, src_lons):
    """Resample 2D array para grid alvo 0.05°."""
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


def generate_source_specific_surfaces(year):
    """
    Gera superfícies PM2.5 fracionadas por fonte para um ano.

    Returns:
        dict com arrays (n_lats, n_lons):
          w_veicular, w_biomassa, w_outros (frações 0-1)
        dict com metadados
    """
    # 1. Carregar componentes MERRA-2
    components = load_merra2_components(year)
    if components is None:
        return None, None

    # 2. Calcular frações brutas
    fractions = compute_source_fractions(components)

    # 3. Calibrar com PMF
    calibrated = calibrate_fractions(fractions)

    # 4. Resample para grid alvo
    src_lats = components["BCSMASS"].coords["lat"].values
    src_lons = components["BCSMASS"].coords["lon"].values

    w_veh = resample_to_target(calibrated["w_veicular"].values, src_lats, src_lons)
    w_bio = resample_to_target(calibrated["w_biomassa"].values, src_lats, src_lons)
    w_out = resample_to_target(calibrated["w_outros"].values, src_lats, src_lons)

    # Renormalizar após interpolação
    total = w_veh + w_bio + w_out
    w_veh = w_veh / total
    w_bio = w_bio / total
    w_out = w_out / total

    # Metadata
    meta = {
        "year": year,
        "bc_oc_rmsp": calibrated["bc_oc_rmsp"],
        "w_veh_mean": round(float(np.nanmean(w_veh)), 3),
        "w_bio_mean": round(float(np.nanmean(w_bio)), 3),
        "w_out_mean": round(float(np.nanmean(w_out)), 3),
        "w_veh_rmsp": round(float(w_veh[
            np.argmin(np.abs(TARGET_LATS - RMSP_LAT)),
            np.argmin(np.abs(TARGET_LONS - RMSP_LON))
        ]), 3),
        "w_bio_rmsp": round(float(w_bio[
            np.argmin(np.abs(TARGET_LATS - RMSP_LAT)),
            np.argmin(np.abs(TARGET_LONS - RMSP_LON))
        ]), 3),
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
                            SOURCE="MERRA-2 calibrated PMF Pereira2025")
        paths[name] = filepath

    return paths


# ============================================================
# LOOKUP SOURCE-SPECIFIC
# ============================================================
def lookup_source_specific(lat, lon, year, pm25_total, output_dir=None):
    """
    Consulta instantânea da exposição fracionada por fonte.

    Args:
        lat, lon: coordenadas
        year: ano
        pm25_total: PM2.5 total (µg/m³) — do lookup_pm25()

    Returns:
        dict {pm25_veicular, pm25_biomassa, pm25_outros, w_veicular, w_biomassa, w_outros}
    """
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

    # Calcular PM2.5 por fonte
    if pm25_total is not None:
        result["pm25_veicular"] = round(pm25_total * result["w_veicular"], 2)
        result["pm25_biomassa"] = round(pm25_total * result["w_biomassa"], 2)
        result["pm25_outros"] = round(pm25_total * result["w_outros"], 2)

    return result


# ============================================================
# VISUALIZAÇÕES
# ============================================================
def plot_fraction_maps(all_fractions, all_metas, output_dir):
    """Painel mostrando evolução das frações veicular e biomassa."""
    years = sorted(all_fractions.keys())
    n = len(years)
    ncols = 4
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows * 2, ncols, figsize=(ncols * 3.2, nrows * 5.5))
    if axes.ndim == 1:
        axes = axes.reshape(-1, ncols)

    # Colormaps
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

        # Veicular
        ax = axes[row_veh, col]
        im_v = ax.imshow(fracs["w_veicular"], extent=extent,
                         cmap=cmap_veh, vmin=0.05, vmax=0.55,
                         interpolation="bilinear")
        ax.set_title(f"{year} — Veicular\nµ={meta['w_veh_mean']:.0%}",
                     fontsize=7, fontweight="bold")
        ax.set_xlim(SP_LON_MIN, SP_LON_MAX)
        ax.set_ylim(SP_LAT_MIN, SP_LAT_MAX)
        ax.tick_params(labelsize=5)

        # Biomassa
        ax = axes[row_bio, col]
        im_b = ax.imshow(fracs["w_biomassa"], extent=extent,
                         cmap=cmap_bio, vmin=0.05, vmax=0.55,
                         interpolation="bilinear")
        ax.set_title(f"{year} — Biomassa\nµ={meta['w_bio_mean']:.0%}",
                     fontsize=7, fontweight="bold")
        ax.set_xlim(SP_LON_MIN, SP_LON_MAX)
        ax.set_ylim(SP_LAT_MIN, SP_LAT_MAX)
        ax.tick_params(labelsize=5)

    # Esconder eixos extras
    for j in range(i + 1, nrows * ncols):
        r_v = (j // ncols) * 2
        r_b = r_v + 1
        c = j % ncols
        if r_v < axes.shape[0]:
            axes[r_v, c].set_visible(False)
        if r_b < axes.shape[0]:
            axes[r_b, c].set_visible(False)

    plt.suptitle("Frações Source-Specific de PM2.5 — São Paulo\n"
                 "Veicular (laranja) e Biomassa (verde) — Calibração PMF Pereira et al. 2025",
                 fontsize=11, fontweight="bold", y=1.01)
    plt.tight_layout()
    path = os.path.join(output_dir, "source_specific_panel.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Painel de frações: {path}")
    return path


def plot_trend_fractions(all_metas, output_dir):
    """Tendência temporal das frações na RMSP e no estado."""
    years = sorted(all_metas.keys())
    w_veh_state = [all_metas[y]["w_veh_mean"] for y in years]
    w_bio_state = [all_metas[y]["w_bio_mean"] for y in years]
    w_veh_rmsp = [all_metas[y]["w_veh_rmsp"] for y in years]
    w_bio_rmsp = [all_metas[y]["w_bio_rmsp"] for y in years]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Estado
    ax = axes[0]
    ax.plot(years, [v * 100 for v in w_veh_state], "o-", color="#e34a33",
            linewidth=2, markersize=5, label="Veicular")
    ax.plot(years, [v * 100 for v in w_bio_state], "s-", color="#238b45",
            linewidth=2, markersize=5, label="Biomassa")
    ax.plot(years, [(1 - v - b) * 100 for v, b in zip(w_veh_state, w_bio_state)],
            "^-", color="#6a51a3", linewidth=2, markersize=5, label="Outros", alpha=0.5)
    ax.axhline(y=41, color="#e34a33", linestyle=":", alpha=0.4, label="PMF veh (41%)")
    ax.axhline(y=25, color="#238b45", linestyle=":", alpha=0.4, label="PMF bio (25%)")
    ax.set_xlabel("Ano")
    ax.set_ylabel("Fração (%)")
    ax.set_title("Frações Médias — Estado de SP", fontweight="bold")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 60)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    # RMSP
    ax = axes[1]
    ax.plot(years, [v * 100 for v in w_veh_rmsp], "o-", color="#e34a33",
            linewidth=2, markersize=5, label="Veicular")
    ax.plot(years, [v * 100 for v in w_bio_rmsp], "s-", color="#238b45",
            linewidth=2, markersize=5, label="Biomassa")
    ax.axhline(y=41, color="#e34a33", linestyle="--", alpha=0.6, label="PMF veh (41%)")
    ax.axhline(y=25, color="#238b45", linestyle="--", alpha=0.6, label="PMF bio (25%)")
    ax.set_xlabel("Ano")
    ax.set_ylabel("Fração (%)")
    ax.set_title("Frações na RMSP (calibração = Pereira 2025)", fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 60)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.suptitle("Evolução Temporal das Frações Source-Specific de PM2.5",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    path = os.path.join(output_dir, "source_specific_trend.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Tendência temporal: {path}")
    return path


# ============================================================
# MAIN
# ============================================================
def main():
    print("\n" + "=" * 60)
    print("🔬 EXPOSIÇÃO SOURCE-SPECIFIC — PM2.5 SP (Etapa 5)")
    print("   Fracionamento: Veicular / Biomassa / Outros")
    print("   Calibração: PMF Pereira et al. (2025)")
    print("=" * 60)

    all_fractions = {}
    all_metas = {}

    for year in range(YEAR_START, YEAR_END + 1):
        print(f"\n{'─' * 40}")
        print(f"📅 {year}")
        print(f"{'─' * 40}")

        fracs, meta = generate_source_specific_surfaces(year)
        if fracs is None:
            print(f"  ❌ Dados indisponíveis para {year}")
            continue

        # Salvar GeoTIFFs
        paths = save_fraction_geotiffs(fracs, year, OUTPUT_DIR)
        all_fractions[year] = fracs
        all_metas[year] = meta

        print(f"  ✅ Frações: veicular={meta['w_veh_mean']:.1%}, "
              f"biomassa={meta['w_bio_mean']:.1%}, "
              f"outros={meta['w_out_mean']:.1%}")
        print(f"     RMSP: veicular={meta['w_veh_rmsp']:.1%}, "
              f"biomassa={meta['w_bio_rmsp']:.1%}")
        for name, path in paths.items():
            print(f"     📄 {name}: {os.path.basename(path)}")

    # Resumo
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

    csv_path = os.path.join(OUTPUT_DIR, "source_specific_summary.csv")
    summary.to_csv(csv_path, index=False)

    # Visualizações
    print("\n" + "=" * 60)
    print("📈 GERANDO VISUALIZAÇÕES")
    print("=" * 60)

    plot_fraction_maps(all_fractions, all_metas, OUTPUT_DIR)
    plot_trend_fractions(all_metas, OUTPUT_DIR)

    # Teste de lookup integrado
    print("\n" + "=" * 60)
    print("🔎 TESTE DE EXPOSIÇÃO SOURCE-SPECIFIC")
    print("=" * 60)

    try:
        from pm25_surface import lookup_pm25
        test_points = [
            ("Av. Paulista, SP", -23.5636, -46.6544),
            ("Campinas", -22.9056, -47.0608),
            ("S.J. Rio Preto", -20.8113, -49.3758),
            ("Itaí (rural)", -23.4167, -49.0919),
            ("Ribeirão Preto", -21.1767, -47.8208),
        ]

        for name, lat, lon in test_points:
            pm25_total = lookup_pm25(lat, lon, 2023)
            if pm25_total is not None:
                ss = lookup_source_specific(lat, lon, 2023, pm25_total)
                if ss:
                    print(f"\n  📍 {name} (PM2.5 total = {pm25_total:.1f} µg/m³)")
                    print(f"     Veicular: {ss['pm25_veicular']:.1f} µg/m³ "
                          f"({ss['w_veicular']:.0%})")
                    print(f"     Biomassa: {ss['pm25_biomassa']:.1f} µg/m³ "
                          f"({ss['w_biomassa']:.0%})")
                    print(f"     Outros:   {ss['pm25_outros']:.1f} µg/m³ "
                          f"({ss['w_outros']:.0%})")

    except ImportError:
        print("  ⚠ pm25_surface.py não disponível — pulando teste integrado")

    print("\n" + "=" * 60)
    print("✅ EXPOSIÇÃO SOURCE-SPECIFIC COMPLETA!")
    print(f"   {len(all_fractions) * 3} GeoTIFFs em: {OUTPUT_DIR}/")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
