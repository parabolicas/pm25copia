#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pm25_surface.py — Geração de superfícies anuais de PM2.5 para SP
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

Metodologia: Bias-Corrected Satellite Fusion
  1. Base: MERRA-2 anual → resample bilinear para 0.05°
  2. Calibração: resíduos CETESB − MERRA-2 interpolados via IDW
  3. Fusão: PM2.5_final = MERRA-2_resampled + campo_resíduos

Funções de consulta (usadas pelo exposure.py / cumulative.py):
  - lookup_pm25(lat, lon, year): valor do pixel mais próximo
  - lookup_pm25_buffer(lat, lon, year, buffer_km): média dos pixels
    dentro do raio buffer_km — fonte primária para exposição individual.

Uso:
    python3 pm25_surface.py
"""
import os
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
from data_loaders import load_cetesb_annual, load_merra2_pm25


# ============================================================
# CONFIGURAÇÃO DO GRID ALVO
# ============================================================
TARGET_RES = 0.05           # graus (~5.5 km)
SP_LAT_MIN, SP_LAT_MAX = -25.5, -19.5
SP_LON_MIN, SP_LON_MAX = -53.5, -44.0

TARGET_LATS = np.arange(SP_LAT_MAX, SP_LAT_MIN - TARGET_RES/2, -TARGET_RES)  # N→S
TARGET_LONS = np.arange(SP_LON_MIN, SP_LON_MAX + TARGET_RES/2, TARGET_RES)    # W→E

YEAR_START, YEAR_END = 2008, 2024

# Bias global médio MERRA-2 (do cross-validation)
MERRA2_GLOBAL_BIAS = -3.1  # MERRA-2 subestima em ~3.1 µg/m³

# Diretórios de saída
SURFACE_DIR = os.path.join(config.OUTPUT_DIR, "surfaces")
os.makedirs(SURFACE_DIR, exist_ok=True)


def _haversine_km(lat1, lon1, lat2, lon2):
    """Distância haversine entre dois pontos em km."""
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return R * 2 * np.arcsin(np.sqrt(a))


def resample_to_target(data_array, src_lat_name, src_lon_name):
    """
    Resample um DataArray para o grid alvo usando interpolação bilinear.
    """
    data_2d = data_array.squeeze(drop=True)
    src_lats = data_2d.coords[src_lat_name].values
    src_lons = data_2d.coords[src_lon_name].values
    src_values = data_2d.values

    if src_lats[0] > src_lats[-1]:
        src_lats = src_lats[::-1]
        src_values = src_values[::-1, :]

    if np.any(np.isnan(src_values)):
        from scipy.ndimage import generic_filter
        def nanmean_filter(x):
            valid = x[~np.isnan(x)]
            return np.mean(valid) if len(valid) > 0 else np.nan
        src_values = generic_filter(src_values, nanmean_filter, size=3)

    interpolator = RegularGridInterpolator(
        (src_lats, src_lons), src_values,
        method="linear", bounds_error=False, fill_value=None
    )

    grid_lats, grid_lons = np.meshgrid(TARGET_LATS, TARGET_LONS, indexing="ij")
    points = np.column_stack([grid_lats.ravel(), grid_lons.ravel()])
    resampled = interpolator(points).reshape(grid_lats.shape)

    return resampled


def compute_residual_field(cetesb_gdf, year, merra2_surface):
    """
    Calcula campo de resíduos (CETESB − MERRA-2) interpolado via IDW.
    """
    year_data = cetesb_gdf[cetesb_gdf["ano"] == year].copy()
    if len(year_data) < 1:
        return None

    residuals = []
    for _, station in year_data.iterrows():
        s_lat, s_lon = station["lat"], station["lon"]
        cetesb_val = station["media_anual_pm25_ugm3"]

        lat_idx = np.argmin(np.abs(TARGET_LATS - s_lat))
        lon_idx = np.argmin(np.abs(TARGET_LONS - s_lon))

        if 0 <= lat_idx < len(TARGET_LATS) and 0 <= lon_idx < len(TARGET_LONS):
            merra2_val = merra2_surface[lat_idx, lon_idx]
            if not np.isnan(merra2_val) and not np.isnan(cetesb_val):
                residual = cetesb_val - merra2_val
                residuals.append({
                    "lat": s_lat, "lon": s_lon,
                    "residual": residual,
                    "station": station["estacao_nome"],
                })

    if len(residuals) < 1:
        return None

    res_df = pd.DataFrame(residuals)

    n_lats = len(TARGET_LATS)
    n_lons = len(TARGET_LONS)
    residual_field = np.zeros((n_lats, n_lons))
    power = 2.0
    max_dist_km = 150.0

    station_lats = res_df["lat"].values
    station_lons = res_df["lon"].values
    station_res = res_df["residual"].values

    for i in range(n_lats):
        for j in range(n_lons):
            g_lat = TARGET_LATS[i]
            g_lon = TARGET_LONS[j]

            dists = np.array([
                _haversine_km(g_lat, g_lon, s_lat, s_lon)
                for s_lat, s_lon in zip(station_lats, station_lons)
            ])

            within = dists <= max_dist_km

            if not within.any():
                residual_field[i, j] = -MERRA2_GLOBAL_BIAS
                continue

            d = dists[within]
            r = station_res[within]

            if d.min() < 0.5:
                residual_field[i, j] = r[d.argmin()]
                continue

            weights = 1.0 / (d ** power)
            residual_field[i, j] = np.sum(weights * r) / np.sum(weights)

    return residual_field


def generate_surface(year, cetesb_gdf):
    """Gera superfície PM2.5 para um ano."""
    merra2_data = load_merra2_pm25(year)
    if merra2_data is None:
        print(f"  ⚠ MERRA-2 não disponível para {year}")
        return None, None

    lat_name = "lat" if "lat" in merra2_data.dims else "latitude"
    lon_name = "lon" if "lon" in merra2_data.dims else "longitude"
    merra2_resampled = resample_to_target(merra2_data, lat_name, lon_name)

    has_cetesb = year >= 2012 and cetesb_gdf is not None
    if has_cetesb:
        residual_field = compute_residual_field(cetesb_gdf, year, merra2_resampled)
    else:
        residual_field = None

    if residual_field is not None:
        surface = merra2_resampled + residual_field
        n_stations = len(cetesb_gdf[cetesb_gdf["ano"] == year])
        method = f"merra2_cetesb_fusion_{n_stations}st"
    else:
        surface = merra2_resampled + (-MERRA2_GLOBAL_BIAS)
        method = "merra2_global_bias_corrected"

    surface = np.clip(surface, 0.5, None)

    meta = {
        "year": year,
        "method": method,
        "mean": round(float(np.nanmean(surface)), 2),
        "median": round(float(np.nanmedian(surface)), 2),
        "min": round(float(np.nanmin(surface)), 2),
        "max": round(float(np.nanmax(surface)), 2),
        "std": round(float(np.nanstd(surface)), 2),
        "shape": surface.shape,
    }

    return surface, meta


def save_geotiff(surface, year, output_dir):
    """Salva superfície como GeoTIFF."""
    filepath = os.path.join(output_dir, f"pm25_sp_{year}.tif")

    n_lats = len(TARGET_LATS)
    n_lons = len(TARGET_LONS)

    west = TARGET_LONS[0] - TARGET_RES / 2
    east = TARGET_LONS[-1] + TARGET_RES / 2
    north = TARGET_LATS[0] + TARGET_RES / 2
    south = TARGET_LATS[-1] - TARGET_RES / 2

    transform = from_bounds(west, south, east, north, n_lons, n_lats)

    with rasterio.open(
        filepath, "w",
        driver="GTiff",
        height=n_lats,
        width=n_lons,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
        nodata=-9999.0,
    ) as dst:
        data = surface.astype(np.float32)
        data[np.isnan(data)] = -9999.0
        dst.write(data, 1)
        dst.update_tags(
            YEAR=str(year),
            UNITS="ug/m3",
            SOURCE="MERRA-2 + CETESB fusion",
            PROJECT="PM2.5 CPNPC FMUSP",
        )

    return filepath


def lookup_pm25(lat, lon, year, surface_dir=None):
    """
    Consulta de PM2.5 no pixel mais próximo da superfície bias-corrected.

    Args:
        lat, lon: coordenadas do ponto
        year: ano de interesse
        surface_dir: diretório com os GeoTIFFs (default: SURFACE_DIR)

    Returns:
        float PM2.5 em µg/m³, ou None se indisponível.
    """
    if surface_dir is None:
        surface_dir = SURFACE_DIR

    filepath = os.path.join(surface_dir, f"pm25_sp_{year}.tif")
    if not os.path.exists(filepath):
        return None

    with rasterio.open(filepath) as src:
        row, col = src.index(lon, lat)
        if 0 <= row < src.height and 0 <= col < src.width:
            val = src.read(1)[row, col]
            if val == src.nodata or np.isnan(val):
                return None
            return float(val)
    return None


def lookup_pm25_buffer(lat, lon, year, buffer_km, surface_dir=None):
    """
    Consulta de PM2.5 em buffer: média dos pixels da superfície
    bias-corrected dentro de um raio (km) ao redor do ponto.

    Esta é a função PREFERIDA para estimativa de exposição individual:
    usa a superfície já calibrada com CETESB (via IDW de resíduos),
    em vez da hierarquia ad-hoc CETESB→CAMS→MERRA-2 bruto.

    Args:
        lat, lon: coordenadas do ponto (paciente)
        year: ano de interesse
        buffer_km: raio em km (5, 10, 25, etc.)
        surface_dir: diretório com os GeoTIFFs (default: SURFACE_DIR)

    Returns:
        dict com:
          - pm25_value (float): média dos pixels no buffer (µg/m³)
          - n_pixels (int): número de pixels usados
          - method (str): "buffer_mean" ou "nearest_pixel"
        ou None se a superfície/ponto estiver indisponível.
    """
    if surface_dir is None:
        surface_dir = SURFACE_DIR

    filepath = os.path.join(surface_dir, f"pm25_sp_{year}.tif")
    if not os.path.exists(filepath):
        return None

    with rasterio.open(filepath) as src:
        data = src.read(1)
        nodata = src.nodata

        # Pixel central
        row_c, col_c = src.index(lon, lat)
        if not (0 <= row_c < src.height and 0 <= col_c < src.width):
            return None

        # Resolução do raster em graus (esperada: 0.05° ≈ 5.5 km)
        res_deg = abs(src.res[0])
        # Janela de pixels suficiente para cobrir o buffer (com folga de 1 px)
        buffer_deg = buffer_km / 111.0
        n_pix_window = int(np.ceil(buffer_deg / res_deg)) + 1

        r0 = max(0, row_c - n_pix_window)
        r1 = min(src.height, row_c + n_pix_window + 1)
        c0 = max(0, col_c - n_pix_window)
        c1 = min(src.width, col_c + n_pix_window + 1)

        values = []
        for r in range(r0, r1):
            for c in range(c0, c1):
                v = data[r, c]
                if v == nodata or np.isnan(v):
                    continue
                # Centro do pixel em coordenadas geográficas
                p_lon, p_lat = src.xy(r, c)
                d_km = _haversine_km(lat, lon, p_lat, p_lon)
                if d_km <= buffer_km:
                    values.append(float(v))

        if len(values) == 0:
            # Buffer menor que 1 pixel — usa o pixel central
            v = data[row_c, col_c]
            if v == nodata or np.isnan(v):
                return None
            return {
                "pm25_value": round(float(v), 2),
                "n_pixels": 1,
                "method": "nearest_pixel",
            }

        return {
            "pm25_value": round(float(np.mean(values)), 2),
            "n_pixels": len(values),
            "method": "buffer_mean",
        }


def generate_annual_maps(surfaces, metas, cetesb_gdf, output_dir):
    """Painel de mapas anuais de PM2.5 para publicação."""
    from matplotlib.gridspec import GridSpec

    years = sorted(surfaces.keys())
    n = len(years)
    ncols = 4
    nrows = (n + ncols - 1) // ncols

    colors_list = ["#1a9850", "#91cf60", "#d9ef8b", "#fee08b",
                   "#fdae61", "#f46d43", "#d73027", "#a50026"]
    cmap = mcolors.LinearSegmentedColormap.from_list("pm25", colors_list, N=256)

    all_vals = np.concatenate([surfaces[y].ravel() for y in years])
    vmin = max(0, np.nanpercentile(all_vals, 2))
    vmax = np.nanpercentile(all_vals, 98)

    fig = plt.figure(figsize=(ncols * 3.5 + 1.2, nrows * 3.2))
    gs = GridSpec(nrows, ncols + 1, figure=fig,
                  width_ratios=[0.3] + [1] * ncols,
                  wspace=0.25, hspace=0.35)

    axes_map = []
    for r in range(nrows):
        for c in range(ncols):
            ax = fig.add_subplot(gs[r, c + 1], aspect="equal")
            axes_map.append(ax)

    im = None
    for i, year in enumerate(years):
        ax = axes_map[i]
        surface = surfaces[year]
        meta = metas[year]

        im = ax.imshow(
            surface,
            extent=[TARGET_LONS[0], TARGET_LONS[-1], TARGET_LATS[-1], TARGET_LATS[0]],
            cmap=cmap, vmin=vmin, vmax=vmax,
            interpolation="bilinear",
        )

        if cetesb_gdf is not None:
            yr_data = cetesb_gdf[cetesb_gdf["ano"] == year]
            if len(yr_data) > 0:
                ax.scatter(yr_data["lon"], yr_data["lat"],
                           c="white", s=8, edgecolors="black", linewidth=0.3,
                           zorder=5, marker="^")

        ax.set_title(f"{year}\nµ={meta['mean']:.1f} (n={meta.get('method','').split('_')[-1] if 'st' in meta.get('method','') else 'global'})",
                     fontsize=8, fontweight="bold")
        ax.set_xlim(SP_LON_MIN, SP_LON_MAX)
        ax.set_ylim(SP_LAT_MIN, SP_LAT_MAX)
        ax.tick_params(labelsize=6)

        if i % ncols != 0:
            ax.set_yticklabels([])
        if i < (nrows - 1) * ncols:
            ax.set_xticklabels([])

    for j in range(i + 1, len(axes_map)):
        axes_map[j].set_visible(False)

    cbar_ax = fig.add_subplot(gs[:, 0])
    cbar_ax.set_facecolor("white")
    cbar = fig.colorbar(im, ax=cbar_ax, location="left", shrink=0.8, aspect=40, pad=0.05)
    cbar.set_label("PM2.5 (µg/m³)", fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    cbar_ax.set_visible(False)

    cbar.ax.axhline(y=5, color="darkred", linewidth=1, linestyle="--")
    cbar.ax.axhline(y=15, color="red", linewidth=1, linestyle="--")

    plt.suptitle("Superfícies Anuais de PM2.5 — Estado de São Paulo\n"
                 "Método: Bias-Corrected MERRA-2 + CETESB Fusion",
                 fontsize=13, fontweight="bold", y=1.01)

    path = os.path.join(output_dir, "pm25_annual_panel.png")
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  📄 Painel anual salvo: {path}")
    return path


def generate_trend_plot(metas, output_dir):
    """Tendência temporal da média de PM2.5 em SP."""
    years = sorted(metas.keys())
    means = [metas[y]["mean"] for y in years]
    stds = [metas[y]["std"] for y in years]

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.fill_between(years,
                    [m - s for m, s in zip(means, stds)],
                    [m + s for m, s in zip(means, stds)],
                    alpha=0.2, color="#3498db")
    ax.plot(years, means, "o-", color="#2c3e50", markersize=6, linewidth=2,
            label="Média estadual PM2.5")

    ax.axhline(y=5, color="darkred", linestyle=":", linewidth=1, alpha=0.7,
               label="OMS Annual (5 µg/m³)")
    ax.axhline(y=15, color="red", linestyle="--", linewidth=1, alpha=0.5,
               label="OMS Interim-4 (15 µg/m³)")

    ax.axvspan(2012, 2024, alpha=0.05, color="green", label="Com dados CETESB")
    ax.axvspan(2008, 2011, alpha=0.05, color="orange", label="Só MERRA-2 (corrigido)")

    ax.set_xlabel("Ano", fontsize=12)
    ax.set_ylabel("PM2.5 (µg/m³)", fontsize=12)
    ax.set_title("Tendência Temporal de PM2.5 — Estado de São Paulo (2008–2024)",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlim(2007.5, 2024.5)

    plt.tight_layout()
    path = os.path.join(output_dir, "pm25_trend.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Tendência temporal salva: {path}")
    return path


def run_leave_one_out(cetesb_gdf, year=2023):
    """Validação leave-one-out para um ano."""
    year_data = cetesb_gdf[cetesb_gdf["ano"] == year].copy()
    if len(year_data) < 5:
        print(f"  ⚠ Insuficientes estações para LOO em {year}")
        return None

    merra2_data = load_merra2_pm25(year)
    if merra2_data is None:
        return None

    lat_name = "lat" if "lat" in merra2_data.dims else "latitude"
    lon_name = "lon" if "lon" in merra2_data.dims else "longitude"
    merra2_resampled = resample_to_target(merra2_data, lat_name, lon_name)

    results = []
    stations = year_data.index.tolist()

    for idx in stations:
        left_out = year_data.loc[idx]
        remaining = cetesb_gdf.drop(index=idx)

        residual_field = compute_residual_field(remaining, year, merra2_resampled)
        if residual_field is None:
            continue

        surface = merra2_resampled + residual_field
        surface = np.clip(surface, 0.5, None)

        s_lat, s_lon = left_out["lat"], left_out["lon"]
        lat_idx = np.argmin(np.abs(TARGET_LATS - s_lat))
        lon_idx = np.argmin(np.abs(TARGET_LONS - s_lon))

        pred_val = surface[lat_idx, lon_idx]
        obs_val = left_out["media_anual_pm25_ugm3"]
        merra2_val = merra2_resampled[lat_idx, lon_idx]

        results.append({
            "station": left_out["estacao_nome"],
            "lat": s_lat,
            "lon": s_lon,
            "observed": obs_val,
            "predicted_fused": pred_val,
            "predicted_merra2": merra2_val,
            "error_fused": pred_val - obs_val,
            "error_merra2": merra2_val - obs_val,
        })

    return pd.DataFrame(results)


# ============================================================
# MAIN
# ============================================================
def main():
    print("\n" + "=" * 60)
    print("🗺️  GERAÇÃO DE SUPERFÍCIES ANUAIS PM2.5 — SP")
    print("   Método: Bias-Corrected MERRA-2 + CETESB Fusion")
    print("=" * 60)

    print("\n📡 Carregando dados CETESB...")
    cetesb = load_cetesb_annual()

    print(f"\n🎯 Grid alvo: {len(TARGET_LATS)}×{len(TARGET_LONS)} pixels "
          f"({TARGET_RES}° ≈ {TARGET_RES*111:.1f} km)")
    print(f"   Extensão: lat [{SP_LAT_MIN}, {SP_LAT_MAX}], "
          f"lon [{SP_LON_MIN}, {SP_LON_MAX}]")

    surfaces = {}
    metas = {}

    for year in range(YEAR_START, YEAR_END + 1):
        print(f"\n{'─'*40}")
        print(f"📅 {year}")
        print(f"{'─'*40}")

        surface, meta = generate_surface(year, cetesb)
        if surface is None:
            print(f"  ❌ Não foi possível gerar superfície para {year}")
            continue

        tif_path = save_geotiff(surface, year, SURFACE_DIR)
        surfaces[year] = surface
        metas[year] = meta

        print(f"  ✅ PM2.5: µ={meta['mean']:.1f}, med={meta['median']:.1f}, "
              f"range=[{meta['min']:.1f}, {meta['max']:.1f}] µg/m³")
        print(f"  📄 GeoTIFF: {tif_path}")
        print(f"     Método: {meta['method']}")

    print("\n" + "=" * 60)
    print("📊 RESUMO DAS SUPERFÍCIES GERADAS")
    print("=" * 60)

    summary_rows = []
    for year in sorted(metas.keys()):
        m = metas[year]
        summary_rows.append({
            "year": year, "mean": m["mean"], "median": m["median"],
            "min": m["min"], "max": m["max"], "std": m["std"],
            "method": m["method"],
        })
    summary_df = pd.DataFrame(summary_rows)
    print(summary_df.to_string(index=False))

    summary_csv = os.path.join(SURFACE_DIR, "surface_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"\n  📄 Resumo: {summary_csv}")

    print("\n" + "=" * 60)
    print("📈 GERANDO VISUALIZAÇÕES")
    print("=" * 60)

    generate_annual_maps(surfaces, metas, cetesb, SURFACE_DIR)
    generate_trend_plot(metas, SURFACE_DIR)

    print("\n" + "=" * 60)
    print("🔬 VALIDAÇÃO LEAVE-ONE-OUT (2023)")
    print("=" * 60)

    loo_df = run_leave_one_out(cetesb, year=2023)
    if loo_df is not None and len(loo_df) > 0:
        rmse_fused = np.sqrt(np.mean(loo_df["error_fused"] ** 2))
        rmse_merra2 = np.sqrt(np.mean(loo_df["error_merra2"] ** 2))
        mae_fused = np.mean(np.abs(loo_df["error_fused"]))
        mae_merra2 = np.mean(np.abs(loo_df["error_merra2"]))
        bias_fused = np.mean(loo_df["error_fused"])
        bias_merra2 = np.mean(loo_df["error_merra2"])
        from scipy import stats
        r_fused, _ = stats.pearsonr(loo_df["observed"], loo_df["predicted_fused"])
        r_merra2, _ = stats.pearsonr(loo_df["observed"], loo_df["predicted_merra2"])

        print(f"\n  {'Métrica':<12} {'Fusão':>10} {'MERRA-2 puro':>14} {'Melhoria':>10}")
        print(f"  {'─'*48}")
        print(f"  {'RMSE':<12} {rmse_fused:>10.2f} {rmse_merra2:>14.2f} "
              f"{(1-rmse_fused/rmse_merra2)*100:>+9.1f}%")
        print(f"  {'MAE':<12} {mae_fused:>10.2f} {mae_merra2:>14.2f} "
              f"{(1-mae_fused/mae_merra2)*100:>+9.1f}%")
        print(f"  {'Bias':<12} {bias_fused:>+10.2f} {bias_merra2:>+14.2f}")
        print(f"  {'R':<12} {r_fused:>10.3f} {r_merra2:>14.3f}")
        print(f"  {'n estações':<12} {len(loo_df):>10}")

        loo_csv = os.path.join(SURFACE_DIR, "leave_one_out_2023.csv")
        loo_df.to_csv(loo_csv, index=False)
        print(f"\n  📄 LOO detalhes: {loo_csv}")

    print("\n" + "=" * 60)
    print("🔎 TESTE DE LOOKUP RÁPIDO")
    print("=" * 60)

    test_points = [
        ("Av. Paulista, SP", -23.5636, -46.6544),
        ("Campinas", -22.9056, -47.0608),
        ("S.J. Rio Preto", -20.8113, -49.3758),
        ("Itaí (rural)", -23.4167, -49.0919),
        ("Santos", -23.9536, -46.3330),
    ]

    for name, lat, lon in test_points:
        print(f"\n  📍 {name} ({lat:.4f}, {lon:.4f}):")
        for year in [2012, 2018, 2024]:
            val = lookup_pm25(lat, lon, year)
            if val is not None:
                print(f"     {year}: {val:.1f} µg/m³ (pixel)")
            else:
                print(f"     {year}: N/D")
            # Demonstrar lookup_pm25_buffer (raio 10 km)
            buf = lookup_pm25_buffer(lat, lon, year, 10)
            if buf is not None:
                print(f"            buffer 10km: {buf['pm25_value']:.1f} µg/m³ "
                      f"(n={buf['n_pixels']} px, {buf['method']})")

    print("\n" + "=" * 60)
    print("✅ SUPERFÍCIES GERADAS COM SUCESSO!")
    print(f"   {len(surfaces)} GeoTIFFs em: {SURFACE_DIR}/")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
