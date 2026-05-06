#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cross_validation.py — Validação cruzada CETESB vs CAMS vs MERRA-2
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

Para cada estação CETESB com dados, extrai o valor correspondente
dos grids CAMS e MERRA-2, gerando métricas e gráficos de comparação.

Uso:
    python3 cross_validation.py
"""
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

import config
from data_loaders import load_cetesb_annual, load_cams_pm25, load_merra2_pm25


def extract_satellite_at_stations(cetesb_gdf, years):
    """
    Para cada estação CETESB e ano, extrai o valor PM2.5 correspondente
    dos grids CAMS e MERRA-2 (pixel mais próximo).

    Returns: DataFrame com colunas:
        estacao, ano, lat, lon, cetesb, cams, merra2
    """
    rows = []
    cams_cache = {}
    merra2_cache = {}

    total = len(cetesb_gdf)
    for i, (_, row) in enumerate(cetesb_gdf.iterrows()):
        year = int(row["ano"])
        if year not in years:
            continue

        lat, lon = row["lat"], row["lon"]
        cetesb_val = row["media_anual_pm25_ugm3"]

        # CAMS
        if year not in cams_cache:
            print(f"  📥 Carregando CAMS {year}...")
            cams_cache[year] = load_cams_pm25(year)

        cams_val = None
        if cams_cache[year] is not None:
            try:
                data = cams_cache[year].squeeze(drop=True)
                nearest = data.sel(latitude=lat, longitude=lon, method="nearest")
                cams_val = float(nearest.values)
            except Exception:
                pass

        # MERRA-2
        if year not in merra2_cache:
            print(f"  📥 Carregando MERRA-2 {year}...")
            merra2_cache[year] = load_merra2_pm25(year)

        merra2_val = None
        if merra2_cache[year] is not None:
            try:
                data = merra2_cache[year].squeeze(drop=True)
                lat_name = "lat" if "lat" in data.dims else "latitude"
                lon_name = "lon" if "lon" in data.dims else "longitude"
                nearest = data.sel(**{lat_name: lat, lon_name: lon}, method="nearest")
                merra2_val = float(nearest.values)
            except Exception:
                pass

        rows.append({
            "estacao": row["estacao_nome"],
            "ano": year,
            "lat": lat,
            "lon": lon,
            "cetesb_ugm3": cetesb_val,
            "cams_ugm3": cams_val,
            "merra2_ugm3": merra2_val,
        })

        if (i + 1) % 50 == 0:
            print(f"    Progresso: {i+1}/{total}")

    return pd.DataFrame(rows)


def compute_metrics(obs, pred, name):
    """Calcula métricas estatísticas entre observado e predito."""
    mask = ~(np.isnan(obs) | np.isnan(pred))
    o, p = obs[mask], pred[mask]

    if len(o) < 3:
        return {"source": name, "n": len(o), "error": "dados insuficientes"}

    bias = np.mean(p - o)
    mae = np.mean(np.abs(p - o))
    rmse = np.sqrt(np.mean((p - o) ** 2))
    nrmse = rmse / np.mean(o) * 100
    r, p_val = stats.pearsonr(o, p)
    r2 = r ** 2
    slope, intercept, _, _, std_err = stats.linregress(o, p)

    return {
        "source": name,
        "n": len(o),
        "mean_obs": round(np.mean(o), 2),
        "mean_pred": round(np.mean(p), 2),
        "bias": round(bias, 2),
        "mae": round(mae, 2),
        "rmse": round(rmse, 2),
        "nrmse_pct": round(nrmse, 1),
        "r": round(r, 3),
        "r2": round(r2, 3),
        "p_value": f"{p_val:.2e}",
        "slope": round(slope, 3),
        "intercept": round(intercept, 2),
    }


def plot_scatter(df, output_dir):
    """Gera scatter plots CETESB vs CAMS e CETESB vs MERRA-2."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sources = [
        ("cams_ugm3", "CAMS EAC4", "#3498db", axes[0]),
        ("merra2_ugm3", "MERRA-2", "#e74c3c", axes[1]),
    ]

    for col, label, color, ax in sources:
        mask = df["cetesb_ugm3"].notna() & df[col].notna()
        subset = df[mask]

        if len(subset) < 3:
            ax.text(0.5, 0.5, "Dados insuficientes", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        obs = subset["cetesb_ugm3"].values
        pred = subset[col].values

        ax.scatter(obs, pred, c=color, alpha=0.4, s=25, edgecolors="white", linewidth=0.3)

        # Linha 1:1
        lims = [min(obs.min(), pred.min()) * 0.8, max(obs.max(), pred.max()) * 1.1]
        ax.plot(lims, lims, "k--", linewidth=1, alpha=0.5, label="1:1")

        # Regressão
        slope, intercept, r, p_val, _ = stats.linregress(obs, pred)
        x_fit = np.linspace(lims[0], lims[1], 100)
        ax.plot(x_fit, slope * x_fit + intercept, color=color, linewidth=2,
                label=f"y = {slope:.2f}x + {intercept:.1f}")

        # Métricas no gráfico
        rmse = np.sqrt(np.mean((pred - obs) ** 2))
        bias = np.mean(pred - obs)
        textstr = f"n = {len(obs)}\nR² = {r**2:.3f}\nRMSE = {rmse:.1f}\nBias = {bias:+.1f}"
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment="top", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

        ax.set_xlabel("CETESB (µg/m³)", fontsize=11)
        ax.set_ylabel(f"{label} (µg/m³)", fontsize=11)
        ax.set_title(f"CETESB vs {label}", fontsize=12, fontweight="bold")
        ax.legend(fontsize=9)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(output_dir, "cross_validation_scatter.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Scatter plots salvos: {path}")


def plot_annual_bias(df, output_dir):
    """Gera gráfico de bias anual por fonte."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for col, label, color, ax in [
        ("cams_ugm3", "CAMS EAC4", "#3498db", axes[0]),
        ("merra2_ugm3", "MERRA-2", "#e74c3c", axes[1]),
    ]:
        mask = df["cetesb_ugm3"].notna() & df[col].notna()
        subset = df[mask].copy()
        if len(subset) == 0:
            continue

        subset["bias"] = subset[col] - subset["cetesb_ugm3"]
        annual = subset.groupby("ano").agg(
            mean_bias=("bias", "mean"),
            std_bias=("bias", "std"),
            n=("bias", "count"),
            mean_cetesb=("cetesb_ugm3", "mean"),
            mean_sat=pd.NamedAgg(column=col, aggfunc="mean"),
        ).reset_index()

        ax.bar(annual["ano"], annual["mean_bias"], color=color, alpha=0.7,
               yerr=annual["std_bias"], capsize=3, label=f"Bias ± DP")
        ax.axhline(y=0, color="black", linestyle="-", linewidth=0.8)
        ax.set_xlabel("Ano", fontsize=11)
        ax.set_ylabel("Bias (µg/m³)", fontsize=11)
        ax.set_title(f"Bias anual: {label} − CETESB", fontsize=12, fontweight="bold")
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

        # Texto com bias médio global
        global_bias = subset["bias"].mean()
        ax.text(0.95, 0.05, f"Bias médio: {global_bias:+.1f} µg/m³",
                transform=ax.transAxes, ha="right", fontsize=9,
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

    plt.tight_layout()
    path = os.path.join(output_dir, "cross_validation_annual_bias.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Bias anual salvo: {path}")


def plot_spatial_bias(df, output_dir):
    """Gera mapa de bias espacial médio por estação."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for col, label, ax in [
        ("cams_ugm3", "CAMS EAC4", axes[0]),
        ("merra2_ugm3", "MERRA-2", axes[1]),
    ]:
        mask = df["cetesb_ugm3"].notna() & df[col].notna()
        subset = df[mask].copy()
        if len(subset) == 0:
            continue

        subset["bias"] = subset[col] - subset["cetesb_ugm3"]
        station_bias = subset.groupby("estacao").agg(
            mean_bias=("bias", "mean"),
            lat=("lat", "first"),
            lon=("lon", "first"),
            n=("bias", "count"),
        ).reset_index()

        # Color map: negative=blue, positive=red
        vmax = max(abs(station_bias["mean_bias"].min()), abs(station_bias["mean_bias"].max()))
        sc = ax.scatter(station_bias["lon"], station_bias["lat"],
                        c=station_bias["mean_bias"], cmap="RdBu_r",
                        vmin=-vmax, vmax=vmax, s=60, edgecolors="black", linewidth=0.5)
        plt.colorbar(sc, ax=ax, label="Bias médio (µg/m³)", shrink=0.8)

        ax.set_xlabel("Longitude", fontsize=10)
        ax.set_ylabel("Latitude", fontsize=10)
        ax.set_title(f"Bias espacial: {label} − CETESB", fontsize=11, fontweight="bold")
        ax.set_xlim(-54, -44)
        ax.set_ylim(-26, -19)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(output_dir, "cross_validation_spatial_bias.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Mapa de bias espacial salvo: {path}")


def plot_time_series_comparison(df, output_dir):
    """Séries temporais comparativas para estações selecionadas."""
    # Selecionar 6 estações com mais dados
    station_counts = df.groupby("estacao")["cetesb_ugm3"].count()
    top_stations = station_counts.nlargest(6).index.tolist()

    fig, axes = plt.subplots(2, 3, figsize=(16, 8), sharey=True)
    axes = axes.flatten()

    for i, station in enumerate(top_stations):
        ax = axes[i]
        sdata = df[df["estacao"] == station].sort_values("ano")

        ax.plot(sdata["ano"], sdata["cetesb_ugm3"], "ko-", markersize=4,
                linewidth=1.5, label="CETESB", zorder=3)
        if sdata["cams_ugm3"].notna().any():
            ax.plot(sdata["ano"], sdata["cams_ugm3"], "b^--", markersize=4,
                    linewidth=1, alpha=0.8, label="CAMS")
        if sdata["merra2_ugm3"].notna().any():
            ax.plot(sdata["ano"], sdata["merra2_ugm3"], "rs:", markersize=4,
                    linewidth=1, alpha=0.8, label="MERRA-2")

        # Linha OMS
        ax.axhline(y=15, color="red", linestyle="-", linewidth=0.8, alpha=0.4)
        ax.axhline(y=5, color="darkred", linestyle=":", linewidth=0.5, alpha=0.3)

        ax.set_title(station[:20], fontsize=9, fontweight="bold")
        ax.set_xlabel("Ano", fontsize=8)
        if i % 3 == 0:
            ax.set_ylabel("PM2.5 (µg/m³)", fontsize=9)
        ax.legend(fontsize=6, loc="upper right")
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis="x", labelsize=7, rotation=45)

    plt.suptitle("Séries Temporais: CETESB vs Satélite (Top 6 Estações)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    path = os.path.join(output_dir, "cross_validation_timeseries.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Séries temporais salvas: {path}")


# ============================================================
# MAIN
# ============================================================
def main():
    print("\n" + "=" * 60)
    print("🔍 VALIDAÇÃO CRUZADA: CETESB vs CAMS vs MERRA-2")
    print("=" * 60)

    output_dir = os.path.join(config.OUTPUT_DIR, "cross_validation")
    os.makedirs(output_dir, exist_ok=True)

    # 1. Carregar CETESB
    print("\n📡 Carregando dados CETESB...")
    cetesb = load_cetesb_annual()
    years = sorted(cetesb["ano"].unique())
    print(f"   Anos disponíveis: {years}")

    # 2. Extrair valores satelitais nas estações
    print("\n🛰️  Extraindo valores CAMS e MERRA-2 nas estações CETESB...")
    comparison_df = extract_satellite_at_stations(cetesb, set(years))

    # Salvar dados brutos
    csv_path = os.path.join(output_dir, "cross_validation_data.csv")
    comparison_df.to_csv(csv_path, index=False)
    print(f"\n  📄 Dados salvos: {csv_path}")
    print(f"  📊 Total registros: {len(comparison_df)}")
    print(f"  📊 Com CAMS: {comparison_df['cams_ugm3'].notna().sum()}")
    print(f"  📊 Com MERRA-2: {comparison_df['merra2_ugm3'].notna().sum()}")

    # 3. Métricas estatísticas
    print("\n" + "=" * 60)
    print("📊 MÉTRICAS ESTATÍSTICAS")
    print("=" * 60)

    metrics = []
    for col, name in [("cams_ugm3", "CAMS EAC4"), ("merra2_ugm3", "MERRA-2")]:
        mask = comparison_df["cetesb_ugm3"].notna() & comparison_df[col].notna()
        subset = comparison_df[mask]
        if len(subset) >= 3:
            m = compute_metrics(subset["cetesb_ugm3"].values, subset[col].values, name)
            metrics.append(m)

    metrics_df = pd.DataFrame(metrics)
    print("\n", metrics_df.to_string(index=False))

    metrics_csv = os.path.join(output_dir, "cross_validation_metrics.csv")
    metrics_df.to_csv(metrics_csv, index=False)
    print(f"\n  📄 Métricas salvas: {metrics_csv}")

    # Métricas por ano
    print("\n--- Métricas por ano ---")
    annual_metrics = []
    for year in years:
        year_data = comparison_df[comparison_df["ano"] == year]
        for col, name in [("cams_ugm3", "CAMS"), ("merra2_ugm3", "MERRA-2")]:
            mask = year_data["cetesb_ugm3"].notna() & year_data[col].notna()
            subset = year_data[mask]
            if len(subset) >= 3:
                m = compute_metrics(subset["cetesb_ugm3"].values, subset[col].values, f"{name}")
                m["year"] = year
                annual_metrics.append(m)

    if annual_metrics:
        annual_df = pd.DataFrame(annual_metrics)
        annual_csv = os.path.join(output_dir, "cross_validation_annual_metrics.csv")
        annual_df.to_csv(annual_csv, index=False)
        # Print summary
        for source in ["CAMS", "MERRA-2"]:
            src_data = annual_df[annual_df["source"] == source]
            if len(src_data) > 0:
                print(f"\n  {source}: R² médio={src_data['r2'].mean():.3f}, "
                      f"RMSE médio={src_data['rmse'].mean():.1f}, "
                      f"Bias médio={src_data['bias'].mean():+.1f}")

    # 4. Visualizações
    print("\n" + "=" * 60)
    print("📈 GERANDO VISUALIZAÇÕES")
    print("=" * 60)

    plot_scatter(comparison_df, output_dir)
    plot_annual_bias(comparison_df, output_dir)
    plot_spatial_bias(comparison_df, output_dir)
    plot_time_series_comparison(comparison_df, output_dir)

    print("\n" + "=" * 60)
    print("✅ VALIDAÇÃO CRUZADA COMPLETA!")
    print(f"   Resultados em: {output_dir}/")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
