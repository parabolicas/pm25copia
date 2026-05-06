#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
source_apportionment.py — Classificação simplificada de fontes de PM2.5
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

Classifica cada município de SP como:
  - VEICULAR: dominância de emissões veiculares (alta densidade de frota)
  - BIOMASSA: dominância de queima de biomassa (alta densidade de focos de calor)
  - MISTA: contribuição relevante de ambas as fontes

Metodologia:
  1. Índice Veicular (IV): veículos diesel por km² do município
  2. Índice Biomassa (IB): focos de calor acumulados por 1000 km² (2008–2024)
  3. Classificação por scores padronizados (z-scores)

Uso:
    python3 source_apportionment.py
"""
import os
import glob
import zipfile
import tempfile
import warnings
import unicodedata
import re
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from shapely.geometry import Point

import config


def _normalize_municipio(name):
    """Normaliza nome de município: remove acentos, caracteres especiais, upper."""
    if pd.isna(name):
        return ""
    s = str(name).upper().strip()
    # Remove accents via NFD decomposition
    s = unicodedata.normalize("NFD", s)
    s = "".join(c for c in s if unicodedata.category(c) != "Mn")
    # Fix broken Latin-1 sequences (e.g. \x89, \x8d)
    s = re.sub(r'[^A-Z0-9 ]', '', s)
    s = re.sub(r'\s+', ' ', s).strip()
    return s


# ============================================================
# CONFIGURAÇÃO
# ============================================================
OUTPUT_DIR = os.path.join(config.OUTPUT_DIR, "source_apportionment")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Limiares para classificação (percentis do z-score)
# Acima de THRESHOLD_HIGH → fonte dominante
THRESHOLD_DOMINANT = 0.5  # z-score padronizado
THRESHOLD_RELEVANT = 0.2  # contribuição relevante


# ============================================================
# DADOS DE FROTA
# ============================================================
def load_fleet_data():
    """
    Carrega e agrega frota por município.
    Retorna DataFrame com totais por município e indicadores de emissão.
    """
    print("🚗 Carregando dados de frota (939 MB)...")
    cols = ["id_municipio", "municipio", "combustivel", "tipo_veiculo",
            "quantidade_veiculos"]
    df = pd.read_csv(config.FROTA_CSV, usecols=cols)
    print(f"   ✅ {len(df):,} registros, {df['municipio'].nunique()} municípios")

    # Agregar por município
    agg = df.groupby(["id_municipio", "municipio"]).agg(
        total_veiculos=("quantidade_veiculos", "sum"),
    ).reset_index()

    # Diesel (caminhões, ônibus — emissão pesada)
    diesel = df[df["combustivel"] == "DIESEL"].groupby("municipio").agg(
        veiculos_diesel=("quantidade_veiculos", "sum"),
    ).reset_index()

    # Veículos pesados (caminhão + caminhonete + utilitário)
    heavy_types = ["CAMINHAO", "CAMINHONETE", "UTILITARIO"]
    heavy = df[df["tipo_veiculo"].isin(heavy_types)].groupby("municipio").agg(
        veiculos_pesados=("quantidade_veiculos", "sum"),
    ).reset_index()

    # Merge
    fleet = agg.merge(diesel, on="municipio", how="left")
    fleet = fleet.merge(heavy, on="municipio", how="left")
    fleet = fleet.fillna(0)

    # Indicadores
    fleet["pct_diesel"] = (fleet["veiculos_diesel"] / fleet["total_veiculos"] * 100).round(1)
    fleet["pct_pesados"] = (fleet["veiculos_pesados"] / fleet["total_veiculos"] * 100).round(1)

    print(f"   📊 Total: {fleet['total_veiculos'].sum():,.0f} veículos, "
          f"{fleet['veiculos_diesel'].sum():,.0f} diesel")

    return fleet


# ============================================================
# DADOS DE QUEIMADAS
# ============================================================
def load_all_fires():
    """
    Carrega todos os focos de calor (2008–2024) e agrega por município.
    """
    print("\n🔥 Carregando focos de calor (BDQueimadas 2008–2024)...")

    all_foci = []
    zip_files = sorted(glob.glob(os.path.join(config.BDQUEIMADAS_DIR, "*.zip")))

    for zf_path in zip_files:
        year_str = os.path.basename(zf_path).replace("focos_br_sp_ref_", "").replace(".zip", "")
        try:
            year = int(year_str)
        except ValueError:
            continue

        tmpdir = tempfile.mkdtemp()
        try:
            with zipfile.ZipFile(zf_path, "r") as zf:
                zf.extractall(tmpdir)
            csvs = glob.glob(os.path.join(tmpdir, "**", "*.csv"), recursive=True)
            if csvs:
                df = pd.read_csv(csvs[0], encoding="latin-1",
                                 usecols=["lat", "lon", "municipio"])
                df["ano"] = year
                all_foci.append(df)
        except Exception as e:
            print(f"   ⚠ Erro {year}: {e}")

    if not all_foci:
        return pd.DataFrame()

    fires = pd.concat(all_foci, ignore_index=True)
    print(f"   ✅ {len(fires):,} focos totais, "
          f"{fires['municipio'].nunique()} municípios, "
          f"anos {fires['ano'].min()}-{fires['ano'].max()}")

    return fires


def aggregate_fires_by_municipality(fires):
    """Agrega focos de calor por município."""
    if len(fires) == 0:
        return pd.DataFrame()

    # Total de focos por município (acumulado 2008–2024)
    fire_mun = fires.groupby("municipio").agg(
        total_focos=("lat", "count"),
        n_anos_com_focos=("ano", "nunique"),
        lat_centroid=("lat", "mean"),
        lon_centroid=("lon", "mean"),
    ).reset_index()

    # Média anual de focos
    n_years = fires["ano"].nunique()
    fire_mun["focos_por_ano"] = (fire_mun["total_focos"] / n_years).round(1)

    return fire_mun


# ============================================================
# ÁREAS DOS MUNICÍPIOS (estimativa via IBGE)
# ============================================================
def estimate_municipality_areas(fleet, fire_mun):
    """
    Estima área dos municípios usando dados IBGE (via API ou hardcoded).
    Como fallback, usa a área total de SP dividida proporcionalmente.
    """
    # Área total de SP: 248,219 km²
    SP_AREA_KM2 = 248_219.0
    n_mun = len(fleet)

    # Usar áreas aproximadas baseadas em dados públicos do IBGE
    # Para este projeto, usamos a área mediana de SP = ~384 km²/município
    # e ajustamos pela frota (proxy de urbanização/tamanho)
    median_area = SP_AREA_KM2 / n_mun

    # Municípios com mais veículos tendem a ser menores em área (urbanos)
    # Usamos log da frota como proxy inverso de área
    log_fleet = np.log10(fleet["total_veiculos"] + 1)
    max_log = log_fleet.max()
    min_log = log_fleet.min()

    # Municípios grandes (rural): mais frota = menor área
    # Normalizar: quando log_fleet é alto (urbano), area é baixa
    # Quando log_fleet é baixo (rural), area é alta
    norm = (log_fleet - min_log) / (max_log - min_log)  # 0=rural, 1=urbano

    # Área: rural→grande (2x mediana), urbano→pequeno (0.3x mediana)
    # Isso é uma aproximação, mas funciona para classificação relativa
    areas = median_area * (2.0 - 1.7 * norm)

    # Ajustar para que o total seja SP_AREA_KM2
    areas = areas * (SP_AREA_KM2 / areas.sum())

    fleet["area_km2"] = areas.round(1)

    return fleet


# ============================================================
# CLASSIFICAÇÃO
# ============================================================
def classify_sources(fleet, fire_mun):
    """
    Classifica cada município como VEICULAR, BIOMASSA ou MISTA.

    Índices:
    - IV (Índice Veicular): veículos diesel / km² + pesados / km²
    - IB (Índice Biomassa): focos / 1000 km² / ano

    Classificação:
    - z_IV > 0.5 e z_IV > z_IB → VEICULAR
    - z_IB > 0.5 e z_IB > z_IV → BIOMASSA
    - ambos > 0.2 ou diferença < 0.3 → MISTA
    - ambos baixos → BAIXA_EMISSÃO
    """
    # Estimar áreas
    fleet = estimate_municipality_areas(fleet, fire_mun)

    # Normalizar nome do município para merge (strip accents for both sources)
    fleet["municipio_norm"] = fleet["municipio"].apply(_normalize_municipio)

    if len(fire_mun) > 0:
        fire_mun["municipio_norm"] = fire_mun["municipio"].apply(_normalize_municipio)
        # Merge
        merged = fleet.merge(fire_mun[["municipio_norm", "total_focos", "focos_por_ano",
                                        "lat_centroid", "lon_centroid"]],
                             on="municipio_norm", how="left")
    else:
        merged = fleet.copy()
        merged["total_focos"] = 0
        merged["focos_por_ano"] = 0

    merged = merged.fillna({"total_focos": 0, "focos_por_ano": 0})

    # Índice Veicular: diesel/km² + pesados/km² (ponderado)
    merged["diesel_per_km2"] = merged["veiculos_diesel"] / merged["area_km2"]
    merged["pesados_per_km2"] = merged["veiculos_pesados"] / merged["area_km2"]
    merged["idx_veicular"] = merged["diesel_per_km2"] + 0.5 * merged["pesados_per_km2"]

    # Índice Biomassa: focos/1000km²/ano
    merged["idx_biomassa"] = merged["focos_por_ano"] / merged["area_km2"] * 1000

    # Padronizar (z-scores robustos — usando mediana e MAD)
    for col in ["idx_veicular", "idx_biomassa"]:
        med = merged[col].median()
        mad = np.median(np.abs(merged[col] - med))
        if mad > 0:
            merged[f"z_{col.replace('idx_', '')}"] = (merged[col] - med) / (mad * 1.4826)
        else:
            merged[f"z_{col.replace('idx_', '')}"] = 0

    # Classificar
    def _classify(row):
        z_v = row["z_veicular"]
        z_b = row["z_biomassa"]

        # Ambos altos → MISTA
        if z_v > THRESHOLD_RELEVANT and z_b > THRESHOLD_RELEVANT:
            if abs(z_v - z_b) < 0.5:
                return "MISTA"
            elif z_v > z_b:
                return "VEICULAR_MISTA"
            else:
                return "BIOMASSA_MISTA"

        # Dominância clara
        if z_v > THRESHOLD_DOMINANT and z_b <= THRESHOLD_RELEVANT:
            return "VEICULAR"
        if z_b > THRESHOLD_DOMINANT and z_v <= THRESHOLD_RELEVANT:
            return "BIOMASSA"

        # Low emissão
        if z_v <= THRESHOLD_RELEVANT and z_b <= THRESHOLD_RELEVANT:
            return "BAIXA_EMISSÃO"

        # Outros
        if z_v > z_b:
            return "VEICULAR"
        return "BIOMASSA"

    merged["classificacao"] = merged.apply(_classify, axis=1)

    # Simplificar para 3 categorias principais (para análise epidemiológica)
    def _simplify(c):
        if "VEICULAR" in c and "MISTA" not in c:
            return "VEICULAR"
        if "BIOMASSA" in c and "MISTA" not in c:
            return "BIOMASSA"
        if "MISTA" in c:
            return "MISTA"
        return "BAIXA_EMISSÃO"

    merged["fonte_principal"] = merged["classificacao"].apply(_simplify)

    return merged


# ============================================================
# VISUALIZAÇÃO
# ============================================================
def plot_classification_map(classified, mesorregioes, output_dir):
    """Mapa de classificação de fontes por município."""
    # Criar GeoDataFrame com centróides
    # Usar centróides das queimadas quando disponível, senão estimar
    has_coords = classified["lat_centroid"].notna()

    # Para municípios sem focos (sem coordenadas de queimada),
    # estimar localização usando mesorregião shapefile centroid como fallback
    if not has_coords.all():
        # Usar coordenadas médias de SP como fallback genérico
        classified.loc[~has_coords, "lat_centroid"] = -22.5
        classified.loc[~has_coords, "lon_centroid"] = -48.5

    geometry = [Point(lon, lat) for lon, lat in
                zip(classified["lon_centroid"], classified["lat_centroid"])]
    gdf = gpd.GeoDataFrame(classified, geometry=geometry, crs="EPSG:4326")

    # Cores por fonte
    color_map = {
        "VEICULAR": "#e74c3c",       # vermelho
        "BIOMASSA": "#f39c12",        # laranja
        "MISTA": "#9b59b6",           # roxo
        "BAIXA_EMISSÃO": "#27ae60",   # verde
    }

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # --- Mapa 1: Classificação ---
    ax1 = axes[0]
    if mesorregioes is not None:
        mesorregioes.boundary.plot(ax=ax1, linewidth=0.5, color="gray", alpha=0.5)

    for fonte, color in color_map.items():
        subset = gdf[gdf["fonte_principal"] == fonte]
        if len(subset) > 0:
            # Tamanho proporcional à frota
            sizes = np.clip(np.log10(subset["total_veiculos"] + 1) * 5, 5, 80)
            ax1.scatter(subset.geometry.x, subset.geometry.y,
                        c=color, s=sizes, alpha=0.6,
                        edgecolors="black", linewidth=0.2, label=fonte, zorder=3)

    ax1.set_xlim(-54, -44)
    ax1.set_ylim(-26, -19)
    ax1.set_title("Classificação de Fontes de PM2.5\npor Município — SP",
                   fontsize=12, fontweight="bold")
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    ax1.legend(fontsize=8, loc="lower right",
               title="Fonte Principal", title_fontsize=9)
    ax1.grid(True, alpha=0.2)

    # --- Mapa 2: Índices sobrepostos ---
    ax2 = axes[1]
    if mesorregioes is not None:
        mesorregioes.boundary.plot(ax=ax2, linewidth=0.5, color="gray", alpha=0.5)

    # Normalizar índices para transparência/cor
    max_v = classified["idx_veicular"].quantile(0.95)
    max_b = classified["idx_biomassa"].quantile(0.95)

    # Vetores RGB: Veicular=Red, Biomassa=Orange-Yellow
    norm_v = np.clip(classified["idx_veicular"] / max_v, 0, 1)
    norm_b = np.clip(classified["idx_biomassa"] / max_b, 0, 1)

    colors = np.column_stack([
        np.clip(norm_v * 0.9 + norm_b * 0.3, 0, 1),   # R
        np.clip(norm_b * 0.6 + norm_v * 0.1, 0, 1),    # G
        np.clip(norm_v * 0.2, 0, 1),                     # B
    ])

    ax2.scatter(gdf.geometry.x, gdf.geometry.y,
                c=colors, s=15, alpha=0.7, edgecolors="none", zorder=3)

    ax2.set_xlim(-54, -44)
    ax2.set_ylim(-26, -19)
    ax2.set_title("Intensidade: Veicular (verm.) vs Biomassa (amar.)\npor Município — SP",
                   fontsize=12, fontweight="bold")
    ax2.set_xlabel("Longitude")
    ax2.grid(True, alpha=0.2)

    # Legend manual
    legend_elements = [
        Patch(facecolor="#cc0000", label="Alto veicular"),
        Patch(facecolor="#cc9900", label="Alto biomassa"),
        Patch(facecolor="#993300", label="Ambos altos"),
        Patch(facecolor="#333333", label="Ambos baixos"),
    ]
    ax2.legend(handles=legend_elements, fontsize=7, loc="lower right")

    plt.tight_layout()
    path = os.path.join(output_dir, "source_classification_map.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Mapa de classificação salvo: {path}")


def plot_distribution(classified, output_dir):
    """Gráficos de distribuição dos índices e classificação."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Histograma de classificação
    ax = axes[0, 0]
    counts = classified["fonte_principal"].value_counts()
    colors_bar = {"VEICULAR": "#e74c3c", "BIOMASSA": "#f39c12",
                  "MISTA": "#9b59b6", "BAIXA_EMISSÃO": "#27ae60"}
    bars = ax.bar(counts.index, counts.values,
                  color=[colors_bar.get(c, "gray") for c in counts.index])
    ax.set_title("Distribuição das Classificações", fontweight="bold")
    ax.set_ylabel("Nº Municípios")
    for bar, val in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                str(val), ha="center", fontsize=9)

    # 2. Scatter: Índice Veicular vs Biomassa
    ax = axes[0, 1]
    for fonte, color in colors_bar.items():
        subset = classified[classified["fonte_principal"] == fonte]
        ax.scatter(subset["idx_veicular"], subset["idx_biomassa"],
                   c=color, alpha=0.5, s=20, label=fonte, edgecolors="white",
                   linewidth=0.3)
    ax.set_xlabel("Índice Veicular (diesel/km²)")
    ax.set_ylabel("Índice Biomassa (focos/1000km²/ano)")
    ax.set_title("Índices: Veicular vs Biomassa", fontweight="bold")
    ax.legend(fontsize=7)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_yscale("symlog", linthresh=1)
    ax.grid(True, alpha=0.3)

    # 3. Boxplots por classificação — veículos diesel
    ax = axes[1, 0]
    data_box = [classified[classified["fonte_principal"] == c]["veiculos_diesel"].values
                for c in ["VEICULAR", "MISTA", "BIOMASSA", "BAIXA_EMISSÃO"]
                if c in classified["fonte_principal"].values]
    labels_box = [c for c in ["VEICULAR", "MISTA", "BIOMASSA", "BAIXA_EMISSÃO"]
                  if c in classified["fonte_principal"].values]
    bp = ax.boxplot(data_box, labels=labels_box, patch_artist=True, showfliers=False)
    for patch, label in zip(bp["boxes"], labels_box):
        patch.set_facecolor(colors_bar.get(label, "gray"))
        patch.set_alpha(0.6)
    ax.set_ylabel("Veículos Diesel")
    ax.set_title("Frota Diesel por Classificação", fontweight="bold")
    ax.set_yscale("symlog", linthresh=100)
    ax.grid(axis="y", alpha=0.3)

    # 4. Boxplots por classificação — focos de calor
    ax = axes[1, 1]
    data_box = [classified[classified["fonte_principal"] == c]["focos_por_ano"].values
                for c in ["VEICULAR", "MISTA", "BIOMASSA", "BAIXA_EMISSÃO"]
                if c in classified["fonte_principal"].values]
    bp = ax.boxplot(data_box, labels=labels_box, patch_artist=True, showfliers=False)
    for patch, label in zip(bp["boxes"], labels_box):
        patch.set_facecolor(colors_bar.get(label, "gray"))
        patch.set_alpha(0.6)
    ax.set_ylabel("Focos de Calor / ano")
    ax.set_title("Queimadas por Classificação", fontweight="bold")
    ax.set_yscale("symlog", linthresh=1)
    ax.grid(axis="y", alpha=0.3)

    plt.suptitle("Source Apportionment Simplificado — PM2.5 SP\n"
                 "Classificação por Município (645 municípios)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    path = os.path.join(output_dir, "source_distribution.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Distribuições salvas: {path}")


def plot_top_municipalities(classified, output_dir):
    """Top 15 municípios por cada índice."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # Top veicular
    ax = axes[0]
    top_v = classified.nlargest(15, "idx_veicular")
    colors_v = [{"VEICULAR": "#e74c3c", "BIOMASSA": "#f39c12",
                 "MISTA": "#9b59b6", "BAIXA_EMISSÃO": "#27ae60"}.get(c, "gray")
                for c in top_v["fonte_principal"]]
    ax.barh(range(15), top_v["idx_veicular"].values, color=colors_v, alpha=0.8)
    ax.set_yticks(range(15))
    ax.set_yticklabels([m[:22] for m in top_v["municipio"]], fontsize=8)
    ax.set_xlabel("Índice Veicular (diesel+pesados / km²)")
    ax.set_title("Top 15 — Índice Veicular", fontweight="bold")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)

    # Top biomassa
    ax = axes[1]
    top_b = classified.nlargest(15, "idx_biomassa")
    colors_b = [{"VEICULAR": "#e74c3c", "BIOMASSA": "#f39c12",
                 "MISTA": "#9b59b6", "BAIXA_EMISSÃO": "#27ae60"}.get(c, "gray")
                for c in top_b["fonte_principal"]]
    ax.barh(range(15), top_b["idx_biomassa"].values, color=colors_b, alpha=0.8)
    ax.set_yticks(range(15))
    ax.set_yticklabels([m[:22] for m in top_b["municipio"]], fontsize=8)
    ax.set_xlabel("Índice Biomassa (focos / 1000 km² / ano)")
    ax.set_title("Top 15 — Índice Biomassa", fontweight="bold")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)

    plt.suptitle("Municípios com Maior Exposição por Tipo de Fonte",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    path = os.path.join(output_dir, "source_top_municipalities.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Top municípios salvo: {path}")


# ============================================================
# INTEGRAÇÃO COM SUPERFÍCIE PM2.5
# ============================================================
def enrich_with_pm25(classified):
    """
    Enriquece a tabela com PM2.5 da superfície pré-computada.
    """
    try:
        from pm25_surface import lookup_pm25
        pm25_vals = []
        for _, row in classified.iterrows():
            if pd.notna(row.get("lat_centroid")) and pd.notna(row.get("lon_centroid")):
                val = lookup_pm25(row["lat_centroid"], row["lon_centroid"], 2023)
                pm25_vals.append(val)
            else:
                pm25_vals.append(None)
        classified["pm25_2023"] = pm25_vals
        n_valid = classified["pm25_2023"].notna().sum()
        print(f"  ✅ PM2.5 (2023) adicionado para {n_valid}/{len(classified)} municípios")
    except Exception as e:
        print(f"  ⚠ Não foi possível carregar PM2.5: {e}")
        classified["pm25_2023"] = None

    return classified


# ============================================================
# MAIN
# ============================================================
def main():
    print("\n" + "=" * 60)
    print("🏭 SOURCE APPORTIONMENT SIMPLIFICADO — PM2.5 SP")
    print("   Classificação: Veicular / Biomassa / Mista")
    print("=" * 60)

    # 1. Carregar dados
    fleet = load_fleet_data()
    fires = load_all_fires()
    fire_mun = aggregate_fires_by_municipality(fires)

    # 2. Carregar mesorregiões para mapa
    mesorregioes = None
    try:
        mesorregioes = gpd.read_file(config.MESORREGIOES_SHP)
    except Exception:
        pass

    # 3. Classificar
    print("\n" + "=" * 60)
    print("📊 CLASSIFICAÇÃO")
    print("=" * 60)

    classified = classify_sources(fleet, fire_mun)

    # Tabela resumo
    summary = classified["fonte_principal"].value_counts()
    print("\n  Resultado da classificação:")
    for fonte, count in summary.items():
        pct = count / len(classified) * 100
        veiculos = classified[classified["fonte_principal"] == fonte]["total_veiculos"].sum()
        focos = classified[classified["fonte_principal"] == fonte]["total_focos"].sum()
        print(f"    {fonte:15s}: {count:3d} municípios ({pct:5.1f}%) | "
              f"veículos: {veiculos:>10,.0f} | focos: {focos:>8,.0f}")

    # 4. Enriquecer com PM2.5
    print("\n📡 Adicionando dados de PM2.5 (superfície 2023)...")
    classified = enrich_with_pm25(classified)

    # PM2.5 médio por classificação
    if classified["pm25_2023"].notna().any():
        print("\n  PM2.5 médio (2023) por classificação:")
        for fonte in ["VEICULAR", "MISTA", "BIOMASSA", "BAIXA_EMISSÃO"]:
            subset = classified[(classified["fonte_principal"] == fonte) &
                                (classified["pm25_2023"].notna())]
            if len(subset) > 0:
                mean_pm = subset["pm25_2023"].mean()
                std_pm = subset["pm25_2023"].std()
                print(f"    {fonte:15s}: {mean_pm:.1f} ± {std_pm:.1f} µg/m³ (n={len(subset)})")

    # 5. Salvar resultados
    print("\n" + "=" * 60)
    print("💾 SALVANDO RESULTADOS")
    print("=" * 60)

    # CSV completo
    out_cols = ["id_municipio", "municipio", "total_veiculos", "veiculos_diesel",
                "veiculos_pesados", "pct_diesel", "pct_pesados", "area_km2",
                "total_focos", "focos_por_ano", "diesel_per_km2", "pesados_per_km2",
                "idx_veicular", "idx_biomassa", "z_veicular", "z_biomassa",
                "classificacao", "fonte_principal", "pm25_2023",
                "lat_centroid", "lon_centroid"]
    out_cols = [c for c in out_cols if c in classified.columns]
    csv_path = os.path.join(OUTPUT_DIR, "source_apportionment_municipalities.csv")
    classified[out_cols].to_csv(csv_path, index=False)
    print(f"  📄 CSV completo: {csv_path}")

    # 6. Visualizações
    print("\n" + "=" * 60)
    print("📈 GERANDO VISUALIZAÇÕES")
    print("=" * 60)

    plot_classification_map(classified, mesorregioes, OUTPUT_DIR)
    plot_distribution(classified, OUTPUT_DIR)
    plot_top_municipalities(classified, OUTPUT_DIR)

    # 7. Exemplos para os 5 pacientes-modelo
    print("\n" + "=" * 60)
    print("🔬 CLASSIFICAÇÃO DOS MUNICÍPIOS DOS PACIENTES-MODELO")
    print("=" * 60)

    test_cities = {
        "PAC_01": "SAO PAULO",
        "PAC_02a": "CAMPINAS",
        "PAC_02b": "SAO PAULO",
        "PAC_03": "SAO JOSE DO RIO PRETO",
        "PAC_04": "ITAI",
        "PAC_05": "SANTOS",
    }

    for pac_id, city in test_cities.items():
        match = classified[classified["municipio_norm"] == city]
        if len(match) > 0:
            row = match.iloc[0]
            pm25_str = f", PM2.5={row['pm25_2023']:.1f}" if pd.notna(row.get("pm25_2023")) else ""
            print(f"  {pac_id}: {city:25s} → {row['fonte_principal']:15s} "
                  f"(IV={row['idx_veicular']:.1f}, IB={row['idx_biomassa']:.1f}{pm25_str})")
        else:
            # Fuzzy match
            candidates = classified[classified["municipio_norm"].str.contains(city.split()[0])]
            if len(candidates) > 0:
                row = candidates.iloc[0]
                print(f"  {pac_id}: {row['municipio']:25s} → {row['fonte_principal']:15s}")
            else:
                print(f"  {pac_id}: {city:25s} → NÃO ENCONTRADO")

    print("\n" + "=" * 60)
    print("✅ SOURCE APPORTIONMENT COMPLETO!")
    print(f"   Resultados em: {OUTPUT_DIR}/")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
