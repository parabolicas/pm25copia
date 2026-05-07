#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
source_apportionment.py — Classificação simplificada de fontes de PM2.5
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

ATUALIZAÇÃO v5: substituição da heurística de áreas (log-frota) por
ÁREAS OFICIAIS DO IBGE (Censo 2022, BR_Municipios_2022.shp). A heurística
antiga produzia erros de até 125× para municípios pequenos da RMSP
(ex: Águas de S. Pedro: 447 km² heurística vs. 3,6 km² IBGE) e ~20×
de subestimação para São Paulo capital (77 km² heurística vs. 1521 km²
IBGE). A correção inverte a classificação de vários municípios em
relação à versão anterior.

Classifica cada município de SP como:
  - VEICULAR: dominância de emissões veiculares (alta diesel/km²)
  - BIOMASSA: dominância de queima de biomassa (alta focos/km²)
  - MISTA: contribuição relevante de ambas as fontes
  - BAIXA_EMISSÃO: ambos os índices abaixo do limiar

Metodologia:
  1. Índice Veicular (IV): (diesel + 0.5*pesados) / km²
  2. Índice Biomassa (IB): focos / 1000 km² / ano
  3. Padronização robusta (z-score via mediana/MAD)
  4. Classificação por z-scores e diferenças relativas

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


# ============================================================
# CONFIGURAÇÃO
# ============================================================
OUTPUT_DIR = os.path.join(config.OUTPUT_DIR, "source_apportionment")
os.makedirs(OUTPUT_DIR, exist_ok=True)

THRESHOLD_DOMINANT = 0.5
THRESHOLD_RELEVANT = 0.2

# Aliases para reconciliar grafias entre DETRAN-SP e IBGE.
# Após normalização estrita (remove acentos + tudo que não é A-Z 0-9),
# os 7 casos abaixo permanecem com diferença genuína:
#   - Apóstrofes tratadas diferente ("D'Oeste" vs "Doeste"/"do Oeste")
#   - Variações Luís/Luiz (S vs Z)
#   - Diferenças de preposição (DE vs DO)
#   - Inclusão/omissão de "DO" antes de Aracanguá
# Todos confirmados via inspeção dos shapefiles oficiais (Mar/2023).
MUNICIPIO_ALIASES = {
    # DETRAN_NORM → IBGE_NORM
    "APARECIDADOOESTE": "APARECIDADOESTE",
    "BOMSUCESSODOITARARE": "BOMSUCESSODEITARARE",
    "GUARANIDOOESTE": "GUARANIDOESTE",
    "LUISIANIA": "LUIZIANIA",
    "SANTABARBARADOOESTE": "SANTABARBARADOESTE",
    "SANTOANTONIOARACANGUA": "SANTOANTONIODOARACANGUA",
    "SAOLUISDOPARAITINGA": "SAOLUIZDOPARAITINGA",
}


def _normalize_municipio(name):
    """
    Normalização estrita: maiúsculas + remove acentos + remove tudo que
    não é A-Z ou 0-9 (apóstrofes, hífens, espaços).
    """
    if pd.isna(name):
        return ""
    s = str(name).upper().strip()
    s = unicodedata.normalize("NFD", s)
    s = "".join(c for c in s if unicodedata.category(c) != "Mn")
    s = re.sub(r"[^A-Z0-9]", "", s)
    return s


def _apply_alias(name_norm):
    """Aplica alias DETRAN→IBGE quando aplicável."""
    return MUNICIPIO_ALIASES.get(name_norm, name_norm)


# ============================================================
# CARREGAMENTO DE DADOS
# ============================================================
def load_ibge_municipalities():
    """
    Carrega o shapefile oficial IBGE (Censo 2022) e filtra UF=SP.

    Returns:
        GeoDataFrame com colunas:
          - CD_MUN, NM_MUN, SIGLA_UF (originais IBGE)
          - AREA_KM2 (área oficial em km², calculada pelo IBGE)
          - municipio_norm_ibge (nome normalizado para merge)
          - geometry (polígono em EPSG:4674 SIRGAS 2000)
          - lon_centroid, lat_centroid (centroides oficiais)
    """
    print(f"\n🗺️  Carregando shapefile IBGE de municípios...")
    print(f"   Arquivo: {os.path.basename(config.MUNICIPIOS_SHP)}")
    gdf = gpd.read_file(config.MUNICIPIOS_SHP)
    print(f"   Total Brasil: {len(gdf)} municípios")

    sp = gdf[gdf["SIGLA_UF"] == "SP"].copy()
    print(f"   ✅ Filtrados {len(sp)} municípios em SP")

    sp["municipio_norm_ibge"] = sp["NM_MUN"].apply(_normalize_municipio)

    # Centroides oficiais (a partir dos polígonos)
    centroids = sp.geometry.centroid
    sp["lon_centroid"] = centroids.x
    sp["lat_centroid"] = centroids.y

    print(f"   Áreas IBGE: min={sp['AREA_KM2'].min():.1f}, "
          f"max={sp['AREA_KM2'].max():.1f}, "
          f"mediana={sp['AREA_KM2'].median():.1f} km²")
    print(f"   Total SP: {sp['AREA_KM2'].sum():.0f} km² "
          f"(referência IBGE: 248.219 km²)")

    return sp


def load_fleet_data():
    """Carrega e agrega frota por município."""
    print("\n🚗 Carregando dados de frota...")
    cols = ["id_municipio", "municipio", "combustivel", "tipo_veiculo",
            "quantidade_veiculos"]
    df = pd.read_csv(config.FROTA_CSV, usecols=cols)
    print(f"   ✅ {len(df):,} registros, {df['municipio'].nunique()} municípios")

    agg = df.groupby(["id_municipio", "municipio"]).agg(
        total_veiculos=("quantidade_veiculos", "sum"),
    ).reset_index()

    diesel = df[df["combustivel"] == "DIESEL"].groupby("municipio").agg(
        veiculos_diesel=("quantidade_veiculos", "sum"),
    ).reset_index()

    heavy_types = ["CAMINHAO", "CAMINHONETE", "UTILITARIO"]
    heavy = df[df["tipo_veiculo"].isin(heavy_types)].groupby("municipio").agg(
        veiculos_pesados=("quantidade_veiculos", "sum"),
    ).reset_index()

    fleet = agg.merge(diesel, on="municipio", how="left")
    fleet = fleet.merge(heavy, on="municipio", how="left")
    fleet = fleet.fillna(0)

    fleet["pct_diesel"] = (fleet["veiculos_diesel"] / fleet["total_veiculos"] * 100).round(1)
    fleet["pct_pesados"] = (fleet["veiculos_pesados"] / fleet["total_veiculos"] * 100).round(1)

    print(f"   📊 Total: {fleet['total_veiculos'].sum():,.0f} veículos, "
          f"{fleet['veiculos_diesel'].sum():,.0f} diesel")

    # Normalização + alias para merge com IBGE
    fleet["municipio_norm"] = fleet["municipio"].apply(_normalize_municipio)
    fleet["municipio_norm_ibge"] = fleet["municipio_norm"].apply(_apply_alias)

    return fleet


def load_all_fires():
    """Carrega todos os focos de calor (BDQueimadas)."""
    print("\n🔥 Carregando focos de calor (BDQueimadas)...")

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
    """Agrega focos de calor por município (com alias para merge IBGE)."""
    if len(fires) == 0:
        return pd.DataFrame()

    fire_mun = fires.groupby("municipio").agg(
        total_focos=("lat", "count"),
        n_anos_com_focos=("ano", "nunique"),
    ).reset_index()

    n_years = fires["ano"].nunique()
    fire_mun["focos_por_ano"] = (fire_mun["total_focos"] / n_years).round(1)
    fire_mun["municipio_norm"] = fire_mun["municipio"].apply(_normalize_municipio)
    fire_mun["municipio_norm_ibge"] = fire_mun["municipio_norm"].apply(_apply_alias)

    return fire_mun


# ============================================================
# CLASSIFICAÇÃO
# ============================================================
def classify_sources(fleet, fire_mun, ibge):
    """
    Classifica cada município como VEICULAR, BIOMASSA, MISTA ou BAIXA_EMISSÃO.

    A diferença em relação à versão anterior é o uso de AREAS OFICIAIS DO
    IBGE (coluna AREA_KM2) em vez da heurística log-frota.

    Args:
        fleet: DataFrame com frota por município (com municipio_norm_ibge)
        fire_mun: DataFrame com focos por município (com municipio_norm_ibge)
        ibge: GeoDataFrame com 645 municípios SP do IBGE

    Returns:
        GeoDataFrame com classificação + áreas + geometrias.
    """
    # Merge frota com IBGE (autoritativo)
    merged = ibge.merge(fleet, on="municipio_norm_ibge", how="left",
                        suffixes=("", "_detran"))

    # Merge focos
    if len(fire_mun) > 0:
        merged = merged.merge(
            fire_mun[["municipio_norm_ibge", "total_focos", "focos_por_ano"]],
            on="municipio_norm_ibge", how="left"
        )
    else:
        merged["total_focos"] = 0
        merged["focos_por_ano"] = 0

    # Preencher NaN (municípios sem focos = 0; sem frota = mantém NaN)
    merged["total_focos"] = merged["total_focos"].fillna(0)
    merged["focos_por_ano"] = merged["focos_por_ano"].fillna(0)
    merged["total_veiculos"] = merged["total_veiculos"].fillna(0)
    merged["veiculos_diesel"] = merged["veiculos_diesel"].fillna(0)
    merged["veiculos_pesados"] = merged["veiculos_pesados"].fillna(0)
    merged["pct_diesel"] = merged["pct_diesel"].fillna(0)
    merged["pct_pesados"] = merged["pct_pesados"].fillna(0)

    # Diagnóstico do match
    n_with_fleet = (merged["total_veiculos"] > 0).sum()
    n_with_fires = (merged["total_focos"] > 0).sum()
    print(f"\n   Match diagnóstico:")
    print(f"     {n_with_fleet}/{len(merged)} municípios IBGE com dados de frota")
    print(f"     {n_with_fires}/{len(merged)} municípios IBGE com focos de calor")

    # ÍNDICES (agora usando AREA_KM2 oficial!)
    merged["diesel_per_km2"] = merged["veiculos_diesel"] / merged["AREA_KM2"]
    merged["pesados_per_km2"] = merged["veiculos_pesados"] / merged["AREA_KM2"]
    merged["idx_veicular"] = merged["diesel_per_km2"] + 0.5 * merged["pesados_per_km2"]
    merged["idx_biomassa"] = merged["focos_por_ano"] / merged["AREA_KM2"] * 1000

    # z-scores robustos (mediana / MAD)
    for col in ["idx_veicular", "idx_biomassa"]:
        med = merged[col].median()
        mad = np.median(np.abs(merged[col] - med))
        if mad > 0:
            merged[f"z_{col.replace('idx_', '')}"] = (merged[col] - med) / (mad * 1.4826)
        else:
            merged[f"z_{col.replace('idx_', '')}"] = 0

    def _classify(row):
        z_v = row["z_veicular"]
        z_b = row["z_biomassa"]

        if z_v > THRESHOLD_RELEVANT and z_b > THRESHOLD_RELEVANT:
            if abs(z_v - z_b) < 0.5:
                return "MISTA"
            elif z_v > z_b:
                return "VEICULAR_MISTA"
            else:
                return "BIOMASSA_MISTA"

        if z_v > THRESHOLD_DOMINANT and z_b <= THRESHOLD_RELEVANT:
            return "VEICULAR"
        if z_b > THRESHOLD_DOMINANT and z_v <= THRESHOLD_RELEVANT:
            return "BIOMASSA"

        if z_v <= THRESHOLD_RELEVANT and z_b <= THRESHOLD_RELEVANT:
            return "BAIXA_EMISSÃO"

        if z_v > z_b:
            return "VEICULAR"
        return "BIOMASSA"

    merged["classificacao"] = merged.apply(_classify, axis=1)

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
    """Mapa de classificação por município (agora usando POLÍGONOS reais)."""
    color_map = {
        "VEICULAR": "#e74c3c",
        "BIOMASSA": "#f39c12",
        "MISTA": "#9b59b6",
        "BAIXA_EMISSÃO": "#27ae60",
    }

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # --- Mapa 1: Polígonos coloridos por classificação ---
    ax1 = axes[0]
    if mesorregioes is not None:
        mesorregioes.boundary.plot(ax=ax1, linewidth=0.8, color="black", alpha=0.6)

    for fonte, color in color_map.items():
        subset = classified[classified["fonte_principal"] == fonte]
        if len(subset) > 0:
            subset.plot(ax=ax1, color=color, edgecolor="white",
                        linewidth=0.2, alpha=0.85, label=fonte)

    ax1.set_xlim(-54, -44)
    ax1.set_ylim(-26, -19)
    ax1.set_title("Classificação de Fontes de PM2.5\n"
                  "por Município — SP (645 municípios, áreas IBGE oficial)",
                  fontsize=11, fontweight="bold")
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    ax1.legend(fontsize=8, loc="lower right",
               title="Fonte Principal", title_fontsize=9)
    ax1.grid(True, alpha=0.2)

    # --- Mapa 2: índice contínuo (escolha: idx_veicular) ---
    ax2 = axes[1]
    if mesorregioes is not None:
        mesorregioes.boundary.plot(ax=ax2, linewidth=0.8, color="black", alpha=0.6)

    # Escala log para idx_veicular (varia 0 a milhares)
    iv = classified["idx_veicular"].copy()
    vmax = float(np.nanpercentile(iv, 95))  # corta extremos para visualização
    classified.plot(
        ax=ax2,
        column="idx_veicular",
        cmap="OrRd",
        vmin=0,
        vmax=vmax,
        edgecolor="white",
        linewidth=0.2,
        legend=True,
        legend_kwds={"label": "Diesel + 0.5·pesados / km²",
                     "orientation": "horizontal",
                     "shrink": 0.7, "pad": 0.05},
    )

    ax2.set_xlim(-54, -44)
    ax2.set_ylim(-26, -19)
    ax2.set_title("Intensidade Veicular (diesel + pesados / km²)\n"
                  "por Município — SP",
                  fontsize=11, fontweight="bold")
    ax2.set_xlabel("Longitude")
    ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    path = os.path.join(output_dir, "source_classification_map.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Mapa de classificação salvo: {path}")


def plot_distribution(classified, output_dir):
    """Gráficos de distribuição."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors_bar = {"VEICULAR": "#e74c3c", "BIOMASSA": "#f39c12",
                  "MISTA": "#9b59b6", "BAIXA_EMISSÃO": "#27ae60"}

    # 1. Histograma de classificação
    ax = axes[0, 0]
    counts = classified["fonte_principal"].value_counts()
    bars = ax.bar(counts.index, counts.values,
                  color=[colors_bar.get(c, "gray") for c in counts.index])
    ax.set_title("Distribuição das Classificações", fontweight="bold")
    ax.set_ylabel("Nº Municípios")
    for bar, val in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                str(val), ha="center", fontsize=9)

    # 2. Scatter
    ax = axes[0, 1]
    for fonte, color in colors_bar.items():
        subset = classified[classified["fonte_principal"] == fonte]
        ax.scatter(subset["idx_veicular"], subset["idx_biomassa"],
                   c=color, alpha=0.5, s=20, label=fonte,
                   edgecolors="white", linewidth=0.3)
    ax.set_xlabel("Índice Veicular (diesel/km²)")
    ax.set_ylabel("Índice Biomassa (focos/1000km²/ano)")
    ax.set_title("Índices: Veicular vs Biomassa", fontweight="bold")
    ax.legend(fontsize=7)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_yscale("symlog", linthresh=1)
    ax.grid(True, alpha=0.3)

    # 3. Boxplots — diesel
    ax = axes[1, 0]
    order = ["VEICULAR", "MISTA", "BIOMASSA", "BAIXA_EMISSÃO"]
    data_box = [classified[classified["fonte_principal"] == c]["veiculos_diesel"].values
                for c in order if c in classified["fonte_principal"].values]
    labels_box = [c for c in order if c in classified["fonte_principal"].values]
    bp = ax.boxplot(data_box, labels=labels_box, patch_artist=True, showfliers=False)
    for patch, label in zip(bp["boxes"], labels_box):
        patch.set_facecolor(colors_bar.get(label, "gray"))
        patch.set_alpha(0.6)
    ax.set_ylabel("Veículos Diesel")
    ax.set_title("Frota Diesel por Classificação", fontweight="bold")
    ax.set_yscale("symlog", linthresh=100)
    ax.grid(axis="y", alpha=0.3)

    # 4. Boxplots — focos
    ax = axes[1, 1]
    data_box = [classified[classified["fonte_principal"] == c]["focos_por_ano"].values
                for c in order if c in classified["fonte_principal"].values]
    bp = ax.boxplot(data_box, labels=labels_box, patch_artist=True, showfliers=False)
    for patch, label in zip(bp["boxes"], labels_box):
        patch.set_facecolor(colors_bar.get(label, "gray"))
        patch.set_alpha(0.6)
    ax.set_ylabel("Focos de Calor / ano")
    ax.set_title("Queimadas por Classificação", fontweight="bold")
    ax.set_yscale("symlog", linthresh=1)
    ax.grid(axis="y", alpha=0.3)

    plt.suptitle("Source Apportionment — PM2.5 SP\n"
                 "Classificação por Município (645 municípios, áreas IBGE)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    path = os.path.join(output_dir, "source_distribution.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Distribuições salvas: {path}")


def plot_top_municipalities(classified, output_dir):
    """Top 15 municípios por cada índice."""
    colors_bar = {"VEICULAR": "#e74c3c", "BIOMASSA": "#f39c12",
                  "MISTA": "#9b59b6", "BAIXA_EMISSÃO": "#27ae60"}

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # Top veicular
    ax = axes[0]
    top_v = classified.nlargest(15, "idx_veicular")
    colors_v = [colors_bar.get(c, "gray") for c in top_v["fonte_principal"]]
    ax.barh(range(15), top_v["idx_veicular"].values, color=colors_v, alpha=0.8)
    ax.set_yticks(range(15))
    ax.set_yticklabels([m[:22] for m in top_v["NM_MUN"]], fontsize=8)
    ax.set_xlabel("Índice Veicular (diesel + 0.5·pesados / km²)")
    ax.set_title("Top 15 — Índice Veicular", fontweight="bold")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)

    # Top biomassa
    ax = axes[1]
    top_b = classified.nlargest(15, "idx_biomassa")
    colors_b = [colors_bar.get(c, "gray") for c in top_b["fonte_principal"]]
    ax.barh(range(15), top_b["idx_biomassa"].values, color=colors_b, alpha=0.8)
    ax.set_yticks(range(15))
    ax.set_yticklabels([m[:22] for m in top_b["NM_MUN"]], fontsize=8)
    ax.set_xlabel("Índice Biomassa (focos / 1000 km² / ano)")
    ax.set_title("Top 15 — Índice Biomassa", fontweight="bold")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)

    plt.suptitle("Municípios com Maior Exposição por Tipo de Fonte (áreas IBGE)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    path = os.path.join(output_dir, "source_top_municipalities.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Top municípios salvo: {path}")


# ============================================================
# INTEGRAÇÃO COM PM2.5
# ============================================================
def enrich_with_pm25(classified):
    """Adiciona PM2.5 (2023) às coordenadas dos centroides oficiais IBGE."""
    try:
        from pm25_surface import lookup_pm25
        pm25_vals = []
        for _, row in classified.iterrows():
            if pd.notna(row["lat_centroid"]) and pd.notna(row["lon_centroid"]):
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
    print("🏭 SOURCE APPORTIONMENT — PM2.5 SP")
    print("   Classificação: Veicular / Biomassa / Mista")
    print("   Áreas: IBGE Censo 2022 (oficial, BR_Municipios_2022.shp)")
    print("=" * 60)

    # 1. Carregar dados
    ibge = load_ibge_municipalities()
    fleet = load_fleet_data()
    fires = load_all_fires()
    fire_mun = aggregate_fires_by_municipality(fires)

    # 2. Mesorregiões para mapa de fundo
    mesorregioes = None
    try:
        mesorregioes = gpd.read_file(config.MESORREGIOES_SHP)
    except Exception:
        pass

    # 3. Classificar
    print("\n" + "=" * 60)
    print("📊 CLASSIFICAÇÃO")
    print("=" * 60)

    classified = classify_sources(fleet, fire_mun, ibge)

    # Tabela resumo
    summary = classified["fonte_principal"].value_counts()
    print("\n  Resultado da classificação:")
    for fonte, count in summary.items():
        pct = count / len(classified) * 100
        veiculos = classified[classified["fonte_principal"] == fonte]["total_veiculos"].sum()
        focos = classified[classified["fonte_principal"] == fonte]["total_focos"].sum()
        area = classified[classified["fonte_principal"] == fonte]["AREA_KM2"].sum()
        print(f"    {fonte:15s}: {count:3d} mun. ({pct:5.1f}%) | "
              f"área={area:>8.0f} km² | "
              f"veíc={veiculos:>10,.0f} | focos={focos:>8,.0f}")

    # 4. Enriquecer com PM2.5
    print("\n📡 Adicionando PM2.5 (superfície 2023)...")
    classified = enrich_with_pm25(classified)

    if classified["pm25_2023"].notna().any():
        print("\n  PM2.5 médio (2023) por classificação:")
        for fonte in ["VEICULAR", "MISTA", "BIOMASSA", "BAIXA_EMISSÃO"]:
            subset = classified[(classified["fonte_principal"] == fonte) &
                                (classified["pm25_2023"].notna())]
            if len(subset) > 0:
                mean_pm = subset["pm25_2023"].mean()
                std_pm = subset["pm25_2023"].std()
                print(f"    {fonte:15s}: {mean_pm:.1f} ± {std_pm:.1f} µg/m³ "
                      f"(n={len(subset)})")

    # 5. Salvar resultados
    print("\n" + "=" * 60)
    print("💾 SALVANDO RESULTADOS")
    print("=" * 60)

    out_cols = ["CD_MUN", "NM_MUN", "SIGLA_UF", "AREA_KM2",
                "id_municipio", "municipio",
                "total_veiculos", "veiculos_diesel", "veiculos_pesados",
                "pct_diesel", "pct_pesados",
                "total_focos", "focos_por_ano",
                "diesel_per_km2", "pesados_per_km2",
                "idx_veicular", "idx_biomassa",
                "z_veicular", "z_biomassa",
                "classificacao", "fonte_principal",
                "pm25_2023",
                "lat_centroid", "lon_centroid"]
    out_cols = [c for c in out_cols if c in classified.columns]
    csv_path = os.path.join(OUTPUT_DIR, "source_apportionment_municipalities.csv")
    classified[out_cols].to_csv(csv_path, index=False)
    print(f"  📄 CSV completo: {csv_path}")

    # GeoPackage com polígonos (formato moderno, abre em QGIS)
    try:
        gpkg_path = os.path.join(OUTPUT_DIR, "source_apportionment_municipalities.gpkg")
        cols_gpkg = out_cols + ["geometry"]
        cols_gpkg = [c for c in cols_gpkg if c in classified.columns]
        classified[cols_gpkg].to_file(gpkg_path, driver="GPKG")
        print(f"  📄 GeoPackage (polígonos): {gpkg_path}")
    except Exception as e:
        print(f"  ⚠ Não foi possível salvar GeoPackage: {e}")

    # 6. Visualizações
    print("\n" + "=" * 60)
    print("📈 GERANDO VISUALIZAÇÕES")
    print("=" * 60)

    plot_classification_map(classified, mesorregioes, OUTPUT_DIR)
    plot_distribution(classified, OUTPUT_DIR)
    plot_top_municipalities(classified, OUTPUT_DIR)

    # 7. Pacientes-modelo
    print("\n" + "=" * 60)
    print("🔬 CLASSIFICAÇÃO DOS MUNICÍPIOS DOS PACIENTES-MODELO")
    print("=" * 60)

    test_cities = {
        "PAC_01": "São Paulo",
        "PAC_02a": "Campinas",
        "PAC_02b": "São Paulo",
        "PAC_03": "São José do Rio Preto",
        "PAC_04": "Itaí",
        "PAC_05": "Santos",
    }

    for pac_id, city in test_cities.items():
        city_norm = _apply_alias(_normalize_municipio(city))
        match = classified[classified["municipio_norm_ibge"] == city_norm]
        if len(match) > 0:
            row = match.iloc[0]
            pm25_str = (f", PM2.5={row['pm25_2023']:.1f}"
                        if pd.notna(row.get("pm25_2023")) else "")
            print(f"  {pac_id}: {row['NM_MUN']:25s} → {row['fonte_principal']:15s} "
                  f"(área={row['AREA_KM2']:.1f} km², "
                  f"IV={row['idx_veicular']:.2f}, IB={row['idx_biomassa']:.2f}{pm25_str})")
        else:
            print(f"  {pac_id}: {city:25s} → NÃO ENCONTRADO")

    print("\n" + "=" * 60)
    print("✅ SOURCE APPORTIONMENT COMPLETO!")
    print(f"   Resultados em: {OUTPUT_DIR}/")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
