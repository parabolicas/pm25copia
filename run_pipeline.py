#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_pipeline.py — Script principal do pipeline de exposição cumulativa ao PM2.5
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

Saídas em output/:
    - exposure_details.csv    — detalhes anuais por paciente
                                 (inclui fraction_of_year)
    - exposure_summary.csv    — resumo cumulativo por paciente × buffer
    - exposure_map.png        — mapa
    - validation_report.txt   — relatório de validação

Pré-requisito: superfícies anuais bias-corrected em output/surfaces/
(geradas por `python3 pm25_surface.py`).
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
from matplotlib.lines import Line2D

import config
from geocoder import geocode_patient
from data_loaders import DataCache
from cumulative import run_all_patients, residence_years, _residence_dates


def geocode_all_patients(patients):
    """Geocodifica todos os pacientes."""
    print("\n" + "=" * 60)
    print("📍 ETAPA 1: GEOCODIFICAÇÃO DOS CEPs")
    print("=" * 60)

    geocoded = []
    for p in patients:
        print(f"\n  {p['id']}: {p['name']}")
        gp = geocode_patient(p)
        for res in gp["residences"]:
            geo = res["geocoded"]
            print(f"    CEP {res['cep']} → ({geo['lat']:.4f}, {geo['lon']:.4f}) "
                  f"[{geo['source']}]")
        geocoded.append(gp)
    return geocoded


def check_surfaces_present():
    """Verifica se as superfícies anuais existem; alerta se faltarem."""
    surface_dir = os.path.join(config.OUTPUT_DIR, "surfaces")
    if not os.path.isdir(surface_dir):
        print("\n  ⚠ Diretório output/surfaces/ não encontrado.")
        print("    Rode primeiro: python3 pm25_surface.py")
        return False

    tifs = [f for f in os.listdir(surface_dir)
            if f.startswith("pm25_sp_") and f.endswith(".tif")]
    if not tifs:
        print("\n  ⚠ Nenhum GeoTIFF de superfície encontrado em output/surfaces/.")
        print("    Rode primeiro: python3 pm25_surface.py")
        return False

    print(f"\n  ✅ {len(tifs)} superfícies anuais encontradas em output/surfaces/")
    return True


def run_validation(summary_df, details_df, patients):
    """
    Executa testes de validação e retorna relatório.
    Inclui novo teste: soma das frações de ano = duração residencial.
    """
    print("\n" + "=" * 60)
    print("✅ ETAPA 4: VALIDAÇÃO")
    print("=" * 60)

    tests = []

    # Teste 1: PM2.5 SP centro coerente com literatura (Pereira 17-24, Vasconcelos 17.7)
    sp_data = details_df[
        (details_df["patient_id"] == "PAC_01") &
        (details_df["buffer_km"] == 5)
    ]
    if len(sp_data) > 0 and sp_data["pm25_ugm3"].notna().any():
        sp_mean = sp_data["pm25_ugm3"].mean()
        passed = 10 <= sp_mean <= 25
        tests.append({
            "test": "PM2.5 SP Centro (PAC_01, buffer 5km, fonte SURFACE)",
            "expected": "PM2.5 entre 10-25 µg/m³ (literatura: 17-24)",
            "actual": f"{sp_mean:.1f} µg/m³",
            "passed": passed,
        })

    # Teste 2: ≥80% SURFACE
    if "source" in details_df.columns and len(details_df) > 0:
        pct_surface = (details_df["source"] == "SURFACE").sum() / len(details_df) * 100
        passed = pct_surface >= 80
        tests.append({
            "test": "Fonte primária = SURFACE (não CAMS/MERRA-2 bruto)",
            "expected": "≥80% dos registros com source=SURFACE",
            "actual": f"{pct_surface:.0f}% SURFACE",
            "passed": passed,
        })

    # Teste 3: Fração de ano nunca > 1.0 (sem sobreposição residencial)
    if "fraction_of_year" in details_df.columns:
        # Soma de frações por paciente×ano×buffer não pode exceder 1.0
        per_year = (details_df.groupby(["patient_id", "year", "buffer_km"])
                    ["fraction_of_year"].sum().reset_index())
        max_overlap = per_year["fraction_of_year"].max()
        passed = max_overlap <= 1.001  # tolerância de arredondamento
        tests.append({
            "test": "Frações de ano não-sobrepostas (residências sequenciais)",
            "expected": "Soma de frações por (paciente, ano, buffer) ≤ 1.000",
            "actual": f"máx observado = {max_overlap:.4f}",
            "passed": passed,
        })

    # Teste 4: Soma de frações = duração residencial total para cada paciente
    for p in patients:
        try:
            duration_days = sum(
                (_residence_dates(r)[1] - _residence_dates(r)[0]).days + 1
                for r in p["residences"]
            )
            duration_years = duration_days / 365.25
        except Exception:
            continue

        # Pega de qualquer buffer (todos devem ter o mesmo total_years)
        p_summary = summary_df[summary_df["patient_id"] == p["id"]]
        if len(p_summary) == 0:
            continue
        observed = p_summary["total_years"].iloc[0]
        # Margem: 0.05 ano (~18 dias) — cobre bissexto e arredondamento
        passed = abs(observed - duration_years) < 0.05
        tests.append({
            "test": f"Soma de frações = duração residencial ({p['id']})",
            "expected": f"≈ {duration_years:.3f} anos (±0.05)",
            "actual": f"{observed:.3f} anos",
            "passed": passed,
        })

    # Teste 5: Variação entre buffers
    for pid in summary_df["patient_id"].unique():
        p_data = summary_df[summary_df["patient_id"] == pid]
        if len(p_data) >= 2:
            values = p_data["cumulative_exposure_ugm3_years"].values
            has_variation = np.std(values) > 0.01 or len(set(values)) > 1
            tests.append({
                "test": f"Variação entre buffers ({pid})",
                "expected": "Valores diferentes para 5/10/25 km",
                "actual": f"std={np.std(values):.2f}, valores={[round(v,1) for v in values]}",
                "passed": has_variation,
            })

    # Teste 6: Coordenadas em SP
    for _, row in details_df.drop_duplicates(subset=["patient_id", "cep"]).iterrows():
        lat, lon = row["lat"], row["lon"]
        in_sp = (-26 <= lat <= -19) and (-54 <= lon <= -44)
        tests.append({
            "test": f"Coordenadas em SP ({row['patient_id']}, {row['cep']})",
            "expected": "lat ∈ [-26,-19], lon ∈ [-54,-44]",
            "actual": f"({lat:.4f}, {lon:.4f})",
            "passed": in_sp,
        })

    # Teste 7: Multi-endereço PAC_02
    pac02 = details_df[
        (details_df["patient_id"] == "PAC_02") &
        (details_df["buffer_km"] == 10)
    ]
    if len(pac02) > 0:
        n_cities = pac02["city"].nunique()
        tests.append({
            "test": "Multi-endereço PAC_02 (Campinas→SP)",
            "expected": "≥2 cidades diferentes",
            "actual": f"{n_cities} cidades: {pac02['city'].unique().tolist()}",
            "passed": n_cities >= 2,
        })

    # Teste 8: Demonstração de fracionamento (PAC_02 deve ter pelo menos um ano
    # com fração < 1.0 — o ano da mudança residencial)
    if "fraction_of_year" in details_df.columns:
        pac02 = details_df[
            (details_df["patient_id"] == "PAC_02") &
            (details_df["buffer_km"] == 5)
        ]
        if len(pac02) > 0:
            partial = pac02[(pac02["fraction_of_year"] > 0) &
                            (pac02["fraction_of_year"] < 0.999)]
            n_partial = len(partial)
            passed = n_partial >= 1
            tests.append({
                "test": "Fracionamento por mudança em meio de ano (PAC_02)",
                "expected": "≥1 registro com 0 < fração < 1.0",
                "actual": (f"{n_partial} registro(s) parcial(is); frações = "
                           f"{[round(f, 3) for f in partial['fraction_of_year'].tolist()]}"),
                "passed": passed,
            })

    # Imprimir e salvar
    n_pass = sum(1 for t in tests if t["passed"])
    n_total = len(tests)

    report_lines = []
    report_lines.append(f"\n  Resultados: {n_pass}/{n_total} testes passaram\n")
    for t in tests:
        icon = "✅" if t["passed"] else "❌"
        report_lines.append(f"  {icon} {t['test']}")
        report_lines.append(f"     Esperado: {t['expected']}")
        report_lines.append(f"     Obtido:   {t['actual']}")

    for line in report_lines:
        print(line)

    report_path = os.path.join(config.OUTPUT_DIR, "validation_report.txt")
    with open(report_path, "w") as f:
        f.write("RELATÓRIO DE VALIDAÇÃO — Pipeline PM2.5\n")
        f.write("=" * 50 + "\n")
        for line in report_lines:
            f.write(line + "\n")
    print(f"\n  📄 Relatório salvo: {report_path}")

    return n_pass, n_total


def generate_map(patients, summary_df, details_df, cache):
    """Mapa com pacientes, buffers e estações CETESB."""
    print("\n" + "=" * 60)
    print("🗺️  ETAPA 5: GERANDO MAPA")
    print("=" * 60)

    fig, ax = plt.subplots(1, 1, figsize=(14, 10))

    try:
        meso = cache.get_mesorregioes()
        meso.plot(ax=ax, color="lightyellow", edgecolor="gray",
                  linewidth=0.5, alpha=0.7)
    except Exception as e:
        print(f"  ⚠ Não foi possível plotar mesorregiões: {e}")

    cetesb = cache.get_cetesb()
    cetesb_unique = cetesb.drop_duplicates(subset=["estacao_nome"])
    ax.scatter(cetesb_unique["lon"], cetesb_unique["lat"],
               c="green", s=20, alpha=0.6, marker="^", zorder=3,
               label="Estações CETESB")

    colors = ["#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6"]
    buffer_styles = {5: "-", 10: "--", 25: ":"}

    all_points = []
    for i, patient in enumerate(patients):
        for res in patient["residences"]:
            geo = res.get("geocoded", res)
            lat = geo.get("lat", res.get("lat"))
            lon = geo.get("lon", res.get("lon"))
            all_points.append((i, patient, res, lat, lon))

    import math
    OVERLAP_THRESHOLD_DEG = 0.05
    OFFSET_DEG = 0.08

    display_coords = {}
    for idx, (i, patient, res, lat, lon) in enumerate(all_points):
        overlap_group = []
        for prev_idx, (pi, pp, pr, plat, plon) in enumerate(all_points[:idx]):
            if (abs(plat - lat) < OVERLAP_THRESHOLD_DEG
                    and abs(plon - lon) < OVERLAP_THRESHOLD_DEG):
                overlap_group.append(prev_idx)

        if overlap_group:
            n = len(overlap_group) + 1
            angle = (2 * math.pi / max(n + 1, 3)) * n
            dlat = OFFSET_DEG * math.sin(angle)
            dlon = OFFSET_DEG * math.cos(angle)
            display_coords[(i, res.get("cep", ""))] = (lat + dlat, lon + dlon)
        else:
            display_coords[(i, res.get("cep", ""))] = (lat, lon)

    for i, patient in enumerate(patients):
        color = colors[i % len(colors)]
        for res in patient["residences"]:
            geo = res.get("geocoded", res)
            lat = geo.get("lat", res.get("lat"))
            lon = geo.get("lon", res.get("lon"))
            disp_lat, disp_lon = display_coords[(i, res.get("cep", ""))]

            ax.scatter(disp_lon, disp_lat, c=color, s=100, zorder=5,
                       edgecolors="black", linewidth=1)
            ax.annotate(f"{patient['id']}\n{res.get('city', '')}",
                        (disp_lon, disp_lat), textcoords="offset points",
                        xytext=(8, 8), fontsize=7, color=color, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.2",
                                  facecolor="white", alpha=0.8))

            for bkm, ls in buffer_styles.items():
                radius_deg = bkm / 111.0
                circle = plt.Circle((lon, lat), radius_deg, fill=False,
                                    color=color, linestyle=ls,
                                    linewidth=0.8, alpha=0.5)
                ax.add_patch(circle)

    ax.set_xlim(-54, -44)
    ax.set_ylim(-26, -19)
    ax.set_xlabel("Longitude", fontsize=11)
    ax.set_ylabel("Latitude", fontsize=11)
    ax.set_title("Pipeline PM2.5 — Pacientes-modelo e Buffers de Análise\n"
                 "Estado de São Paulo", fontsize=13, fontweight="bold")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    legend_elements = [
        Line2D([0], [0], marker="^", color="w", markerfacecolor="green",
               markersize=8, label="Estações CETESB"),
        Line2D([0], [0], color="gray", linestyle="-", label="Buffer 5 km"),
        Line2D([0], [0], color="gray", linestyle="--", label="Buffer 10 km"),
        Line2D([0], [0], color="gray", linestyle=":", label="Buffer 25 km"),
    ]
    for i, patient in enumerate(patients):
        legend_elements.append(
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=colors[i % len(colors)],
                   markersize=8, label=f"{patient['id']}: {patient['name'][:25]}")
        )
    ax.legend(handles=legend_elements, loc="lower left", fontsize=7,
              framealpha=0.9, ncol=1)

    plt.tight_layout()
    map_path = os.path.join(config.OUTPUT_DIR, "exposure_map.png")
    plt.savefig(map_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Mapa salvo: {map_path}")


def generate_comparison_chart(summary_df):
    """Gráfico comparativo."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    ax1 = axes[0]
    patients = summary_df["patient_id"].unique()
    buffers = sorted(summary_df["buffer_km"].unique())
    x = np.arange(len(patients))
    width = 0.25

    for i, bkm in enumerate(buffers):
        data = summary_df[summary_df["buffer_km"] == bkm]
        values = [data[data["patient_id"] == p]["cumulative_exposure_ugm3_years"].values[0]
                  if len(data[data["patient_id"] == p]) > 0 else 0
                  for p in patients]
        ax1.bar(x + i * width, values, width, label=f"{bkm} km", alpha=0.85)

    ax1.set_xlabel("Paciente", fontsize=11)
    ax1.set_ylabel("Exposição Cumulativa (µg/m³-anos)", fontsize=11)
    ax1.set_title("Exposição Cumulativa ao PM2.5\npor Paciente e Raio de Buffer",
                  fontweight="bold")
    ax1.set_xticks(x + width)
    ax1.set_xticklabels(patients, fontsize=9)
    ax1.legend(title="Buffer")
    ax1.grid(axis="y", alpha=0.3)

    ax2 = axes[1]
    for i, bkm in enumerate(buffers):
        data = summary_df[summary_df["buffer_km"] == bkm]
        values = [data[data["patient_id"] == p]["mean_annual_pm25_ugm3"].values[0]
                  if (len(data[data["patient_id"] == p]) > 0
                      and pd.notna(data[data["patient_id"] == p]["mean_annual_pm25_ugm3"].values[0]))
                  else 0
                  for p in patients]
        ax2.bar(x + i * width, values, width, label=f"{bkm} km", alpha=0.85)

    ax2.axhline(y=15, color="red", linestyle="--", linewidth=1.5,
                label="OMS (15 µg/m³)")
    ax2.axhline(y=5, color="darkred", linestyle=":", linewidth=1,
                label="OMS ideal (5 µg/m³)")

    ax2.set_xlabel("Paciente", fontsize=11)
    ax2.set_ylabel("PM2.5 Médio Anual (µg/m³)", fontsize=11)
    ax2.set_title("PM2.5 Médio Anual\npor Paciente e Raio de Buffer",
                  fontweight="bold")
    ax2.set_xticks(x + width)
    ax2.set_xticklabels(patients, fontsize=9)
    ax2.legend(title="Buffer / Referência", fontsize=8)
    ax2.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    chart_path = os.path.join(config.OUTPUT_DIR, "exposure_comparison.png")
    plt.savefig(chart_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📄 Gráfico comparativo salvo: {chart_path}")


def print_summary_table(summary_df):
    """Tabela formatada — agora com total_years e valid_years como floats."""
    print("\n" + "=" * 100)
    print("📊 TABELA RESUMO — EXPOSIÇÃO CUMULATIVA AO PM2.5 (fonte: SURFACE bias-corrected)")
    print("=" * 100)
    print(f"{'Paciente':<10} {'Buffer':>8} {'Expo.Cum':>12} {'PM2.5 Méd':>11} "
          f"{'Anos Vál.':>10} {'Anos Tot.':>10} {'Queim.':>8} "
          f"{'%Surf.':>8} {'%Miss':>7}")
    print("-" * 100)

    for _, row in summary_df.iterrows():
        mean_str = (f"{row['mean_annual_pm25_ugm3']:.1f}"
                    if pd.notna(row['mean_annual_pm25_ugm3']) else "  N/D ")
        print(f"{row['patient_id']:<10} {row['buffer_km']:>6} km "
              f"{row['cumulative_exposure_ugm3_years']:>10.1f} "
              f"{mean_str:>11} "
              f"{row['valid_years']:>10.2f} "
              f"{row['total_years']:>10.2f} "
              f"{row['total_fire_foci']:>8} "
              f"{row['pct_surface']:>7.0f}% "
              f"{row['pct_missing']:>6.0f}%")
    print("=" * 100)


# ============================================================
# MAIN
# ============================================================
def main():
    print("\n" + "🔬" * 30)
    print("  PIPELINE DE EXPOSIÇÃO CUMULATIVA AO PM2.5")
    print("  Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP")
    print("  Fonte: superfície bias-corrected (MERRA-2 + CETESB)")
    print("  Histórico residencial: ponderação por fração-de-ano")
    print("🔬" * 30)

    cache = DataCache()

    if not check_surfaces_present():
        print("\n❌ ABORTANDO: rode primeiro `python3 pm25_surface.py`.\n")
        sys.exit(1)

    patients = geocode_all_patients(config.TEST_PATIENTS)

    print("\n" + "=" * 60)
    print("📦 ETAPA 2: CARREGANDO DADOS AMBIENTAIS")
    print("=" * 60)

    cache.get_cetesb()

    # Coletar todos os anos cobertos por qualquer residência (suporta ambos
    # date_start/date_end e year_start/year_end via residence_years)
    all_years = set()
    for p in patients:
        for res in p["residences"]:
            for y in residence_years(res):
                all_years.add(y)

    print(f"\n  Anos necessários: {sorted(all_years)}")
    print(f"  Total: {len(all_years)} anos calendário cobertos")

    print("\n" + "=" * 60)
    print("🧮 ETAPA 3: CALCULANDO EXPOSIÇÃO CUMULATIVA")
    print("=" * 60)

    details_df, summary_df = run_all_patients(
        patients, config.BUFFER_KM, cache
    )

    details_path = os.path.join(config.OUTPUT_DIR, "exposure_details.csv")
    summary_path = os.path.join(config.OUTPUT_DIR, "exposure_summary.csv")
    details_df.to_csv(details_path, index=False)
    summary_df.to_csv(summary_path, index=False)
    print(f"\n  📄 Detalhes salvos: {details_path}")
    print(f"  📄 Resumo salvo: {summary_path}")

    print_summary_table(summary_df)

    n_pass, n_total = run_validation(summary_df, details_df, patients)

    generate_map(patients, summary_df, details_df, cache)
    generate_comparison_chart(summary_df)

    print("\n" + "🏁" * 30)
    print(f"  PIPELINE COMPLETO!")
    print(f"  Validação: {n_pass}/{n_total} testes passaram")
    print(f"  Resultados em: {config.OUTPUT_DIR}/")
    print("🏁" * 30 + "\n")


if __name__ == "__main__":
    main()
