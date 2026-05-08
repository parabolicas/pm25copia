#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sensibility_h.py — Análise de Sensibilidade (h):
                   Incerteza Temporal da Data de Mudança Residencial
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

OBJETIVO:
  Quantificar o viés residual na exposição cumulativa devido à imprecisão da
  data de mudança residencial. Em pacientes reais com histórico residencial
  imperfeito (apenas mês ou ano de mudança), três cenários podem ser
  contemplados:

    h1 (extremo precoce): mudança em 01-jan do ano declarado
                          → todo o ano em residência NOVA
    h2 (central, default): mudança em 01-jul (meio do ano)
                          → divisão proporcional ~0.5 / 0.5
    h3 (extremo tardio): mudança em 31-dez do ano declarado
                          → todo o ano em residência ANTERIOR

  Demonstrado com PAC_02 (Campinas → São Paulo em 2017).

METODOLOGIA:
  Recalcula a exposição cumulativa analiticamente a partir de
  output/exposure_details.csv (resultados do pipeline em cenário central),
  substituindo apenas a fração-de-ano do PAC_02 no ano 2017. Como o PM2.5
  anual de cada localização é invariante à data de mudança, basta combinar
  os valores existentes com as três frações alternativas:

    Cenário h1 (jan): fração_Camp=0.000, fração_SP=1.000
    Cenário h2 (jul): fração_Camp=0.496, fração_SP=0.504  (default atual)
    Cenário h3 (dez): fração_Camp=1.000, fração_SP=0.000

USO:
    python3 sensibility_h.py

SAÍDA:
    output/sensibility_h_results.csv  — tabela detalhada
    Console: tabela formatada para inclusão no protocolo
"""
import os
import pandas as pd
import config

DETAILS_CSV = os.path.join(config.OUTPUT_DIR, "exposure_details.csv")
OUTPUT_CSV = os.path.join(config.OUTPUT_DIR, "sensibility_h_results.csv")


def main():
    print("\n" + "=" * 70)
    print("🔬 ANÁLISE DE SENSIBILIDADE *h* — Incerteza Temporal de Mudança")
    print("   Cenários: 01-jan / 01-jul (default) / 31-dez")
    print("   Demonstrado com PAC_02 (Campinas → São Paulo, mudança em 2017)")
    print("=" * 70)

    if not os.path.exists(DETAILS_CSV):
        print(f"\n❌ {DETAILS_CSV} não encontrado.")
        print("   Rode primeiro: python3 run_pipeline.py")
        return

    df = pd.read_csv(DETAILS_CSV)
    print(f"\n📂 Carregado: {len(df)} registros de exposure_details.csv")

    # Filtrar PAC_02
    pac02 = df[df["patient_id"] == "PAC_02"].copy()
    if pac02.empty:
        print("❌ PAC_02 não encontrado no CSV.")
        return

    # ============================================================
    # Para cada buffer, calcular as 3 versões da exposição cumulativa
    # ============================================================
    results = []
    for buffer_km in sorted(pac02["buffer_km"].unique()):
        sub = pac02[pac02["buffer_km"] == buffer_km].copy()

        # Anos em Campinas (2012-2017): cidade = Campinas
        camp = sub[sub["city"].str.contains("Campinas", na=False)].copy()
        # Anos em SP (2017-2024): cidade = São Paulo
        sp = sub[sub["city"].str.contains("São Paulo|SP", na=False, regex=True)].copy()

        if camp.empty or sp.empty:
            print(f"  ⚠ Buffer {buffer_km}km: dados de Campinas ou SP ausentes.")
            continue

        # PM2.5 anual em cada localização nos anos relevantes
        camp_2017 = camp[camp["year"] == 2017]["pm25_ugm3"].iloc[0] if (camp["year"] == 2017).any() else None
        sp_2017 = sp[sp["year"] == 2017]["pm25_ugm3"].iloc[0] if (sp["year"] == 2017).any() else None

        if camp_2017 is None or sp_2017 is None:
            print(f"  ⚠ Buffer {buffer_km}km: PM2.5 2017 não disponível.")
            continue

        # Soma das contribuições FORA do ano 2017 (não muda entre cenários)
        # Campinas 2012-2016: fração=1.0 cada
        camp_outros_anos = camp[camp["year"] != 2017]["pm25_ugm3"].sum()
        # SP 2018-2024: fração=1.0 cada
        sp_outros_anos = sp[sp["year"] != 2017]["pm25_ugm3"].sum()
        base_exposure = camp_outros_anos + sp_outros_anos

        # Os 3 cenários
        scenarios = {
            "h1_jan_01": {
                "frac_camp": 0.000, "frac_sp": 1.000,
                "descricao": "Mudança 01-jan-2017 (todo 2017 em SP)",
            },
            "h2_jul_01_default": {
                "frac_camp": 181.0 / 365, "frac_sp": 184.0 / 365,
                "descricao": "Mudança 01-jul-2017 (atual; meio do ano)",
            },
            "h3_dec_31": {
                "frac_camp": 1.000, "frac_sp": 0.000,
                "descricao": "Mudança 31-dez-2017 (todo 2017 em Campinas)",
            },
        }

        for scenario_id, params in scenarios.items():
            contrib_2017 = (camp_2017 * params["frac_camp"]
                            + sp_2017 * params["frac_sp"])
            cum_exposure = round(base_exposure + contrib_2017, 2)
            mean_pm25 = round(cum_exposure / 13.0, 2)  # PAC_02 cobre 13 anos

            results.append({
                "scenario": scenario_id,
                "buffer_km": buffer_km,
                "frac_campinas_2017": round(params["frac_camp"], 3),
                "frac_sp_2017": round(params["frac_sp"], 3),
                "pm25_camp_2017": camp_2017,
                "pm25_sp_2017": sp_2017,
                "contrib_2017": round(contrib_2017, 2),
                "cumulative_exposure": cum_exposure,
                "mean_annual_pm25": mean_pm25,
                "descricao": params["descricao"],
            })

    if not results:
        print("\n❌ Sem dados suficientes para análise.")
        return

    rdf = pd.DataFrame(results)

    # Calcular diferenças relativas vs. cenário central (h2)
    central = rdf[rdf["scenario"] == "h2_jul_01_default"].set_index("buffer_km")[
        "cumulative_exposure"]
    rdf["delta_abs_vs_central"] = rdf.apply(
        lambda r: round(r["cumulative_exposure"] - central[r["buffer_km"]], 2),
        axis=1
    )
    rdf["delta_pct_vs_central"] = rdf.apply(
        lambda r: round(
            (r["cumulative_exposure"] - central[r["buffer_km"]])
            / central[r["buffer_km"]] * 100, 2
        ),
        axis=1
    )

    # ============================================================
    # IMPRIMIR TABELA RESUMO PARA O PROTOCOLO
    # ============================================================
    print("\n" + "=" * 80)
    print("📊 TABELA — Análise de Sensibilidade *h* (PAC_02, mudança 2017)")
    print("=" * 80)
    print(f"\n{'Cenário':<22} {'Buffer':>8} {'Cum.Exp':>10} {'PM2.5 méd':>12} "
          f"{'Δ vs central':>14} {'Δ %':>8}")
    print("-" * 80)

    for buffer_km in sorted(rdf["buffer_km"].unique()):
        for scenario_id in ["h1_jan_01", "h2_jul_01_default", "h3_dec_31"]:
            row = rdf[
                (rdf["buffer_km"] == buffer_km) & (rdf["scenario"] == scenario_id)
            ].iloc[0]
            scenario_short = {
                "h1_jan_01": "h1: 01-jan (extremo)",
                "h2_jul_01_default": "h2: 01-jul (default)",
                "h3_dec_31": "h3: 31-dez (extremo)",
            }[scenario_id]
            delta_str = (f"{row['delta_abs_vs_central']:+.2f}"
                         if row['delta_abs_vs_central'] != 0 else "0.00")
            pct_str = (f"{row['delta_pct_vs_central']:+.2f}%"
                       if row['delta_pct_vs_central'] != 0 else "0.00%")
            print(f"{scenario_short:<22} {buffer_km:>6} km "
                  f"{row['cumulative_exposure']:>10.2f} "
                  f"{row['mean_annual_pm25']:>12.2f} "
                  f"{delta_str:>14} "
                  f"{pct_str:>8}")
        print("-" * 80)

    # ============================================================
    # SUMÁRIO QUANTITATIVO
    # ============================================================
    extreme_deltas = rdf[rdf["scenario"] != "h2_jul_01_default"]["delta_abs_vs_central"]
    pct_extreme = rdf[rdf["scenario"] != "h2_jul_01_default"]["delta_pct_vs_central"]

    print("\n" + "=" * 80)
    print("🔍 INTERPRETAÇÃO QUANTITATIVA")
    print("=" * 80)
    print(f"  Range absoluto entre extremos: "
          f"{extreme_deltas.abs().max():.2f} µg/m³·anos")
    print(f"  Range relativo entre extremos: ±{pct_extreme.abs().max():.2f}%")
    print(f"  Magnitude do viés residual (vs. exposição cum. ~218 µg/m³·anos): "
          f"≤{pct_extreme.abs().max():.1f}%")
    print(f"\n  Conclusão: a incerteza temporal de mudança em UM único ano "
          f"em pacientes com 13 anos de exposição cumulativa produz viés "
          f"residual <{pct_extreme.abs().max():.0f}% — pequeno em relação "
          f"aos efeitos esperados da literatura (HR 1,10-1,20 por 10 µg/m³).")

    # ============================================================
    # SALVAR
    # ============================================================
    rdf.to_csv(OUTPUT_CSV, index=False)
    print(f"\n  📄 Tabela detalhada salva: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
