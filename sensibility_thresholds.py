#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sensibility_thresholds.py — Análise de Sensibilidade dos Thresholds
                             do Source Apportionment (item 11)
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

OBJETIVO:
  Os limiares THRESHOLD_DOMINANT=0,5 e THRESHOLD_RELEVANT=0,2 utilizados
  em source_apportionment.classify_sources são valores convencionais
  derivados de z-scores robustos (mediana/MAD), porém podem ser
  questionados como arbitrários. Esta análise re-classifica os 645
  municípios sob 4 combinações de thresholds, quantifica a estabilidade
  da classificação e identifica os municípios "sensíveis" (que mudam de
  categoria entre cenários).

PROCEDIMENTO:
  Como source_apportionment_municipalities.csv já contém z_veicular e
  z_biomassa pré-computados, aplicamos a lógica de classificação
  diretamente sobre o CSV existente — sem necessidade de re-rodar
  o pipeline pesado.

Cenários testados:
    1. (0,4 ; 0,1) — mais permissivo (mais classificados como VEIC/BIO/MISTA)
    2. (0,5 ; 0,2) — DEFAULT (referência)
    3. (0,6 ; 0,3) — mais restritivo (mais classificados como BAIXA_EMISSÃO)
    4. (0,7 ; 0,1) — assimétrico (dominância estrita, relevância frouxa)

USO:
    python3 sensibility_thresholds.py

SAÍDA:
    output/sensibility_thresholds_summary.csv     — distribuição por cenário
    output/sensibility_thresholds_concordance.csv — matriz % concordância
    output/sensibility_thresholds_changers.csv    — municípios que mudam
    Console: tabelas formatadas para inclusão no protocolo §4.6
"""
import os
import pandas as pd
import config


SA_CSV = os.path.join(
    config.OUTPUT_DIR, "source_apportionment",
    "source_apportionment_municipalities.csv"
)
OUT_DIR = config.OUTPUT_DIR
SUM_CSV = os.path.join(OUT_DIR, "sensibility_thresholds_summary.csv")
CONC_CSV = os.path.join(OUT_DIR, "sensibility_thresholds_concordance.csv")
CHG_CSV = os.path.join(OUT_DIR, "sensibility_thresholds_changers.csv")


def classify(z_v, z_b, t_dom, t_rel):
    """
    Replica a lógica de source_apportionment.classify_sources, mas com
    thresholds parametrizáveis. Retorna: VEICULAR, BIOMASSA, MISTA ou
    BAIXA_EMISSÃO.
    """
    # Ambos altos → MISTA (com possível dominância)
    if z_v > t_rel and z_b > t_rel:
        if abs(z_v - z_b) < 0.5:
            return "MISTA"
        elif z_v > z_b:
            return "VEICULAR_MISTA"
        else:
            return "BIOMASSA_MISTA"
    # Dominância clara
    if z_v > t_dom and z_b <= t_rel:
        return "VEICULAR"
    if z_b > t_dom and z_v <= t_rel:
        return "BIOMASSA"
    # Baixa emissão
    if z_v <= t_rel and z_b <= t_rel:
        return "BAIXA_EMISSÃO"
    # Outros: vence o maior
    if z_v > z_b:
        return "VEICULAR"
    return "BIOMASSA"


def simplify(c):
    if "VEICULAR" in c and "MISTA" not in c:
        return "VEICULAR"
    if "BIOMASSA" in c and "MISTA" not in c:
        return "BIOMASSA"
    if "MISTA" in c:
        return "MISTA"
    return "BAIXA_EMISSÃO"


def main():
    print("\n" + "=" * 72)
    print("🔬 SENSIBILIDADE DOS THRESHOLDS DO SOURCE APPORTIONMENT (item 11)")
    print("=" * 72)

    if not os.path.exists(SA_CSV):
        print(f"\n❌ Arquivo não encontrado: {SA_CSV}")
        print("   Rode primeiro: python3 source_apportionment.py")
        return

    df = pd.read_csv(SA_CSV)
    print(f"\n📂 Carregado: {len(df)} municípios de SP")
    print(f"   Colunas-chave: z_veicular, z_biomassa")

    if not {"z_veicular", "z_biomassa", "NM_MUN"}.issubset(df.columns):
        print("❌ Colunas necessárias ausentes do CSV.")
        return

    # 4 cenários
    scenarios = [
        ("permissivo", 0.4, 0.1),
        ("default",    0.5, 0.2),
        ("restritivo", 0.6, 0.3),
        ("assimétrico", 0.7, 0.1),
    ]

    # Aplicar cada cenário
    for name, t_dom, t_rel in scenarios:
        col = f"fonte_{name}"
        df[col] = df.apply(
            lambda r: simplify(classify(r["z_veicular"], r["z_biomassa"],
                                         t_dom, t_rel)),
            axis=1
        )

    # ============================================================
    # Tabela 1: Distribuição por cenário
    # ============================================================
    print("\n" + "=" * 72)
    print("📊 DISTRIBUIÇÃO POR CENÁRIO (% de 645 municípios)")
    print("=" * 72)
    print(f"\n{'Cenário':<14} {'(t_dom; t_rel)':<18} "
          f"{'VEIC':>8} {'BIO':>8} {'MISTA':>8} {'BAIXA':>8}")
    print("-" * 72)

    summary_rows = []
    for name, t_dom, t_rel in scenarios:
        col = f"fonte_{name}"
        counts = df[col].value_counts()
        pct = (counts / len(df) * 100).round(1)
        veic = pct.get("VEICULAR", 0)
        bio = pct.get("BIOMASSA", 0)
        mista = pct.get("MISTA", 0)
        baixa = pct.get("BAIXA_EMISSÃO", 0)
        marker = " ← DEFAULT" if name == "default" else ""
        print(f"{name:<14} ({t_dom:.1f} ; {t_rel:.1f}){'':<5} "
              f"{veic:>7.1f}% {bio:>7.1f}% {mista:>7.1f}% {baixa:>7.1f}%{marker}")
        summary_rows.append({
            "scenario": name,
            "t_dom": t_dom,
            "t_rel": t_rel,
            "n_VEICULAR": int(counts.get("VEICULAR", 0)),
            "n_BIOMASSA": int(counts.get("BIOMASSA", 0)),
            "n_MISTA": int(counts.get("MISTA", 0)),
            "n_BAIXA": int(counts.get("BAIXA_EMISSÃO", 0)),
            "pct_VEICULAR": float(veic),
            "pct_BIOMASSA": float(bio),
            "pct_MISTA": float(mista),
            "pct_BAIXA": float(baixa),
        })

    pd.DataFrame(summary_rows).to_csv(SUM_CSV, index=False)
    print(f"\n  📄 Distribuição salva: {SUM_CSV}")

    # ============================================================
    # Tabela 2: Concordância pareada (% municípios na mesma categoria)
    # ============================================================
    print("\n" + "=" * 72)
    print("📊 CONCORDÂNCIA PAREADA (% de municípios com classificação idêntica)")
    print("=" * 72)
    print(f"\n{'':<14}", end="")
    for name, _, _ in scenarios:
        print(f"{name:>14}", end="")
    print()
    print("-" * 72)

    concordance = {}
    for name1, _, _ in scenarios:
        col1 = f"fonte_{name1}"
        line = []
        for name2, _, _ in scenarios:
            col2 = f"fonte_{name2}"
            if name1 == name2:
                pct = 100.0
            else:
                pct = (df[col1] == df[col2]).sum() / len(df) * 100
            line.append(pct)
            concordance[(name1, name2)] = pct
        print(f"{name1:<14}", end="")
        for v in line:
            print(f"{v:>13.1f}%", end="")
        print()

    # Salvar matriz de concordância
    conc_df = pd.DataFrame(
        index=[s[0] for s in scenarios],
        columns=[s[0] for s in scenarios],
        data=[[concordance[(s1[0], s2[0])] for s2 in scenarios]
              for s1 in scenarios]
    )
    conc_df.to_csv(CONC_CSV)
    print(f"\n  📄 Matriz de concordância: {CONC_CSV}")

    # ============================================================
    # Tabela 3: Municípios sensíveis (mudam de categoria)
    # ============================================================
    print("\n" + "=" * 72)
    print("📊 MUNICÍPIOS SENSÍVEIS (categoria muda em ≥1 cenário vs. default)")
    print("=" * 72)

    df["unique_classifications"] = df[
        [f"fonte_{n}" for n, _, _ in scenarios]
    ].nunique(axis=1)
    sensitive = df[df["unique_classifications"] > 1].copy()
    n_sens = len(sensitive)
    print(f"\n  Total de municípios sensíveis: {n_sens}/645 "
          f"({100*n_sens/645:.1f}%)")

    if n_sens > 0:
        # Top 10 mais sensíveis (mais classificações distintas)
        top = sensitive.nlargest(15, "unique_classifications")
        print(f"\n  Top 15 mais sensíveis:")
        print(f"  {'Município':<28} {'IV':>7} {'IB':>7}", end="")
        for name, _, _ in scenarios:
            print(f"  {name[:9]:>10}", end="")
        print()
        print("  " + "-" * 100)
        for _, row in top.iterrows():
            print(f"  {row['NM_MUN']:<28} {row['idx_veicular']:>7.1f} "
                  f"{row['idx_biomassa']:>7.1f}", end="")
            for name, _, _ in scenarios:
                print(f"  {row[f'fonte_{name}']:>10}", end="")
            print()

        # Salvar tabela completa
        cols_save = ["CD_MUN", "NM_MUN", "idx_veicular", "idx_biomassa",
                     "z_veicular", "z_biomassa"] + \
                    [f"fonte_{n}" for n, _, _ in scenarios]
        cols_save = [c for c in cols_save if c in sensitive.columns]
        sensitive[cols_save].to_csv(CHG_CSV, index=False)
        print(f"\n  📄 Municípios sensíveis: {CHG_CSV}")

    # ============================================================
    # INTERPRETAÇÃO
    # ============================================================
    pct_concord = concordance[("default", "permissivo")]
    pct_concord_restr = concordance[("default", "restritivo")]
    pct_concord_min = min(pct_concord, pct_concord_restr)

    print("\n" + "=" * 72)
    print("🔍 INTERPRETAÇÃO QUANTITATIVA")
    print("=" * 72)
    print(f"  Concordância default vs. permissivo (±0,1): {pct_concord:.1f}%")
    print(f"  Concordância default vs. restritivo (±0,1): {pct_concord_restr:.1f}%")
    print(f"  Mínima entre cenários simétricos: {pct_concord_min:.1f}%")
    print(f"  Municípios sensíveis: {n_sens}/645 ({100*n_sens/645:.1f}%)")

    if pct_concord_min >= 90:
        print("\n  ✅ Concordância ≥90% nos cenários simétricos — ")
        print("     classificação é ROBUSTA aos thresholds escolhidos.")
        print("     Pode-se afirmar no protocolo que a classificação")
        print("     municipal não depende criticamente dos limiares 0,5/0,2.")
    elif pct_concord_min >= 80:
        print("\n  ⚠ Concordância 80–90% — robustez moderada.")
        print("     Considere reportar análise primária com default e")
        print("     análise de sensibilidade em ambos os extremos.")
    else:
        print("\n  ⚠ Concordância <80% — classificação SENSÍVEL aos thresholds.")
        print("     Reportar resultados de TODOS os cenários no Artigo 3.")

    print("\n" + "=" * 72)
    print(f"✅ ANÁLISE CONCLUÍDA")
    print(f"   Resultados em: {OUT_DIR}/sensibility_thresholds_*.csv")
    print("=" * 72 + "\n")


if __name__ == "__main__":
    main()
