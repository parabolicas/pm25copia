# -*- coding: utf-8 -*-
"""
cumulative.py — Cálculo da exposição cumulativa ao PM2.5
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

Fórmula (alinhada com o protocolo §4.4):
    Expo_cumulativa (µg/m³-anos) = Σ_residências Σ_anos
        PM2.5(ano, lat, lon) × fração_do_ano_residência

A "fração_do_ano_residência" é calculada com precisão diária (anos bissextos
tratados corretamente), permitindo modelar mudanças residenciais em qualquer
mês — não apenas em janeiro.

Estruturas aceitas para uma residência (em config.TEST_PATIENTS ou em
banco de dados externo):

    a) date_start / date_end — formato "AAAA-MM-DD" (preferido):
       {"cep": "...", "lat": ..., "lon": ...,
        "date_start": "2015-05-15", "date_end": "2024-12-31"}

    b) year_start / year_end — anos inteiros (compatibilidade retroativa):
       {"cep": "...", "lat": ..., "lon": ...,
        "year_start": 2015, "year_end": 2024}
       Equivale a date_start="2015-01-01", date_end="2024-12-31".

Importante: assumimos residências SEQUENCIAIS (sem sobreposição temporal).
Caso uma coorte real tenha pacientes com endereços simultâneos, a soma das
frações em algum ano pode exceder 1,0 — o relatório de validação alerta.

Fonte primária do PM2.5: superfície bias-corrected via exposure.estimate_pm25
(que internamente usa pm25_surface.lookup_pm25_buffer).
"""
from datetime import date, datetime
from typing import Tuple, Union

import pandas as pd
from exposure import estimate_pm25


# ============================================================
# UTILITÁRIOS DE DATA
# ============================================================
def _parse_date(value: Union[str, int, date]) -> date:
    """
    Converte string/int para datetime.date.
    Aceita: 'AAAA-MM-DD', 'AAAA-MM' (→dia 01), 'AAAA' (→01/01), int ano.
    """
    if isinstance(value, date):
        return value
    if isinstance(value, int):
        return date(int(value), 1, 1)
    if isinstance(value, str):
        s = value.strip()
        for fmt, normalizer in [
            ("%Y-%m-%d", lambda x: x),
            ("%Y-%m", lambda x: x + "-01"),
            ("%Y", lambda x: x + "-01-01"),
        ]:
            try:
                return datetime.strptime(normalizer(s), "%Y-%m-%d").date()
            except ValueError:
                continue
    raise ValueError(f"Não foi possível interpretar data: {value!r}")


def _residence_dates(residence: dict) -> Tuple[date, date]:
    """
    Retorna (date_start, date_end) para uma residência.
    Aceita date_start/date_end (preferido) ou year_start/year_end (compat).
    """
    if "date_start" in residence and "date_end" in residence:
        d_start = _parse_date(residence["date_start"])
        d_end = _parse_date(residence["date_end"])
    elif "year_start" in residence and "year_end" in residence:
        # Compatibilidade retroativa: ano inteiro = de 01/jan a 31/dez
        d_start = date(int(residence["year_start"]), 1, 1)
        d_end = date(int(residence["year_end"]), 12, 31)
    else:
        raise KeyError(
            "Residência precisa ter 'date_start'/'date_end' "
            "ou 'year_start'/'year_end'."
        )
    if d_start > d_end:
        raise ValueError(
            f"date_start ({d_start}) é posterior a date_end ({d_end}) — "
            f"residência inválida."
        )
    return d_start, d_end


def _fraction_of_year(year: int, d_start: date, d_end: date) -> float:
    """
    Fração do ano calendário 'year' coberta pelo intervalo [d_start, d_end].

    Trata anos bissextos: o denominador é 366 quando year é bissexto, 365 caso
    contrário. O numerador é o número de dias da interseção entre o ano e o
    intervalo da residência (inclusivo nas duas pontas).

    Returns: float entre 0.0 e 1.0.
    """
    yr_start = date(year, 1, 1)
    yr_end = date(year, 12, 31)

    overlap_start = max(d_start, yr_start)
    overlap_end = min(d_end, yr_end)

    if overlap_start > overlap_end:
        return 0.0

    days_in_year = (yr_end - yr_start).days + 1  # 365 ou 366
    days_overlap = (overlap_end - overlap_start).days + 1
    return days_overlap / days_in_year


def residence_years(residence: dict) -> range:
    """
    Iterador dos anos calendário cobertos pela residência (inclusive
    anos onde a fração é < 1.0). Útil para o run_pipeline determinar
    quais anos de dados ambientais carregar.
    """
    d_start, d_end = _residence_dates(residence)
    return range(d_start.year, d_end.year + 1)


# ============================================================
# CÁLCULO PRINCIPAL
# ============================================================
def calc_cumulative_exposure(patient: dict, buffer_km: float, cache) -> dict:
    """
    Calcula a exposição cumulativa de um paciente para um dado raio de buffer,
    com ponderação por fração do ano em cada residência.
    """
    yearly_details = []
    source_counts = {"SURFACE": 0, "N/A": 0}
    total_exposure = 0.0          # µg/m³ × anos (somatório das contribuições)
    total_year_fraction = 0.0     # soma das frações de anos cobertos
    valid_year_fraction = 0.0     # soma das frações com pm25 válido
    total_fire_foci = 0
    years_missing = []

    for residence in patient["residences"]:
        geo = residence.get("geocoded", residence)
        lat = geo.get("lat", residence.get("lat"))
        lon = geo.get("lon", residence.get("lon"))
        d_start, d_end = _residence_dates(residence)
        city = residence.get("city", "N/A")
        cep = residence.get("cep", "N/A")

        for year in range(d_start.year, d_end.year + 1):
            fraction = _fraction_of_year(year, d_start, d_end)
            if fraction <= 0.0:
                continue

            result = estimate_pm25(lat, lon, year, buffer_km, cache)
            pm25_val = result["pm25_value"]

            if pm25_val is not None:
                exposure_year = pm25_val * fraction  # µg/m³ × fração-de-ano
                total_exposure += exposure_year
                valid_year_fraction += fraction
                source = result["primary_source"] or "SURFACE"
            else:
                exposure_year = None
                source = "N/A"
                years_missing.append((cep, year))

            source_counts[source] = source_counts.get(source, 0) + 1
            total_year_fraction += fraction

            fire_info = result.get("fires", {}) or {}
            n_foci = fire_info.get("n_foci", 0)
            total_fire_foci += n_foci

            surface_info = result.get("surface") or {}
            cetesb_info = result.get("cetesb") or {}
            cams_info = result.get("cams") or {}
            merra2_info = result.get("merra2") or {}

            yearly_details.append({
                "patient_id": patient["id"],
                "year": year,
                "cep": cep,
                "city": city,
                "lat": lat,
                "lon": lon,
                "buffer_km": buffer_km,
                "fraction_of_year": round(fraction, 4),
                # Fonte primária
                "pm25_ugm3": pm25_val,
                "source": source,
                "exposure_ugm3_year": (round(exposure_year, 4)
                                       if exposure_year is not None else None),
                # Detalhes da superfície
                "surface_pm25": surface_info.get("pm25_value"),
                "surface_n_pixels": surface_info.get("n_pixels"),
                "surface_method": surface_info.get("method"),
                # Focos de calor
                "fire_foci": n_foci,
                "fire_density_km2": fire_info.get("density_per_km2", 0),
                # Diagnóstico
                "cetesb_pm25": cetesb_info.get("pm25_value"),
                "cetesb_n_stations": cetesb_info.get("n_stations", 0),
                "cams_pm25_raw": cams_info.get("pm25_value"),
                "merra2_pm25_raw": merra2_info.get("pm25_value"),
            })

    mean_pm25 = (round(total_exposure / valid_year_fraction, 2)
                 if valid_year_fraction > 0 else None)

    return {
        "patient_id": patient["id"],
        "patient_name": patient["name"],
        "buffer_km": buffer_km,
        "cumulative_exposure": round(total_exposure, 2),
        # Agora total_years e valid_years são FLOATS (soma de frações)
        "total_years": round(total_year_fraction, 4),
        "valid_years": round(valid_year_fraction, 4),
        "mean_annual_pm25": mean_pm25,
        "total_fire_foci": total_fire_foci,
        "source_summary": source_counts,
        "years_missing": years_missing,
        "yearly_details": yearly_details,
    }


def run_all_patients(patients: list, buffer_km_list: list, cache) -> pd.DataFrame:
    """Executa o pipeline para todos os pacientes × buffers."""
    all_results = []
    summary_rows = []

    for patient in patients:
        print(f"\n{'='*60}")
        print(f"🔬 {patient['id']}: {patient['name']}")
        print(f"{'='*60}")

        # Mostrar período residencial em formato legível
        for res in patient["residences"]:
            try:
                d_start, d_end = _residence_dates(res)
                print(f"   📅 {res.get('city', '?')}: "
                      f"{d_start.isoformat()} → {d_end.isoformat()} "
                      f"({(d_end - d_start).days + 1} dias)")
            except Exception as e:
                print(f"   ⚠ Erro ao parsear datas de residência: {e}")

        for buffer_km in buffer_km_list:
            print(f"\n  📐 Buffer: {buffer_km} km")
            result = calc_cumulative_exposure(patient, buffer_km, cache)

            mean_str = (f"{result['mean_annual_pm25']:.1f}"
                        if result['mean_annual_pm25'] is not None else "N/D")
            print(f"     Exposição cumulativa: {result['cumulative_exposure']:.1f} µg/m³-anos")
            print(f"     PM2.5 médio anual:    {mean_str} µg/m³")
            print(f"     Anos válidos / total: {result['valid_years']:.2f} / "
                  f"{result['total_years']:.2f}")
            print(f"     Focos de calor total: {result['total_fire_foci']}")
            print(f"     Fontes: {result['source_summary']}")
            if result['years_missing']:
                print(f"     ⚠ Anos sem superfície: {result['years_missing']}")

            all_results.extend(result["yearly_details"])

            n_total_records = max(
                result["source_summary"].get("SURFACE", 0)
                + result["source_summary"].get("N/A", 0),
                1
            )
            summary_rows.append({
                "patient_id": result["patient_id"],
                "patient_name": result["patient_name"],
                "buffer_km": buffer_km,
                "cumulative_exposure_ugm3_years": result["cumulative_exposure"],
                "mean_annual_pm25_ugm3": result["mean_annual_pm25"],
                "total_years": result["total_years"],
                "valid_years": result["valid_years"],
                "total_fire_foci": result["total_fire_foci"],
                "pct_surface": (result["source_summary"].get("SURFACE", 0)
                                / n_total_records * 100),
                "pct_missing": (result["source_summary"].get("N/A", 0)
                                / n_total_records * 100),
            })

    details_df = pd.DataFrame(all_results)
    summary_df = pd.DataFrame(summary_rows)

    return details_df, summary_df
