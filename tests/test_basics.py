# -*- coding: utf-8 -*-
"""
test_basics.py — Testes mínimos do pipeline PM2.5
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

Cobertura (3 grupos temáticos, ~14 asserções):
  1. Fórmula NASA MERRA-2 PM2.5 = DUS + SS + BC + 1.4·OC + 1.375·SO4
     Garante que coeficientes não foram alterados acidentalmente.
  2. Normalização de CEPs brasileiros (geocoder._normalize_cep)
     Aceita "01310-100", "01310100", "01310 100"; rejeita inválidos.
  3. Fração de ano residencial (cumulative._fraction_of_year)
     Datas precisas, anos bissextos (2020, 2024), sem sobreposição.

Execução:
    pytest tests/                  # roda tudo (verbose por padrão via pytest.ini)
    pytest tests/test_basics.py::test_merra2_formula_factors  # um único teste
"""
from datetime import date
import pytest


# ============================================================
# 1. FÓRMULA NASA MERRA-2
# ============================================================
class TestMERRA2Formula:
    """
    PM2.5 = DUSMASS25 + SSSMASS25 + BCSMASS + 1.4·OCSMASS + 1.375·SO4SMASS
    Coeficientes:
      - 1.4 converte OC em massa de matéria orgânica (Aiken et al. 2008)
      - 1.375 converte SO4 em sulfato de amônio (estequiometria)
    """

    def test_dust_factor_is_one(self):
        import config
        assert config.MERRA2_FORMULA["DUSMASS25"] == 1.0

    def test_seasalt_factor_is_one(self):
        import config
        assert config.MERRA2_FORMULA["SSSMASS25"] == 1.0

    def test_black_carbon_factor_is_one(self):
        import config
        assert config.MERRA2_FORMULA["BCSMASS"] == 1.0

    def test_organic_carbon_factor_is_1_4(self):
        """OC → matéria orgânica: fator 1.4 (Aiken 2008, Pang 2020)."""
        import config
        assert config.MERRA2_FORMULA["OCSMASS"] == 1.4

    def test_sulfate_factor_is_1_375(self):
        """SO4 → sulfato de amônio: fator 1.375 = 132/96."""
        import config
        assert config.MERRA2_FORMULA["SO4SMASS"] == 1.375

    def test_all_five_components_present(self):
        """As 5 espécies da fórmula MERRA-2 devem existir."""
        import config
        required = {"DUSMASS25", "SSSMASS25", "BCSMASS",
                    "OCSMASS", "SO4SMASS"}
        assert required.issubset(set(config.MERRA2_FORMULA.keys()))


# ============================================================
# 2. NORMALIZAÇÃO DE CEPs
# ============================================================
class TestCepNormalization:
    """
    geocoder._normalize_cep deve:
      - aceitar formatos com hífen, sem hífen, com espaços;
      - retornar 8 dígitos puros quando válido;
      - retornar string vazia quando inválido.
    """

    def test_cep_with_dash(self):
        from geocoder import _normalize_cep
        assert _normalize_cep("01310-100") == "01310100"

    def test_cep_already_clean(self):
        from geocoder import _normalize_cep
        assert _normalize_cep("01310100") == "01310100"

    def test_cep_with_spaces(self):
        from geocoder import _normalize_cep
        assert _normalize_cep("01310 100") == "01310100"

    def test_cep_with_letters_strips_them(self):
        """Letras são removidas; resultado precisa ter 8 dígitos para passar."""
        from geocoder import _normalize_cep
        assert _normalize_cep("CEP: 01310-100") == "01310100"

    def test_cep_too_short_returns_empty(self):
        from geocoder import _normalize_cep
        assert _normalize_cep("123") == ""

    def test_cep_too_long_returns_empty(self):
        """Mais de 8 dígitos é inválido."""
        from geocoder import _normalize_cep
        assert _normalize_cep("0131010012") == ""

    def test_cep_empty_returns_empty(self):
        from geocoder import _normalize_cep
        assert _normalize_cep("") == ""
        assert _normalize_cep(None) == ""


# ============================================================
# 3. FRAÇÃO DE ANO RESIDENCIAL (com bissextos)
# ============================================================
class TestFractionOfYear:
    """
    cumulative._fraction_of_year(year, d_start, d_end) deve:
      - retornar 1.0 para ano completo coberto;
      - retornar 0.0 quando o intervalo não cruza o ano;
      - usar 366 dias em anos bissextos (2020, 2024);
      - bater com cálculo manual em ponto-a-ponto.
    """

    def test_full_year_is_1(self):
        from cumulative import _fraction_of_year
        f = _fraction_of_year(2023, date(2020, 1, 1), date(2025, 12, 31))
        assert f == 1.0

    def test_first_half_of_2017(self):
        """01/jan a 30/jun de 2017: 181 dias / 365 ≈ 0.4959"""
        from cumulative import _fraction_of_year
        f = _fraction_of_year(2017, date(2017, 1, 1), date(2017, 6, 30))
        assert abs(f - (181 / 365)) < 1e-9

    def test_second_half_of_2017(self):
        """01/jul a 31/dez de 2017: 184 dias / 365 ≈ 0.5041"""
        from cumulative import _fraction_of_year
        f = _fraction_of_year(2017, date(2017, 7, 1), date(2017, 12, 31))
        assert abs(f - (184 / 365)) < 1e-9

    def test_two_halves_of_2017_sum_to_one(self):
        """Soma das duas metades = exatamente 1.0 (sem dia perdido)."""
        from cumulative import _fraction_of_year
        f1 = _fraction_of_year(2017, date(2017, 1, 1), date(2017, 6, 30))
        f2 = _fraction_of_year(2017, date(2017, 7, 1), date(2017, 12, 31))
        assert abs((f1 + f2) - 1.0) < 1e-9

    def test_leap_year_2020_has_366_days(self):
        """2020 é bissexto: janeiro (31 dias) divide por 366, não 365."""
        from cumulative import _fraction_of_year
        f = _fraction_of_year(2020, date(2020, 1, 1), date(2020, 1, 31))
        assert abs(f - (31 / 366)) < 1e-9

    def test_leap_year_2024_full_is_1(self):
        """2024 é bissexto. Ano completo ainda retorna 1.0."""
        from cumulative import _fraction_of_year
        f = _fraction_of_year(2024, date(2024, 1, 1), date(2024, 12, 31))
        assert f == 1.0

    def test_no_overlap_returns_zero(self):
        from cumulative import _fraction_of_year
        # Residência terminou em 2024, ano consultado é 2025
        f = _fraction_of_year(2025, date(2020, 1, 1), date(2024, 12, 31))
        assert f == 0.0

    def test_mid_year_move_not_january(self):
        """Mudança em 15/maio: cobertura 01/jan a 14/maio = 134 dias/365."""
        from cumulative import _fraction_of_year
        f = _fraction_of_year(2019, date(2019, 1, 1), date(2019, 5, 14))
        assert abs(f - (134 / 365)) < 1e-9
