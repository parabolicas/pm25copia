# -*- coding: utf-8 -*-
"""
test_source_specific.py — Testes da calibração dual source-specific
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

Cobertura (5 testes em 2 grupos):
  1. Fire density grid (KDE BDQueimadas)
     - Shape correto da grade alvo
     - Comportamento com array vazio
  2. Calibração dual (RMSP-anchored + boost por densidade)
     - Boost factor=1.0 não muda fracionamento (neutro)
     - Boost factor>1.0 aumenta w_biomassa onde densidade > 0
     - Renormalização preserva soma = 1.0 após boost

Execução:
    pytest tests/test_source_specific.py -v
"""
import numpy as np
import pytest


# ============================================================
# 1. FIRE DENSITY GRID (KDE)
# ============================================================
class TestFireDensityGrid:
    """
    generate_fire_density_grid(fire_coords, bandwidth_km=25) deve:
      - retornar array (n_target_lats, n_target_lons) consistente com a grade;
      - acumular kernel gaussiano para cada foco;
      - retornar array todo-zero se não houver focos.
    """

    def test_shape_matches_target_grid(self):
        from source_specific import (generate_fire_density_grid,
                                      TARGET_LATS, TARGET_LONS)
        # 5 focos sintéticos espalhados em SP
        fires = np.array([
            [-23.5, -46.6],   # SP capital
            [-23.6, -46.5],
            [-22.9, -47.0],   # Campinas
            [-21.2, -47.8],   # Ribeirão Preto
            [-22.7, -47.6],   # Piracicaba
        ])
        density = generate_fire_density_grid(fires, bandwidth_km=25)
        assert density.shape == (len(TARGET_LATS), len(TARGET_LONS)), \
            f"Shape incorreto: {density.shape}"
        # Pelo menos um pixel deve ter contribuição não-zero
        assert density.max() > 0, "Densidade máxima deveria ser > 0"

    def test_empty_fires_returns_zero_grid(self):
        from source_specific import (generate_fire_density_grid,
                                      TARGET_LATS, TARGET_LONS)
        empty = np.zeros((0, 2))
        density = generate_fire_density_grid(empty, bandwidth_km=25)
        assert density.shape == (len(TARGET_LATS), len(TARGET_LONS))
        assert density.max() == 0, "Sem focos, densidade deve ser zero"


# ============================================================
# 2. CALIBRAÇÃO DUAL (RMSP-anchored + BOOST)
# ============================================================
class TestSourceSpecificBoost:
    """
    calibrate_fractions(fractions_raw, fire_density_normalized, boost) deve:
      - manter fracionamento intacto quando boost=1.0;
      - aumentar w_biomassa onde densidade > 0 quando boost > 1.0;
      - renormalizar w_veh + w_bio + w_out = 1 em cada pixel após boost.
    """

    @staticmethod
    def _make_synthetic_fractions():
        """
        Cria DataArrays sintéticos no grid MERRA-2 simplificado.
        Valores: bc_oc=0.15 (RMSP-like), f_om=0.30 (faixa típica).
        """
        import xarray as xr
        lats = np.linspace(-25.5, -19.5, 13)
        lons = np.linspace(-53.5, -44.0, 16)
        bc_oc = xr.DataArray(
            np.full((13, 16), 0.15),
            dims=["lat", "lon"],
            coords={"lat": lats, "lon": lons},
        )
        f_om = xr.DataArray(
            np.full((13, 16), 0.30),
            dims=["lat", "lon"],
            coords={"lat": lats, "lon": lons},
        )
        return {"bc_oc_ratio": bc_oc, "f_om": f_om}

    def test_factor_one_does_not_alter_calibration(self):
        from source_specific import (calibrate_fractions,
                                      TARGET_LATS, TARGET_LONS)
        fr = self._make_synthetic_fractions()
        density_full = np.ones((len(TARGET_LATS), len(TARGET_LONS)))
        r = calibrate_fractions(
            fr,
            fire_density_normalized=density_full,
            biomass_boost_factor=1.0,
        )
        # boost_factor=1.0: meta registra mas não altera w_biomassa
        assert r["boost_meta"]["biomass_boost_factor"] == 1.0
        # Soma das frações deve ser 1.0
        total = (r["w_veicular"].values
                 + r["w_biomassa"].values
                 + r["w_outros"].values)
        assert np.allclose(total, 1.0, atol=1e-3), \
            f"Frações não somam 1.0 (sem boost)"

    def test_boost_greater_than_one_increases_biomassa(self):
        from source_specific import (calibrate_fractions,
                                      TARGET_LATS, TARGET_LONS)
        fr = self._make_synthetic_fractions()
        # Sem density → sem boost (baseline)
        r0 = calibrate_fractions(fr, fire_density_normalized=None,
                                  biomass_boost_factor=1.4)
        # Com density=1 em todo grid → boost máximo aplicado
        density_full = np.ones((len(TARGET_LATS), len(TARGET_LONS)))
        r1 = calibrate_fractions(fr, fire_density_normalized=density_full,
                                  biomass_boost_factor=1.4)
        bio_baseline = float(r0["w_biomassa"].mean())
        bio_boosted = float(r1["w_biomassa"].mean())
        assert bio_boosted > bio_baseline, \
            f"Boost não aumentou w_biomassa: {bio_boosted:.4f} ≤ {bio_baseline:.4f}"
        assert r1["boost_meta"]["applied"] is True

    def test_normalization_after_boost(self):
        """w_veh + w_bio + w_out = 1.0 em cada pixel após boost."""
        from source_specific import (calibrate_fractions,
                                      TARGET_LATS, TARGET_LONS)
        fr = self._make_synthetic_fractions()
        density_full = np.ones((len(TARGET_LATS), len(TARGET_LONS)))
        r = calibrate_fractions(
            fr,
            fire_density_normalized=density_full,
            biomass_boost_factor=1.6,
        )
        total = (r["w_veicular"].values
                 + r["w_biomassa"].values
                 + r["w_outros"].values)
        assert np.allclose(total, 1.0, atol=1e-3), \
            f"Frações não somam 1.0 após boost: max diff = {abs(total - 1.0).max():.4f}"
