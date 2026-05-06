# -*- coding: utf-8 -*-
"""
exposure.py — Estimativa espacial de PM2.5 por múltiplas fontes
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

ARQUITETURA (versão corrigida — alinhada com o protocolo v2.2 §4.4):
  • Fonte PRIMÁRIA: superfície anual bias-corrected (MERRA-2 + CETESB IDW
    de resíduos) gerada por pm25_surface.py — consulta via lookup_pm25_buffer.
  • Fontes auxiliares (apenas DIAGNÓSTICO): CETESB IDW puro, CAMS bruto e
    MERRA-2 bruto. Esses valores são SEMPRE computados e armazenados, mas
    NUNCA entram no pm25_value final. Servem para análise de sensibilidade
    e validação cruzada (cross_validation.py).

Mudança em relação à versão anterior:
  Antes: hierarquia ad-hoc CETESB→CAMS→MERRA-2 (CAMS bruto tem bias
  documentado de +16 µg/m³ → produzia valores irrealistas em pacientes
  fora do raio de qualquer estação CETESB).
  Agora: fonte única e calibrada — superfície bias-corrected.
"""
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from scipy.spatial.distance import cdist


def _haversine_km(lat1, lon1, lat2, lon2):
    """Distância haversine entre dois pontos em km."""
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return R * 2 * np.arcsin(np.sqrt(a))


# ============================================================
# FONTES AUXILIARES (DIAGNÓSTICO) — não entram no pm25_value
# ============================================================
def cetesb_idw(lat: float, lon: float, year: int, buffer_km: float,
               cetesb_gdf: gpd.GeoDataFrame, power: float = 2.0) -> dict:
    """
    Estimativa PM2.5 por IDW usando estações CETESB DENTRO do buffer.
    Mantida apenas para fins de diagnóstico/comparação.
    """
    year_data = cetesb_gdf[cetesb_gdf["ano"] == year].copy()
    if len(year_data) == 0:
        return None

    distances = year_data.apply(
        lambda row: _haversine_km(lat, lon, row["lat"], row["lon"]), axis=1
    )
    year_data = year_data.copy()
    year_data["dist_km"] = distances

    within_buffer = year_data[year_data["dist_km"] <= buffer_km].copy()

    if len(within_buffer) < 1:
        return None

    if within_buffer["dist_km"].min() < 0.1:
        closest = within_buffer.loc[within_buffer["dist_km"].idxmin()]
        return {
            "pm25_value": float(closest["media_anual_pm25_ugm3"]),
            "n_stations": 1,
            "stations_used": [closest["estacao_nome"]],
            "distances_km": [float(closest["dist_km"])],
            "method": "exact_station",
        }

    weights = 1.0 / (within_buffer["dist_km"].values ** power)
    weights_norm = weights / weights.sum()
    pm25_idw = np.sum(weights_norm * within_buffer["media_anual_pm25_ugm3"].values)

    return {
        "pm25_value": float(round(pm25_idw, 2)),
        "n_stations": len(within_buffer),
        "stations_used": within_buffer["estacao_nome"].tolist(),
        "distances_km": within_buffer["dist_km"].round(1).tolist(),
        "method": "idw",
    }


def cams_buffer_mean(lat: float, lon: float, year: int, buffer_km: float,
                     cams_data) -> dict:
    """
    Média PM2.5 CAMS bruto nos pixels dentro do buffer.

    ATENÇÃO: CAMS apresenta bias sistemático de +16 µg/m³ em SP (validação
    cruzada interna, n=282). Este valor é mantido apenas como referência de
    DIAGNÓSTICO e nunca deve ser usado como exposição individual.
    """
    if cams_data is None:
        return None

    try:
        data_2d = cams_data.squeeze(drop=True)

        lats = data_2d.coords["latitude"].values
        lons = data_2d.coords["longitude"].values

        buffer_deg = buffer_km / 111.0

        lat_mask = (lats >= lat - buffer_deg) & (lats <= lat + buffer_deg)
        lon_mask = (lons >= lon - buffer_deg) & (lons <= lon + buffer_deg)

        sub_lats = lats[lat_mask]
        sub_lons = lons[lon_mask]

        if len(sub_lats) == 0 or len(sub_lons) == 0:
            nearest = data_2d.sel(latitude=lat, longitude=lon, method="nearest")
            val = float(nearest.values)
            return {
                "pm25_value": round(val, 2),
                "n_pixels": 1,
                "method": "nearest_pixel",
            }

        subset = data_2d.sel(latitude=sub_lats, longitude=sub_lons)
        subset_vals = subset.values

        grid_lats, grid_lons = np.meshgrid(sub_lats, sub_lons, indexing="ij")
        distances = np.vectorize(_haversine_km)(lat, lon, grid_lats, grid_lons)
        mask = distances <= buffer_km

        if not mask.any():
            nearest = data_2d.sel(latitude=lat, longitude=lon, method="nearest")
            val = float(nearest.values)
            return {
                "pm25_value": round(val, 2),
                "n_pixels": 1,
                "method": "nearest_pixel",
            }

        values = subset_vals[mask]
        valid = values[~np.isnan(values)]
        if len(valid) == 0:
            return None

        return {
            "pm25_value": round(float(np.mean(valid)), 2),
            "n_pixels": int(len(valid)),
            "method": "buffer_mean",
        }
    except Exception:
        try:
            data_2d = cams_data.squeeze(drop=True)
            nearest = data_2d.sel(latitude=lat, longitude=lon, method="nearest")
            val = float(nearest.values)
            return {
                "pm25_value": round(val, 2),
                "n_pixels": 1,
                "method": "nearest_pixel_fallback",
            }
        except Exception:
            return None


def merra2_buffer_mean(lat: float, lon: float, year: int, buffer_km: float,
                       merra2_data) -> dict:
    """
    Média PM2.5 MERRA-2 BRUTO nos pixels dentro do buffer.

    ATENÇÃO: MERRA-2 bruto subestima PM2.5 em ~3 µg/m³ em SP. Mantido como
    DIAGNÓSTICO. A versão corrigida (com bias-correction CETESB) está
    disponível via pm25_surface.lookup_pm25_buffer e é a fonte primária.
    """
    if merra2_data is None:
        return None

    try:
        data_2d = merra2_data.squeeze(drop=True)

        lat_name = "lat" if "lat" in data_2d.dims else "latitude"
        lon_name = "lon" if "lon" in data_2d.dims else "longitude"

        lats = data_2d.coords[lat_name].values
        lons = data_2d.coords[lon_name].values

        buffer_deg = buffer_km / 111.0

        lat_mask = (lats >= lat - buffer_deg) & (lats <= lat + buffer_deg)
        lon_mask = (lons >= lon - buffer_deg) & (lons <= lon + buffer_deg)

        sub_lats = lats[lat_mask]
        sub_lons = lons[lon_mask]

        if len(sub_lats) == 0 or len(sub_lons) == 0:
            nearest = data_2d.sel(**{lat_name: lat, lon_name: lon}, method="nearest")
            val = float(nearest.values)
            return {
                "pm25_value": round(val, 2),
                "n_pixels": 1,
                "method": "nearest_pixel",
            }

        subset = data_2d.sel(**{lat_name: sub_lats, lon_name: sub_lons})
        subset_vals = subset.values

        grid_lats, grid_lons = np.meshgrid(sub_lats, sub_lons, indexing="ij")
        distances = np.vectorize(_haversine_km)(lat, lon, grid_lats, grid_lons)
        mask = distances <= buffer_km

        if not mask.any():
            nearest = data_2d.sel(**{lat_name: lat, lon_name: lon}, method="nearest")
            val = float(nearest.values)
            return {
                "pm25_value": round(val, 2),
                "n_pixels": 1,
                "method": "nearest_pixel",
            }

        values = subset_vals[mask]
        valid = values[~np.isnan(values)]
        if len(valid) == 0:
            return None

        return {
            "pm25_value": round(float(np.mean(valid)), 2),
            "n_pixels": int(len(valid)),
            "method": "buffer_mean",
        }
    except Exception:
        try:
            data_2d = merra2_data.squeeze(drop=True)
            lat_name = "lat" if "lat" in data_2d.dims else "latitude"
            lon_name = "lon" if "lon" in data_2d.dims else "longitude"
            nearest = data_2d.sel(**{lat_name: lat, lon_name: lon}, method="nearest")
            val = float(nearest.values)
            return {
                "pm25_value": round(val, 2),
                "n_pixels": 1,
                "method": "nearest_pixel_fallback",
            }
        except Exception:
            return None


def fire_density(lat: float, lon: float, year: int, buffer_km: float,
                 fire_gdf: gpd.GeoDataFrame) -> dict:
    """Conta focos de calor dentro do buffer para um ano."""
    if fire_gdf is None or len(fire_gdf) == 0:
        return {"n_foci": 0, "density_per_km2": 0.0}

    distances = fire_gdf.geometry.apply(
        lambda pt: _haversine_km(lat, lon, pt.y, pt.x)
    )

    n_within = int((distances <= buffer_km).sum())
    area_km2 = np.pi * buffer_km ** 2

    return {
        "n_foci": n_within,
        "density_per_km2": round(n_within / area_km2, 4),
    }


# ============================================================
# FONTE PRIMÁRIA: SUPERFÍCIE BIAS-CORRECTED
# ============================================================
def surface_buffer_mean(lat: float, lon: float, year: int, buffer_km: float) -> dict:
    """
    Estimativa PM2.5 a partir da superfície anual bias-corrected
    (MERRA-2 + CETESB IDW de resíduos).

    Esta é a FONTE PRIMÁRIA para exposição individual.

    Returns:
        dict {pm25_value, n_pixels, method} ou None.
    """
    # Import tardio para evitar ciclo (pm25_surface importa data_loaders)
    try:
        from pm25_surface import lookup_pm25_buffer
    except ImportError:
        return None

    return lookup_pm25_buffer(lat, lon, year, buffer_km)


# ============================================================
# FUNÇÃO PRINCIPAL — usada por cumulative.py
# ============================================================
def estimate_pm25(lat: float, lon: float, year: int, buffer_km: float,
                  cache) -> dict:
    """
    Estimativa PM2.5 individual.

    NOVA HIERARQUIA (alinhada com o protocolo v2.2 §4.4):
      • PRIMÁRIO: superfície bias-corrected (lookup_pm25_buffer).
      • DIAGNÓSTICO (sempre calculados, nunca usados como pm25_value):
          - CETESB IDW puro (mostra contribuição direta das estações)
          - CAMS bruto (referência satelital independente, com bias +16)
          - MERRA-2 bruto (referência satelital independente, com bias -3)

    Se a superfície bias-corrected não estiver disponível para o ano, o
    pm25_value retornará None — e o registro deve ser tratado como missing
    pela exposição cumulativa, em vez de cair silenciosamente para CAMS.

    Args:
        cache: DataCache instance

    Returns:
        dict com:
            pm25_value (float|None), primary_source ("SURFACE"|None),
            cetesb (dict|None), cams (dict|None), merra2 (dict|None),
            surface (dict|None), fires (dict)
    """
    result = {
        "lat": lat,
        "lon": lon,
        "year": year,
        "buffer_km": buffer_km,
        "pm25_value": None,
        "primary_source": None,
        "surface": None,   # NOVO: fonte primária
        "cetesb": None,    # diagnóstico
        "cams": None,      # diagnóstico
        "merra2": None,    # diagnóstico
        "fires": None,
    }

    # 1. SUPERFÍCIE BIAS-CORRECTED (fonte primária — única usada para pm25_value)
    surface_result = surface_buffer_mean(lat, lon, year, buffer_km)
    result["surface"] = surface_result

    if surface_result is not None and surface_result.get("pm25_value") is not None:
        result["pm25_value"] = surface_result["pm25_value"]
        result["primary_source"] = "SURFACE"

    # 2. CETESB IDW puro (diagnóstico)
    cetesb_gdf = cache.get_cetesb()
    result["cetesb"] = cetesb_idw(lat, lon, year, buffer_km, cetesb_gdf)

    # 3. CAMS bruto (diagnóstico)
    cams_data = cache.get_cams(year)
    result["cams"] = cams_buffer_mean(lat, lon, year, buffer_km, cams_data)

    # 4. MERRA-2 bruto (diagnóstico)
    merra2_data = cache.get_merra2(year)
    result["merra2"] = merra2_buffer_mean(lat, lon, year, buffer_km, merra2_data)

    # 5. Focos de calor (já era informação independente do PM2.5)
    fire_gdf = cache.get_fires(year)
    result["fires"] = fire_density(lat, lon, year, buffer_km, fire_gdf)

    return result
