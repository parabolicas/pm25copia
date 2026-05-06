# -*- coding: utf-8 -*-
"""
geocoder.py — Geocodificação de CEPs brasileiros para coordenadas
Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP

ARQUITETURA (versão corrigida — ViaCEP + Nominatim + cache):

  Estratégia em 5 etapas hierárquicas:
    1. CACHE local (.geocoder_cache.json) — evita repetir queries.
    2. ViaCEP (api.viacep.com.br) — retorna endereço estruturado
       (logradouro, bairro, localidade/cidade, UF) com cobertura próxima
       de 100% dos CEPs brasileiros ativos.
    3. Nominatim (OpenStreetMap) — geocoda o endereço estruturado em
       três tentativas, do mais preciso ao mais geral:
         a) logradouro + bairro + cidade + UF
         b) logradouro + cidade + UF
         c) cidade + UF
    4. VALIDAÇÃO de cada tentativa: dentro do estado de SP, e (se houver
       fallback) distância haversine ≤ 50 km do fallback fornecido.
    5. FALLBACK: usa coordenadas hardcoded se todas as tentativas
       falharem. Resultado marcado como validated=False (não cacheado).

  Vantagens vs. versão anterior (Nominatim direto com CEP):
    - ViaCEP é nacional, oficial, gratuito, sem rate limit pesado.
    - Endereço estruturado dá precisão de logradouro, não só município.
    - Cache persistente acelera re-execuções em coortes grandes.
    - Resultado guarda também o endereço ViaCEP — facilita auditoria.

  Performance: na PRIMEIRA execução, ~3-4 segundos por CEP (rate limit
  Nominatim de 1 req/s, com até 3 tentativas). Em re-execuções, leitura
  imediata do cache (<1 ms por CEP).

  Dependências:
    - requests (geralmente já instalado): pip install requests
    - geopy:                              pip install geopy
"""
import json
import os
import re
import time
import math
from typing import Optional

import requests
from geopy.geocoders import Nominatim
from geopy.exc import GeocoderTimedOut, GeocoderServiceError


# ============================================================
# CONFIGURAÇÃO
# ============================================================
# Limites do Estado de São Paulo (validação de plausibilidade)
SP_LAT_MIN, SP_LAT_MAX = -26.0, -19.0
SP_LON_MIN, SP_LON_MAX = -54.0, -44.0

# Distância máxima aceitável entre o resultado online e o fallback
# (impede aceitar geocodificações espúrias do Nominatim)
MAX_DEVIATION_KM = 50.0

# Cache persistente (na pasta do projeto)
CACHE_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    ".geocoder_cache.json"
)

# Timeouts
VIACEP_TIMEOUT = 5
NOMINATIM_TIMEOUT = 10

# User-Agent obrigatório para Nominatim
USER_AGENT = "pm25_cpnpc_research_fmusp"

# Pausa entre queries Nominatim (rate limit oficial: 1 req/s)
NOMINATIM_SLEEP = 1.0


# ============================================================
# CACHE
# ============================================================
def _load_cache() -> dict:
    if not os.path.exists(CACHE_PATH):
        return {}
    try:
        with open(CACHE_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def _save_cache(cache: dict) -> None:
    try:
        with open(CACHE_PATH, "w", encoding="utf-8") as f:
            json.dump(cache, f, indent=2, ensure_ascii=False)
    except Exception as e:
        print(f"  ⚠ Não foi possível salvar cache: {e}")


# Cache em memória (carregado uma única vez por execução)
_cache: dict = _load_cache()


# ============================================================
# UTILITIES
# ============================================================
def _normalize_cep(cep: str) -> str:
    """Remove hífen/espaços; retorna 8 dígitos puros, ou '' se inválido."""
    if not cep:
        return ""
    digits = re.sub(r"\D", "", str(cep))
    return digits if len(digits) == 8 else ""


def _is_in_sp(lat: float, lon: float) -> bool:
    return SP_LAT_MIN <= lat <= SP_LAT_MAX and SP_LON_MIN <= lon <= SP_LON_MAX


def _haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) ** 2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


# ============================================================
# VIACEP
# ============================================================
def query_viacep(cep_digits: str) -> Optional[dict]:
    """
    Consulta a API ViaCEP. Retorna dict com endereço estruturado
    {cep, logradouro, bairro, localidade, uf} ou None.
    """
    if not cep_digits or len(cep_digits) != 8:
        return None
    try:
        url = f"https://viacep.com.br/ws/{cep_digits}/json/"
        r = requests.get(url, timeout=VIACEP_TIMEOUT)
        if r.status_code != 200:
            return None
        data = r.json()
        if data.get("erro"):
            return None
        return {
            "cep": data.get("cep", "") or "",
            "logradouro": (data.get("logradouro") or "").strip(),
            "bairro": (data.get("bairro") or "").strip(),
            "localidade": (data.get("localidade") or "").strip(),
            "uf": (data.get("uf") or "").strip(),
        }
    except requests.exceptions.RequestException as e:
        print(f"  ⚠ ViaCEP falhou para {cep_digits}: {e}")
        return None
    except Exception as e:
        print(f"  ⚠ ViaCEP erro inesperado para {cep_digits}: {e}")
        return None


# ============================================================
# NOMINATIM
# ============================================================
_geolocator = None


def _get_geolocator():
    global _geolocator
    if _geolocator is None:
        _geolocator = Nominatim(user_agent=USER_AGENT)
    return _geolocator


def query_nominatim(query: str) -> Optional[dict]:
    """
    Geocoda uma string no Nominatim. Retorna {lat, lon, address} ou None.
    Respeita rate limit oficial via NOMINATIM_SLEEP antes de cada chamada.
    """
    try:
        time.sleep(NOMINATIM_SLEEP)
        geocoder = _get_geolocator()
        location = geocoder.geocode(
            query, timeout=NOMINATIM_TIMEOUT, country_codes="br"
        )
        if location is None:
            return None
        return {
            "lat": float(location.latitude),
            "lon": float(location.longitude),
            "address": location.address,
        }
    except (GeocoderTimedOut, GeocoderServiceError) as e:
        print(f"  ⚠ Nominatim timeout/serviço: {e}")
        return None
    except Exception as e:
        print(f"  ⚠ Nominatim erro inesperado: {e}")
        return None


# ============================================================
# GEOCODE PRINCIPAL
# ============================================================
def geocode_cep(cep: str,
                fallback_lat: Optional[float] = None,
                fallback_lon: Optional[float] = None,
                city: Optional[str] = None) -> dict:
    """
    Geocodifica um CEP brasileiro com estratégia hierárquica.

    Args:
        cep: CEP em qualquer formato ("01310-100", "01310100", etc.)
        fallback_lat, fallback_lon: coordenadas de fallback opcionais.
            Se fornecidas, são usadas para validar resultados online
            (descartados se distância > MAX_DEVIATION_KM) e como último
            recurso se todas as tentativas online falharem.
        city: nome da cidade (usado se ViaCEP falhar ou não retornar
            localidade — ainda permite tentativa de Nominatim por cidade).

    Returns:
        dict com:
          - lat (float), lon (float)
          - source (str): "cache:..." | "nominatim_full" |
                          "nominatim_street_city" | "nominatim_city" |
                          "fallback"
          - address (str): endereço retornado pelo Nominatim, ou descrição
          - validated (bool): True se passou em todas as checagens
          - viacep (dict|None): endereço estruturado retornado pelo ViaCEP

    Raises:
        ValueError se não houver fallback e nenhuma tentativa online
        retornou resultado válido em SP.
    """
    cep_digits = _normalize_cep(cep)

    # 1. CACHE
    if cep_digits and cep_digits in _cache:
        cached = dict(_cache[cep_digits])
        original_source = cached.get("source", "unknown")
        cached["source"] = f"cache:{original_source}"
        return cached

    # 2. VIACEP — endereço estruturado
    viacep = query_viacep(cep_digits) if cep_digits else None

    # Determinar campos efetivos (preferindo ViaCEP, com fallback ao parâmetro city)
    cidade = ((viacep["localidade"] if viacep else None) or (city or "")).strip()
    uf = ((viacep["uf"] if viacep else None) or "SP").strip()
    logradouro = ((viacep["logradouro"] if viacep else None) or "").strip()
    bairro = ((viacep["bairro"] if viacep else None) or "").strip()

    # 3. NOMINATIM — tentativas em ordem decrescente de precisão
    queries: list = []
    if logradouro and bairro and cidade:
        queries.append(("nominatim_full",
                        f"{logradouro}, {bairro}, {cidade}, {uf}, Brasil"))
    if logradouro and cidade:
        queries.append(("nominatim_street_city",
                        f"{logradouro}, {cidade}, {uf}, Brasil"))
    if cidade:
        queries.append(("nominatim_city",
                        f"{cidade}, {uf}, Brasil"))

    for tag, query in queries:
        nom = query_nominatim(query)
        if nom is None:
            continue
        # 4a. Validação geográfica: dentro do estado de SP
        if not _is_in_sp(nom["lat"], nom["lon"]):
            continue
        # 4b. Validação por proximidade ao fallback (se houver)
        if fallback_lat is not None and fallback_lon is not None:
            dist = _haversine_km(nom["lat"], nom["lon"],
                                 fallback_lat, fallback_lon)
            if dist > MAX_DEVIATION_KM:
                print(f"  ⚠ {tag}: {dist:.0f} km do fallback "
                      f"(máx {MAX_DEVIATION_KM} km). Descartado.")
                continue
        # SUCESSO
        result = {
            "lat": nom["lat"],
            "lon": nom["lon"],
            "source": tag,
            "address": nom["address"],
            "validated": True,
            "viacep": ({
                "logradouro": logradouro,
                "bairro": bairro,
                "localidade": cidade,
                "uf": uf,
            } if viacep else None),
        }
        if cep_digits:
            _cache[cep_digits] = result
            _save_cache(_cache)
        return result

    # 5. FALLBACK
    if fallback_lat is not None and fallback_lon is not None:
        addr = f"Fallback para CEP {cep}"
        if city:
            addr += f" ({city})"
        result = {
            "lat": float(fallback_lat),
            "lon": float(fallback_lon),
            "source": "fallback",
            "address": addr,
            "validated": False,
            "viacep": ({
                "logradouro": logradouro,
                "bairro": bairro,
                "localidade": cidade,
                "uf": uf,
            } if viacep else None),
        }
        # NÃO cacheamos fallbacks — para que sejam retentados em execuções futuras
        return result

    raise ValueError(
        f"Não foi possível geocodificar CEP {cep} e nenhum fallback fornecido."
    )


def geocode_patient(patient: dict) -> dict:
    """
    Geocodifica todas as residências de um paciente.

    Args:
        patient: dict com 'id', 'name', 'residences' (lista de residências)

    Returns:
        paciente com campo 'geocoded' adicionado a cada residência.
    """
    result = patient.copy()
    result["residences"] = []

    for res in patient["residences"]:
        geo = geocode_cep(
            cep=res["cep"],
            fallback_lat=res.get("lat"),
            fallback_lon=res.get("lon"),
            city=res.get("city"),
        )
        res_copy = res.copy()
        res_copy["geocoded"] = geo
        result["residences"].append(res_copy)

    return result


# ============================================================
# TESTE STANDALONE
# ============================================================
def _test():
    """Roda alguns CEPs de teste — use com `python3 geocoder.py`."""
    test_ceps = [
        ("01310-100", "São Paulo (Av. Paulista)", -23.5636, -46.6544),
        ("13083-970", "Campinas", -22.8168, -47.0688),
        ("15015-100", "S.J. Rio Preto", -20.8113, -49.3758),
        ("18730-003", "Itaí", -23.4183, -49.0917),
        ("11013-100", "Santos", -23.9618, -46.3322),
    ]
    print("\n🧪 TESTE DO GEOCODER\n" + "=" * 60)
    for cep, name, fb_lat, fb_lon in test_ceps:
        print(f"\n📍 CEP {cep} ({name})")
        try:
            r = geocode_cep(cep, fallback_lat=fb_lat, fallback_lon=fb_lon, city=name)
            print(f"   → ({r['lat']:.4f}, {r['lon']:.4f})")
            print(f"   Fonte: {r['source']}, validado={r['validated']}")
            if r.get("viacep"):
                v = r["viacep"]
                print(f"   ViaCEP: {v['logradouro']}, {v['bairro']}, "
                      f"{v['localidade']}/{v['uf']}")
            print(f"   Endereço: {r['address'][:80]}")
        except Exception as e:
            print(f"   ❌ ERRO: {e}")
    print("\n" + "=" * 60)
    print(f"Cache atual: {len(_cache)} CEPs em {CACHE_PATH}")


if __name__ == "__main__":
    _test()
