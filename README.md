# pm25-cpnpc

[![tests](https://github.com/parabolicas/pm25copia/actions/workflows/tests.yml/badge.svg)](https://github.com/parabolicas/pm25copia/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Pipeline geoespacial para estimativa individualizada de exposição cumulativa ao PM₂.₅ no Estado de São Paulo, desenvolvido para o estudo de doutorado direto **"A Assinatura Ambiental do Câncer de Pulmão em Não-Tabagistas: Associação Dose-Resposta entre Exposição Cumulativa ao PM₂.₅ e o Perfil Mutacional do CPNPC no Estado de São Paulo"** (Faculdade de Medicina, USP).

## Visão geral

O pipeline integra dados terrestres (CETESB), satelitais (NASA MERRA-2, ECMWF CAMS), focos de calor (INPE BDQueimadas), dados de frota (DETRAN-SP) e polígonos oficiais (IBGE Censo 2022) para gerar superfícies anuais bias-corrected de PM₂.₅ na resolução 0,05° (~5,5 km). Para cada paciente, a residência é geocodificada (CEP → coordenadas via ViaCEP + Nominatim) e a exposição cumulativa é calculada como Σ (PM₂.₅ anual × fração-de-ano residencial), com precisão diária (suporte a anos bissextos e mudanças residenciais em qualquer época do ano).

## Arquitetura

```
[CEP do paciente]
      │
      ▼
geocoder.py  ──►  ViaCEP (logradouro estruturado)
                      │
                      ▼
                 Nominatim (lat/lon, validação ≤50 km do município)
                      │
                      ▼
              cumulative.py  ──►  exposure.py  ──►  pm25_surface.lookup_pm25_buffer
                      │                                       │
                      │                       ┌───────────────┴───────────────┐
                      │                       │                               │
                      │              Superfície anual                  Diagnóstico
                      │              (GeoTIFF bias-corrected)          (CETESB IDW,
                      │              MERRA-2 + IDW resíduos             CAMS, MERRA-2
                      │              CETESB)                             brutos)
                      ▼
              Σ (PM₂.₅ × fração-de-ano)
                      │
                      ▼
              exposure_summary.csv + validation_report.txt
```

Os componentes principais ficam em arquivos Python independentes:

| Módulo | Responsabilidade |
|---|---|
| `config.py` | Caminhos, parâmetros (buffers 5/10/25 km), pacientes-modelo, fórmula NASA MERRA-2 |
| `data_loaders.py` | Carregamento de CETESB, MERRA-2, CAMS, BDQueimadas, mesorregiões |
| `geocoder.py` | CEP → coordenadas via ViaCEP + Nominatim, com cache local |
| `pm25_surface.py` | Geração das 17 superfícies anuais bias-corrected (2008–2024) e funções `lookup_pm25` / `lookup_pm25_buffer` |
| `exposure.py` | Estimativa de PM₂.₅ pontual (fonte primária = superfície bias-corrected) |
| `cumulative.py` | Exposição cumulativa com ponderação por fração-de-ano residencial |
| `cross_validation.py` | CETESB vs. CAMS/MERRA-2 brutos (282 pares estação-ano) |
| `source_apportionment.py` | Classificação municipal: VEICULAR / BIOMASSA / MISTA / BAIXA EMISSÃO (645 municípios, áreas IBGE oficial) |
| `source_specific.py` | Frações veicular/biomassa/outros calibradas por PMF (Pereira et al. 2025) + KDE focos BDQueimadas |
| `sensibility_h.py` | Análise de sensibilidade da incerteza temporal de mudança residencial |
| `run_pipeline.py` | Script principal: geocodificação → exposição → validação → mapas |

## Como reproduzir

### 1. Pré-requisitos

- Python 3.10+ (testado em 3.10, 3.11, 3.12 via CI; também roda em 3.9)
- ~2 GB de espaço em disco para dados de entrada e saídas
- Conexão à internet para a primeira execução (ViaCEP + Nominatim)

### 2. Instalar dependências

```bash
pip install -r requirements.txt
```

### 3. Estrutura de dados de entrada esperada

```
pm25/
├── dados_pm25_cetesb/
│   ├── RESUMO_MEDIAS_ANUAIS_PM25.csv     # 282 registros estação-ano
│   └── cetesb_estacoes_coordenadas.csv   # 41 estações com lat/lon
├── MERRA2/
│   └── MERRA2_*.tavgM_2d_aer_Nx.YYYYMM.nc4.nc4   # mensais 1980–2025
├── CAMS/
│   └── data_sfcYY.nc                      # anuais 2008–2025
├── BDQueimadas/
│   └── focos_br_sp_ref_YYYY.zip           # focos por ano
├── BR_Municipios_2022/
│   └── BR_Municipios_2022.shp             # IBGE Censo 2022 (645 mun. SP)
├── SP_Mesorregioes_2022/
│   └── SP_Mesorregioes_2022.shp
└── frota_ativa_sp.csv                     # DETRAN-SP municipal
```

### 4. Geração das superfícies anuais (uma vez)

```bash
python3 pm25_surface.py
```

Saídas em `output/surfaces/`:
- 17 GeoTIFFs anuais (`pm25_sp_2008.tif` … `pm25_sp_2024.tif`)
- `surface_summary.csv` (estatísticas por ano)
- `pm25_annual_panel.png` (painel multi-ano)
- `pm25_trend.png` (tendência temporal)
- `leave_one_out_2023.csv` (validação LOO: MAE 2,13 µg/m³, R = 0,64)

### 5. Pipeline de exposição individual

```bash
python3 run_pipeline.py
```

Saídas em `output/`:
- `exposure_details.csv` — uma linha por (paciente × ano × buffer × fração-de-ano)
- `exposure_summary.csv` — uma linha por (paciente × buffer)
- `validation_report.txt` — 21 testes automatizados de plausibilidade
- `exposure_map.png` — mapa com pacientes, buffers e estações CETESB
- `exposure_comparison.png` — gráfico comparativo

### 6. Source apportionment municipal

```bash
python3 source_apportionment.py    # 645 municípios → áreas IBGE oficial
python3 source_specific.py         # 51 GeoTIFFs com calibração dual
```

Análise de sensibilidade da calibração dual (boost factor):

```bash
python3 source_specific.py --boost 1.2 --output-suffix _boost12
python3 source_specific.py --boost 1.6 --output-suffix _boost16
```

### 7. Análise de sensibilidade temporal

```bash
python3 sensibility_h.py    # cenários 01-jan / 01-jul / 31-dez
```

### 8. Validação cruzada com satélite

```bash
python3 cross_validation.py        # 282 pares estação-ano: CETESB vs. CAMS/MERRA-2
```

## Testes automatizados

```bash
python3 -m pytest -v
```

Saída esperada (26 testes):

```
tests/test_basics.py::TestMERRA2Formula  ............................  6 PASSED
tests/test_basics.py::TestCepNormalization  .........................  7 PASSED
tests/test_basics.py::TestFractionOfYear  ...........................  8 PASSED
tests/test_source_specific.py::TestFireDensityGrid  .................  2 PASSED
tests/test_source_specific.py::TestSourceSpecificBoost  .............  3 PASSED
==================== 26 passed in ~3s ====================
```

Os testes cobrem cinco áreas-chave:

1. **Fórmula NASA MERRA-2** (6 testes) — confirma que os fatores estequiométricos (1,0; 1,0; 1,0; 1,4; 1,375) não foram alterados.
2. **Normalização de CEPs** (7 testes) — aceita formatos brasileiros variados ("01310-100", "01310100", com espaços), rejeita inválidos.
3. **Fração de ano residencial** (8 testes) — valida cálculo dia-a-dia, incluindo anos bissextos (366 dias em 2020 e 2024) e mudanças em qualquer mês.
4. **Fire density grid (KDE)** (2 testes) — shape consistente com a grade alvo, comportamento com array vazio.
5. **Source-specific boost dual** (3 testes) — boost factor=1.0 mantém calibração neutra; factor>1.0 aumenta `w_biomassa`; renormalização preserva soma=1.0.

## Integração contínua (CI/CD)

A cada `git push` ou pull-request em `main`/`master`, o GitHub Actions executa o conjunto de testes em três versões do Python (3.10, 3.11, 3.12) sobre Ubuntu Linux. A *badge* "tests" no topo deste README reflete o status atual: verde quando todos os 26 testes passam, vermelho quando algum falha.

Configuração: [`.github/workflows/tests.yml`](.github/workflows/tests.yml). Tempo médio de execução: ~60-90 segundos por job.

## Validação científica

### Cross-validation CETESB × satélite (282 estação-ano, 2012–2024)

| Fonte | Bias (µg/m³) | MAE | R² | Conclusão |
|---|---|---|---|---|
| CAMS EAC4 (bruto) | **+16,3** | 16,6 | 0,19 | Superestimação sistemática |
| MERRA-2 (bruto) | −3,1 | 3,6 | 0,12 | Subestimação leve |
| **MERRA-2 + CETESB IDW (fusão)** | **+0,3** | **2,1** | **0,64** | Adotado como base |

### Leave-one-out 2023 (35 estações)

A superfície bias-corrected reduz o MAE em **14,8%** em relação ao MERRA-2 puro.

### Validação cruzada interna source-specific (Spearman)

A calibração dual de fração biomassa (PMF Pereira 2025 + KDE focos BDQueimadas) é validada por correlação de Spearman entre `w_biomassa` estimado e densidade de focos por município. Resultados em três cenários de boost factor β:

| β | ρ Spearman médio | w_bio_RMSP | w_bio_SP estado | Status |
|---|---|---|---|---|
| 1,2 | 0,451 | 27% | 35% | ⚠ Abaixo do limiar |
| **1,4 (default)** | **0,590** | 31% | 37% | ✅ Equilíbrio |
| 1,6 | 0,669 | 35% | 38% | ✅ Reproduz melhor o padrão espacial |

### Pipeline de validação automatizado

`run_pipeline.py` executa 21 testes a cada execução, cobrindo: PM₂.₅ médio em SP centro vs. literatura (Pereira et al. 2025; Vasconcelos & Miranda 2025), fonte primária ≥ 80% SURFACE, frações de ano não-sobrepostas, soma de frações = duração residencial, variação espacial entre buffers e plausibilidade de coordenadas.

## Referências

- Hill W et al. *Lung adenocarcinoma promotion by air pollutants.* Nature 2023;616:159–167.
- Pereira GM et al. *Source apportionment and ecotoxicity of PM2.5 pollution events in a major Southern Hemisphere megacity.* Atmos Chem Phys 2025;25:4587–4612.
- Vasconcelos GM, Miranda RM. *Characterization and source identification of fine particulate matter in the Eastern Zone of São Paulo, Brazil.* Atmos Pollut Res 2025;16(5):102426.
- Cortot A et al. *Fine Particulate Matter Affects EGFR-Mutated Lung Cancers — KPB-2020.* Ann Oncol 2024;35(S2):S937.

## Limitações

- **Resolução espacial.** A grade de 0,05° (~5,5 km) não captura gradientes intra-urbanos finos (ex.: proximidade a vias de alto tráfego). Análises de sensibilidade com buffers de 5/10/25 km mitigam parcialmente.
- **Anos pré-CETESB (2008–2011).** Para anos sem estações CETESB, aplica-se correção global de bias (+3,1 µg/m³) em vez de IDW de resíduos. As superfícies desses anos têm incerteza maior.
- **Ausência de PMF dedicado ao interior agrícola.** A calibração source-specific da fração biomassa em áreas não-RMSP utiliza modulador empírico baseado em densidade de focos BDQueimadas (KDE), validada por correlação interna ρ Spearman ≥ 0,5. Limitação reconhecida e mitigada por análise de sensibilidade β ∈ {1,2; 1,4; 1,6}.
- **Geocodificação.** Depende da cobertura ViaCEP (~99% dos CEPs ativos) e Nominatim. Casos de CEP extinto ou inconsistente caem em fallback (sinalizado por `validated=False` no log).

## Citação

Quando este pipeline for usado, cite:

```bibtex
@misc{pm25_cpnpc_2026,
  author       = {Felix, Vitor and colaboradores},
  title        = {pm25-cpnpc: pipeline geoespacial para exposição cumulativa
                  ao PM2.5 em estudos de epidemiologia oncológica},
  year         = {2026},
  howpublished = {\url{https://github.com/parabolicas/pm25copia}},
}
```

## Contato

Vitor Felix — vitor.felixx@gmail.com
Faculdade de Medicina, Universidade de São Paulo

## Licença

A definir. Sugerido: MIT (código) + CC-BY-4.0 (dados derivados).
