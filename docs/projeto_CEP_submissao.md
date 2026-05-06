# PROJETO DE PESQUISA — SUBMISSÃO AO COMITÊ DE ÉTICA EM PESQUISA

---

## 1. IDENTIFICAÇÃO DO PROJETO

**Título**: Exposição cumulativa ao material particulado fino (PM2.5) e risco de câncer de pulmão de não pequenas células em pacientes não tabagistas: estudo caso-controle no Estado de São Paulo

**Título em inglês**: Cumulative fine particulate matter (PM2.5) exposure and risk of non-small cell lung cancer in never-smokers: a case-control study in São Paulo State, Brazil

**Área temática**: Epidemiologia ambiental / Oncologia

**Instituição proponente**: Faculdade de Medicina da Universidade de São Paulo (FMUSP)

**Centro coparticipante**: Instituto do Câncer do Estado de São Paulo (ICESP)

**Pesquisador responsável**: [Nome do orientador]

**Pesquisador executante**: Vitor Hugo Felix — Doutorado Direto FMUSP

**Duração prevista**: 48 meses

---

## 2. RESUMO

Este estudo caso-controle investiga a associação entre exposição cumulativa ao PM2.5 (material particulado fino) e o risco de câncer de pulmão de não pequenas células (CPNPC) em indivíduos não tabagistas no Estado de São Paulo. Utilizando dados do prontuário eletrônico do ICESP e bases ambientais públicas (CETESB, MERRA-2, CAMS, BDQueimadas/INPE e DETRAN-SP), foi desenvolvido um pipeline computacional que estima a exposição individual ao PM2.5 a partir do histórico residencial dos pacientes. O estudo calculará a exposição cumulativa (µg/m³·anos) e avaliará se diferentes fontes de poluição (veicular vs. biomassa) modulam o risco de CPNPC, contribuindo para a compreensão da epidemiologia do câncer de pulmão em não fumantes e para políticas públicas de qualidade do ar.

**Palavras-chave**: PM2.5; câncer de pulmão; não tabagistas; exposição ambiental; caso-controle; São Paulo

---

## 3. INTRODUÇÃO E JUSTIFICATIVA

### 3.1. Contexto

O câncer de pulmão é a principal causa de morte por câncer no mundo, com 1,8 milhão de óbitos anuais (GLOBOCAN, 2022). Embora o tabagismo seja o principal fator de risco, aproximadamente 15–25% dos casos ocorrem em indivíduos que nunca fumaram (never-smokers), proporção que vem aumentando globalmente.

A poluição do ar por material particulado fino (PM2.5 — partículas ≤ 2,5 µm) foi classificada como carcinógeno do Grupo 1 pela IARC/OMS em 2013. Estudos recentes demonstraram mecanismo direto: o PM2.5 promove a expansão clonal de células epiteliais pulmonares portadoras de mutações oncogênicas pré-existentes (especialmente *EGFR* e *KRAS*), sem necessariamente causar novas mutações (Hill et al., *Nature*, 2023).

### 3.2. Lacuna de conhecimento

- **Poucos estudos em países em desenvolvimento**: A maioria das evidências provém de coortes norte-americanas (ACS-CPS II) e europeias (ESCAPE), com níveis de PM2.5 inferiores aos observados no Brasil.
- **Exposição cumulativa pouco estudada**: A maioria dos estudos usa a média de PM2.5 no momento do diagnóstico, sem considerar o histórico residencial.
- **Fontes de PM2.5 negligenciadas**: Em SP, o PM2.5 tem composição variável — veicular na RMSP, queima de biomassa no interior — e a toxicidade pode diferir conforme a fonte.
- **Nenhum estudo brasileiro caso-controle** avaliou especificamente a associação PM2.5 × CPNPC em não tabagistas com dados georeferenciados.

### 3.3. Justificativa

O Estado de São Paulo oferece um cenário único: possui a maior rede de monitoramento da qualidade do ar da América Latina (CETESB, 41 estações), diversidade de fontes emissoras (metrópole + cinturão agrícola), e concentra o maior centro oncológico público do país (ICESP), com registro eletrônico estruturado. Dados ambientais de reanálise global (MERRA-2, CAMS) permitem cobertura espacial e temporal completa do estado, e os dados de queimadas (INPE) e frota veicular (DETRAN) viabilizam a identificação das fontes predominantes.

---

## 4. HIPÓTESE

A exposição cumulativa ao PM2.5 está associada a maior risco de CPNPC em indivíduos não tabagistas residentes no Estado de São Paulo, e a magnitude desta associação pode diferir conforme a fonte predominante de poluição (veicular vs. queima de biomassa).

---

## 5. OBJETIVOS

### 5.1. Objetivo primário

Avaliar a associação entre a exposição cumulativa ao PM2.5 e o risco de CPNPC em pacientes não tabagistas atendidos no ICESP.

### 5.2. Objetivos secundários

1. Estimar a exposição individual ao PM2.5 a partir do histórico residencial e de dados ambientais integrados (CETESB, MERRA-2, CAMS);
2. Classificar os municípios de exposição quanto à fonte predominante de PM2.5 (veicular, biomassa, mista);
3. Avaliar se a fonte predominante de PM2.5 modifica a associação com o risco de CPNPC;
4. Gerar superfícies anuais de PM2.5 para o Estado de São Paulo (2008–2024) como produto cartográfico.

---

## 6. METODOLOGIA

### 6.1. Desenho do estudo

Estudo observacional, analítico, do tipo **caso-controle**, retrospectivo, com base hospitalar.

### 6.2. Local

Instituto do Câncer do Estado de São Paulo Octavio Frias de Oliveira (ICESP-FMUSP), São Paulo, SP.

### 6.3. População

#### Casos
- Pacientes com diagnóstico histologicamente confirmado de CPNPC (adenocarcinoma, carcinoma de células escamosas ou carcinoma de grandes células);
- Atendidos no ICESP entre 2008 e 2024;
- Classificados como **não tabagistas** (nunca fumou) no prontuário eletrônico.

#### Controles
- Pacientes atendidos no ICESP no mesmo período, com outros diagnósticos oncológicos (excluindo neoplasias do trato respiratório e cânceres tabaco-relacionados);
- Classificados como **não tabagistas** no prontuário eletrônico;
- Pareados por sexo, faixa etária (±5 anos) e período de atendimento;
- Proporção: 1 caso : 2–3 controles.

### 6.4. Critérios de inclusão

| Critério | Casos | Controles |
|----------|:-----:|:---------:|
| Diagnóstico anatomopatológico de CPNPC | ✅ | ❌ |
| Outro diagnóstico oncológico (não respiratório) | ❌ | ✅ |
| Registro de "não tabagista" / "nunca fumou" no prontuário | ✅ | ✅ |
| Endereço residencial (CEP, endereço ou bairro/cidade) no prontuário | ✅ | ✅ |
| Residência no Estado de São Paulo | ✅ | ✅ |
| Idade ≥ 18 anos ao diagnóstico | ✅ | ✅ |

### 6.5. Critérios de exclusão

- Pacientes com registro de tabagismo atual ou pregresso;
- Pacientes sem informação de tabagismo no prontuário;
- Pacientes sem qualquer informação de endereço residencial no prontuário;
- Pacientes com residência fora do Estado de São Paulo;
- Diagnósticos histológicos indeterminados ou tumor primário desconhecido;
- Para controles: cânceres da cabeça e pescoço, esôfago, bexiga e pâncreas (por associação com tabagismo passivo ou ambiental).

### 6.6. Estimativa do tamanho amostral

Com base na literatura (HR ≈ 1,08 por 10 µg/m³ de PM2.5 — Huang et al., 2021), considerando:
- Nível de significância α = 0,05 (bicaudal)
- Poder estatístico (1−β) = 0,80
- Razão caso:controle = 1:2
- Prevalência de exposição elevada (> mediana) entre controles ≈ 50%
- OR esperado = 1,5–2,0

Estima-se necessidade de **150–200 casos** e **300–400 controles** (total: 450–600 participantes).

### 6.7. Coleta de dados

#### Fase 1 — Dados de prontuário (exclusivamente)
Extração dos seguintes dados do prontuário eletrônico do ICESP:
- Dados demográficos: sexo, data de nascimento, raça/cor
- Diagnóstico oncológico: tipo histológico, estadiamento (TNM), data do diagnóstico
- Status de tabagismo: conforme registro em anamnese
- Endereço(s) residencial(is): CEP, logradouro, bairro, cidade, com período de residência quando disponível
- Exposição ocupacional: ocupação registrada
- Comorbidades relevantes: DPOC, doença pulmonar prévia

> **Não haverá contato direto com pacientes na Fase 1.** Todos os dados são provenientes do prontuário eletrônico.

#### Fase 2 — Complementação por contato telefônico (quando necessário)
Para pacientes com dados incompletos no prontuário (especialmente endereço residencial e confirmação de tabagismo):
- Contato telefônico com o paciente ou, em caso de óbito, com familiar/representante legal (respondente substituto — *proxy respondent*);
- Aplicação de questionário estruturado breve (Apêndice A).

### 6.8. Exposição ao PM2.5 — Pipeline computacional

O desfecho primário de exposição será estimado por pipeline computacional desenvolvido especificamente para este projeto, utilizando:

**Fontes de dados ambientais (todos públicos):**

| Fonte | Tipo | Resolução | Período |
|-------|------|-----------|---------|
| CETESB | Monitoramento terrestre | Pontual (41 estações) | 2012–2024 |
| MERRA-2 (NASA) | Reanálise global | 0,5° × 0,625° (~55 km) | 1980–2025 |
| CAMS EAC4 (ECMWF) | Reanálise global | 0,75° (~83 km) | 2008–2025 |
| BDQueimadas (INPE) | Focos de calor | Pontual | 2008–2024 |
| DETRAN-SP | Frota veicular ativa | Município | 2025 |

**Método de fusão de dados:**
1. Geocodificação do endereço residencial (CEP → coordenadas lat/lon)
2. Consulta a superfícies anuais de PM2.5 pré-computadas (resolução 0,05° ≈ 5,5 km), geradas por fusão de MERRA-2 (base) corrigido com resíduos CETESB via interpolação IDW (*Inverse Distance Weighting*)
3. Cálculo da exposição cumulativa: **Expo_cum = Σ (PM2.5_anual × tempo_residência)**
4. Classificação do município quanto à fonte predominante: veicular, biomassa, mista ou baixa emissão

**Validação do pipeline:** Validação cruzada (*leave-one-out*) do método de fusão nas 41 estações CETESB demonstrou RMSE = 2,86 µg/m³ e melhoria de 14,7% no MAE em relação ao MERRA-2 puro.

### 6.9. Variáveis do estudo

| Variável | Tipo | Descrição |
|----------|------|-----------|
| **Desfecho**: CPNPC | Categórica | Caso vs. Controle |
| **Exposição principal**: Expo_cum | Contínua | µg/m³·anos |
| PM2.5 médio anual | Contínua | µg/m³ |
| Fonte predominante | Categórica | Veicular / Biomassa / Mista / Baixa |
| Sexo | Categórica | Masculino / Feminino |
| Idade ao diagnóstico | Contínua | Anos |
| Raça/cor | Categórica | Autodeclarada |
| Tipo histológico | Categórica | Adenoca / CEC / Grandes cél. |
| Estadiamento | Ordinal | I–IV |
| Exposição ocupacional | Categórica | Sim/Não (registro) |
| Focos de calor no raio | Contínua | Nº acumulado |

---

## 7. ANÁLISE ESTATÍSTICA

1. **Análise descritiva**: frequências, médias, medianas e distribuições por grupo
2. **Comparação univariada**: teste t/Mann-Whitney para variáveis contínuas; qui-quadrado/Fisher para categóricas
3. **Regressão logística condicional**: OR (IC 95%) da associação PM2.5 → CPNPC, ajustado por covariáveis
4. **Análise dose-resposta**: OR por incremento de 10 µg/m³·anos; modelagem por splines cúbicos restritos para avaliar linearidade
5. **Análise de subgrupos**: estratificação por fonte predominante de PM2.5 (veicular vs. biomassa vs. mista)
6. **Análise de sensibilidade**: excluir pacientes geocodificados apenas a nível de cidade; variar definição de "não tabagista"

**Software**: Python 3.x (scipy, statsmodels, scikit-learn), R 4.x.

**Nível de significância**: p < 0,05 (bicaudal).

---

## 8. RISCOS E BENEFÍCIOS

### 8.1. Riscos

**Fase 1 (prontuário)**: Risco mínimo, limitado a potencial quebra de sigilo dos dados. Mitigação:
- Dados armazenados em computador com acesso restrito e senha
- Anonimização no momento da extração (código alfanumérico)
- Nenhum dado identificador será incluído em publicações ou apresentações

**Fase 2 (telefone)**: Risco de desconforto emocional ao abordar diagnóstico oncológico e/ou óbito de familiar. Mitigação:
- Abordagem por profissional de saúde treinado
- Esclarecimento prévio sobre objetivo da pesquisa
- Liberdade para recusar participação a qualquer momento
- Encaminhamento para apoio psicossocial se necessário

### 8.2. Benefícios

- **Diretos ao participante**: Não há benefício direto individual
- **Indiretos/coletivos**:
  - Contribuição para compreensão da epidemiologia do câncer de pulmão em não fumantes no Brasil
  - Subsídio para políticas públicas de controle da qualidade do ar em SP
  - Produção de mapas anuais de PM2.5 para SP como ferramenta de vigilância em saúde ambiental

---

## 9. ASPECTOS ÉTICOS

### 9.1. Base legal

Este projeto será conduzido em conformidade com:
- Resolução CNS nº 466/2012 e Resolução CNS nº 510/2016
- Declaração de Helsinki (versão 2013)
- Lei Geral de Proteção de Dados (LGPD — Lei nº 13.709/2018)

### 9.2. Termo de Consentimento Livre e Esclarecido (TCLE)

**Fase 1**: Solicita-se **dispensa do TCLE** com fundamento no Art. 1º, parágrafo único, inciso V da Resolução CNS 510/2016, considerando que:
- Trata-se de pesquisa com dados secundários, obtidos exclusivamente de prontuário eletrônico
- Não haverá contato com os participantes
- A obtenção do TCLE tornaria a pesquisa inviável, dado o grande número de participantes e a elevada taxa de óbito na população-alvo
- Todos os dados serão anonimizados no momento da extração
- Não há risco adicional ao participante além do já existente na assistência

**Fase 2**: Para pacientes contatados por telefone, será obtido **TCLE verbal** (gravação autorizada) ou TCLE enviado por meio digital (e-mail/WhatsApp), conforme Carta Circular nº 1/2021 da CONEP. Para familiares de pacientes falecidos (*proxy respondent*), será obtido TCLE específico.

### 9.3. Respondente substituto (proxy respondent)

Justificativa metodológica: Dada a alta mortalidade do CPNPC (sobrevida em 5 anos < 20%), uma proporção significativa dos pacientes-caso terá evoluído a óbito. A exclusão sistemática desses pacientes introduziria viés de sobrevivência (*survival bias*), comprometendo a validade interna do estudo. O uso de respondente substituto (preferencialmente cônjuge ou filho/a) é prática consolidada em estudos epidemiológicos de caso-controle em oncologia (McLaughlin et al., *Am J Epidemiol*, 1990; Percy et al., *Am J Epidemiol*, 1981).

**Variáveis coletadas com proxy**: endereço residencial e período de residência. A informação sobre tabagismo será obtida exclusivamente do prontuário eletrônico.

### 9.4. Armazenamento e proteção de dados

- Dados armazenados em repositório institucional com acesso restrito
- Anonimização por código alfanumérico; chave de identificação em arquivo separado, criptografado
- Dados ambientais (PM2.5, queimadas, frota) são públicos e não contêm informação pessoal
- Período de guarda: 5 anos após término da pesquisa, conforme Resolução CNS 466/2012
- Dados destruídos ao término do período de guarda

---

## 10. CRONOGRAMA

| Etapa | Semestre 1–2 | Semestre 3–4 | Semestre 5–6 | Semestre 7–8 |
|-------|:---:|:---:|:---:|:---:|
| Revisão da literatura | ██ | █ | | |
| Desenvolvimento do pipeline PM2.5 | ██ | | | |
| Submissão e aprovação CEP | ██ | | | |
| Coleta de dados — Fase 1 (prontuário) | | ██ | | |
| Coleta de dados — Fase 2 (telefone) | | █ | █ | |
| Geocodificação e cálculo de exposição | | | ██ | |
| Análise estatística | | | █ | ██ |
| Redação do Artigo 1 (metodológico) | | ██ | | |
| Redação do Artigo 2 (caso-controle) | | | | ██ |
| Redação e defesa da tese | | | | ██ |

---

## 11. ORÇAMENTO

| Item | Custo estimado |
|------|:-:|
| Computador para processamento de dados | Já disponível |
| Software (Python, R) | Gratuito (open-source) |
| API de geocodificação (Nominatim) | Gratuito |
| Dados ambientais (CETESB, MERRA-2, CAMS, INPE) | Gratuito (acesso público) |
| Telefonia para Fase 2 | R$ 500,00 |
| Material de escritório | R$ 200,00 |
| **Total** | **R$ 700,00** |

> Não há financiamento externo. Os custos serão absorvidos pelo pesquisador.

---

## 12. REFERÊNCIAS BIBLIOGRÁFICAS PRINCIPAIS

1. Bray F, et al. Global cancer statistics 2022: GLOBOCAN. *CA Cancer J Clin*. 2024;74(3):229-263.
2. Hill W, et al. Lung adenocarcinoma promotion by air pollutants. *Nature*. 2023;616(7955):159-167.
3. Huang Y, et al. Air pollution and lung cancer incidence in a fine-scale spatial analysis. *Environ Health Perspect*. 2021;129(9):097002.
4. Loomis D, et al. The carcinogenicity of outdoor air pollution. *Lancet Oncol*. 2013;14(13):1262-1263.
5. Turner MC, et al. Long-term ambient fine particulate matter air pollution and lung cancer in a large cohort of never-smokers. *Am J Respir Crit Care Med*. 2011;184(12):1374-1381.
6. McLaughlin JK, et al. Use of proxy respondents in environmental epidemiology. *Am J Epidemiol*. 1990;131(5):761-765.
7. CETESB. Qualidade do ar no Estado de São Paulo — Relatório Anual. 2024.
8. Gelaro R, et al. The Modern-Era Retrospective Analysis for Research and Applications, version 2 (MERRA-2). *J Climate*. 2017;30(14):5419-5454.

---

## APÊNDICE A — QUESTIONÁRIO ESTRUTURADO (Fase 2 — Contato Telefônico)

**Identificação**: Código do participante: _____ | Data: ___/___/_____ | Respondente: ( ) Paciente ( ) Familiar

---

**Bloco 1 — Confirmação de tabagismo** *(apenas se houver dúvida no prontuário)*

1. O(A) Sr(a) / [nome do paciente] já fumou cigarros, charutos, cachimbo ou cigarros eletrônicos em algum momento da vida?
   - ( ) Nunca fumou → Prosseguir
   - ( ) Já fumou → Critério de exclusão — Encerrar

---

**Bloco 2 — Histórico residencial**

2. Nos últimos 20 anos, em quais endereços o(a) Sr(a) / [nome do paciente] morou? *(Começar pelo mais recente)*

| Nº | Endereço (rua, nº, bairro) | Cidade | CEP (se souber) | De (ano) | Até (ano) |
|:--:|---|---|---|---|---|
| 1 | | | | | |
| 2 | | | | | |
| 3 | | | | | |

3. Algum desses endereços era próximo a:
   - ( ) Rodovia ou via de tráfego pesado (< 200m)
   - ( ) Indústria / fábrica
   - ( ) Área de queimadas frequentes (canaviais, pastagens)
   - ( ) Nenhuma das anteriores

---

**Bloco 3 — Exposição ocupacional** *(se não registrada no prontuário)*

4. Qual a principal ocupação do(a) Sr(a) / [nome do paciente] ao longo da vida?
   - Ocupação: _________________
   - Teve contato frequente com poeira, fumaça, químicos ou solventes no trabalho? ( ) Sim ( ) Não

---

*Duração estimada: 5–10 minutos*

---

## APÊNDICE B — TERMO DE CONSENTIMENTO LIVRE E ESCLARECIDO (Fase 2)

### TCLE — Para o próprio paciente

**Título da pesquisa**: Exposição cumulativa ao material particulado fino (PM2.5) e risco de câncer de pulmão de não pequenas células em pacientes não tabagistas

**Pesquisador responsável**: [Nome] — Telefone: [XX] XXXX-XXXX — E-mail: [email]

Você está sendo convidado(a) a participar de uma pesquisa que busca entender se a poluição do ar pode estar relacionada ao câncer de pulmão em pessoas que nunca fumaram. Sua participação é voluntária.

**O que será feito**: Faremos algumas perguntas por telefone sobre seus endereços de moradia nos últimos 20 anos e sua ocupação profissional, com duração de aproximadamente 10 minutos. Não serão realizados exames ou procedimentos.

**Riscos**: O risco é mínimo, limitado a possível desconforto ao falar sobre seu diagnóstico. Você pode interromper a qualquer momento.

**Benefícios**: Não há benefício direto para você, mas os resultados podem ajudar a entender melhor as causas do câncer de pulmão e contribuir para políticas de saúde.

**Sigilo**: Suas informações serão mantidas em sigilo. Seu nome não aparecerá em nenhuma publicação.

**Liberdade**: Sua participação é voluntária. Você pode desistir a qualquer momento, sem prejuízo ao seu atendimento no ICESP.

( ) Concordo em participar da pesquisa e autorizo a gravação desta ligação como registro do meu consentimento.

---

### TCLE — Para familiar/representante legal (respondente substituto)

**Título da pesquisa**: [idem acima]

O(A) Sr(a) está sendo convidado(a) a participar como representante de [nome do paciente] em uma pesquisa sobre poluição do ar e câncer de pulmão. Sua participação consiste em fornecer informações sobre os endereços de moradia de [nome do paciente] nos últimos anos.

O uso de respondente substituto (familiar) é uma prática consolidada em pesquisas médicas e é necessário neste estudo porque muitos pacientes com câncer de pulmão já faleceram. Sua contribuição é importante para que possamos compreender melhor esta doença.

[Demais itens idênticos ao TCLE do paciente]

( ) Concordo em participar como respondente substituto e autorizo a gravação desta ligação como registro do meu consentimento.

**Grau de parentesco com o paciente**: _______________

---

## APÊNDICE C — JUSTIFICATIVA PARA DISPENSA DE TCLE (Fase 1)

Ao Comitê de Ética em Pesquisa,

Solicito a dispensa da obtenção do Termo de Consentimento Livre e Esclarecido para a Fase 1 deste projeto, que consiste exclusivamente na coleta de dados secundários de prontuários eletrônicos, pelos seguintes motivos:

1. **Natureza retrospectiva e não intervencionista**: não haverá qualquer contato com pacientes, familiares ou profissionais de saúde durante a Fase 1;

2. **Inviabilidade prática**: a população estudada apresenta elevada taxa de mortalidade (sobrevida em 5 anos do CPNPC < 20%), tornando impossível a obtenção de consentimento de parte significativa dos participantes;

3. **Anonimização integral**: os dados serão anonimizados no momento da extração, com substituição de identificadores diretos (nome, RGH, CPF) por código alfanumérico sequencial;

4. **Risco mínimo**: não há risco adicional ao participante além do inerente ao armazenamento de dados em prontuário eletrônico;

5. **Interesse público**: a pesquisa aborda questão de saúde pública relevante (relação poluição do ar e câncer), e a exigência de TCLE para todos os participantes inviabilizaria o estudo;

6. **Fundamentação legal**: Resolução CNS nº 510/2016, Art. 1º, parágrafo único, inciso V; Resolução CNS nº 466/2012, item IV.8.

Comprometo-me a manter o sigilo absoluto sobre a identidade dos participantes e a utilizar os dados exclusivamente para os fins descritos neste protocolo.

Atenciosamente,

[Pesquisador executante] | [Pesquisador responsável]
