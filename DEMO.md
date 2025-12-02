# Demo Completo - Rosetta Alanine Scanning

Este guia mostra como rodar o protocolo completo com um exemplo prático.

## Pré-requisitos

Abra o terminal no VSCode e execute:

```bash
# 1. Entre no diretório do projeto
cd /Users/madsonluna/Documents/alaninescanning4eif4e

# 2. Crie e ative o ambiente virtual
python3 -m venv venv
source venv/bin/activate

# 3. Instale as dependências
pip install -e .
```

## Demo Rápido (SEM Rosetta)

### Passo 1: Executar o Script de Demonstração

```bash
# Execute o demo completo
python3 examples/demo_run.py
```

### O que você verá:

```
================================================================================
ROSETTA ALANINE SCANNING - COMPLETE DEMO
================================================================================

Working directory: examples/demo_output
Input structure: example_protein.pdb

--------------------------------------------------------------------------------
> STEP 1: Structure Analysis
--------------------------------------------------------------------------------

Loading PDB structure...

OK Structure loaded successfully!

  Chains found: ['A', 'B']
    • Chain A: 5 residues
    • Chain B: 4 residues

--------------------------------------------------------------------------------
> STEP 2: Generate Alanine Mutations
--------------------------------------------------------------------------------

Generating alanine scanning mutations...

OK Generated 6 mutations

Mutation breakdown by chain:
  • Chain A: 4 mutations
  • Chain B: 2 mutations

First 10 mutations to be tested:
  #    Mutation     Chain    Position   Original→Ala
  ------------------------------------------------------------
  1    A1M          A        1          M→A
  2    A2L          A        2          L→A
  3    A3V          A        3          V→A
  ...
```

### Passo 2: Explorar os Resultados

```bash
# Veja a estrutura de arquivos gerada
ls -la examples/demo_output/

# Visualize o relatório de mutações
cat examples/demo_output/mutations.txt

# Veja os hotspots identificados
cat examples/demo_output/hotspots.csv

# Leia o relatório completo
cat examples/demo_output/analysis_report.txt
```

##  Executando com CLI (Passo a Passo)

### Exemplo 1: Gerar Mutações

```bash
# Scan de todas as cadeias
rosetta-scan scan examples/example_protein.pdb \
    --output examples/scan_output/
```

**Saída esperada:**
```
Alanine Scanning Mutation Generator

Loading structure... --------------------------------- 100%
Generating mutations... ------------------------------ 100%

----------------------------------------┓
| Property           | Value            |
┡--------------------╇------------------┩
| Total Mutations    | 6                |
| Structure File     | example_protein.pdb |
| Chain A            | 4                |
| Chain B            | 2                |
----------------------------------------┘

Mutations saved to: examples/scan_output/mutations.txt
```

### Exemplo 2: Scan Apenas da Interface

```bash
# Identificar apenas resíduos de interface
rosetta-scan scan examples/example_protein.pdb \
    --chains A B \
    --interface-only \
    --interface-cutoff 8.0 \
    --output examples/interface_scan/
```

**Saída esperada:**
```
Alanine Scanning Mutation Generator

Loading structure... OK
Identifying interface residues... OK
Generating mutations... OK

----------------------------------------┓
| Property           | Value            |
┡--------------------╇------------------┩
| Total Mutations    | 3                |
| Structure File     | example_protein.pdb |
| Chain A            | 2                |
| Chain B            | 1                |
----------------------------------------┘

Mutations saved to: examples/interface_scan/mutations.txt
```

### Exemplo 3: Analisar Resultados (Demo com Dados Simulados)

```bash
# Primeiro, vamos criar resultados simulados para análise
python3 -c "
import pandas as pd
import numpy as np
from pathlib import Path

# Simula resultados do Rosetta
np.random.seed(42)
mutations = ['A2L', 'A3V', 'A4K', 'A5E', 'B2W', 'B3F']
results = []

for mut in mutations:
    chain = mut[0]
    pos = int(mut[1])

    # Simula ddG (alguns são hotspots)
    if pos in [2, 4]:
        ddg = np.random.uniform(1.5, 3.0)
    else:
        ddg = np.random.uniform(0.2, 1.2)

    results.append({
        'mutation': mut,
        'chain': chain,
        'position': pos,
        'ddg': ddg
    })

df = pd.DataFrame(results)
Path('examples/demo_results').mkdir(exist_ok=True)
df.to_csv('examples/demo_results/results.csv', index=False)
print('OK Demo results created')
"

# Agora podemos usar o visualizador
python3 -c "
import sys
sys.path.insert(0, 'src')
import pandas as pd
from rosetta_scan.analysis.visualizer import ResultVisualizer

# Carrega resultados
df = pd.read_csv('examples/demo_results/results.csv')

print('\n=== Análise de Resultados ===\n')
print(f'Total de mutações: {len(df)}')
print(f'ΔΔG médio: {df[\"ddg\"].mean():.2f} kcal/mol')
print(f'ΔΔG máximo: {df[\"ddg\"].max():.2f} kcal/mol')

# Identifica hotspots
hotspots = df[df['ddg'] > 1.5].sort_values('ddg', ascending=False)
print(f'\nHotspots identificados (ΔΔG > 1.5): {len(hotspots)}')
print('\nTop hotspots:')
for i, row in hotspots.iterrows():
    print(f'  • {row[\"mutation\"]}: {row[\"ddg\"]:.2f} kcal/mol')
"
```

## Estrutura dos Arquivos de Saída

### 1. Arquivo de Mutações (mutations.txt)

```
Alanine Scanning Mutation Report
==================================================

Structure: example_protein.pdb
Total mutations: 6

Chain    Position  Original  Target
--------------------------------------------------
A        1         M         A
A        2         L         A
A        3         V         A
A        4         K         A
A        5         E         A
B        2         W         A
B        3         F         A
```

### 2. Formato Rosetta (mutations_rosetta.txt)

```
total 1
1
A1M
1
A2L
1
A3V
1
A4K
1
A5E
1
B2W
1
B3F
```

### 3. Formato CSV (mutations.csv)

```csv
chain,position,original_aa,target_aa,mutation
A,1,M,A,M1A (chain A)
A,2,L,A,L2A (chain A)
A,3,V,A,V3A (chain A)
A,4,K,A,K4A (chain A)
A,5,E,A,E5A (chain A)
B,2,W,A,W2A (chain B)
B,3,F,A,F3A (chain B)
```

### 4. Resultados ddG (ddg_results.csv)

```csv
mutation,chain,position,original_aa,ddg,total_score
A2L,A,2,L,2.35,-47.65
A4K,A,4,K,1.87,-48.13
B2W,B,2,W,2.12,-47.88
A3V,A,3,V,0.45,-49.55
A5E,A,5,E,0.78,-49.22
B3F,B,3,F,0.62,-49.38
```

### 5. Hotspots (hotspots.csv)

```csv
mutation,chain,position,original_aa,ddg
A2L,A,2,L,2.35
B2W,B,2,W,2.12
A4K,A,4,K,1.87
```

### 6. Relatório de Análise (analysis_report.txt)

```
================================================================================
ROSETTA ALANINE SCANNING - ANALYSIS REPORT
================================================================================

Structure: example_protein.pdb
Date: 2025-12-02 16:30:00

STRUCTURE INFORMATION
--------------------------------------------------------------------------------
Chain A: 5 residues
Chain B: 4 residues

MUTATION SUMMARY
--------------------------------------------------------------------------------
Total mutations tested: 6
Interface mutations: 3

DDG STATISTICS
--------------------------------------------------------------------------------
Mean ΔΔG: 1.36 kcal/mol
Std ΔΔG: 0.78 kcal/mol
Min ΔΔG: 0.45 kcal/mol
Max ΔΔG: 2.35 kcal/mol

HOTSPOTS (threshold = 1.5 kcal/mol)
--------------------------------------------------------------------------------
Number of hotspots: 3

Rank   Mutation     Chain    Position   ΔΔG
--------------------------------------------------------------------------------
1      A2L          A        2          2.35
2      B2W          B        2          2.12
3      A4K          A        4          1.87

================================================================================
```

### 7. Script PyMOL (visualize_hotspots.pml)

```python
# PyMOL script to visualize alanine scanning hotspots
# Generated from: example_protein.pdb

# Load structure
load /path/to/example_protein.pdb

# Basic visualization
hide everything
show cartoon
color grey80, all
util.cbc

# Highlight hotspots
# Hotspot 1: A2L (ΔΔG = 2.35)
select hot_1, chain A and resi 2
show sticks, hot_1
color red, hot_1
label hot_1, 'A2L\n2.4'

# Hotspot 2: B2W (ΔΔG = 2.12)
select hot_2, chain B and resi 2
show sticks, hot_2
color orange, hot_2
label hot_2, 'B2W\n2.1'

# Hotspot 3: A4K (ΔΔG = 1.87)
select hot_3, chain A and resi 4
show sticks, hot_3
color yellow, hot_3
label hot_3, 'A4K\n1.9'

# Zoom to hotspots
zoom hot_*
```

##  Interpretação dos Resultados

### Valores de ΔΔG:

- **ΔΔG > 2.5 kcal/mol**:  **Hotspot crítico** - Essencial para binding/estabilidade
- **ΔΔG > 2.0 kcal/mol**:  **Hotspot alto** - Muito importante
- **ΔΔG > 1.5 kcal/mol**:  **Hotspot médio** - Contribuição significativa
- **ΔΔG < 1.0 kcal/mol**:  **Não-hotspot** - Contribuição menor

### Exemplo de Análise:

Para o resíduo **A2L** com **ΔΔG = 2.35 kcal/mol**:

- **Interpretação**: Mutação L→A causa perda de 2.35 kcal/mol
- **Significado**: Este resíduo é um hotspot importante
- **Impacto**: Leucina na posição 2 da cadeia A é crítica para a função
- **Recomendação**: Conservar este resíduo, evitar mutações

## Rodando com Rosetta Real

### Pré-requisitos Adicionais:

```bash
# 1. Instale o Rosetta (acadêmico gratuito)
# Baixe de: https://www.rosettacommons.org/software/license-and-download

# 2. Configure a variável de ambiente
export ROSETTA=/path/to/rosetta
export PATH=$ROSETTA/main/source/bin:$PATH

# 3. Verifique a instalação
ls $ROSETTA/main/source/bin/flex_ddg*
```

### Executando o Protocolo Real:

```bash
# 1. Gere as mutações
rosetta-scan scan examples/example_protein.pdb \
    --chains A B \
    --interface-only \
    --output real_run/

# 2. Execute Flex ddG (isso vai demorar!)
rosetta-scan run examples/example_protein.pdb \
    real_run/mutations.txt \
    --nstruct 35 \
    --iterations 3 \
    --interface \
    --output real_run/ddg_results/

# 3. Analise os resultados
rosetta-scan analyze real_run/ddg_results/ \
    --plot \
    --threshold 1.5 \
    --output real_run/analysis.csv
```

### Tempo Estimado:

- **Test run** (nstruct=5): ~5-10 minutos
- **Production run** (nstruct=35): ~30-60 minutos por mutação
- **Para 6 mutações**: ~3-6 horas total

## Visualizações Geradas

O comando `analyze --plot` gera automaticamente:

1. **ddg_distribution.png** - Histograma da distribuição de ΔΔG
2. **hotspot_heatmap.png** - Mapa de calor por posição e cadeia
3. **chain_analysis.png** - Box plots por cadeia
4. **top_hotspots.png** - Gráfico de barras dos top hotspots
5. **position_scan_chain_*.png** - Scan ao longo da sequência

## Casos de Uso Comuns

### Caso 1: Análise de Interface Proteína-Proteína

```bash
rosetta-scan pipeline complex.pdb \
    --chains A B \
    --interface-only \
    --nstruct 35 \
    --output interface_analysis/
```

### Caso 2: Estabilidade de Dobramento

```bash
rosetta-scan pipeline protein.pdb \
    --chains A \
    --nstruct 35 \
    --output stability_analysis/
```

### Caso 3: Região Específica

```bash
rosetta-scan scan protein.pdb \
    -r "A:100-150" \
    --output region_scan/
```

## Troubleshooting

### Erro: "Rosetta not found"

```bash
# Verifique se ROSETTA está configurado
echo $ROSETTA

# Configure se necessário
export ROSETTA=/path/to/rosetta
```

### Erro: "ModuleNotFoundError"

```bash
# Reinstale o pacote
pip install -e .
```

### Erro: Matplotlib backend

```bash
# Use backend não-interativo
export MPLBACKEND=Agg
```

## Próximos Passos

1. Execute o demo: `python3 examples/demo_run.py`
2. Leia os arquivos gerados em `examples/demo_output/`
3. Teste com sua própria estrutura
4. Visualize no PyMOL: `pymol examples/demo_output/visualize_hotspots.pml`
5. Execute com Rosetta para resultados reais

---

**Pronto para começar?** Execute: `python3 examples/demo_run.py`
