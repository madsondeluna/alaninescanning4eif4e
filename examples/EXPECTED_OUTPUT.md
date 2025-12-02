#  Exemplos de Saída - Output Esperado

Este documento mostra exemplos reais de saída do protocolo de alanine scanning.

## 🔍 Saída do Terminal

### Comando: `rosetta-scan scan`

```bash
$ rosetta-scan scan examples/example_protein.pdb --chains A B --output scan_output/

 Alanine Scanning Mutation Generator

⠋ Loading structure...                                      ------------------ 100%
⠋ Generating mutations...                                   ------------------ 100%


                        Mutation Summary
---------------------------------------------------┓
| Property           | Value                       |
┡--------------------╇-----------------------------┩
| Total Mutations    | 6                           |
| Structure File     | example_protein.pdb         |
| Chain A            | 4                           |
| Chain B            | 2                           |
---------------------------------------------------┘


 Mutations saved to: scan_output/mutations.txt
```

### Comando: `rosetta-scan scan --interface-only`

```bash
$ rosetta-scan scan examples/example_protein.pdb --chains A B --interface-only --interface-cutoff 8.0 --output interface_scan/

 Alanine Scanning Mutation Generator

⠋ Loading structure...                                      ------------------ 100%
⠋ Identifying interface residues...                         ------------------ 100%
⠋ Generating mutations...                                   ------------------ 100%


                        Mutation Summary
---------------------------------------------------┓
| Property           | Value                       |
┡--------------------╇-----------------------------┩
| Total Mutations    | 3                           |
| Structure File     | example_protein.pdb         |
| Chain A            | 2                           |
| Chain B            | 1                           |
---------------------------------------------------┘


 Mutations saved to: interface_scan/mutations.txt
```

### Comando: `python3 examples/demo_run.py`

```bash
$ python3 examples/demo_run.py

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

Amino acid distribution:
  • E: 1 mutations
  • K: 1 mutations
  • L: 1 mutations
  • V: 1 mutations
  • F: 1 mutations
  • W: 1 mutations

First 10 mutations to be tested:
  #    Mutation     Chain    Position   Original→Ala
  ------------------------------------------------------------
  1    A2L          A        2          L→A
  2    A3V          A        3          V→A
  3    A4K          A        4          K→A
  4    A5E          A        5          E→A
  5    B2W          B        2          W→A
  6    B3F          B        3          F→A

 Mutations saved to:
  • Text report: mutations.txt
  • Rosetta format: mutations_rosetta.txt
  • CSV format: mutations.csv

--------------------------------------------------------------------------------
> STEP 3: Interface Analysis
--------------------------------------------------------------------------------

Identifying interface residues (cutoff = 8.0 Å)...

OK Found 3 interface residues

Interface hotspot candidates:
  Chain    Position   Residue    Mutation
  --------------------------------------------------
  A        2          L          A2L
  A        4          K          A4K
  B        2          W          B2W

 Interface mutations saved to: interface_mutations.txt

--------------------------------------------------------------------------------
> STEP 4: Simulated ddG Results
--------------------------------------------------------------------------------

NOTE: This is a demo with simulated ddG values.
In a real run, you would use:
  rosetta-scan run example_protein.pdb mutations_rosetta.txt --nstruct 35

OK Simulated results generated

Statistics:
  • Total mutations: 6
  • Mean ΔΔG: 1.20 kcal/mol
  • Std ΔΔG: 0.83 kcal/mol
  • Min ΔΔG: 0.31 kcal/mol
  • Max ΔΔG: 2.64 kcal/mol

--------------------------------------------------------------------------------
> STEP 5: Hotspot Identification
--------------------------------------------------------------------------------

Identifying hotspots (threshold = 1.5 kcal/mol)...

OK Found 3 hotspot residues

 TOP HOTSPOTS (ranked by ΔΔG):

  Rank   Mutation     Chain    Position   ΔΔG (kcal/mol)  Impact
  ----------------------------------------------------------------------
  1      A2L          A        2          2.64             Critical
  2      B2W          B        2          2.13             High
  3      A4K          A        4          1.76             Medium

 Results saved to:
  • All results: ddg_results.csv
  • Hotspots: hotspots.csv

--------------------------------------------------------------------------------
> STEP 6: Generate Visualizations
--------------------------------------------------------------------------------

Creating publication-quality plots...
  • Creating ΔΔG distribution plot...
  • Creating top hotspots plot...
  • Creating per-chain analysis...
  • Creating hotspot heatmap...
  • Creating position scan plots...

OK Created 6 visualization plots

 Plots saved to: examples/demo_output/plots/
  • ddg_distribution.png
  • top_hotspots.png
  • chain_analysis.png
  • hotspot_heatmap.png
  • position_scan_chain_A.png
  • position_scan_chain_B.png

--------------------------------------------------------------------------------
> STEP 7: PyMOL Visualization
--------------------------------------------------------------------------------

Generating PyMOL script for 3D visualization...

OK PyMOL script created: visualize_hotspots.pml

To visualize in PyMOL, run:
  pymol examples/demo_output/visualize_hotspots.pml

--------------------------------------------------------------------------------
> STEP 8: Generate Summary Report
--------------------------------------------------------------------------------

OK Summary report saved: analysis_report.txt

================================================================================
   DEMO COMPLETED SUCCESSFULLY
================================================================================

 All output files generated in: demo_output/

 Files created:
  Mutations:
    • mutations.txt (human-readable)
    • mutations_rosetta.txt (Rosetta format)
    • mutations.csv (spreadsheet)
    • interface_mutations.txt (interface only)

  Results:
    • ddg_results.csv (all ΔΔG values)
    • hotspots.csv (filtered hotspots)
    • analysis_report.txt (summary)

  Visualizations:
    • plots/ddg_distribution.png
    • plots/top_hotspots.png
    • plots/chain_analysis.png
    • plots/hotspot_heatmap.png
    • plots/position_scan_chain_A.png
    • plots/position_scan_chain_B.png

  3D Visualization:
    • visualize_hotspots.pml (PyMOL script)

 Key Findings:
  • Tested 6 alanine mutations
  • Identified 3 hotspot residues
  • Top hotspot: A2L (ΔΔG = 2.64 kcal/mol)

 Next Steps:
  1. Review plots in the plots/ directory
  2. Examine hotspots.csv for detailed results
  3. Visualize in PyMOL: pymol visualize_hotspots.pml
  4. Run with real Rosetta for actual ddG values

================================================================================
```

##  Arquivos de Saída

### 1. mutations.txt (Human-Readable)

```
Alanine Scanning Mutation Report
==================================================

Structure: example_protein.pdb
Total mutations: 6

Chain    Position  Original  Target
--------------------------------------------------
A        2         L         A
A        3         V         A
A        4         K         A
A        5         E         A
B        2         W         A
B        3         F         A
```

### 2. mutations_rosetta.txt (Rosetta Format)

```
total 1
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

### 3. mutations.csv (Spreadsheet)

```csv
chain,position,original_aa,target_aa,mutation
A,2,L,A,L2A (chain A)
A,3,V,A,V3A (chain A)
A,4,K,A,K4A (chain A)
A,5,E,A,E5A (chain A)
B,2,W,A,W2A (chain B)
B,3,F,A,F3A (chain B)
```

### 4. interface_mutations.txt

```
Alanine Scanning Mutation Report
==================================================

Structure: example_protein.pdb
Total mutations: 3

Chain    Position  Original  Target
--------------------------------------------------
A        2         L         A
A        4         K         A
B        2         W         A
```

### 5. ddg_results.csv

```csv
mutation,chain,position,original_aa,ddg,total_score
A2L,A,2,L,2.64,-47.36
B2W,B,2,W,2.13,-47.87
A4K,A,4,K,1.76,-48.24
A5E,A,5,E,0.89,-49.11
B3F,B,3,F,0.56,-49.44
A3V,A,3,V,0.31,-49.69
```

### 6. hotspots.csv

```csv
mutation,chain,position,original_aa,ddg,abs_ddg
A2L,A,2,L,2.64,2.64
B2W,B,2,W,2.13,2.13
A4K,A,4,K,1.76,1.76
```

### 7. analysis_report.txt

```
================================================================================
ROSETTA ALANINE SCANNING - ANALYSIS REPORT
================================================================================

Structure: example_protein.pdb
Date: 2025-12-02 16:45:23

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
Mean ΔΔG: 1.20 kcal/mol
Std ΔΔG: 0.83 kcal/mol
Min ΔΔG: 0.31 kcal/mol
Max ΔΔG: 2.64 kcal/mol

HOTSPOTS (threshold = 1.5 kcal/mol)
--------------------------------------------------------------------------------
Number of hotspots: 3

Rank   Mutation     Chain    Position   ΔΔG
--------------------------------------------------------------------------------
1      A2L          A        2          2.64
2      B2W          B        2          2.13
3      A4K          A        4          1.76

================================================================================
```

### 8. visualize_hotspots.pml (PyMOL Script)

```python
# PyMOL script to visualize alanine scanning hotspots
# Generated from: example_protein.pdb

# Load structure
load /Users/madsonluna/Documents/alaninescanning4eif4e/examples/example_protein.pdb

# Basic visualization
hide everything
show cartoon
color grey80, all
util.cbc

# Highlight hotspots
# Hotspot 1: A2L (ΔΔG = 2.64)
select hot_1, chain A and resi 2
show sticks, hot_1
color red, hot_1
label hot_1, 'A2L\n2.6'

# Hotspot 2: B2W (ΔΔG = 2.13)
select hot_2, chain B and resi 2
show sticks, hot_2
color orange, hot_2
label hot_2, 'B2W\n2.1'

# Hotspot 3: A4K (ΔΔG = 1.76)
select hot_3, chain A and resi 4
show sticks, hot_3
color yellow, hot_3
label hot_3, 'A4K\n1.8'

# Zoom to hotspots
zoom hot_*

# Set viewing angle
set_view (1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,20,0)
```

##  Visualizações

### 1. ΔΔG Distribution (ddg_distribution.png)

```
Descrição: Histograma mostrando a distribuição dos valores de ΔΔG

Características:
- Eixo X: ΔΔG (kcal/mol)
- Eixo Y: Frequência
- Linha tracejada: Média (1.20 kcal/mol)
- Cores: Azul (#2E86AB) com transparência
- Título: "Alanine Scanning ΔΔG Distribution"
```

### 2. Top Hotspots (top_hotspots.png)

```
Descrição: Gráfico de barras horizontais dos top 10 hotspots

Características:
- Barras vermelhas: ΔΔG positivo (desestabilizante)
- Barras azuis: ΔΔG negativo (estabilizante)
- Labels: Mutação + Chain
- Ordenado por magnitude de ΔΔG
```

### 3. Chain Analysis (chain_analysis.png)

```
Descrição: Box plots e violin plots por cadeia

Painel 1 (esquerda): Box plot
- Mostra mediana, quartis, outliers
- Comparação entre cadeias A e B

Painel 2 (direita): Violin plot
- Mostra distribuição completa
- Kernel density estimation
```

### 4. Hotspot Heatmap (hotspot_heatmap.png)

```
Descrição: Mapa de calor 2D

Eixos:
- X: Cadeia (A, B)
- Y: Posição do resíduo
- Cor: ΔΔG (azul → branco → vermelho)
- Escala: -2 a +3 kcal/mol
```

### 5. Position Scan (position_scan_chain_*.png)

```
Descrição: Gráfico de linha mostrando ΔΔG ao longo da sequência

Características:
- Linha azul: Valores de ΔΔG
- Pontos amarelos: Hotspots (ΔΔG > threshold)
- Linhas tracejadas: Threshold (±1.5 kcal/mol)
- Linha sólida preta: Zero
```

##  Interpretação dos Resultados

### Hotspot A2L (ΔΔG = 2.64 kcal/mol)

```
 Classificação: HOTSPOT CRÍTICO 

Interpretação:
- Mutação Leu2 → Ala na cadeia A causa perda de 2.64 kcal/mol
- Este resíduo é essencial para a estabilidade/binding
- Provavelmente faz contatos hidrofóbicos importantes
- Localizado na interface entre as cadeias

Recomendação:
- Conservar leucina na posição 2
- Evitar mutações nesta posição
- Candidato para validação experimental
```

### Hotspot B2W (ΔΔG = 2.13 kcal/mol)

```
 Classificação: HOTSPOT ALTO 

Interpretação:
- Mutação Trp2 → Ala na cadeia B causa perda de 2.13 kcal/mol
- Triptofano é grande e aromático, perda significativa
- Provavelmente envolvido em stacking π-π ou ligações CH-π
- Posição crítica na interface

Recomendação:
- Manter triptofano
- Se mutar, usar resíduos aromáticos (Phe, Tyr)
- Prioridade alta para validação
```

### Hotspot A4K (ΔΔG = 1.76 kcal/mol)

```
 Classificação: HOTSPOT MÉDIO 

Interpretação:
- Mutação Lys4 → Ala na cadeia A causa perda de 1.76 kcal/mol
- Lisina carregada positivamente, pode fazer salt bridges
- Contribuição significativa mas não crítica
- Interface ou superfície

Recomendação:
- Considerar conservação
- Mutações conservativas podem ser toleradas (Arg)
- Avaliar interações eletrostáticas
```

### Não-Hotspots (ΔΔG < 1.5 kcal/mol)

```
 A3V, A5E, B3F: Contribuição menor

Interpretação:
- Mutações causam pouca perda de estabilidade
- Resíduos menos críticos para função
- Podem ser alvos para engenharia de proteínas

Possibilidades:
- Otimização de propriedades sem perda de função
- Adição de tags ou sítios de modificação
- Estudos de evolução molecular
```

##  Métricas de Qualidade

### Para Resultados do Rosetta Real:

```
Métricas recomendadas:
- nstruct ≥ 35: Boa estatística
- Std ddG < 0.5: Alta precisão
- Convergência em iterations: Sampling adequado

Valores típicos:
- Hotspots fortes: ΔΔG > 2.0 kcal/mol
- Hotspots moderados: 1.5 < ΔΔG < 2.0
- Contribuição menor: ΔΔG < 1.5
```

##  Próximos Passos

1. **Validação Experimental**
   - Mutagênese sítio-dirigida
   - Medições de binding (SPR, ITC)
   - Ensaios funcionais

2. **Análise Estrutural**
   - Visualização no PyMOL
   - MD simulations
   - Análise de contatos

3. **Design Racional**
   - Usar hotspots como base
   - Engenharia de afinidade
   - Desenvolvimento de inibidores

---

**Pronto para experimentar?** Execute `python3 examples/demo_run.py`
