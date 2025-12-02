# Gerador de Heatmap ΔΔG

Script standalone para gerar visualização em heatmap dos valores de ΔΔG com normalização e colormap seismic.

## Como Usar

### Execução Simples

```bash
cd /Users/madsonluna/Documents/alaninescanning4eif4e/gamma-example
python3 generate_ddg_heatmap.py
```

### Ou como executável

```bash
./generate_ddg_heatmap.py
```

## O Que o Script Faz

1. Carrega os resultados de ΔΔG do arquivo `analysis_results/ddg_results.csv`
2. Normaliza os valores de ΔΔG para o intervalo [0, 1]
3. Gera um heatmap usando o colormap "seismic" do matplotlib
4. Anota cada célula com o valor real de ΔΔG
5. Salva a visualização em `analysis_results/plots/ddg_heatmap_normalized.png`
6. Exibe estatísticas dos resultados

## Características do Heatmap

### Colormap Seismic
- Azul: Valores baixos de ΔΔG (menor impacto)
- Branco: Valores intermediários
- Vermelho: Valores altos de ΔΔG (maior impacto)

### Normalização
- Valores de ΔΔG são normalizados para [0, 1]
- A barra de cores mostra a correspondência com valores reais de ΔΔG
- Permite comparação visual mesmo com ranges diferentes

### Anotações
- Cada célula mostra o valor real de ΔΔG
- Cor do texto ajustada automaticamente para contraste
- Labels mostram: Aminoácido original + Posição

## Requisitos

O script requer que a análise de alanine scanning já tenha sido executada e que exista o arquivo:
```
analysis_results/ddg_results.csv
```

## Saída

### Arquivo Gerado
```
analysis_results/plots/ddg_heatmap_normalized.png
```

### Estatísticas Exibidas
- Total de mutações
- Média, mediana, desvio padrão de ΔΔG
- Valores mínimo e máximo
- Top 5 mutações por impacto
- Bottom 5 mutações por impacto

## Exemplo de Saída no Terminal

```
============================================================
ΔΔG HEATMAP GENERATOR
Normalized colormap: seismic
============================================================

Loading data from: analysis_results/ddg_results.csv
Loaded 9 mutations from analysis_results/ddg_results.csv

============================================================
ΔΔG STATISTICS
============================================================

Total mutations: 9
Mean ΔΔG: 0.600 kcal/mol
Median ΔΔG: 0.619 kcal/mol
Std ΔΔG: 0.318 kcal/mol
Min ΔΔG: 0.180 kcal/mol
Max ΔΔG: 1.284 kcal/mol

Top 5 mutations by ΔΔG:
------------------------------------------------------------
  A7F      | Chain A | F7   | ΔΔG =  1.284 kcal/mol
  A10R     | Chain A | R10  | ΔΔG =  0.715 kcal/mol
  A4C      | Chain A | C4   | ΔΔG =  0.713 kcal/mol
  A3R      | Chain A | R3   | ΔΔG =  0.684 kcal/mol
  A9R      | Chain A | R9   | ΔΔG =  0.619 kcal/mol

Bottom 5 mutations by ΔΔG:
------------------------------------------------------------
  A11C     | Chain A | C11  | ΔΔG =  0.180 kcal/mol
  A12F     | Chain A | F12  | ΔΔG =  0.327 kcal/mol
  A8R      | Chain A | R8   | ΔΔG =  0.354 kcal/mol
  A5R      | Chain A | R5   | ΔΔG =  0.503 kcal/mol
  A9R      | Chain A | R9   | ΔΔG =  0.619 kcal/mol

============================================================

Generating heatmap...
Heatmap saved to: analysis_results/plots/ddg_heatmap_normalized.png

Visualization complete!
Output saved to: analysis_results/plots/ddg_heatmap_normalized.png

============================================================
```

## Customização

Para modificar o script, você pode ajustar:

### Colormap
Linha 107: `cmap='seismic'`

Opções alternativas:
- `'coolwarm'` - Azul para vermelho
- `'RdYlGn_r'` - Vermelho-amarelo-verde invertido
- `'viridis'` - Amarelo-verde-azul
- `'plasma'` - Roxo-laranja

### Tamanho da Figura
Linha 99: `figsize=(14, max(4, len(row_labels) * 1.5))`

### Resolução DPI
Linha 166: `dpi=300`

### Threshold para Cor do Texto
Linha 141: `text_color = 'white' if matrix_normalized[i, j] < 0.5 else 'black'`

## Estrutura do Código

```python
Funções principais:
- load_ddg_results()          # Carrega CSV
- normalize_ddg()             # Normaliza valores
- create_heatmap_matrix()     # Cria matriz para heatmap
- generate_heatmap()          # Gera visualização
- generate_statistics_table() # Exibe estatísticas
- main()                      # Função principal
```

## Troubleshooting

### Erro: Arquivo não encontrado
```
Error: Input file not found: analysis_results/ddg_results.csv
```
Solução: Execute primeiro a análise de alanine scanning com `python3 run_gamma_analysis.py`

### Erro: Módulo não encontrado
```
ModuleNotFoundError: No module named 'matplotlib'
```
Solução: Instale as dependências com `pip install -e .` no diretório raiz do projeto

## Integração com Pipeline

Este script pode ser integrado ao pipeline principal ou executado separadamente após a análise.

### Uso Standalone
```bash
python3 generate_ddg_heatmap.py
```

### Uso no Pipeline Python
```python
from pathlib import Path
import subprocess

# Executar análise
subprocess.run(['python3', 'run_gamma_analysis.py'])

# Gerar heatmap
subprocess.run(['python3', 'generate_ddg_heatmap.py'])
```

## Análise Visual do Heatmap

### Interpretação das Cores

**Vermelho Intenso (ΔΔG alto):**
- Resíduos críticos para estabilidade
- Mutação para alanina causa grande perda de energia
- Candidatos a hotspots

**Azul Intenso (ΔΔG baixo):**
- Resíduos com menor contribuição energética
- Mutação para alanina tem impacto mínimo
- Posições tolerantes a mutação

**Branco/Neutro (ΔΔG intermediário):**
- Contribuição moderada
- Importância média para estrutura/função

### Para o Peptídeo Gamma

No heatmap gerado, você verá:
- Phe7 em vermelho (ΔΔG = 1.28 kcal/mol) - Resíduo mais importante
- Cys11 em azul (ΔΔG = 0.18 kcal/mol) - Menor impacto
- Gradient de cores mostrando a importância relativa de cada posição

## Formato do Arquivo de Entrada

O script espera um CSV com as seguintes colunas:
```
mutation,chain,position,original_aa,ddg,total_score
A7F,A,7,F,1.2844406193379763,-48.51844286972222
A10R,A,10,R,0.7154929429483751,-48.84362653757732
...
```

## Arquivos Relacionados

- `run_gamma_analysis.py` - Script principal de análise
- `INSTRUCTIONS.md` - Instruções gerais
- `GAMMA_ANALYSIS_SUMMARY.md` - Resumo executivo da análise
- `analysis_results/ANALYSIS_REPORT.md` - Relatório técnico completo
