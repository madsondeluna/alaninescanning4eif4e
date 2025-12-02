#  Estrutura do Projeto

Visão completa da organização do framework Rosetta Alanine Scanning.

## 🌲 Árvore de Diretórios

```
alaninescanning4eif4e/
|
---  README.md                      # Documentação principal
---  INSTALL.md                     # Guia de instalação detalhado
---  DEMO.md                        # Tutorial de demonstração
---  PROJECT_STRUCTURE.md           # Este arquivo
---  LICENSE                        # Licença MIT
---  .gitignore                     # Arquivos ignorados pelo Git
|
--- 📦 setup.py                       # Configuração do pacote Python
--- 📦 requirements.txt               # Dependências Python
|
--- 📂 src/rosetta_scan/              #  Código-fonte principal
|   --- __init__.py
|   --- cli.py                        # Interface CLI (Click)
|   |
|   --- 📂 protocols/                 # Protocolos do Rosetta
|   |   --- __init__.py
|   |   --- flex_ddg.py              # Protocolo Flex ddG
|   |   --- alanine_scanner.py       # Gerador de mutações
|   |
|   --- 📂 analysis/                  # Análise e visualização
|   |   --- __init__.py
|   |   --- parser.py                # Parser de resultados
|   |   --- visualizer.py            # Gerador de plots
|   |
|   --- 📂 utils/                     # Utilidades gerais
|       --- __init__.py
|
--- 📂 config/                        # Configurações
|   --- example_config.yaml          # Configuração exemplo
|
--- 📂 examples/                      #  Exemplos e demos
|   --- example_protein.pdb          # Estrutura PDB exemplo
|   --- example_workflow.py          # Workflow Python exemplo
|   --- demo_run.py                  # Demo interativo completo
|   --- quick_test.sh                # Script de teste rápido
|   --- EXPECTED_OUTPUT.md           # Exemplos de saída
|   |
|   --- 📂 demo_output/              # Saída do demo (gerado)
|       --- mutations.txt
|       --- mutations_rosetta.txt
|       --- mutations.csv
|       --- interface_mutations.txt
|       --- ddg_results.csv
|       --- hotspots.csv
|       --- analysis_report.txt
|       --- visualize_hotspots.pml
|       --- 📂 plots/
|           --- ddg_distribution.png
|           --- top_hotspots.png
|           --- chain_analysis.png
|           --- hotspot_heatmap.png
|           --- position_scan_chain_A.png
|           --- position_scan_chain_B.png
|
--- 📂 tests/                         # Testes unitários (futuro)
    --- __init__.py
```

##  Componentes Principais

### 1.  Core Protocols (`src/rosetta_scan/protocols/`)

#### `flex_ddg.py`
```python
Classes:
- FlexDdGConfig     # Configuração do protocolo
- FlexDdGProtocol   # Executor do Flex ddG

Principais métodos:
- run_flex_ddg()              # Executa cálculos ddG
- prepare_input_structure()   # Prepara estrutura
- parse_results()             # Parseia resultados
- generate_mutation_file()    # Gera arquivo de mutações
```

#### `alanine_scanner.py`
```python
Classes:
- MutationSite      # Representa uma mutação
- AlanineScan       # Gerador de mutações

Principais métodos:
- generate_mutations()              # Gera mutações sistemáticas
- _identify_interface_residues()    # Encontra interface
- filter_by_solvent_accessibility() # Filtra por SASA
- get_rosetta_mutation_list()       # Lista para Rosetta
- save_mutation_report()            # Salva relatório
```

### 2.  Analysis & Visualization (`src/rosetta_scan/analysis/`)

#### `parser.py`
```python
Classes:
- ResultParser      # Parser de resultados Rosetta

Principais métodos:
- parse_results()               # Parse score files
- identify_hotspots()           # Identifica hotspots
- get_per_chain_statistics()    # Stats por cadeia
- export_pymol_script()         # Script PyMOL
- export_summary_report()       # Relatório texto
```

#### `visualizer.py`
```python
Classes:
- ResultVisualizer  # Gerador de visualizações

Principais métodos:
- plot_ddg_distribution()       # Histograma
- plot_hotspot_heatmap()        # Heatmap 2D
- plot_per_chain_analysis()     # Box/violin plots
- plot_top_hotspots()           # Top N hotspots
- plot_position_scan()          # Scan por posição
- create_dashboard()            # Todas visualizações
```

### 3.  CLI Interface (`src/rosetta_scan/cli.py`)

```bash
Comandos disponíveis:

rosetta-scan scan           # Gerar mutações
rosetta-scan run            # Executar Flex ddG
rosetta-scan analyze        # Analisar resultados
rosetta-scan pipeline       # Pipeline completo
rosetta-scan init-config    # Gerar config exemplo
```

##  Tipos de Arquivos

### Arquivos de Entrada

```
Input Files:
--- *.pdb               # Estrutura PDB
--- *.yaml              # Configuração
--- mutations*.txt      # Lista de mutações (Rosetta format)
```

### Arquivos de Saída

```
Output Files:
---  Mutações
|   --- mutations.txt              # Formato texto
|   --- mutations_rosetta.txt      # Formato Rosetta
|   --- mutations.csv              # Formato CSV
|   --- interface_mutations.txt    # Apenas interface
|
---  Resultados ddG
|   --- ddg_results.csv            # Todos resultados
|   --- hotspots.csv               # Apenas hotspots
|   --- *.sc                       # Score files Rosetta
|
---  Relatórios
|   --- analysis_report.txt        # Relatório completo
|   --- rosetta.log                # Log do Rosetta
|
---  Visualizações
|   --- ddg_distribution.png       # Distribuição
|   --- hotspot_heatmap.png        # Heatmap
|   --- chain_analysis.png         # Análise por cadeia
|   --- top_hotspots.png           # Top hotspots
|   --- position_scan_*.png        # Scan por posição
|
---  Scripts 3D
    --- visualize_hotspots.pml     # PyMOL script
```

## 🔄 Fluxo de Dados

```
------------------------------------------------------------------┐
|                        INPUT                                    |
|                     protein.pdb                                 |
------------------------------------------------------------------┘
                           |
                           ▼
------------------------------------------------------------------┐
|                    STEP 1: SCAN                                 |
|                  AlanineScan                                    |
|  • Load PDB structure (BioPython)                              |
|  • Identify scannable residues                                 |
|  • Generate mutation list                                      |
|  • Optional: identify interface                                |
------------------------------------------------------------------┘
                           |
                           ▼
                    mutations.txt
                           |
                           ▼
------------------------------------------------------------------┐
|                  STEP 2: FLEX DDG                               |
|                 FlexDdGProtocol                                 |
|  • Prepare structure (relax)                                   |
|  • Run Flex ddG calculations (Rosetta)                         |
|  • For each mutation: backrub + repack + score                |
|  • Generate *.sc score files                                   |
------------------------------------------------------------------┘
                           |
                           ▼
                   ddg_results/*.sc
                           |
                           ▼
------------------------------------------------------------------┐
|                  STEP 3: PARSE                                  |
|                   ResultParser                                  |
|  • Parse Rosetta score files                                   |
|  • Extract ddG values                                          |
|  • Identify hotspots                                           |
|  • Calculate statistics                                        |
------------------------------------------------------------------┘
                           |
                           ▼
              results.csv + hotspots.csv
                           |
                           ▼
------------------------------------------------------------------┐
|               STEP 4: VISUALIZE                                 |
|                ResultVisualizer                                 |
|  • Create distribution plots                                   |
|  • Generate heatmaps                                           |
|  • Per-chain analysis                                          |
|  • Position scans                                              |
------------------------------------------------------------------┘
                           |
                           ▼
------------------------------------------------------------------┐
|                       OUTPUT                                    |
|  • CSV files (data)                                            |
|  • PNG plots (visualizations)                                  |
|  • TXT reports (summaries)                                     |
|  • PML scripts (PyMOL)                                         |
------------------------------------------------------------------┘
```

## 🔌 Dependências Externas

### Python Packages

```
Core:
--- click           # CLI framework
--- biopython       # PDB parsing
--- pandas          # Data manipulation
--- numpy           # Numerical operations

Visualization:
--- matplotlib      # Plotting
--- seaborn         # Statistical plots

Configuration:
--- pyyaml          # YAML parsing
--- rich            # Terminal formatting
```

### External Software

```
Required:
--- Rosetta         # Flex ddG calculations

Optional:
--- PyMOL           # 3D visualization
--- PyRosetta       # Python Rosetta bindings
```

##  Pontos de Entrada

### CLI (Command Line)

```bash
# Instalado como comando do sistema
$ rosetta-scan --help

# Acessa src/rosetta_scan/cli.py:main()
```

### Python API

```python
# Importação programática
from rosetta_scan import FlexDdGProtocol, AlanineScan

# Acessa src/rosetta_scan/__init__.py
```

### Demo Script

```bash
# Script de demonstração
$ python3 examples/demo_run.py

# Acessa examples/demo_run.py:main()
```

##  Formato de Dados

### Estrutura DataFrame Principal

```python
# results_df
columns = [
    'mutation',      # str:  "A2L"
    'chain',         # str:  "A"
    'position',      # int:  2
    'original_aa',   # str:  "L"
    'ddg',           # float: 2.35
    'total_score',   # float: -47.65
]
```

### Configuração YAML

```yaml
# config.yaml
nstruct: 35
iterations: 3
repack_radius: 8.0
use_backrub: true
interface_ddg: false
num_processors: 1
memory_gb: 4
rosetta_path: /path/to/rosetta
rosetta_database: /path/to/rosetta/main/database
```

##  Customização

### Adicionar Novo Protocolo

```python
# 1. Criar novo arquivo em src/rosetta_scan/protocols/
# src/rosetta_scan/protocols/my_protocol.py

class MyProtocol:
    def __init__(self, config):
        self.config = config

    def run(self, pdb_path):
        # Implementação
        pass

# 2. Adicionar ao __init__.py
# src/rosetta_scan/__init__.py
from .protocols.my_protocol import MyProtocol
__all__.append('MyProtocol')

# 3. Adicionar comando CLI (opcional)
# src/rosetta_scan/cli.py
@main.command()
def my_command():
    protocol = MyProtocol(config)
    protocol.run(pdb_path)
```

### Adicionar Nova Visualização

```python
# 1. Adicionar método em src/rosetta_scan/analysis/visualizer.py

class ResultVisualizer:
    def plot_my_visualization(self, output_path):
        # Implementação usando matplotlib/seaborn
        pass

# 2. Usar no CLI ou API
visualizer = ResultVisualizer(results_df)
visualizer.plot_my_visualization('output.png')
```

## 🧪 Testing

```python
# Estrutura para testes (futuro)
tests/
--- __init__.py
--- test_alanine_scanner.py
--- test_flex_ddg.py
--- test_parser.py
--- test_visualizer.py
--- test_cli.py

# Executar testes
$ pytest tests/
```

##  Documentação

```
Arquivos de documentação:
--- README.md              # Overview e quick start
--- INSTALL.md             # Instalação detalhada
--- DEMO.md                # Tutorial completo
--- PROJECT_STRUCTURE.md   # Este arquivo
--- examples/
    --- EXPECTED_OUTPUT.md # Exemplos de saída
```

##  Comandos Úteis

```bash
# Desenvolvimento
pip install -e .                    # Instalar modo dev
pip install -e ".[dev]"             # Com deps de dev
python -m pytest                    # Rodar testes

# Uso
rosetta-scan --help                 # Ver comandos
rosetta-scan scan --help            # Help específico
python examples/demo_run.py         # Demo completo

# Manutenção
pip list --outdated                 # Verificar updates
pip install --upgrade -r requirements.txt
```

##  Estatísticas do Código

```
Linhas de código (aproximado):
--- flex_ddg.py:         ~350 linhas
--- alanine_scanner.py:  ~400 linhas
--- parser.py:           ~350 linhas
--- visualizer.py:       ~350 linhas
--- cli.py:              ~400 linhas
--- demo_run.py:         ~450 linhas
------------------------------------
Total:                   ~2300 linhas

Arquivos Python:         14
Arquivos Markdown:       5
Arquivos Config:         2
Arquivos Exemplo:        4
```

---

**Navegação Rápida:**
- [🏠 README](README.md)
- [📦 Installation](INSTALL.md)
- [ Demo](DEMO.md)
- [ Expected Output](examples/EXPECTED_OUTPUT.md)
