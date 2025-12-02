# COMECE AQUI - Execute no Terminal do VSCode

## Setup Rápido (2 minutos)

Abra o terminal integrado do VSCode e execute os comandos abaixo:

### Opção 1: Script Automático (Recomendado)

```bash
# Entre no diretório do projeto
cd /Users/madsonluna/Documents/alaninescanning4eif4e

# Execute o script de setup
./QUICK_START.sh
```

### Opção 2: Comandos Manuais

```bash
# 1. Entre no diretório
cd /Users/madsonluna/Documents/alaninescanning4eif4e

# 2. Crie ambiente virtual
python3 -m venv venv

# 3. Ative o ambiente
source venv/bin/activate

# 4. Instale o pacote
pip install -e .

# 5. Execute o demo
python3 examples/demo_run.py
```

## O Que Você Vai Ver

Após executar, você verá uma saída similar a:

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

[... continua ...]
```

## Arquivos Gerados

Após a execução, verifique os arquivos criados:

```bash
# Liste todos os arquivos gerados
ls -la examples/demo_output/

# Veja o relatório principal
cat examples/demo_output/analysis_report.txt

# Veja os hotspots identificados
cat examples/demo_output/hotspots.csv

# Veja as mutações geradas
cat examples/demo_output/mutations.txt
```

### Estrutura Esperada:

```
examples/demo_output/
--- mutations.txt              # Lista de mutações (texto)
--- mutations_rosetta.txt      # Formato Rosetta
--- mutations.csv              # Formato CSV
--- interface_mutations.txt    # Apenas interface
--- ddg_results.csv            # Resultados completos
--- hotspots.csv               # Hotspots identificados
--- analysis_report.txt        # Relatório completo
--- visualize_hotspots.pml     # Script PyMOL
--- plots/                     # Visualizações
    --- ddg_distribution.png
    --- top_hotspots.png
    --- chain_analysis.png
    --- hotspot_heatmap.png
    --- position_scan_chain_A.png
    --- position_scan_chain_B.png
```

## Visualizar os Plots

```bash
# macOS
open examples/demo_output/plots/

# Linux
xdg-open examples/demo_output/plots/

# Ou abra no VSCode
code examples/demo_output/plots/
```

## Testar os Comandos CLI

```bash
# Ver ajuda geral
rosetta-scan --help

# Ver ajuda do comando scan
rosetta-scan scan --help

# Testar comando scan
rosetta-scan scan examples/example_protein.pdb
    --chains A B
    --output test_scan/

# Ver resultados
cat test_scan/mutations.txt
```

## Próximos Passos

### 1. Explorar a Documentação

```bash
# Ler o README principal
cat README.md

# Ler o tutorial completo
cat DEMO.md

# Ver estrutura do projeto
cat PROJECT_STRUCTURE.md

# Ver exemplos de saída
cat examples/EXPECTED_OUTPUT.md
```

### 2. Testar com Sua Própria Estrutura

```bash
# Substitua 'sua_proteina.pdb' pelo seu arquivo PDB
rosetta-scan scan sua_proteina.pdb --chains A B --interface-only --output meu_scan/
```

### 3. Executar Pipeline Completo (Requer Rosetta)

```bash
# Primeiro configure o Rosetta
export ROSETTA=/caminho/para/rosetta

# Execute pipeline completo
rosetta-scan pipeline sua_proteina.pdb --chains A B --interface-only --nstruct 35 --output meus_resultados/
```

## Troubleshooting

### Erro: "command not found: python3"

```bash
# Tente com python ao invés de python3
python --version
python examples/demo_run.py
```

### Erro: "No module named 'rosetta_scan'"

```bash
# Certifique-se que está no ambiente virtual
source venv/bin/activate

# Reinstale o pacote
pip install -e .
```

### Erro: "Permission denied"

```bash
# Torne o script executável
chmod +x QUICK_START.sh
chmod +x examples/*.sh
```

### Plots não estão sendo gerados

```bash
# Configure backend não-interativo
export MPLBACKEND=Agg

# Execute novamente
python3 examples/demo_run.py
```

## Exemplos de Uso Comuns

### Exemplo 1: Scan Rápido

```bash
# Gerar mutações de uma estrutura
rosetta-scan scan protein.pdb --output mutacoes/
```

### Exemplo 2: Apenas Interface

```bash
# Identificar apenas resíduos de interface
rosetta-scan scan complex.pdb --chains A B --interface-only --interface-cutoff 8.0 --output interface/
```

### Exemplo 3: Região Específica

```bash
# Scan de região específica da sequência
rosetta-scan scan protein.pdb --chains A --range "A:100-200" --output regiao/
```

### Exemplo 4: Pipeline Completo (com Rosetta)

```bash
# ATENÇÃO: Requer Rosetta instalado!
export ROSETTA=/path/to/rosetta

rosetta-scan pipeline protein.pdb --chains A B --interface-only --nstruct 35 --output resultados_completos/
```

## Interpretando os Resultados

### Valores de ΔΔG:

- **ΔΔG > 2.5**:  Hotspot CRÍTICO
- **ΔΔG > 2.0**:  Hotspot ALTO
- **ΔΔG > 1.5**:  Hotspot MÉDIO
- **ΔΔG < 1.5**:  Contribuição menor

### Exemplo de Análise:

```csv
# hotspots.csv
mutation,chain,position,original_aa,ddg
A2L,A,2,L,2.64
B2W,B,2,W,2.13
A4K,A,4,K,1.76
```

**Interpretação:**
- **A2L**: Leucina na posição 2 é CRÍTICA (ΔΔG = 2.64)
- **B2W**: Triptofano na posição 2 é IMPORTANTE (ΔΔG = 2.13)
- **A4K**: Lisina na posição 4 é SIGNIFICATIVA (ΔΔG = 1.76)

## Visualização 3D (PyMOL)

Se você tem PyMOL instalado:

```bash
# Carregar visualização dos hotspots
pymol examples/demo_output/visualize_hotspots.pml
```

Isso abrirá o PyMOL com:
- Estrutura em cartoon
- Hotspots destacados em sticks
- Cores: vermelho (crítico), laranja (alto), amarelo (médio)
- Labels com mutação e ΔΔG

## Documentos Disponíveis

| Arquivo | Descrição |
|---------|-----------|
| [README.md](README.md) | Overview geral e quick start |
| [INSTALL.md](INSTALL.md) | Instalação detalhada e troubleshooting |
| [DEMO.md](DEMO.md) | Tutorial passo a passo completo |
| [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) | Estrutura do código |
| [examples/EXPECTED_OUTPUT.md](examples/EXPECTED_OUTPUT.md) | Exemplos de saída |

## Precisa de Ajuda?

1. **Verifique a documentação**: Leia os arquivos .md
2. **Execute o demo**: `python3 examples/demo_run.py`
3. **Teste os exemplos**: Veja `examples/`
4. **Abra uma issue**: GitHub issues

## Checklist Rápido

Marque conforme você completa:

- [ ] Clonou/baixou o repositório
- [ ] Criou ambiente virtual (`python3 -m venv venv`)
- [ ] Ativou ambiente (`source venv/bin/activate`)
- [ ] Instalou dependências (`pip install -e .`)
- [ ] Executou o demo (`python3 examples/demo_run.py`)
- [ ] Explorou os arquivos gerados (`ls examples/demo_output/`)
- [ ] Visualizou os plots (`open examples/demo_output/plots/`)
- [ ] Testou comando CLI (`rosetta-scan --help`)
- [ ] Leu a documentação (`cat README.md`)
- [ ] Pronto para usar! 

---

## Comando Único para Setup Completo

Se quiser fazer tudo de uma vez:

```bash
cd /Users/madsonluna/Documents/alaninescanning4eif4e &&
python3 -m venv venv &&
source venv/bin/activate &&
pip install -e . &&
python3 examples/demo_run.py
```

**Tempo estimado:** 2-3 minutos

---

**Pronto para começar?** Execute no terminal: `./QUICK_START.sh`
