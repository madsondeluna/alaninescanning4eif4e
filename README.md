# Alanine Scanning com PyRosetta - Protocolo Flex ddG

Framework Python para análise computacional de alanine scanning usando PyRosetta e o protocolo Flex ddG.

## Descrição

Este framework implementa o protocolo Flex ddG para calcular mudanças na energia livre (ΔΔG) ao introduzir mutações pontuais para alanina em proteínas. O método considera flexibilidade conformacional da cadeia principal através de movimentos backrub, proporcionando estimativas mais acuradas que métodos de estrutura fixa.

## Fundamentos Teóricos

### Função de Energia REF2015

O Rosetta utiliza a função de energia REF2015 (Alford et al., 2017), que combina termos físicos e estatísticos:

```
E_total = Σ w_i × E_i
```

**Componentes principais:**

1. `fa_atr` - van der Waals atrativo (Lennard-Jones)
   - Peso: ~1.0
   - Captura interações hidrofóbicas e empacotamento

2. `fa_rep` - van der Waals repulsivo
   - Peso: ~0.55
   - Previne overlaps estéricos

3. `fa_sol` - Solvatação (modelo Lazaridis-Karplus)
   - Peso: ~1.0
   - Efeito hidrofóbico e dessolvatação polar

4. `fa_elec` - Eletrostático (coulômbico)
   - Peso: ~1.0
   - Interações carga-carga, pontes salinas

5. `hbond_sc/hbond_bb` - Ligações de hidrogênio
   - Peso: ~1.0-1.5
   - Termos direcionais para cadeia lateral e principal

6. Termos entrópicos:
   - `pro_close`: Fechamento de prolina
   - `dslf_fa53`: Pontes dissulfeto
   - `rama_prepro`: Mapa de Ramachandran

**Unidade de energia:**
- Rosetta Energy Units (REU)
- Aproximadamente 1 REU ≈ 1 kcal/mol (métrica relativa)

### Protocolo Flex ddG

O protocolo Flex ddG (Barlow et al., 2018) é um método ensemble-based para calcular ΔΔG:

**Workflow:**

1. **Minimização da estrutura wild-type**
   - Relaxamento com packing e minimização
   - Gera estrutura de referência otimizada

2. **Backrub ensemble generation**
   - Aplica movimentos "backrub" na cadeia principal
   - Simula flexibilidade local em torno da mutação
   - Gera diversidade conformacional realista

3. **Mutação e repacking**
   - Introduz mutação (ex: Trp → Ala)
   - Reempacota cadeias laterais na vizinhança (raio ~8 Å)

4. **Cálculo de ΔΔG**
   ```
   ΔΔG = E_mutante - E_wild-type
   ```
   - ΔΔG > 0: Mutação desestabiliza (hotspot)
   - ΔΔG < 0: Mutação estabiliza
   - ΔΔG ≈ 0: Mutação neutra

5. **Estatística sobre ensemble**
   - Calcula média e desvio padrão sobre n estruturas
   - Erro típico: ±0.5-1.0 kcal/mol

### Critério de Hotspots (±2 REU)

O critério **ΔΔG > +2.0 kcal/mol** para definir hotspots vem de estudos experimentais e computacionais:

**Justificativa experimental:**
- ΔΔG = +2.0 kcal/mol corresponde a redução de ~30x na afinidade a 25°C:
  ```
  K_d,mut / K_d,wt = exp(ΔΔG / RT) = exp(2.0 / 0.592) ≈ 30
  ```
- Efeito biologicamente significativo
- Consenso na literatura (Clackson & Wells, 1995; Kortemme & Baker, 2002)

**Justificativa computacional:**
- Erro médio do Flex ddG: ~1.0-1.3 kcal/mol (RMSE)
- ΔΔG > 2.0 kcal/mol está acima do ruído computacional
- Correlação experimental: r = 0.69 (Barlow et al., 2018)

**Interpretação:**
- **ΔΔG > +2.0 kcal/mol**: Hotspot confirmado (alta confiança)
- **ΔΔG = +1.0 a +2.0 kcal/mol**: Contribuição moderada
- **ΔΔG < +1.0 kcal/mol**: Contribuição fraca ou neutra

## Instalação

### Requisitos

- Python 3.8+
- PyRosetta (instalado via conda ou pip)

### Passo 1: Criar ambiente conda

```bash
conda create -n pyrosetta python=3.10
conda activate pyrosetta
```

### Passo 2: Instalar PyRosetta

```bash
# Via conda (recomendado)
conda install -c conda-forge pyrosetta

# OU via pip (requer licença acadêmica)
pip install pyrosetta-installer
python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"
```

### Passo 3: Instalar framework

```bash
git clone https://github.com/madsondeluna/alaninescanning4eif4e.git
cd alaninescanning4eif4e
pip install -e .
```

### Testar instalação

```python
import pyrosetta
pyrosetta.init()
print("PyRosetta instalado corretamente")
```

## Uso

### Análise completa do eIF4E

```bash
cd eif4e-test
python run_analysis.py
```

Este script:
1. Carrega a estrutura PDB
2. Identifica todos os resíduos mutáveis
3. Executa Flex ddG para todas as mutações
4. Identifica hotspots (ΔΔG > +2.0 kcal/mol)
5. Gera relatórios e visualizações

### Uso programático

```python
from rosetta_scan.protocols.flex_ddg_pyrosetta import (
    FlexDdGPyRosetta,
    FlexDdGConfig
)

# Configurar protocolo
config = FlexDdGConfig(
    nstruct=35,      # Número de estruturas no ensemble
    iterations=3,     # Iterações backrub
    temperature=0.6,  # Temperatura Monte Carlo
    repack_radius=8.0 # Raio de repacking (Å)
)

protocol = FlexDdGPyRosetta(config)

# Definir mutações
mutations = [
    {'mutation': 'W39A', 'position': 39, 'original_aa': 'W'},
    {'mutation': 'Y22A', 'position': 22, 'original_aa': 'Y'},
]

# Executar
results = protocol.run_flex_ddg(
    'protein.pdb',
    mutations,
    output_dir='results'
)

# Identificar hotspots
hotspots = results[results['ddg'] > 2.0]
print(hotspots)
```

## Parâmetros do Protocolo

### `nstruct` (padrão: 35)
Número de estruturas no ensemble por mutação.

- **Aumentar (50-100)**: Maior precisão estatística, menor desvio padrão, mas tempo aumenta linearmente
- **Diminuir (10-20)**: Teste rápido, menor confiabilidade
- **Recomendado**: 35 (padrão do paper original)

### `iterations` (padrão: 3)
Número de iterações do movimento backrub.

- **Aumentar (5-10)**: Maior exploração conformacional, melhor modelagem de flexibilidade, mas muito mais lento
- **Diminuir (1-2)**: Mais rápido, menor flexibilidade, pode perder efeitos alostéricos
- **Recomendado**: 3 (balanço entre precisão e custo)

### `temperature` (padrão: 0.6 K)
Temperatura para aceitação Monte Carlo dos movimentos backrub.

### `repack_radius` (padrão: 8.0 Å)
Raio ao redor da mutação para reempacotamento de cadeias laterais.

## Tempo de Execução

Para 132 mutações (caso eIF4E):

- **Teste (nstruct=5)**: ~10-20 minutos
- **Padrão (nstruct=35)**: ~2-4 horas
- **Alta precisão (nstruct=100)**: ~6-10 horas

Tempo por mutação: ~1-2 minutos com nstruct=35

## Estrutura do Projeto

```
alaninescanning4eif4e/
├── README.md                      # Este arquivo
├── requirements.txt               # Dependências Python
├── setup.py                       # Instalação do pacote
│
├── src/
│   └── rosetta_scan/
│       ├── protocols/
│       │   └── flex_ddg_pyrosetta.py  # Implementação Flex ddG
│       └── analysis/
│           ├── parser.py          # Parser de resultados
│           └── visualizer.py      # Visualizações
│
└── eif4e-test/
    ├── README.md                  # Documentação do exemplo
    ├── model.pdb                  # Estrutura eIF4E
    ├── run_analysis.py            # Script principal
    └── analysis_results/          # Resultados gerados
```

## Exemplo: eIF4E

O diretório `eif4e-test/` contém análise completa da proteína eIF4E (Eukaryotic translation initiation factor 4E).

**Principais descobertas:**
- 132 mutações para alanina
- ~28 hotspots identificados (ΔΔG > +2.0 kcal/mol)
- Cluster aromático crítico (resíduos 38-42) para ligação ao cap m7G
- Interface C-terminal para interação com eIF4G

Ver `eif4e-test/README.md` para detalhes.

## Referências

**Protocolo Flex ddG:**
> Barlow KA, et al. (2018) Flex ddG: Rosetta Ensemble-Based Estimation of Changes in Protein-Protein Binding Affinity upon Mutation. J Phys Chem B. 122(21):5389-5399.

**Função de energia REF2015:**
> Alford RF, et al. (2017) The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design. J Chem Theory Comput. 13(6):3031-3048.

**Hotspots em alanine scanning:**
> Clackson T, Wells JA. (1995) A hot spot of binding energy in a hormone-receptor interface. Science. 267(5196):383-6.

**Recursos:**
- [Rosetta Commons](https://www.rosettacommons.org/)
- [PyRosetta Tutorials](http://www.pyrosetta.org/tutorials)
- [Flex ddG Documentation](https://www.rosettacommons.org/docs/latest/application_documentation/analysis/flex-ddg)

## Licença

MIT License

## Contato

Madson Luna - madsondeluna@gmail.com

GitHub: https://github.com/madsondeluna/alaninescanning4eif4e
