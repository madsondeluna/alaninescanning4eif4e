# Análise de Alanine Scanning - eIF4E

Exemplo completo de alanine scanning da proteína eIF4E (Eukaryotic translation initiation factor 4E) usando PyRosetta.

## Descrição

Este diretório contém análise computacional da proteína eIF4E para identificar resíduos críticos (hotspots) através de alanine scanning sistemático. O método calcula ΔΔG usando o protocolo Flex ddG implementado com PyRosetta.

## Sobre a Proteína eIF4E

**Função biológica:**
- Reconhece e liga-se ao cap 5' m7-guanosina do mRNA
- Componente essencial do complexo eIF4F (iniciação da tradução)
- Alvo terapêutico em câncer (superexpressão em tumores)

**Estrutura:**
- 157 resíduos de aminoácidos
- Domínio cap-binding com cluster aromático conservado
- Interface C-terminal para ligação com eIF4G

## Arquivos

- [model.pdb](model.pdb) - Estrutura da eIF4E (157 resíduos)
- [run_analysis.py](run_analysis.py) - Script principal de análise
- **analysis_results/** - Diretório com resultados gerados

## Como Executar

### Requisitos

Certifique-se de ter o ambiente configurado:

```bash
conda activate pyrosetta
```

### Execução

```bash
cd /Volumes/promethion/alaninescanning4eif4e/eif4e-test
python run_analysis.py
```

### O Que o Script Faz

1. **Carrega estrutura PDB**
   - Lê model.pdb usando BioPython

2. **Identifica resíduos mutáveis**
   - Exclui GLY, ALA, PRO (não podem virar alanina)
   - Total esperado: ~132 mutações

3. **Executa protocolo Flex ddG**
   ```
   Para cada mutação X→A:
     - Gera ensemble backrub (nstruct conformações)
     - Introduz mutação para alanina
     - Reempacota cadeias laterais (raio 8 Å)
     - Calcula: ΔΔG = E_mutante - E_wild-type
   ```

4. **Identifica hotspots**
   - Critério: ΔΔG > +2.0 kcal/mol
   - Hotspots = resíduos críticos para estabilidade/função

5. **Gera relatórios**
   - ddg_results.csv: Todos os valores de ΔΔG
   - hotspots.csv: Apenas hotspots identificados
   - resumo.txt: Estatísticas e top 10 hotspots

## Tempo de Execução

Com configuração padrão (nstruct=5 para teste):
- **Teste rápido**: ~10-20 minutos (nstruct=5)
- **Produção**: ~2-4 horas (nstruct=35)
- **Alta precisão**: ~6-10 horas (nstruct=100)

Tempo por mutação: ~1-2 minutos com nstruct=35

## Configuração do Protocolo

Em [run_analysis.py](run_analysis.py) (linhas 97-102):

```python
config = FlexDdGConfig(
    nstruct=5,        # Alterar para 35 em produção
    iterations=3,
    temperature=0.6,
    repack_radius=8.0
)
```

**Parâmetros:**

- **nstruct=5**: Número de estruturas no ensemble
  - Teste rápido: 5
  - Produção: 35 (recomendado)
  - Alta precisão: 50-100

- **iterations=3**: Iterações do movimento backrub
  - Menor flexibilidade: 1-2
  - Balanceado: 3 (recomendado)
  - Maior exploração: 5-10

- **temperature=0.6**: Temperatura Monte Carlo (K)

- **repack_radius=8.0**: Raio de repacking em Å

## Resultados Esperados

### Estatísticas Gerais

```
Total de mutações: 132
Hotspots identificados: ~20-30 (ΔΔG > +2.0 kcal/mol)
ΔΔG médio: ~1.0-1.5 kcal/mol
```

### Hotspots Esperados

Baseado na literatura sobre eIF4E:

**1. Sítio de ligação ao cap (resíduos 22-42)**
- Trp39: Empilhamento π-π com m7-guanosina (hotspot crítico)
- Tyr22, Phe38, Tyr42: Cluster aromático conservado
- Esperado: ΔΔG > +2.5 kcal/mol para Trp39

**2. Core estrutural**
- Trp68, Trp126: Triptofanos enterrados
- Contribuição para estabilidade do fold
- Esperado: ΔΔG > +2.0 kcal/mol

**3. Interface eIF4G (resíduos 145-157)**
- Resíduos carregados (Arg, Lys) e aromáticos (Tyr, Phe)
- Mediação da interação eIF4E-eIF4G
- Esperado: ΔΔG moderado (+1.5 a +2.5 kcal/mol)

## Arquivos Gerados

```
analysis_results/
├── ddg_results.csv      # Todos os valores de ΔΔG
├── hotspots.csv         # Apenas hotspots (ΔΔG > +2.0)
└── resumo.txt           # Resumo executivo
```

### Formato ddg_results.csv

```csv
mutation,position,chain,original_aa,ddg,ddg_std,E_wt,E_mut
W39A,39,A,W,3.25,0.45,-450.2,-446.95
Y22A,22,A,Y,2.67,0.38,-450.2,-447.53
...
```

**Colunas:**
- `mutation`: Mutação (ex: W39A)
- `position`: Posição no PDB
- `chain`: Cadeia (geralmente A)
- `original_aa`: Aminoácido original
- `ddg`: ΔΔG médio (kcal/mol)
- `ddg_std`: Desvio padrão
- `E_wt`: Energia wild-type (REU)
- `E_mut`: Energia mutante (REU)

## Interpretação dos Resultados

### Critério de Hotspots: ΔΔG > +2.0 kcal/mol

**Fundamentação:**

1. **Experimental:**
   ```
   K_d,mut / K_d,wt = exp(ΔΔG / RT)

   Para ΔΔG = +2.0 kcal/mol a 25°C:
   K_d,mut / K_d,wt = exp(2.0 / 0.592) ≈ 30
   ```
   Redução de ~30x na afinidade é biologicamente significativa.

2. **Computacional:**
   - Erro médio do Flex ddG: ~1.0-1.3 kcal/mol (RMSE)
   - ΔΔG > 2.0 está acima do ruído computacional
   - Correlação com dados experimentais: r = 0.69

**Classificação:**
- **ΔΔG > +2.0 kcal/mol**: Hotspot (alta confiança)
- **ΔΔG = +1.0 a +2.0 kcal/mol**: Contribuição moderada
- **ΔΔG < +1.0 kcal/mol**: Contribuição fraca/neutra
- **ΔΔG < 0**: Mutação estabilizadora (raro em alanine scanning)

## Validação Experimental

Hotspots identificados computacionalmente podem ser validados por:

**1. Mutagênese sítio-dirigida**
- Construir mutantes W39A, Y22A, etc.
- Medir afinidade ao cap por fluorescence anisotropy
- Comparar K_d,mut / K_d,wt experimental vs. computacional

**2. Ensaios funcionais**
- Pull-down assays para interface eIF4G
- Translation initiation assays in vitro
- Western blot para formação do complexo eIF4F

**3. Estruturais**
- Cristalografia dos mutantes
- NMR chemical shift perturbations
- Hydrogen-deuterium exchange (HDX-MS)

## Aplicações Terapêuticas

**Contexto:** eIF4E está superexpressa em diversos cânceres (mama, próstata, pulmão).

**Estratégias de inibição baseadas em hotspots:**

1. **Inibidores do sítio de ligação ao cap**
   - Explorar Trp39 e cluster aromático 38-42
   - Análogos do m7GTP (ex: 4EGI-1)
   - Design racional baseado em hotspots

2. **Disruptores da interface eIF4E-eIF4G**
   - Explorar região C-terminal (145-157)
   - Peptídeos miméticos
   - Small molecules que bloqueiam PPI

3. **Moduladores alostéricos**
   - Explorar core estrutural
   - Desestabilizar conformação ativa

## Limitações

1. **Modelagem computacional**
   - Função de energia não é perfeita
   - Erro médio de ~1.0-1.3 kcal/mol
   - Alguns falsos positivos/negativos esperados

2. **Estrutura estática**
   - Baseado em cristalografia (estrutura congelada)
   - Não captura totalmente dinâmica em solução

3. **Contexto celular**
   - Resultados in vitro podem diferir in vivo
   - Interações com outros fatores não modeladas

4. **Apenas mutações para alanina**
   - Não testa outras substituições
   - Alanina pode não ser a pior mutação possível

## Próximos Passos

1. **Ampliar análise**
   - Testar substituições além de alanina
   - Mutantes duplos/triplos para cooperatividade

2. **Validação experimental**
   - Focar nos top 10 hotspots
   - Priorizar Trp39 e cluster 38-42

3. **Dinâmica molecular**
   - Simular mutantes em solução
   - Analisar impacto na flexibilidade

4. **Drug design**
   - Virtual screening usando estrutura 3D
   - Explorar hotspots para docking

## Referências

**Protocolo Flex ddG:**
> Barlow KA, et al. (2018) Flex ddG: Rosetta Ensemble-Based Estimation of Changes in Protein-Protein Binding Affinity upon Mutation. J Phys Chem B. 122(21):5389-5399.

**eIF4E e câncer:**
> Pelletier J, et al. (2015) Targeting the eIF4F translation initiation complex: a critical nexus for cancer development. Cancer Res. 75(2):250-263.

**Estrutura eIF4E:**
> Marcotrigiano J, et al. (1997) A conserved HEAT domain within eIF4G directs assembly of the translation initiation machinery. Mol Cell. 7(1):193-203.

## Contato

Madson Luna - madsondeluna@gmail.com

GitHub: https://github.com/madsondeluna/alaninescanning4eif4e
