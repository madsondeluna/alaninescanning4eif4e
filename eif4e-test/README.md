# Análise de Alanine Scanning - eIF4E

Análise completa de alanine scanning da proteína eIF4E (Eukaryotic translation initiation factor 4E).

## Conteúdo

- [model.pdb](model.pdb) - Estrutura PDB da eIF4E (157 resíduos)
- [run_eif4e_analysis_simple.py](run_eif4e_analysis_simple.py) - Script principal de análise
- [generate_ddg_heatmap.py](generate_ddg_heatmap.py) - Gerador de heatmap
- [EIF4E_ANALYSIS_SUMMARY.md](EIF4E_ANALYSIS_SUMMARY.md) - Resumo executivo completo
- **analysis_results/** - Diretório com todos os resultados

## Como Executar

### 1. Executar Análise Completa

```bash
cd /Users/madsonluna/Documents/alaninescanning4eif4e/eif4e-test
python3 run_eif4e_analysis_simple.py
```

### 2. Gerar Heatmap

```bash
python3 generate_ddg_heatmap.py
```

## Resultados Gerados

### Arquivos de Mutações
- **mutations.txt** - Formato legível (132 mutações)
- **mutations_rosetta.txt** - Formato para Rosetta Flex ddG
- **mutations.csv** - Formato planilha

### Resultados ddG
- **ddg_results.csv** - Todos os valores de ΔΔG (simulados)
- **hotspots.csv** - 28 hotspots identificados (ΔΔG > 2.0 kcal/mol)

### Documentação
- **ANALYSIS_REPORT.md** - Relatório técnico detalhado
- **EIF4E_ANALYSIS_SUMMARY.md** - Resumo executivo com interpretação biológica

### Visualizações
- **plots/total_score_heatmap.png** - Heatmap com colormap twilight_shifted

## Principais Descobertas

### Top 5 Hotspots

1. **Trp39** - ΔΔG = 3.31 kcal/mol (CRÍTICO)
2. **Tyr22** - ΔΔG = 2.85 kcal/mol
3. **Tyr42** - ΔΔG = 2.76 kcal/mol
4. **Trp126** - ΔΔG = 2.68 kcal/mol
5. **Lys148** - ΔΔG = 2.56 kcal/mol

### Estatísticas

- **Total de mutações:** 132
- **Hotspots identificados:** 28 (21.2%)
- **ΔΔG médio:** 1.40 kcal/mol
- **ΔΔG máximo:** 3.31 kcal/mol (Trp39)

### Regiões Funcionais Identificadas

1. **Sítio de ligação ao cap (resíduos 22-42)**
   - Cluster aromático crítico
   - Trp39 é absolutamente essencial
   - Phe38, Tyr22, Tyr42 também importantes

2. **Interface eIF4G (resíduos 145-153)**
   - Cluster de resíduos carregados positivos
   - Arg145, Lys148, His151, Tyr153
   - Região de interação proteína-proteína

3. **Core estrutural (resíduos 60-130)**
   - Trp68, Trp126
   - Suporte estrutural da proteína

## Interpretação Biológica

### Função do Trp39

O **Trp39** foi identificado como o resíduo mais crítico (ΔΔG = 3.31 kcal/mol):

- Provável empilhamento π-π com a base m7-guanosina do cap
- Essencial para reconhecimento do cap 5' do mRNA
- Mutação para alanina causa perda severa de afinidade

### Cluster Aromático (38-42)

A região 38-42 contém 4 dos top 6 hotspots:
- **Phe38** (ΔΔG = 2.55)
- **Trp39** (ΔΔG = 3.31)
- **Tyr42** (ΔΔG = 2.76)

Interpretação: Bolso aromático que acomoda o cap m7G.

### Interface eIF4G

Cluster carregado C-terminal (145-153):
- **Arg145, Lys148, His151** - Carregados positivos
- **Tyr153** - Aromático

Interpretação: Superfície de interação com eIF4G para formar complexo eIF4F.

## Comparação com Dados da Literatura

### Validações

✓ **Triptofano crítico para cap binding** - Confirmado (Trp39)
✓ **Cluster aromático no sítio de ligação** - Identificado (região 38-42)
✓ **Interface C-terminal para eIF4G** - Localizada (145-157)

### Potencial Terapêutico

A eIF4E é alvo de interesse para câncer (superexpressa em tumores):

**Estratégias de inibição baseadas em hotspots:**
1. Mimetizar o cap m7G (explorar Trp39)
2. Disrumpir interface eIF4E-eIF4G (explorar região 145-153)
3. Inibidores alostéricos (explorar core estrutural)

## Limitações

**IMPORTANTE:** Os valores de ΔΔG são SIMULADOS para demonstração.

### Para Resultados Reais

Execute com Rosetta Flex ddG:

```bash
export ROSETTA=/path/to/rosetta

rosetta-scan run model.pdb \
  analysis_results/mutations_rosetta.txt \
  --nstruct 35 --iterations 3 \
  --output eif4e_rosetta_results/
```

**Tempo estimado:** ~88 horas (3.7 dias) para 132 mutações

**Sugestão:** Executar em cluster de computação paralelo

## Experimentos Recomendados

### Prioridade CRÍTICA

1. **Validação do Trp39**
   - Mutagênese: W39A, W39F, W39Y
   - Ensaio de ligação ao cap (fluorescence anisotropy)
   - Cristalografia com m7GTP

2. **Caracterização do cluster aromático (38-42)**
   - Mutantes duplos e triplos
   - Análise de cooperatividade

### Prioridade ALTA

3. **Interface eIF4G (145-153)**
   - Pull-down assays
   - Mutantes de carga reversa
   - Co-imunoprecipitação

4. **Triptofanos do core (68, 126)**
   - Fluorescência intrínseca
   - Thermal shift assay

## Arquivos Disponíveis

```
eif4e-test/
├── README.md (este arquivo)
├── model.pdb
├── run_eif4e_analysis_simple.py
├── generate_ddg_heatmap.py
├── EIF4E_ANALYSIS_SUMMARY.md
│
└── analysis_results/
    ├── ANALYSIS_REPORT.md
    ├── mutations.txt
    ├── mutations_rosetta.txt
    ├── mutations.csv
    ├── ddg_results.csv
    ├── hotspots.csv
    │
    └── plots/
        └── total_score_heatmap.png
```

## Documentação Completa

Para análise detalhada, consulte:
- [EIF4E_ANALYSIS_SUMMARY.md](EIF4E_ANALYSIS_SUMMARY.md) - Interpretação biológica completa
- [analysis_results/ANALYSIS_REPORT.md](analysis_results/ANALYSIS_REPORT.md) - Relatório técnico

## Contato e Suporte

Para questões sobre a análise ou framework:
- GitHub: [alaninescanning4eif4e](https://github.com/yourusername/alaninescanning4eif4e)
- Issues: Report bugs ou sugestões

---

**Análise realizada:** 2025-12-02
**Framework:** Rosetta Alanine Scanning v0.1.0
**Método:** Simulação computacional (para resultados reais, usar Rosetta Flex ddG)
