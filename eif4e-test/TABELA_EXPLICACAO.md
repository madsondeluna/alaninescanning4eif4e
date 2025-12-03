# Explicação da Tabela de Resultados

## Arquivo: `analysis_results/ddg_results.csv`

### Estrutura do CSV

A tabela contém **7 colunas** com informações sobre cada mutação:

| Coluna | Tipo | Descrição | Exemplo |
|--------|------|-----------|---------|
| `mutation` | string | Identificador da mutação no formato RosettaScript | `A39W` |
| `chain` | string | Cadeia da proteína | `A` |
| `position` | int | Posição do resíduo na sequência | `39` |
| `original_aa` | string | Aminoácido original (código de 1 letra) | `W` |
| `ddg` | float | **ΔΔG bruto (positivo)** em kcal/mol | `3.305` |
| `ddg_normalized` | float | **ΔΔG normalizado (negativo)** em kcal/mol | `-3.305` |
| `total_score` | float | Score total de energia (REU) | `-47.66` |

---

## Diferença entre ΔΔG Bruto e Normalizado

### ΔΔG Bruto (Coluna `ddg`)
**Valores POSITIVOS** - Convenção computacional

```
ΔΔG = E(mutante) - E(wild-type)
```

- **Valor positivo** = Mutante é MENOS estável que wild-type
- **Quanto MAIOR o valor** = MAIOR o impacto da mutação
- **Interpretação:**
  - `ddg = 3.31` → Mutação causa grande desestabilização
  - `ddg = 0.13` → Mutação causa pequena desestabilização

**Usado para:**
- Classificar mutações por impacto
- Identificar hotspots (ddg > 2.0)
- Análises estatísticas

---

### ΔΔG Normalizado (Coluna `ddg_normalized`)
**Valores NEGATIVOS** - Convenção bioquímica tradicional

```
ΔΔG_norm = - ΔΔG_bruto
```

- **Valor negativo** = Perda de energia livre na mutação
- **Quanto MAIS NEGATIVO** = MAIOR a perda energética
- **Interpretação:**
  - `ddg_normalized = -3.31` → Grande perda de energia
  - `ddg_normalized = -0.13` → Pequena perda de energia

**Usado para:**
- Visualizações (gráficos lollipop)
- Comunicação com biólogos experimentais
- Consistência com literatura de termodinâmica

---

## Por que duas colunas?

### Vantagens de ter ambas:

1. **Compatibilidade computacional** (`ddg` positivo)
   - Usado por Rosetta e outros softwares
   - Ordenação direta por impacto (maior = pior)
   - Thresholds intuitivos (>2.0 = hotspot)

2. **Interpretação bioquímica** (`ddg_normalized` negativo)
   - Convenção da termodinâmica clássica
   - ΔG negativo = processo espontâneo/favorável
   - Publicações científicas usam valores negativos

3. **Flexibilidade de análise**
   - Escolha o formato que preferir
   - Conversão simples: `ddg_normalized = -ddg`
   - Ambos contêm a mesma informação

---

## Exemplos Práticos

### Exemplo 1: Trp39 (Hotspot Crítico)

| Coluna | Valor | Interpretação |
|--------|-------|---------------|
| `ddg` | `+3.305` | Mutação causa 3.31 kcal/mol de desestabilização |
| `ddg_normalized` | `-3.305` | Perda de 3.31 kcal/mol de energia livre |

**Conclusão:** Trp39 é CRÍTICO para estabilidade

---

### Exemplo 2: Ser88 (Baixo Impacto)

| Coluna | Valor | Interpretação |
|--------|-------|---------------|
| `ddg` | `+0.129` | Mutação causa 0.13 kcal/mol de desestabilização |
| `ddg_normalized` | `-0.129` | Perda de apenas 0.13 kcal/mol |

**Conclusão:** Ser88 tem impacto mínimo na estabilidade

---

## Classificação de Impacto

Usando a coluna `ddg` (bruto, positivo):

| Categoria | ΔΔG Bruto | ΔΔG Normalizado | Cor | Descrição |
|-----------|-----------|-----------------|-----|-----------|
| **Alto impacto** | ≥ 2.5 | ≤ -2.5 | 🔴 Vermelho | Hotspot crítico |
| **Médio impacto** | 1.5 - 2.5 | -2.5 a -1.5 | 🟠 Laranja | Contribuição moderada |
| **Baixo impacto** | < 1.5 | > -1.5 | 🔵 Azul | Impacto mínimo |

---

## Estatísticas - eIF4E

Baseado em 132 mutações:

| Métrica | ΔΔG Bruto | ΔΔG Normalizado |
|---------|-----------|-----------------|
| **Mínimo** | 0.129 | -3.305 |
| **Máximo** | 3.305 | -0.129 |
| **Média** | 1.399 | -1.399 |
| **Mediana** | 1.290 | -1.290 |
| **Desvio Padrão** | 0.643 | 0.643 |

**Observação:** As estatísticas são idênticas (apenas com sinal invertido)

---

## Uso nos Gráficos

### Gráfico Lollipop
- **Eixo Y:** Usa `ddg_normalized` (valores negativos)
- **Eixo X:** Posição do resíduo
- **Cores:** Baseadas em `ddg` (bruto)
- **Linhas tracejadas:**
  - Laranja: -1.5 kcal/mol
  - Vermelha: -2.5 kcal/mol

### Heatmap
- **Valores das células:** `total_score` (REU)
- **Colormap:** twilight_shifted
- **Anotações:** Valores em cada célula

---

## Como Usar no Excel/Planilhas

### Ordenar por Impacto
```excel
=SORT(Tabela, "ddg", DESCENDENTE)
```

### Filtrar Hotspots
```excel
=FILTER(Tabela, ddg >= 2.5)
```

### Converter Bruto ↔ Normalizado
```excel
ddg_normalized = -ddg
ddg = -ddg_normalized
```

---

## Exportando Dados

### Para Publicação
Use `ddg_normalized` (valores negativos):
- Consistente com literatura
- Convenção termodinâmica
- ΔG negativo = desfavorável para mutante

### Para Análise Computacional
Use `ddg` (valores positivos):
- Input para Rosetta
- Thresholds diretos
- Maior = pior

---

## Referências Teóricas

### Definição de ΔΔG

```
ΔΔG = ΔG_mutante - ΔG_wild-type

Onde:
  ΔG = energia livre de Gibbs

Se ΔΔG > 0: Mutante é menos estável
Se ΔΔG < 0: Mutante é mais estável
Se ΔΔG = 0: Sem mudança
```

### Para Alanine Scanning

Em alanine scanning, **sempre esperamos ΔΔG > 0** porque:
- Alanina remove interações específicas
- Perda de contatos = desestabilização
- ΔΔG mede a **contribuição energética** do resíduo original

---

## Top 10 Hotspots - Comparação

| Rank | Mutação | ΔΔG Bruto | ΔΔG Normalizado | Interpretação |
|------|---------|-----------|-----------------|---------------|
| 1 | A39W | +3.31 | -3.31 | Perda crítica de 3.31 kcal/mol |
| 2 | A22Y | +2.85 | -2.85 | Perda de 2.85 kcal/mol |
| 3 | A42Y | +2.76 | -2.76 | Perda de 2.76 kcal/mol |
| 4 | A126W | +2.68 | -2.68 | Perda de 2.68 kcal/mol |
| 5 | A148K | +2.56 | -2.56 | Perda de 2.56 kcal/mol |
| 6 | A38F | +2.55 | -2.55 | Perda de 2.55 kcal/mol |
| 7 | A12R | +2.54 | -2.54 | Perda de 2.54 kcal/mol |
| 8 | A151H | +2.44 | -2.44 | Perda de 2.44 kcal/mol |
| 9 | A145R | +2.40 | -2.40 | Perda de 2.40 kcal/mol |
| 10 | A153Y | +2.39 | -2.39 | Perda de 2.39 kcal/mol |

**Ambas as colunas mostram a mesma informação, apenas com sinal diferente!**

---

## Perguntas Frequentes

### Q1: Qual coluna devo usar?
**R:** Depende do contexto:
- Análise computacional → `ddg` (positivo)
- Visualizações → `ddg_normalized` (negativo)
- Publicações → `ddg_normalized` (negativo)

### Q2: Por que valores negativos são mais intuitivos?
**R:** Na termodinâmica clássica:
- ΔG negativo = processo favorável (liberação de energia)
- ΔG positivo = processo desfavorável (requer energia)
- Bioquímicos estão acostumados com essa convenção

### Q3: Como converter entre as colunas?
**R:** Simples multiplicação por -1:
```python
ddg_normalized = -ddg
ddg = -ddg_normalized
```

### Q4: Os gráficos usam qual coluna?
**R:** Lollipop plot usa `ddg_normalized` (negativo) no eixo Y, mas cores baseadas em `ddg` (positivo)

### Q5: Posso deletar uma das colunas?
**R:** Sim, mas mantenha ambas para flexibilidade. Elas ocupam pouco espaço e são úteis para diferentes análises.

---

## Arquivo CSV Completo

**Localização:** `analysis_results/ddg_results.csv`

**Total de linhas:** 132 mutações + 1 cabeçalho = 133 linhas

**Tamanho:** ~7 KB

**Formato:** UTF-8, separador vírgula

**Pode ser aberto em:**
- Excel
- Google Sheets
- Python (pandas)
- R (read.csv)
- Qualquer editor de texto

---

## Comandos Python Úteis

### Carregar dados
```python
import pandas as pd
df = pd.read_csv('analysis_results/ddg_results.csv')
```

### Ver estatísticas
```python
df[['ddg', 'ddg_normalized']].describe()
```

### Filtrar hotspots
```python
hotspots = df[df['ddg'] >= 2.5]
```

### Plotar distribuição
```python
import matplotlib.pyplot as plt
plt.hist(df['ddg_normalized'], bins=30)
plt.xlabel('ΔΔG Normalizado (kcal/mol)')
plt.ylabel('Frequência')
plt.show()
```

---

**Data:** 2025-12-03
**Análise:** eIF4E Alanine Scanning
**Framework:** Rosetta-based Protocol
