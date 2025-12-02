# Instruções para Análise do Peptídeo Gamma

## Passo 1: Preparar o Ambiente

Abra o terminal do VSCode e execute:

```bash
cd /Users/madsonluna/Documents/alaninescanning4eif4e
```

## Passo 2: Criar Ambiente Virtual (se ainda não criou)

```bash
python3 -m venv venv
source venv/bin/activate
pip install -e .
```

## Passo 3: Executar a Análise

```bash
python3 gamma-example/run_gamma_analysis.py
```

## O que o script faz:

1. Carrega a estrutura gamma.pdb
2. Gera todas as mutações de alanina possíveis
3. Simula valores de ddG (demonstração)
4. Identifica resíduos hotspot
5. Cria visualizações (gráficos PNG)
6. Gera relatório completo em Markdown

## Resultados esperados:

Todos os arquivos serão salvos em `gamma-example/analysis_results/`:

- **mutations.txt** - Lista de mutações
- **mutations_rosetta.txt** - Formato para Rosetta
- **mutations.csv** - Planilha
- **ddg_results.csv** - Resultados de ddG
- **hotspots.csv** - Resíduos críticos
- **ANALYSIS_REPORT.md** - Relatório completo
- **plots/** - Gráficos PNG

## Tempo estimado:

- 1-2 minutos para análise completa

## Próximos passos após execução:

1. Leia `ANALYSIS_REPORT.md` para resultados detalhados
2. Visualize os gráficos em `plots/`
3. Examine `hotspots.csv` para resíduos críticos
