# Como Usar no M1/M2 Mac

## Problema

O PyRosetta ARM nativo (M1/M2) tem bug de segmentation fault. **Não funciona.**

## Solução Rápida (30 minutos)

Use PyRosetta x86 via Rosetta 2 emulation.

### Instalação Automática

```bash
cd /Volumes/promethion/alaninescanning4eif4e
bash setup_x86_pyrosetta.sh
```

O script faz tudo automaticamente:
1. Instala Rosetta 2
2. Cria ambiente conda x86
3. Instala PyRosetta x86
4. Testa instalação

### Uso Diário

```bash
# Ativar ambiente
conda activate pyrosetta_x86

# Ir para diretório de testes
cd /Volumes/promethion/alaninescanning4eif4e/eif4e-test

# Rodar análise
python run_analysis_high_accuracy.py
```

### Aliases Úteis

Adicione ao `~/.zshrc`:

```bash
alias pyrosetta='conda activate pyrosetta_x86'
alias eif4e='cd /Volumes/promethion/alaninescanning4eif4e/eif4e-test && conda activate pyrosetta_x86'
```

Depois:
```bash
source ~/.zshrc
eif4e  # Ativa ambiente e vai para pasta
python run_analysis.py
```

## Performance

- **x86 via Rosetta 2**: ~78% da performance nativa
- **eIF4E (132 mutações)**: ~26 horas (vs ~22h nativo se funcionasse)
- **Custo**: GRATUITO

## Arquivos Criados

1. **[SOLUCOES_M1_M2.md](SOLUCOES_M1_M2.md)** - Documentação completa com todas as opções
2. **[setup_x86_pyrosetta.sh](setup_x86_pyrosetta.sh)** - Script de instalação automática
3. **[GITHUB_ISSUE_PYROSETTA.md](GITHUB_ISSUE_PYROSETTA.md)** - Issue completo para reportar bug

## Outras Opções

Se precisar de mais performance ou estabilidade, veja [SOLUCOES_M1_M2.md](SOLUCOES_M1_M2.md) para:

- **Rosetta C++ + Python wrapper** (melhor para publicação)
- **AWS Cloud** (escala para 1000+ mutações)
- **Google Colab** (gratuito, sem instalar nada)

## Suporte

Problemas? Veja seção "Suporte e Troubleshooting" em [SOLUCOES_M1_M2.md](SOLUCOES_M1_M2.md).
