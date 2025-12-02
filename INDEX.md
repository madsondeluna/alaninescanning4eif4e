#  Índice de Navegação - Rosetta Alanine Scanning

Guia rápido para encontrar o que você precisa.

---

##  COMEÇANDO (Leia Primeiro!)

###  Se você quer rodar AGORA:
1. **[COMANDOS_TERMINAL.txt](COMANDOS_TERMINAL.txt)**  - Copie/cole comandos prontos no terminal
2. **[RUN_ME_FIRST.md](RUN_ME_FIRST.md)** - Guia de setup passo a passo
3. **[QUICK_START.sh](QUICK_START.sh)** - Script automático de instalação

###  Se você quer entender antes:
1. **[README.md](README.md)** - Overview geral do projeto
2. **[SUMMARY.md](SUMMARY.md)** - Sumário completo de tudo que foi criado

---

##  DOCUMENTAÇÃO COMPLETA

###  Tutoriais e Guias

| Arquivo | Conteúdo | Quando Usar |
|---------|----------|-------------|
| **[DEMO.md](DEMO.md)** | Tutorial completo com exemplos | Para aprender a usar o framework |
| **[INSTALL.md](INSTALL.md)** | Guia detalhado de instalação | Se tiver problemas de instalação |
| **[examples/EXPECTED_OUTPUT.md](examples/EXPECTED_OUTPUT.md)** | Exemplos de saída do programa | Para ver como ficam os resultados |

###  Arquitetura e Código

| Arquivo | Conteúdo | Quando Usar |
|---------|----------|-------------|
| **[PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)** | Estrutura completa do código | Para entender a organização |
| **[SUMMARY.md](SUMMARY.md)** | Sumário técnico detalhado | Para visão geral técnica |

---

##  NAVEGAÇÃO RÁPIDA POR TAREFA

### Tarefa 1: "Quero INSTALAR"
```
1. README.md (seção Installation)
2. INSTALL.md (guia completo)
3. QUICK_START.sh (script automático)
```

### Tarefa 2: "Quero RODAR um demo"
```
1. COMANDOS_TERMINAL.txt (comandos prontos)
2. RUN_ME_FIRST.md (passo a passo)
3. examples/demo_run.py (script demo)
```

### Tarefa 3: "Quero APRENDER a usar"
```
1. DEMO.md (tutorial completo)
2. README.md (quick start)
3. examples/EXPECTED_OUTPUT.md (ver exemplos)
```

### Tarefa 4: "Quero ENTENDER o código"
```
1. PROJECT_STRUCTURE.md (arquitetura)
2. SUMMARY.md (visão técnica)
3. Código fonte em src/rosetta_scan/
```

### Tarefa 5: "Tenho um PROBLEMA"
```
1. INSTALL.md (troubleshooting)
2. RUN_ME_FIRST.md (FAQ)
3. GitHub issues
```

### Tarefa 6: "Quero usar com MINHA proteína"
```
1. DEMO.md (exemplos práticos)
2. COMANDOS_TERMINAL.txt (comandos prontos)
3. config/example_config.yaml (configuração)
```

---

##  ESTRUTURA DE ARQUIVOS

###  Documentação (9 arquivos)

```
Essenciais:
--- README.md                     Overview principal
--- COMANDOS_TERMINAL.txt         Comandos prontos
--- RUN_ME_FIRST.md               Quick start
--- INDEX.md                     ← Você está aqui

Guias Detalhados:
--- DEMO.md                      Tutorial passo a passo
--- INSTALL.md                   Instalação completa
--- PROJECT_STRUCTURE.md         Arquitetura do código
--- SUMMARY.md                   Sumário técnico
--- examples/EXPECTED_OUTPUT.md  Exemplos de saída
```

###  Código Fonte (7 arquivos principais)

```
src/rosetta_scan/
--- cli.py                       CLI com Click
--- protocols/
|   --- flex_ddg.py             Protocolo Flex ddG
|   --- alanine_scanner.py      Gerador de mutações
--- analysis/
    --- parser.py               Parser de resultados
    --- visualizer.py           Visualizações
```

###  Exemplos (4 arquivos)

```
examples/
--- demo_run.py                   Demo interativo completo
--- example_workflow.py          Workflow Python
--- example_protein.pdb          Estrutura PDB teste
--- quick_test.sh               Script de teste
```

###  Configuração (3 arquivos)

```
--- setup.py                    Setup do pacote
--- requirements.txt            Dependências Python
--- config/example_config.yaml  Configuração exemplo
```

---

##  FLUXO RECOMENDADO

### Para Iniciantes:

```
1. Leia: README.md
   ↓
2. Execute: QUICK_START.sh
   ↓
3. Rode: python3 examples/demo_run.py
   ↓
4. Explore: examples/demo_output/
   ↓
5. Leia: DEMO.md
   ↓
6. Use com sua proteína!
```

### Para Usuários Experientes:

```
1. Leia: SUMMARY.md (visão técnica)
   ↓
2. Execute: pip install -e .
   ↓
3. Use: rosetta-scan pipeline sua_proteina.pdb
   ↓
4. Customize: config/example_config.yaml
```

### Para Desenvolvedores:

```
1. Leia: PROJECT_STRUCTURE.md
   ↓
2. Estude: src/rosetta_scan/
   ↓
3. Modifique e teste
   ↓
4. Contribua!
```

---

##  TAMANHO DOS ARQUIVOS

```
Documentação:
--- README.md                  8.4 KB   Comece aqui
--- COMANDOS_TERMINAL.txt      9.6 KB   Comandos prontos
--- RUN_ME_FIRST.md            8.4 KB   Quick start
--- DEMO.md                   13.0 KB  Tutorial completo
--- INSTALL.md                 7.0 KB  Instalação
--- PROJECT_STRUCTURE.md      16.0 KB  Arquitetura
--- SUMMARY.md                12.0 KB  Sumário técnico
--- EXPECTED_OUTPUT.md        ~12 KB   Exemplos

Total Documentação: ~86 KB

Código Python: ~2,300 linhas (~70 KB)
```

---

##  LINKS RÁPIDOS

### Documentação Essencial
- [README.md](README.md) - Overview
- [COMANDOS_TERMINAL.txt](COMANDOS_TERMINAL.txt) - Comandos
- [RUN_ME_FIRST.md](RUN_ME_FIRST.md) - Quick start

### Tutoriais
- [DEMO.md](DEMO.md) - Tutorial completo
- [INSTALL.md](INSTALL.md) - Instalação
- [examples/EXPECTED_OUTPUT.md](examples/EXPECTED_OUTPUT.md) - Exemplos

### Referência Técnica
- [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) - Arquitetura
- [SUMMARY.md](SUMMARY.md) - Sumário técnico
- [src/rosetta_scan/](src/rosetta_scan/) - Código fonte

### Exemplos Práticos
- [examples/demo_run.py](examples/demo_run.py) - Demo completo
- [examples/example_workflow.py](examples/example_workflow.py) - Workflow
- [config/example_config.yaml](config/example_config.yaml) - Config

### Scripts de Setup
- [QUICK_START.sh](QUICK_START.sh) - Setup automático
- [examples/quick_test.sh](examples/quick_test.sh) - Teste rápido

---

##  GLOSSÁRIO DE TERMOS

| Termo | Significado |
|-------|-------------|
| **Alanine Scanning** | Mutagênese sistemática para alanina |
| **Flex ddG** | Protocolo Rosetta para calcular ΔΔG |
| **Hotspot** | Resíduo crítico (alto ΔΔG) |
| **ΔΔG** | Variação de energia livre (kcal/mol) |
| **Interface** | Região de contato entre cadeias |
| **SASA** | Solvent Accessible Surface Area |
| **PDB** | Protein Data Bank (formato de estrutura) |

---

##  SUPORTE

### Problemas Comuns:
1. Erro de instalação → [INSTALL.md](INSTALL.md) seção Troubleshooting
2. Rosetta não encontrado → [INSTALL.md](INSTALL.md) seção Rosetta
3. Módulo não encontrado → Reinstale: `pip install -e .`
4. Plots não aparecem → `export MPLBACKEND=Agg`

### Recursos Adicionais:
- Issues no GitHub
- Documentação do Rosetta: https://www.rosettacommons.org/
- Paper do Flex ddG: Barlow et al. (2018)

---

##  CHECKLIST DE USO

### Primeira Vez:
- [ ] Li o README.md
- [ ] Executei QUICK_START.sh ou setup manual
- [ ] Rodei demo: `python3 examples/demo_run.py`
- [ ] Explorei os outputs gerados
- [ ] Li pelo menos um tutorial (DEMO.md ou RUN_ME_FIRST.md)
- [ ] Testei com estrutura exemplo
- [ ] Pronto para usar com minha proteína!

### Para Cada Projeto:
- [ ] Tenho arquivo PDB da proteína
- [ ] Identifiquei cadeias de interesse
- [ ] Decidi se analiso toda proteína ou apenas interface
- [ ] Configurei parâmetros (se usando Rosetta)
- [ ] Executei scan ou pipeline
- [ ] Analisei resultados
- [ ] Visualizei hotspots
- [ ] Interpretei ΔΔG values

---

##  ONDE ENCONTRAR CADA INFORMAÇÃO

| Quero saber sobre... | Arquivo |
|----------------------|---------|
| Como instalar | INSTALL.md |
| Como usar CLI | README.md, COMANDOS_TERMINAL.txt |
| Como interpretar resultados | DEMO.md, EXPECTED_OUTPUT.md |
| Estrutura do código | PROJECT_STRUCTURE.md |
| Exemplos práticos | examples/*.py |
| Configuração | config/example_config.yaml |
| Troubleshooting | INSTALL.md |
| Visão geral técnica | SUMMARY.md |
| API Python | src/rosetta_scan/*.py |

---

##  COMEÇAR AGORA - 3 OPÇÕES

### 1 Mais Rápido (1 comando)
```bash
./QUICK_START.sh
```

### 2 Manual (4 comandos)
```bash
python3 -m venv venv
source venv/bin/activate
pip install -e .
python3 examples/demo_run.py
```

### 3 Explorar Primeiro
```bash
cat README.md
cat RUN_ME_FIRST.md
```

---

**Escolha seu caminho e boa análise! **

**Dúvidas?** Consulte [COMANDOS_TERMINAL.txt](COMANDOS_TERMINAL.txt) para comandos prontos.
