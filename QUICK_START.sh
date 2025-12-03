#!/bin/bash
# ============================================================================
# ROSETTA ALANINE SCANNING - QUICK START
# ============================================================================
#
# Execute este script no terminal do VSCode para configurar e testar
# o framework completo de alanine scanning.
#
# Uso:
#   chmod +x QUICK_START.sh
#   ./QUICK_START.sh
#
# Ou execute os comandos manualmente, linha por linha.
# ============================================================================

set -e  # Para na primeira erro

echo "============================================================================"
echo "ROSETTA ALANINE SCANNING - Setup e Demo"
echo "============================================================================"
echo ""

# ============================================================================
# PASSO 1: Verificar Localização
# ============================================================================

echo "PASSO 1: Verificando diretório..."
echo ""

if [ ! -f "setup.py" ]; then
    echo "Erro: Execute este script no diretório raiz do projeto"
    echo "   cd /Volumes/promethion/alaninescanning4eif4e"
    exit 1
fi

echo "Diretório correto!"
echo "   $(pwd)"
echo ""

# ============================================================================
# PASSO 2: Criar Ambiente Virtual
# ============================================================================

echo "PASSO 2: Configurando ambiente virtual..."
echo ""

if [ ! -d "venv" ]; then
    echo "Criando ambiente virtual..."
    python3 -m venv venv
    echo "Ambiente virtual criado!"
else
    echo "Ambiente virtual já existe!"
fi

echo ""
echo "Ativando ambiente virtual..."
source venv/bin/activate

echo "Ambiente ativado: $VIRTUAL_ENV"
echo ""

# ============================================================================
# PASSO 3: Instalar Dependências
# ============================================================================

echo "PASSO 3: Instalando dependências..."
echo ""

# Upgrade pip
echo "Atualizando pip..."
pip install --upgrade pip -q

# Instalar pacote
echo "Instalando rosetta-scan..."
pip install -e . -q

echo ""
echo "Dependências instaladas!"
echo ""

# Verificar instalação
echo "Verificando instalação..."
pip show rosetta-alanine-scanning 2>/dev/null || echo "Pacote instalado em modo dev"
echo ""

# ============================================================================
# PASSO 4: Verificar Comandos
# ============================================================================

echo "PASSO 4: Verificando comandos CLI..."
echo ""

if command -v rosetta-scan &> /dev/null; then
    echo "Comando 'rosetta-scan' disponível!"
    rosetta-scan --version 2>/dev/null || rosetta-scan --help | head -5
else
    echo "Comando 'rosetta-scan' não encontrado no PATH"
    echo "Você pode usar: python -m rosetta_scan.cli"
fi

echo ""

# ============================================================================
# PASSO 5: Executar Demo
# ============================================================================

echo "PASSO 5: Executando demonstração completa..."
echo ""
echo "Isso irá:"
echo "  1. Analisar estrutura PDB exemplo"
echo "  2. Gerar mutações de alanina"
echo "  3. Simular resultados ddG"
echo "  4. Identificar hotspots"
echo "  5. Criar visualizações"
echo ""

read -p "Pressione ENTER para continuar ou Ctrl+C para cancelar..."
echo ""

# Executar demo
python3 examples/demo_run.py

echo ""
echo "============================================================================"
echo "SETUP COMPLETO!"
echo "============================================================================"
echo ""
echo "Arquivos gerados em: examples/demo_output/"
echo ""
echo "Próximos passos:"
echo ""
echo "1. Explorar resultados:"
echo "   ls -la examples/demo_output/"
echo "   cat examples/demo_output/analysis_report.txt"
echo ""
echo "2. Visualizar plots:"
echo "   open examples/demo_output/plots/"
echo ""
echo "3. Ver hotspots:"
echo "   cat examples/demo_output/hotspots.csv"
echo ""
echo "4. Testar CLI:"
echo "   rosetta-scan --help"
echo "   rosetta-scan scan examples/example_protein.pdb --help"
echo ""
echo "5. Ler documentação:"
echo "   cat README.md"
echo "   cat DEMO.md"
echo ""
echo "============================================================================"
echo "Documentação disponível:"
echo "============================================================================"
echo ""
echo "  • README.md              - Overview e quick start"
echo "  • INSTALL.md             - Guia de instalação"
echo "  • DEMO.md                - Tutorial passo a passo"
echo "  • PROJECT_STRUCTURE.md   - Estrutura do código"
echo "  • examples/EXPECTED_OUTPUT.md - Exemplos de saída"
echo ""
echo "============================================================================"
echo "Pronto para usar!"
echo "============================================================================"
echo ""
