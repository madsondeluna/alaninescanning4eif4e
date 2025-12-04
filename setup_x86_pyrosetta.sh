#!/bin/bash
################################################################################
# Setup x86 PyRosetta em M1/M2 Mac via Rosetta 2
#
# Este script automatiza a instalação do PyRosetta x86 que funciona via emulação
# Rosetta 2, contornando o bug de segmentation fault do PyRosetta ARM nativo.
#
# Uso: bash setup_x86_pyrosetta.sh
################################################################################

set -e  # Exit on error

echo "================================================================================"
echo "Setup x86 PyRosetta para M1/M2 Mac"
echo "================================================================================"
echo ""

# Cores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Verificar se está em Mac
if [[ "$OSTYPE" != "darwin"* ]]; then
    echo -e "${RED}ERRO: Este script é apenas para macOS${NC}"
    exit 1
fi

# Verificar se é M1/M2
ARCH=$(uname -m)
if [[ "$ARCH" != "arm64" ]]; then
    echo -e "${YELLOW}AVISO: Sistema não é ARM64 (M1/M2). Você está em: $ARCH${NC}"
    echo "Este script é otimizado para M1/M2, mas pode funcionar."
    read -p "Continuar? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo -e "${GREEN}✓${NC} Sistema: macOS ARM64 (M1/M2)"
echo ""

# Passo 1: Verificar/Instalar Rosetta 2
echo "================================================================================
Passo 1/5: Verificar Rosetta 2
================================================================================"

if /usr/bin/pgrep -q oahd; then
    echo -e "${GREEN}✓${NC} Rosetta 2 já instalado"
else
    echo -e "${YELLOW}!${NC} Instalando Rosetta 2..."
    /usr/sbin/softwareupdate --install-rosetta --agree-to-license

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓${NC} Rosetta 2 instalado com sucesso"
    else
        echo -e "${RED}✗${NC} Falha ao instalar Rosetta 2"
        exit 1
    fi
fi
echo ""

# Passo 2: Verificar conda
echo "================================================================================"
echo "Passo 2/5: Verificar Conda"
echo "================================================================================"

if command -v conda &> /dev/null; then
    CONDA_PATH=$(which conda)
    echo -e "${GREEN}✓${NC} Conda encontrado: $CONDA_PATH"
else
    echo -e "${RED}✗${NC} Conda não encontrado"
    echo ""
    echo "Por favor, instale Miniconda ou Anaconda primeiro:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Source conda para poder usar conda activate
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
echo ""

# Passo 3: Criar ambiente x86
echo "================================================================================"
echo "Passo 3/5: Criar ambiente conda x86"
echo "================================================================================"

ENV_NAME="pyrosetta_x86"

if conda env list | grep -q "^$ENV_NAME "; then
    echo -e "${YELLOW}!${NC} Ambiente $ENV_NAME já existe"
    read -p "Remover e recriar? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removendo ambiente existente..."
        conda env remove -n $ENV_NAME -y
    else
        echo "Usando ambiente existente"
        conda activate $ENV_NAME
        echo ""
        echo "================================================================================"
        echo "Passo 4/5: Instalar PyRosetta x86"
        echo "================================================================================"
        # Pular para instalação de pacotes
    fi
fi

if ! conda env list | grep -q "^$ENV_NAME "; then
    echo "Criando ambiente conda x86..."
    CONDA_SUBDIR=osx-64 conda create -n $ENV_NAME python=3.10 -y

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓${NC} Ambiente criado com sucesso"
    else
        echo -e "${RED}✗${NC} Falha ao criar ambiente"
        exit 1
    fi
fi

# Ativar ambiente
conda activate $ENV_NAME

# Configurar para sempre usar x86
conda config --env --set subdir osx-64

echo -e "${GREEN}✓${NC} Ambiente x86 configurado"
echo ""

# Passo 4: Instalar PyRosetta
echo "================================================================================"
echo "Passo 4/5: Instalar PyRosetta x86"
echo "================================================================================"

# Verificar se PyRosetta já está instalado
if python -c "import pyrosetta" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} PyRosetta já instalado"
    PYROSETTA_VERSION=$(python -c "import pyrosetta; print(pyrosetta.version())" 2>/dev/null || echo "unknown")
    echo "   Versão: $PYROSETTA_VERSION"
else
    echo "Instalando PyRosetta x86..."
    echo "Isso pode demorar alguns minutos..."

    # Tentar canal JHU primeiro
    conda install -c https://conda.graylab.jhu.edu pyrosetta -y

    if [ $? -ne 0 ]; then
        echo -e "${YELLOW}!${NC} Tentando canal alternativo..."
        conda install -c https://conda.rosettacommons.org pyrosetta -y
    fi

    # Verificar instalação
    if python -c "import pyrosetta" 2>/dev/null; then
        echo -e "${GREEN}✓${NC} PyRosetta instalado com sucesso"
    else
        echo -e "${RED}✗${NC} Falha ao instalar PyRosetta"
        echo ""
        echo "Tente instalar manualmente:"
        echo "  conda activate $ENV_NAME"
        echo "  conda install -c https://conda.graylab.jhu.edu pyrosetta"
        exit 1
    fi
fi
echo ""

# Passo 5: Instalar dependências
echo "================================================================================"
echo "Passo 5/5: Instalar dependências Python"
echo "================================================================================"

echo "Instalando pandas, numpy, matplotlib, seaborn, biopython..."
pip install pandas numpy matplotlib seaborn biopython -q

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓${NC} Dependências instaladas"
else
    echo -e "${YELLOW}!${NC} Algumas dependências falharam (não crítico)"
fi
echo ""

# Verificação final
echo "================================================================================"
echo "VERIFICAÇÃO FINAL"
echo "================================================================================"

echo -n "Arquitetura Python: "
PYTHON_ARCH=$(python -c "import platform; print(platform.machine())")
if [ "$PYTHON_ARCH" == "x86_64" ]; then
    echo -e "${GREEN}$PYTHON_ARCH ✓${NC}"
else
    echo -e "${RED}$PYTHON_ARCH ✗${NC}"
    echo -e "${RED}ERRO: Python não está rodando em x86_64!${NC}"
    exit 1
fi

echo -n "PyRosetta importável: "
if python -c "import pyrosetta" 2>/dev/null; then
    echo -e "${GREEN}SIM ✓${NC}"
else
    echo -e "${RED}NÃO ✗${NC}"
    exit 1
fi

echo -n "Versão PyRosetta: "
PYROSETTA_VERSION=$(python -c "import pyrosetta; print(pyrosetta.version())" 2>/dev/null)
echo "$PYROSETTA_VERSION"

echo ""
echo -n "Teste de inicialização: "
if python -c "import pyrosetta; pyrosetta.init('-mute all')" 2>/dev/null; then
    echo -e "${GREEN}OK ✓${NC}"
else
    echo -e "${RED}FALHOU ✗${NC}"
    exit 1
fi

echo ""
echo "================================================================================"
echo -e "${GREEN}INSTALAÇÃO COMPLETA!${NC}"
echo "================================================================================"
echo ""
echo "Para usar PyRosetta x86, sempre ative o ambiente primeiro:"
echo ""
echo "  ${GREEN}conda activate $ENV_NAME${NC}"
echo ""
echo "Aliases úteis (adicione ao ~/.zshrc):"
echo ""
echo "  alias pyrosetta='conda activate pyrosetta_x86'"
echo "  alias eif4e='cd $(pwd)/eif4e-test && conda activate pyrosetta_x86'"
echo ""
echo "Próximos passos:"
echo ""
echo "  1. Ativar ambiente:"
echo "     ${GREEN}conda activate $ENV_NAME${NC}"
echo ""
echo "  2. Testar com gamma.pdb:"
echo "     ${GREEN}cd $(pwd)/eif4e-test${NC}"
echo "     ${GREEN}python test_gamma.py${NC}"
echo ""
echo "  3. Rodar análise completa:"
echo "     ${GREEN}python run_analysis_high_accuracy.py${NC}"
echo ""
echo "================================================================================"
echo ""

# Criar arquivo de teste rápido
cat > /tmp/test_pyrosetta_x86.py << 'EOF'
#!/usr/bin/env python3
"""Teste rápido de PyRosetta x86"""
import sys
import platform

print("=" * 80)
print("TESTE PyRosetta x86")
print("=" * 80)

# Verificar arquitetura
arch = platform.machine()
print(f"Arquitetura Python: {arch}")
if arch != "x86_64":
    print("❌ ERRO: Não está em x86_64!")
    sys.exit(1)
print("✅ Arquitetura correta")

# Importar PyRosetta
try:
    import pyrosetta
    print("✅ PyRosetta importado")
except ImportError as e:
    print(f"❌ ERRO ao importar PyRosetta: {e}")
    sys.exit(1)

# Versão
try:
    version = pyrosetta.version()
    print(f"✅ Versão: {version}")
except Exception as e:
    print(f"⚠️  Não foi possível obter versão: {e}")

# Inicializar
try:
    pyrosetta.init("-ignore_unrecognized_res -mute all")
    print("✅ PyRosetta inicializado")
except Exception as e:
    print(f"❌ ERRO ao inicializar: {e}")
    sys.exit(1)

# Criar pose simples
try:
    pose = pyrosetta.pose_from_sequence("AAAAAA")
    print(f"✅ Pose criada: {pose.total_residue()} residues")
except Exception as e:
    print(f"❌ ERRO ao criar pose: {e}")
    sys.exit(1)

# Testar mutação (crítico!)
try:
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    mutate = MutateResidue(target=3, new_res='W')
    mutate.apply(pose)
    print(f"✅ Mutação aplicada: {pose.residue(3).name()}")
except Exception as e:
    print(f"❌ ERRO ao mutar: {e}")
    sys.exit(1)

print("=" * 80)
print("✅ TODOS OS TESTES PASSARAM!")
print("=" * 80)
print("\nPyRosetta x86 está funcionando corretamente via Rosetta 2.")
print("Agora você pode rodar suas análises de alanine scanning.")
EOF

echo "Deseja rodar um teste completo agora? (y/n)"
read -p "" -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo ""
    echo "Rodando teste completo..."
    echo ""
    python /tmp/test_pyrosetta_x86.py

    if [ $? -eq 0 ]; then
        echo ""
        echo -e "${GREEN}✓ Teste completo passou!${NC}"
        echo ""
        echo "PyRosetta x86 está pronto para uso!"
    else
        echo ""
        echo -e "${RED}✗ Teste falhou${NC}"
        echo "Verifique os erros acima"
    fi
fi

echo ""
echo "Script finalizado."
echo ""
