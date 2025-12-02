# Installation Guide

Complete installation instructions for Rosetta Alanine Scanning.

## Prerequisites

### 1. Rosetta Software Suite

The Rosetta software suite is required for running Flex ddG calculations.

#### Download and Installation

1. **Get Rosetta License** (free for academic use):
   - Visit: <https://www.rosettacommons.org/software/license-and-download>
   - Register and download the latest version

2. **Install Rosetta**:

   ```bash
   # Extract archive
   tar -xvzf rosetta_[version].tar.gz
   cd rosetta_[version]

   # Compile (this may take 30-60 minutes)
   cd main/source
   ./scons.py -j8 mode=release bin
   ```

3. **Set Environment Variable**:

   ```bash
   # Add to ~/.bashrc or ~/.zshrc
   export ROSETTA=/path/to/rosetta_[version]
   export PATH=$ROSETTA/main/source/bin:$PATH

   # Reload shell configuration
   source ~/.bashrc  # or source ~/.zshrc
   ```

4. **Verify Installation**:

   ```bash
   # Check if flex_ddg is available
   ls $ROSETTA/main/source/bin/flex_ddg*
   ```

### 2. Python Environment

**Required**: Python 3.8 or higher

Check your Python version:

```bash
python --version
# or
python3 --version
```

If you need to install Python:

- **macOS**: `brew install python@3.11`
- **Ubuntu/Debian**: `sudo apt-get install python3.11`
- **Windows**: Download from <https://www.python.org/downloads/>

## Installation Options

### Option 1: Install from Source (Recommended for Development)

```bash
# Clone repository
git clone https://github.com/yourusername/alaninescanning4eif4e.git
cd alaninescanning4eif4e

# Create virtual environment
python -m venv venv

# Activate virtual environment
source venv/bin/activate  # macOS/Linux
# or
venv\Scripts\activate  # Windows

# Install package in editable mode
pip install -e .

# Verify installation
rosetta-scan --version
```

### Option 2: Install from PyPI (When Available)

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install package
pip install rosetta-alanine-scanning

# Verify installation
rosetta-scan --version
```

## Dependency Installation

### Required Dependencies

All Python dependencies are automatically installed with the package:

- `click>=8.1.7` - CLI framework
- `biopython>=1.81` - PDB parsing
- `pandas>=2.0.0` - Data analysis
- `numpy>=1.24.0` - Numerical operations
- `matplotlib>=3.7.0` - Plotting
- `seaborn>=0.12.0` - Statistical visualizations
- `pyyaml>=6.0` - Configuration files
- `rich>=13.0.0` - Terminal formatting

### Optional Dependencies

For PyRosetta integration (advanced users):

```bash
# PyRosetta requires separate installation and licensing
pip install pyrosetta
```

## Verification

### Quick Test

```bash
# Check CLI is working
rosetta-scan --help

# Generate example configuration
rosetta-scan init-config -o test_config.yaml

# Verify Rosetta integration (requires PDB file)
# Download test structure
wget https://files.rcsb.org/download/1ABC.pdb -O test.pdb

# Run scan
rosetta-scan scan test.pdb -c A --output test_mutations/
```

### Run Example Workflow

```bash
cd examples/
chmod +x quick_test.sh
./quick_test.sh
```

## Troubleshooting

### Common Issues

#### 1. Rosetta Not Found

**Error**: `Rosetta path not found. Set ROSETTA environment variable`

**Solution**:

```bash
# Set ROSETTA variable
export ROSETTA=/path/to/rosetta
echo $ROSETTA  # Verify it's set

# Make it permanent by adding to ~/.bashrc
echo 'export ROSETTA=/path/to/rosetta' >> ~/.bashrc
source ~/.bashrc
```

#### 2. Flex ddG Executable Not Found

**Error**: `Rosetta executable not found`

**Solution**:

```bash
# Check if flex_ddg was compiled
ls $ROSETTA/main/source/bin/flex_ddg*

# If not found, compile it
cd $ROSETTA/main/source
./scons.py -j4 mode=release bin/flex_ddg
```

#### 3. Permission Denied

**Error**: `Permission denied` when running executables

**Solution**:

```bash
# Make Rosetta binaries executable
chmod +x $ROSETTA/main/source/bin/*

# Make test scripts executable
chmod +x examples/*.sh
```

#### 4. Import Errors

**Error**: `ModuleNotFoundError: No module named 'rosetta_scan'`

**Solution**:

```bash
# Ensure virtual environment is activated
source venv/bin/activate

# Reinstall package
pip install -e .

# Verify installation
pip show rosetta-alanine-scanning
```

#### 5. BioPython PDB Parsing Warnings

**Warning**: `PDBConstructionWarning: Ignoring unrecognized record 'HETATM'`

**Solution**: These warnings are normal and can be safely ignored. To suppress:

```python
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
```

#### 6. Matplotlib Backend Issues

**Error**: `No module named '_tkinter'`

**Solution**:

```bash
# macOS
brew install python-tk

# Ubuntu/Debian
sudo apt-get install python3-tk

# Or use non-interactive backend
export MPLBACKEND=Agg
```

## Platform-Specific Notes

### macOS

```bash
# Install dependencies via Homebrew
brew install python@3.11

# If using Apple Silicon (M1/M2)
# Rosetta may need to run under Rosetta 2
arch -x86_64 ./scons.py -j4 mode=release bin
```

### Linux (Ubuntu/Debian)

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install -y python3 python3-pip python3-venv
    build-essential zlib1g-dev

# For plotting
sudo apt-get install -y python3-tk
```

### Windows

- Install via WSL (Windows Subsystem for Linux) recommended
- Native Windows support for Rosetta is limited

```powershell
# Install WSL
wsl --install

# Follow Linux installation instructions within WSL
```

## Advanced Configuration

### Using Custom Rosetta Build

If you have a custom Rosetta installation:

```yaml
# In config.yaml
rosetta_path: /path/to/custom/rosetta
rosetta_database: /path/to/custom/rosetta/main/database
```

### MPI Support (Parallel Runs)

For large-scale runs with MPI:

```bash
# Compile Rosetta with MPI
cd $ROSETTA/main/source
./scons.py -j8 mode=release extras=mpi bin

# Configure package to use MPI
# In config.yaml:
num_processors: 8
```

## Updating

### Update Package

```bash
# If installed from source
cd alaninescanning4eif4e
git pull
pip install -e . --upgrade

# If installed from PyPI
pip install --upgrade rosetta-alanine-scanning
```

### Update Dependencies

```bash
pip install --upgrade -r requirements.txt
```

## Uninstallation

```bash
# Deactivate virtual environment
deactivate

# Remove virtual environment
rm -rf venv/

# Uninstall package
pip uninstall rosetta-alanine-scanning

# Remove cloned repository
rm -rf alaninescanning4eif4e/
```

## Getting Help

If you encounter issues:

1. Check the [troubleshooting](#troubleshooting) section above
2. Review the [README](README.md) for usage examples
3. Open an issue on GitHub with:
   - Error message
   - Python version (`python --version`)
   - Operating system
   - Steps to reproduce

## Next Steps

After successful installation:

1. Read the [README](README.md) for usage examples
2. Try the example workflow: `python examples/example_workflow.py`
3. Run a test analysis on your own structure
4. Customize configuration in `config/example_config.yaml`

---

**Need help?** Open an issue on GitHub or email support.
