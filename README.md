# Rosetta Alanine Scanning - Flex ddG Protocol

A Python framework for computational alanine scanning and hotspot identification using Rosetta's Flex ddG protocol.

## Features

- **Automated Alanine Scanning**: Generate systematic alanine mutations from PDB structures
- **Flex ddG Protocol**: Run Rosetta's state-of-the-art Flex ddG calculations
- **Interface Analysis**: Identify interface residues and hotspots in protein-protein interactions
- **Visualizations**: Generate publication-quality plots and heatmaps
- **CLI Interface**: Intuitive command-line interface with progress tracking
- **Complete Pipeline**: End-to-end workflow from structure to results

## Installation

### Prerequisites

1. **Rosetta Software Suite** (required)
   - Download from: https://www.rosettacommons.org/software/license-and-download
   - Compile with `flex_ddg` application
   - Set `ROSETTA` environment variable:
     ```bash
     export ROSETTA=/path/to/rosetta
     ```

2. **Python 3.8+**

### Install Package

```bash
# Clone repository
git clone https://github.com/yourusername/alaninescanning4eif4e.git
cd alaninescanning4eif4e

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install package and dependencies
pip install -e .
```

## Quick Demo

Try the complete workflow with example data (no Rosetta required):

```bash
# Run the interactive demo
python3 examples/demo_run.py

# View the generated outputs
ls examples/demo_output/
```

**[See detailed demo guide](DEMO.md)** with step-by-step instructions and expected outputs.

## Quick Start

### 1. Complete Pipeline (Recommended)

Run the entire workflow with a single command:

```bash
rosetta-scan pipeline protein.pdb --chains A B --interface-only --nstruct 35 --output results/
```

This will:
1. Generate alanine mutations
2. Run Flex ddG calculations
3. Analyze results and identify hotspots
4. Create visualizations

### 2. Step-by-Step Workflow

#### Generate Mutations

```bash
rosetta-scan scan protein.pdb --chains A B --interface-only --output mutations/
```

#### Run Flex ddG

```bash
rosetta-scan run protein.pdb mutations/mutations.txt --nstruct 35 --iterations 3 --output ddg_results/
```

#### Analyze Results

```bash
rosetta-scan analyze ddg_results/ --plot --threshold 1.5 --output analysis.csv
```

## Detailed Usage

### Mutation Generation

```bash
rosetta-scan scan protein.pdb [OPTIONS]

Options:
  -c, --chains TEXT              Chain IDs to scan (e.g., -c A -c B)
  -r, --range TEXT               Residue range (e.g., "A:1-100,B:50-150")
  --interface-only               Only scan interface residues
  --interface-cutoff FLOAT       Distance cutoff for interface (default: 8.0 Å)
  -o, --output PATH              Output directory
  --format [txt|csv|rosetta]     Output format
```

**Examples:**

```bash
# Scan specific chains
rosetta-scan scan protein.pdb -c A -c B

# Scan interface residues only
rosetta-scan scan protein.pdb --interface-only

# Scan specific residue range
rosetta-scan scan protein.pdb -r "A:1-100,B:50-150"
```

### Flex ddG Protocol

```bash
rosetta-scan run pdb_file mutation_file [OPTIONS]

Options:
  --config PATH              YAML configuration file
  --nstruct INTEGER          Number of structures per mutation (default: 35)
  --iterations INTEGER       Backrub iterations (default: 3)
  --interface                Enable interface ddG mode
  -o, --output PATH          Output directory
  --rosetta-path PATH        Path to Rosetta installation
```

**Examples:**

```bash
# Production run
rosetta-scan run protein.pdb mutations.txt --nstruct 35

# Interface ddG
rosetta-scan run complex.pdb mutations.txt --interface --nstruct 35

# Use configuration file
rosetta-scan run protein.pdb mutations.txt --config config.yaml
```

### Result Analysis

```bash
rosetta-scan analyze results_dir [OPTIONS]

Options:
  -o, --output PATH          Output file (CSV/JSON)
  --plot                     Generate visualizations
  --threshold FLOAT          Hotspot threshold (kcal/mol, default: 1.0)
```

**Examples:**

```bash
# Basic analysis
rosetta-scan analyze ddg_results/

# Generate plots and export results
rosetta-scan analyze ddg_results/ --plot --output results.csv

# Custom hotspot threshold
rosetta-scan analyze ddg_results/ --threshold 2.0 --plot
```

## Configuration

Generate example configuration:

```bash
rosetta-scan init-config -o my_config.yaml
```

Edit `my_config.yaml`:

```yaml
# Protocol parameters
nstruct: 35          # Structures per mutation
iterations: 3        # Backrub iterations
repack_radius: 8.0   # Repacking shell (Å)

# Options
use_backrub: true
interface_ddg: false

# Resources
num_processors: 4
memory_gb: 8

# Rosetta paths
rosetta_path: /path/to/rosetta
rosetta_database: /path/to/rosetta/main/database
```

## Output Files

```
results/
--- mutations.txt             # Rosetta mutation file
--- ddg_results/              # Flex ddG output
|   --- score.sc              # Score file
|   --- rosetta.log           # Rosetta log
--- results.csv               # Parsed results
--- analysis_report.txt       # Summary report
--- plots/                    # Visualizations
    --- ddg_distribution.png
    --- hotspot_heatmap.png
    --- chain_analysis.png
    --- top_hotspots.png
```

## Python API

Use the package programmatically:

```python
from rosetta_scan.protocols.flex_ddg import FlexDdGProtocol, FlexDdGConfig
from rosetta_scan.protocols.alanine_scanner import AlanineScan
from rosetta_scan.analysis.parser import ResultParser
from rosetta_scan.analysis.visualizer import ResultVisualizer

# Generate mutations
scanner = AlanineScan("protein.pdb")
mutations = scanner.generate_mutations(
    chains=['A', 'B'],
    interface_only=True
)

# Run Flex ddG
config = FlexDdGConfig(nstruct=35, interface_ddg=True)
protocol = FlexDdGProtocol(config)
protocol.run_flex_ddg(
    pdb_path="protein.pdb",
    mutations=scanner.get_rosetta_mutation_list(),
    output_dir="results/"
)

# Analyze results
parser = ResultParser("results/")
results_df = parser.parse_results()
hotspots = parser.identify_hotspots(threshold=1.5)

# Visualize
visualizer = ResultVisualizer(results_df)
visualizer.create_dashboard("plots/", threshold=1.5)
```

## Examples

### Example 1: Protein-Protein Interface Hotspots

```bash
# Identify hotspots at protein-protein interface
rosetta-scan pipeline complex.pdb --chains A B --interface-only --interface-cutoff 8.0 --nstruct 35 --output interface_hotspots/
```

### Example 2: Folding Stability Analysis

```bash
# Scan entire protein for folding stability
rosetta-scan pipeline protein.pdb --chains A --nstruct 35 --output stability_analysis/
```

### Example 3: Custom Residue Range

```bash
# Scan specific region
rosetta-scan scan protein.pdb -c A -r "A:100-200" --output region_scan/

rosetta-scan run protein.pdb region_scan/mutations.txt --nstruct 35 --output region_ddg/
```

## Visualization Examples

The analysis generates several types of plots:

1. **ΔΔG Distribution**: Histogram of all ΔΔG values
2. **Hotspot Heatmap**: Position-based heatmap of mutations
3. **Chain Analysis**: Box/violin plots per chain
4. **Top Hotspots**: Bar chart of strongest hotspots
5. **Position Scan**: ΔΔG along sequence

## Citation

If you use this tool, please cite:

**Flex ddG Protocol:**
> Barlow KA, et al. (2018) Flex ddG: Rosetta Ensemble-Based Estimation of Changes in Protein-Protein Binding Affinity upon Mutation. J Phys Chem B. 122(21):5389-5399.

**Rosetta:**
> Alford RF, et al. (2017) The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design. J Chem Theory Comput. 13(6):3031-3048.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file.

## cknowledgments

- Rosetta Commons for the Rosetta software suite
- Barlow et al. for the Flex ddG protocol

## Support

For questions or issues:
- Open an issue on GitHub
- Email: madsondeluna@example.com

## Resources

- [Rosetta Commons](https://www.rosettacommons.org/)
- [Flex ddG Documentation](https://www.rosettacommons.org/docs/latest/application_documentation/analysis/ddg-monomer)
- [PyRosetta](http://www.pyrosetta.org/)

---

Made with  for computational protein engineering
