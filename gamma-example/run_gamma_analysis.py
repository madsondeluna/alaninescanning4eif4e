#!/usr/bin/env python3
"""
Alanine Scanning Analysis for Gamma Peptide
=============================================

This script performs a complete alanine scanning analysis on the gamma peptide.
It will:
1. Load the PDB structure
2. Generate all possible alanine mutations
3. Simulate ddG values (since Rosetta is not required for demo)
4. Identify hotspot residues
5. Generate visualizations
6. Create a comprehensive analysis report
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from rosetta_scan.protocols.alanine_scanner import AlanineScan
from rosetta_scan.analysis.visualizer import ResultVisualizer
import pandas as pd
import numpy as np
import csv
from datetime import datetime

# Set random seed for reproducibility
np.random.seed(42)

def main():
    print('='*80)
    print('ALANINE SCANNING ANALYSIS - GAMMA PEPTIDE')
    print('='*80)
    print(f'Analysis started: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    print()

    # Setup paths
    pdb_path = Path(__file__).parent / 'gamma.pdb'
    output_dir = Path(__file__).parent / 'analysis_results'
    output_dir.mkdir(exist_ok=True)
    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)

    # ========================================================================
    # STEP 1: Structure Analysis
    # ========================================================================
    print('STEP 1: Structure Analysis')
    print('-' * 80)

    scanner = AlanineScan(str(pdb_path))

    print(f'✓ PDB file loaded: gamma.pdb')
    print(f'✓ Chains found: {[chain.id for chain in scanner.structure[0]]}')

    # Get sequence info
    chain_a = scanner.structure[0]['A']
    residues = [res for res in chain_a if res.id[0] == ' ']
    sequence = ''.join([res.resname for res in residues])

    print(f'✓ Chain A contains {len(residues)} residues')
    print(f'✓ Sequence: {" -> ".join([res.resname for res in residues])}')
    print()

    # ========================================================================
    # STEP 2: Generate Alanine Mutations
    # ========================================================================
    print('STEP 2: Generate Alanine Mutations')
    print('-' * 80)
    print('Parameters:')
    print('  • Target: Complete peptide (all residues)')
    print('  • Chain: A')
    print('  • Excluded residues: GLY (already glycine), ALA (already alanine)')
    print()

    mutations = scanner.generate_mutations(
        chains=['A'],
        interface_only=False  # Scan ALL residues
    )

    print(f'✓ Generated {len(mutations)} mutations')
    print()

    # Mutation breakdown
    print('Mutation breakdown by amino acid:')
    aa_count = {}
    for mut in mutations:
        aa = mut.original_aa
        aa_count[aa] = aa_count.get(aa, 0) + 1

    for aa, count in sorted(aa_count.items()):
        aa_full = {
            'ARG': 'Arginine',
            'CYS': 'Cysteine',
            'PHE': 'Phenylalanine',
            'GLY': 'Glycine'
        }.get(aa, aa)
        print(f'  • {aa} ({aa_full}): {count} mutation(s)')
    print()

    # Show all mutations
    print('Complete mutation list:')
    print(f'{'#':<4} {'Mutation':<12} {'Position':<10} {'Original':<12} {'Target'}')
    print('-' * 60)
    for i, mut in enumerate(mutations, 1):
        print(f'{i:<4} {mut.chain}{mut.position}{mut.original_aa:<10} {mut.position:<10} {mut.original_aa:<12} A')
    print()

    # ========================================================================
    # STEP 3: Save Mutation Files
    # ========================================================================
    print('STEP 3: Save Mutation Files')
    print('-' * 80)

    # Text format
    scanner.save_mutation_report(output_dir / 'mutations.txt', format='txt')
    print(f'✓ Saved: mutations.txt (human-readable)')

    # Rosetta format
    rosetta_muts = scanner.get_rosetta_mutation_list()
    with open(output_dir / 'mutations_rosetta.txt', 'w') as f:
        f.write('total 1\n')
        for mut in rosetta_muts:
            f.write(f'1\n{mut}\n')
    print(f'✓ Saved: mutations_rosetta.txt (Rosetta format)')

    # CSV format - use scanner method
    scanner.save_mutation_report(output_dir / 'mutations.csv', format='csv')
    print(f'✓ Saved: mutations.csv (spreadsheet format)')
    print()

    # ========================================================================
    # STEP 4: Simulate ddG Results
    # ========================================================================
    print('STEP 4: Simulate ddG Results')
    print('-' * 80)
    print('NOTE: This is a SIMULATION for demonstration purposes.')
    print('      For real ddG values, run with Rosetta Flex ddG protocol.')
    print()

    # Simulate realistic ddG values based on amino acid properties
    results = []
    for mut in mutations:
        # Different amino acids have different contributions
        # ARG: charged, often important -> higher ddG
        # CYS: disulfide bonds, structural -> medium-high ddG
        # PHE: aromatic, hydrophobic -> medium ddG

        base_ddg = {
            'ARG': np.random.uniform(1.2, 2.8),  # Charged residues often critical
            'CYS': np.random.uniform(0.8, 2.2),  # Cysteine can form disulfides
            'PHE': np.random.uniform(0.5, 1.8),  # Aromatic, hydrophobic
        }.get(mut.original_aa, np.random.uniform(0.3, 1.0))

        # Add some position-dependent variation
        # Residues in the middle might be more critical
        if 4 <= mut.position <= 9:
            base_ddg *= np.random.uniform(1.0, 1.3)

        # Add noise
        ddg = base_ddg + np.random.normal(0, 0.15)

        # Total score (more negative is more stable)
        total_score = -50.0 + ddg + np.random.normal(0, 0.5)

        results.append({
            'mutation': f'{mut.chain}{mut.position}{mut.original_aa}',
            'chain': mut.chain,
            'position': mut.position,
            'original_aa': mut.original_aa,
            'ddg': ddg,
            'total_score': total_score
        })

    # Create DataFrame
    df = pd.DataFrame(results)
    df = df.sort_values('ddg', ascending=False)

    print(f'✓ Simulated ddG values for {len(df)} mutations')
    print()
    print('Statistics:')
    print(f'  • Mean ΔΔG: {df["ddg"].mean():.2f} kcal/mol')
    print(f'  • Std ΔΔG: {df["ddg"].std():.2f} kcal/mol')
    print(f'  • Min ΔΔG: {df["ddg"].min():.2f} kcal/mol')
    print(f'  • Max ΔΔG: {df["ddg"].max():.2f} kcal/mol')
    print()

    # Save results
    df.to_csv(output_dir / 'ddg_results.csv', index=False)
    print(f'✓ Saved: ddg_results.csv')
    print()

    # ========================================================================
    # STEP 5: Identify Hotspots
    # ========================================================================
    print('STEP 5: Identify Hotspot Residues')
    print('-' * 80)

    hotspot_threshold = 1.5  # kcal/mol
    print(f'Hotspot threshold: ΔΔG > {hotspot_threshold} kcal/mol')
    print()

    hotspots = df[df['ddg'] > hotspot_threshold].sort_values('ddg', ascending=False)

    print(f'✓ Identified {len(hotspots)} hotspot residues')
    print()

    if len(hotspots) > 0:
        print('HOTSPOT RANKING:')
        print(f'{'Rank':<6} {'Mutation':<12} {'Position':<10} {'Residue':<10} {'ΔΔG (kcal/mol)':<16} {'Classification'}')
        print('-' * 80)

        for i, (idx, row) in enumerate(hotspots.iterrows(), 1):
            classification = ('CRITICAL' if row['ddg'] > 2.5 else
                            'HIGH' if row['ddg'] > 2.0 else
                            'MEDIUM')
            print(f'{i:<6} {row["mutation"]:<12} {row["position"]:<10} {row["original_aa"]:<10} {row["ddg"]:<16.2f} {classification}')
        print()

        # Save hotspots
        hotspots.to_csv(output_dir / 'hotspots.csv', index=False)
        print(f'✓ Saved: hotspots.csv')
    else:
        print('No hotspots identified above threshold.')
    print()

    # ========================================================================
    # STEP 6: Generate Visualizations
    # ========================================================================
    print('STEP 6: Generate Visualizations')
    print('-' * 80)

    try:
        visualizer = ResultVisualizer(df)

        # 1. ddG distribution
        visualizer.plot_ddg_distribution(
            output_path=str(plots_dir / 'ddg_distribution.png'),
            threshold=hotspot_threshold
        )
        print('✓ Created: ddg_distribution.png')

        # 2. Top hotspots
        visualizer.plot_top_hotspots(
            output_path=str(plots_dir / 'top_hotspots.png'),
            top_n=10
        )
        print('✓ Created: top_hotspots.png')

        # 3. Position scan
        visualizer.plot_position_scan(
            chain='A',
            output_path=str(plots_dir / 'position_scan_chain_A.png'),
            threshold=hotspot_threshold
        )
        print('✓ Created: position_scan_chain_A.png')

        # 4. Hotspot heatmap
        visualizer.plot_hotspot_heatmap(
            output_path=str(plots_dir / 'hotspot_heatmap.png')
        )
        print('✓ Created: hotspot_heatmap.png')

        print()
        print(f'✓ All plots saved to: {plots_dir}/')

    except Exception as e:
        print(f'⚠ Warning: Could not generate some visualizations: {e}')

    print()

    # ========================================================================
    # STEP 7: Generate Analysis Report
    # ========================================================================
    print('STEP 7: Generate Analysis Report')
    print('-' * 80)

    report_path = output_dir / 'ANALYSIS_REPORT.md'

    with open(report_path, 'w') as f:
        f.write(f"""# Alanine Scanning Analysis Report - Gamma Peptide

**Date:** {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
**PDB File:** gamma.pdb
**Analysis Type:** Complete peptide alanine scanning

---

## 1. Structure Information

- **Chain:** A
- **Total Residues:** {len(residues)}
- **Sequence:** {' → '.join([res.resname for res in residues])}

### Peptide Composition

""")

        # Amino acid composition
        aa_composition = {}
        for res in residues:
            aa_composition[res.resname] = aa_composition.get(res.resname, 0) + 1

        for aa, count in sorted(aa_composition.items()):
            f.write(f"- **{aa}**: {count} residue(s)\n")

        f.write(f"""

---

## 2. Alanine Scanning Parameters

- **Target:** Complete peptide (all scannable residues)
- **Chain Analyzed:** A
- **Total Mutations Generated:** {len(mutations)}
- **Excluded Residues:** GLY (already glycine)

### Mutation Breakdown

""")

        for aa, count in sorted(aa_count.items()):
            f.write(f"- **{aa}**: {count} mutation(s)\n")

        f.write(f"""

### Complete Mutation List

| # | Mutation | Position | Original AA | Target AA |
|---|----------|----------|-------------|-----------|
""")

        for i, mut in enumerate(mutations, 1):
            f.write(f"| {i} | {mut.chain}{mut.position}{mut.original_aa} | {mut.position} | {mut.original_aa} | A |\n")

        f.write(f"""

---

## 3. ddG Results (Simulated)

**NOTE:** These are SIMULATED values for demonstration. For real results, run with Rosetta Flex ddG.

### Statistics

- **Number of Mutations:** {len(df)}
- **Mean ΔΔG:** {df['ddg'].mean():.2f} kcal/mol
- **Std ΔΔG:** {df['ddg'].std():.2f} kcal/mol
- **Min ΔΔG:** {df['ddg'].min():.2f} kcal/mol
- **Max ΔΔG:** {df['ddg'].max():.2f} kcal/mol

### All Results

| Rank | Mutation | Position | Original AA | ΔΔG (kcal/mol) |
|------|----------|----------|-------------|----------------|
""")

        for i, (idx, row) in enumerate(df.iterrows(), 1):
            f.write(f"| {i} | {row['mutation']} | {row['position']} | {row['original_aa']} | {row['ddg']:.2f} |\n")

        f.write(f"""

---

## 4. Hotspot Analysis

**Hotspot Threshold:** ΔΔG > {hotspot_threshold} kcal/mol

**Number of Hotspots:** {len(hotspots)}

""")

        if len(hotspots) > 0:
            f.write("""
### Hotspot Ranking

| Rank | Mutation | Position | Residue | ΔΔG (kcal/mol) | Classification |
|------|----------|----------|---------|----------------|----------------|
""")

            for i, (idx, row) in enumerate(hotspots.iterrows(), 1):
                classification = ('🔴 CRITICAL' if row['ddg'] > 2.5 else
                                '🟠 HIGH' if row['ddg'] > 2.0 else
                                '🟡 MEDIUM')
                f.write(f"| {i} | {row['mutation']} | {row['position']} | {row['original_aa']} | {row['ddg']:.2f} | {classification} |\n")

            f.write("""

### Hotspot Interpretation

""")

            for i, (idx, row) in enumerate(hotspots.iterrows(), 1):
                classification = ('CRITICAL' if row['ddg'] > 2.5 else
                                'HIGH' if row['ddg'] > 2.0 else
                                'MEDIUM')

                f.write(f"""
#### Hotspot {i}: {row['mutation']} (ΔΔG = {row['ddg']:.2f} kcal/mol)

- **Classification:** {classification}
- **Position:** {row['position']}
- **Original Residue:** {row['original_aa']}
- **Interpretation:** Mutation of {row['original_aa']} to Alanine at position {row['position']} causes a significant loss of {row['ddg']:.2f} kcal/mol in stability/binding energy.
- **Recommendation:** This residue is important for peptide function and should be conserved.

""")
        else:
            f.write("No significant hotspots identified.\n")

        f.write(f"""

---

## 5. Biological Insights

### Peptide Characteristics

The gamma peptide is a **{len(residues)}-residue peptide** with the following notable features:

1. **High Arginine Content ({aa_composition.get('ARG', 0)} residues):** Suggests strong positive charge, potentially important for:
   - Electrostatic interactions
   - DNA/RNA binding
   - Membrane interactions

2. **Cysteine Residues ({aa_composition.get('CYS', 0)} residues):** May form disulfide bonds:
   - Positions: {', '.join([str(mut.position) for mut in mutations if mut.original_aa == 'CYS'])}
   - Potentially critical for structural stability

3. **Aromatic Residues ({aa_composition.get('PHE', 0)} PHE):** Contribute to:
   - Hydrophobic core
   - π-π stacking interactions
   - Binding interfaces

### Functional Predictions

Based on the alanine scanning results, the most critical residues for peptide function appear to be:

""")

        if len(hotspots) > 0:
            for i, (idx, row) in enumerate(hotspots.head(3).iterrows(), 1):
                f.write(f"{i}. **Position {row['position']} ({row['original_aa']})** - ΔΔG = {row['ddg']:.2f} kcal/mol\n")

        f.write(f"""

---

## 6. Recommendations

### For Experimental Validation

1. **Top Priority Mutations:**
""")

        if len(hotspots) > 0:
            for i, (idx, row) in enumerate(hotspots.head(5).iterrows(), 1):
                f.write(f"   - {row['mutation']} (ΔΔG = {row['ddg']:.2f} kcal/mol)\n")

        f.write(f"""

2. **Mutagenesis Strategy:**
   - Use site-directed mutagenesis
   - Test single mutations first
   - Measure binding affinity (if applicable)
   - Assess structural stability (CD, thermal shift)

3. **Controls:**
   - Wild-type peptide
   - Non-hotspot mutations (negative controls)

### For Computational Follow-up

1. **Molecular Dynamics Simulations:**
   - Simulate wild-type and top mutants
   - Assess conformational changes
   - Calculate binding free energies

2. **Rosetta Flex ddG:**
   - Run protocol with actual Rosetta software
   - Use nstruct=35 for statistical rigor
   - Compare with these simulated results

---

## 7. Output Files

All analysis files are saved in `gamma-example/analysis_results/`:

### Mutation Files
- `mutations.txt` - Human-readable mutation list
- `mutations_rosetta.txt` - Rosetta input format
- `mutations.csv` - Spreadsheet format

### Results Files
- `ddg_results.csv` - Complete ddG values
- `hotspots.csv` - Filtered hotspot residues
- `ANALYSIS_REPORT.md` - This report

### Visualizations (`plots/` subdirectory)
- `ddg_distribution.png` - ΔΔG histogram
- `top_hotspots.png` - Bar chart of top mutations
- `position_scan_chain_A.png` - ΔΔG along sequence
- `hotspot_heatmap.png` - 2D heatmap

---

## 8. Notes

- This analysis used **SIMULATED** ddG values for demonstration
- For publication-quality results, use Rosetta Flex ddG protocol
- Recommended Rosetta parameters:
  - `nstruct`: 35 (minimum for statistical significance)
  - `iterations`: 3 (backrub iterations)
  - `score function`: REF2015 (recommended)

---

**Analysis completed:** {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

""")

    print(f'✓ Comprehensive report saved: ANALYSIS_REPORT.md')
    print()

    # ========================================================================
    # Summary
    # ========================================================================
    print('='*80)
    print('ANALYSIS COMPLETE')
    print('='*80)
    print()
    print(f'✓ Total mutations analyzed: {len(mutations)}')
    print(f'✓ Hotspots identified: {len(hotspots)}')
    print(f'✓ Output directory: gamma-example/analysis_results/')
    print()
    print('Files created:')
    print('  Mutations:')
    print('    • mutations.txt')
    print('    • mutations_rosetta.txt')
    print('    • mutations.csv')
    print('  Results:')
    print('    • ddg_results.csv')
    print('    • hotspots.csv')
    print('  Visualizations:')
    print('    • plots/ddg_distribution.png')
    print('    • plots/top_hotspots.png')
    print('    • plots/position_scan_chain_A.png')
    print('    • plots/hotspot_heatmap.png')
    print('  Documentation:')
    print('    • ANALYSIS_REPORT.md')
    print()
    print('='*80)
    print()
    print('Next steps:')
    print('  1. Review ANALYSIS_REPORT.md for detailed results')
    print('  2. Examine plots/ directory for visualizations')
    print('  3. Check hotspots.csv for critical residues')
    print('  4. For real ddG values, run with Rosetta:')
    print('     rosetta-scan run gamma-example/gamma.pdb \\')
    print('       gamma-example/analysis_results/mutations_rosetta.txt \\')
    print('       --nstruct 35 --output gamma_ddg_results/')
    print()
    print('='*80)

if __name__ == '__main__':
    main()
