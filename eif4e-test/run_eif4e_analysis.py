#!/usr/bin/env python3
"""
Complete Alanine Scanning Analysis for eIF4E Protein
Analyzes all scannable residues and generates comprehensive results
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Add parent directory to path for imports
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))
sys.path.insert(0, str(parent_dir / 'src'))

# Import Bio before rosetta_scan to ensure it's available
from Bio.PDB import PDBParser, is_aa
try:
    from Bio.PDB.Polypeptide import three_to_one
except ImportError:
    from Bio.Data.PDBData import protein_letters_3to1
    def three_to_one(three_letter_code):
        return protein_letters_3to1.get(three_letter_code, 'X')

# Now import without triggering full __init__
import importlib.util

# Load alanine_scanner module directly
scanner_path = parent_dir / 'src' / 'rosetta_scan' / 'protocols' / 'alanine_scanner.py'
spec = importlib.util.spec_from_file_location("alanine_scanner", scanner_path)
alanine_scanner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(alanine_scanner)
AlanineScan = alanine_scanner.AlanineScan

# Load visualizer module directly
visualizer_path = parent_dir / 'src' / 'rosetta_scan' / 'analysis' / 'visualizer.py'
spec = importlib.util.spec_from_file_location("visualizer", visualizer_path)
visualizer_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(visualizer_module)
ResultVisualizer = visualizer_module.ResultVisualizer


def main():
    """
    Main analysis pipeline for eIF4E protein
    """
    print("\n" + "="*80)
    print("ALANINE SCANNING ANALYSIS - eIF4E PROTEIN")
    print("Complete Protein Analysis (All Scannable Residues)")
    print("="*80 + "\n")

    # Define paths
    script_dir = Path(__file__).parent
    pdb_path = script_dir / 'model.pdb'
    output_dir = script_dir / 'analysis_results'
    output_dir.mkdir(exist_ok=True)

    plots_dir = output_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)

    # ========================================================================
    # STEP 1: Load and Analyze Structure
    # ========================================================================
    print('STEP 1: Load and Analyze Structure')
    print('-' * 80)

    scanner = AlanineScan(str(pdb_path))
    print(f'PDB file loaded: {pdb_path.name}')

    # Get structure info
    structure = scanner.structure
    model = structure[0]
    chains = list(model.get_chains())

    print(f'Number of chains: {len(chains)}')
    for chain in chains:
        residues = list(chain.get_residues())
        print(f'  Chain {chain.id}: {len(residues)} residues')

    print()

    # ========================================================================
    # STEP 2: Generate Mutations for Complete Protein
    # ========================================================================
    print('STEP 2: Generate Mutations for Complete Protein')
    print('-' * 80)

    # Scan all residues in chain A
    mutations = scanner.generate_mutations(
        chains=['A'],
        interface_only=False,  # Scan complete protein
        interface_cutoff=8.0
    )

    print(f'Total mutations generated: {len(mutations)}')

    # Get summary
    summary = scanner.get_mutation_summary()
    print(f'\nMutation breakdown by amino acid:')
    for aa, count in sorted(summary['by_aa'].items()):
        print(f'  {aa}: {count} mutations')

    print()

    # ========================================================================
    # STEP 3: Save Mutation Files
    # ========================================================================
    print('STEP 3: Save Mutation Files')
    print('-' * 80)

    # Text format
    scanner.save_mutation_report(output_dir / 'mutations.txt', format='txt')
    print(f'Saved: mutations.txt (human-readable)')

    # Rosetta format
    rosetta_muts = scanner.get_rosetta_mutation_list()
    with open(output_dir / 'mutations_rosetta.txt', 'w') as f:
        f.write('total 1\n')
        for mut in rosetta_muts:
            f.write(f'1\n{mut}\n')
    print(f'Saved: mutations_rosetta.txt (Rosetta format)')

    # CSV format
    scanner.save_mutation_report(output_dir / 'mutations.csv', format='csv')
    print(f'Saved: mutations.csv (spreadsheet format)')
    print()

    # ========================================================================
    # STEP 4: Simulate ddG Results
    # ========================================================================
    print('STEP 4: Simulate ddG Results')
    print('-' * 80)
    print('NOTE: Generating simulated ddG values for demonstration')
    print('For real results, run with Rosetta Flex ddG protocol')
    print()

    # Simulate ddG values based on amino acid properties
    np.random.seed(42)

    ddg_data = []
    for mut in mutations:
        # Simulate ddG based on amino acid properties
        aa = mut.original_aa

        # Different amino acids have different contributions
        base_ddg = {
            'W': 2.5, 'Y': 2.2, 'F': 2.0,  # Aromatic (high)
            'R': 1.8, 'K': 1.6, 'E': 1.5, 'D': 1.4,  # Charged (medium-high)
            'L': 1.2, 'I': 1.1, 'V': 1.0, 'M': 1.3,  # Hydrophobic (medium)
            'H': 1.4, 'N': 0.9, 'Q': 0.9,  # Polar (medium-low)
            'S': 0.7, 'T': 0.8, 'C': 1.0,  # Small polar (low)
        }.get(aa, 0.5)

        # Add position-based variation
        position_factor = 0.3 * np.sin(mut.position / 20.0) + 1.0

        # Add random noise
        noise = np.random.normal(0, 0.3)

        ddg = max(0.1, base_ddg * position_factor + noise)

        # Simulate total score (more negative = better)
        total_score = -50.0 + ddg + np.random.normal(0, 1.5)

        ddg_data.append({
            'mutation': f'A{mut.position}{mut.original_aa}',
            'chain': mut.chain,
            'position': mut.position,
            'original_aa': mut.original_aa,
            'ddg': ddg,
            'total_score': total_score
        })

    # Create DataFrame
    ddg_df = pd.DataFrame(ddg_data)
    ddg_df = ddg_df.sort_values('ddg', ascending=False)

    # Save results
    ddg_df.to_csv(output_dir / 'ddg_results.csv', index=False)
    print(f'Saved: ddg_results.csv ({len(ddg_df)} mutations)')

    # Print statistics
    print(f'\nddG Statistics:')
    print(f'  Mean:   {ddg_df["ddg"].mean():.3f} kcal/mol')
    print(f'  Median: {ddg_df["ddg"].median():.3f} kcal/mol')
    print(f'  Std:    {ddg_df["ddg"].std():.3f} kcal/mol')
    print(f'  Min:    {ddg_df["ddg"].min():.3f} kcal/mol')
    print(f'  Max:    {ddg_df["ddg"].max():.3f} kcal/mol')
    print()

    # ========================================================================
    # STEP 5: Identify Hotspots
    # ========================================================================
    print('STEP 5: Identify Hotspots')
    print('-' * 80)

    threshold = 2.0  # kcal/mol
    hotspots = ddg_df[ddg_df['ddg'] > threshold].copy()

    print(f'Hotspot threshold: {threshold} kcal/mol')
    print(f'Number of hotspots identified: {len(hotspots)}')

    if len(hotspots) > 0:
        print(f'\nTop 10 hotspots:')
        for idx, row in hotspots.head(10).iterrows():
            print(f'  {row["mutation"]:8s} | {row["original_aa"]}{row["position"]:<4d} | '
                  f'ddG = {row["ddg"]:5.2f} kcal/mol | Score = {row["total_score"]:6.2f} REU')

        # Save hotspots
        hotspots.to_csv(output_dir / 'hotspots.csv', index=False)
        print(f'\nSaved: hotspots.csv ({len(hotspots)} hotspots)')
    else:
        print('No hotspots above threshold')
        # Save empty file
        hotspots.to_csv(output_dir / 'hotspots.csv', index=False)

    print()

    # ========================================================================
    # STEP 6: Generate Visualizations
    # ========================================================================
    print('STEP 6: Generate Visualizations')
    print('-' * 80)

    visualizer = ResultVisualizer(ddg_df)

    # 1. ddG distribution
    visualizer.plot_ddg_distribution(
        output_path=plots_dir / 'ddg_distribution.png'
    )
    print('Saved: ddg_distribution.png')

    # 2. Top hotspots bar chart
    visualizer.plot_top_hotspots(
        output_path=plots_dir / 'top_hotspots.png',
        top_n=20
    )
    print('Saved: top_hotspots.png (top 20)')

    # 3. Position scan (entire protein)
    visualizer.plot_position_scan(
        output_path=plots_dir / 'position_scan_chain_A.png',
        chain='A'
    )
    print('Saved: position_scan_chain_A.png')

    # 4. Heatmap
    visualizer.plot_hotspot_heatmap(
        output_path=plots_dir / 'hotspot_heatmap.png'
    )
    print('Saved: hotspot_heatmap.png')

    print()

    # ========================================================================
    # STEP 7: Generate Analysis Report
    # ========================================================================
    print('STEP 7: Generate Analysis Report')
    print('-' * 80)

    report_path = output_dir / 'ANALYSIS_REPORT.md'

    with open(report_path, 'w') as f:
        f.write('# Alanine Scanning Analysis Report - eIF4E Protein\n\n')
        f.write(f'**Date:** {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write(f'**PDB File:** {pdb_path.name}\n')
        f.write(f'**Analysis Type:** Complete protein alanine scanning\n\n')
        f.write('---\n\n')

        f.write('## 1. Structure Information\n\n')
        f.write(f'- **Chain:** A\n')
        f.write(f'- **Total Residues:** 157\n')
        f.write(f'- **Total Mutations Generated:** {len(mutations)}\n\n')

        f.write('### Mutation Breakdown\n\n')
        for aa, count in sorted(summary['by_aa'].items()):
            f.write(f'- **{aa}**: {count} mutation(s)\n')

        f.write('\n---\n\n')
        f.write('## 2. ddG Results (Simulated)\n\n')
        f.write('**NOTE:** These are SIMULATED values for demonstration. ')
        f.write('For real results, run with Rosetta Flex ddG.\n\n')

        f.write('### Statistics\n\n')
        f.write(f'- **Number of Mutations:** {len(ddg_df)}\n')
        f.write(f'- **Mean ΔΔG:** {ddg_df["ddg"].mean():.2f} kcal/mol\n')
        f.write(f'- **Std ΔΔG:** {ddg_df["ddg"].std():.2f} kcal/mol\n')
        f.write(f'- **Min ΔΔG:** {ddg_df["ddg"].min():.2f} kcal/mol\n')
        f.write(f'- **Max ΔΔG:** {ddg_df["ddg"].max():.2f} kcal/mol\n\n')

        f.write('### All Results\n\n')
        f.write('| Rank | Mutation | Position | Original AA | ΔΔG (kcal/mol) |\n')
        f.write('|------|----------|----------|-------------|----------------|\n')
        for i, (idx, row) in enumerate(ddg_df.head(50).iterrows(), 1):
            f.write(f'| {i} | {row["mutation"]} | {row["position"]} | '
                   f'{row["original_aa"]} | {row["ddg"]:.2f} |\n')

        if len(ddg_df) > 50:
            f.write(f'\n*Showing top 50 of {len(ddg_df)} mutations*\n')

        f.write('\n---\n\n')
        f.write('## 3. Hotspot Analysis\n\n')
        f.write(f'**Hotspot Threshold:** ΔΔG > {threshold} kcal/mol\n\n')
        f.write(f'**Number of Hotspots:** {len(hotspots)}\n\n')

        if len(hotspots) > 0:
            f.write('### Identified Hotspots\n\n')
            f.write('| Rank | Mutation | Position | AA | ΔΔG (kcal/mol) | Total Score (REU) |\n')
            f.write('|------|----------|----------|----|----------------|-------------------|\n')
            for i, (idx, row) in enumerate(hotspots.iterrows(), 1):
                f.write(f'| {i} | {row["mutation"]} | {row["position"]} | '
                       f'{row["original_aa"]} | {row["ddg"]:.2f} | {row["total_score"]:.2f} |\n')
        else:
            f.write('No significant hotspots identified.\n')

        f.write('\n---\n\n')
        f.write('## 4. Output Files\n\n')
        f.write('### Mutation Files\n')
        f.write('- `mutations.txt` - Human-readable mutation list\n')
        f.write('- `mutations_rosetta.txt` - Rosetta input format\n')
        f.write('- `mutations.csv` - Spreadsheet format\n\n')
        f.write('### Results Files\n')
        f.write('- `ddg_results.csv` - Complete ddG values\n')
        f.write('- `hotspots.csv` - Filtered hotspot residues\n')
        f.write('- `ANALYSIS_REPORT.md` - This report\n\n')
        f.write('### Visualizations (`plots/` subdirectory)\n')
        f.write('- `ddg_distribution.png` - ΔΔG histogram\n')
        f.write('- `top_hotspots.png` - Bar chart of top mutations\n')
        f.write('- `position_scan_chain_A.png` - ΔΔG along sequence\n')
        f.write('- `hotspot_heatmap.png` - 2D heatmap\n\n')
        f.write('---\n\n')
        f.write(f'**Analysis completed:** {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')

    print(f'Saved: ANALYSIS_REPORT.md')
    print()

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print('='*80)
    print('ANALYSIS COMPLETE!')
    print('='*80)
    print(f'\nResults saved to: {output_dir}')
    print(f'\nKey findings:')
    print(f'  - Total mutations: {len(ddg_df)}')
    print(f'  - Hotspots (ΔΔG > {threshold}): {len(hotspots)}')
    print(f'  - Mean ΔΔG: {ddg_df["ddg"].mean():.2f} kcal/mol')

    if len(hotspots) > 0:
        top = hotspots.iloc[0]
        print(f'  - Top hotspot: {top["mutation"]} (ΔΔG = {top["ddg"]:.2f} kcal/mol)')

    print(f'\nNext steps:')
    print(f'  1. Review analysis results in {output_dir}')
    print(f'  2. Check visualizations in {plots_dir}')
    print(f'  3. Generate heatmap: python3 generate_ddg_heatmap.py')
    print()


if __name__ == '__main__':
    main()
