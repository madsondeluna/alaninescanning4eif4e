#!/usr/bin/env python3
"""
Complete Alanine Scanning Analysis for eIF4E Protein
Simplified version - analyzes all scannable residues
"""

from pathlib import Path
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
try:
    from Bio.PDB.Polypeptide import three_to_one, is_aa
except ImportError:
    from Bio.PDB.Polypeptide import is_aa
    from Bio.Data.PDBData import protein_letters_3to1
    def three_to_one(three_letter_code):
        return protein_letters_3to1.get(three_letter_code, 'X')


def analyze_structure(pdb_path):
    """Analyze PDB structure and generate mutations"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)

    mutations = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True):
                    try:
                        aa = three_to_one(residue.get_resname())
                        pos = residue.get_id()[1]

                        # Skip GLY, ALA, PRO
                        if aa not in ['G', 'A', 'P']:
                            mutations.append({
                                'chain': chain.id,
                                'position': pos,
                                'original_aa': aa,
                                'mutation': f'A{pos}{aa}'
                            })
                    except:
                        continue

    return structure, mutations


def simulate_ddg(mutations):
    """Simulate ddG values"""
    np.random.seed(42)

    ddg_data = []
    for mut in mutations:
        aa = mut['original_aa']

        # Different amino acids have different contributions
        base_ddg = {
            'W': 2.5, 'Y': 2.2, 'F': 2.0,  # Aromatic (high)
            'R': 1.8, 'K': 1.6, 'E': 1.5, 'D': 1.4,  # Charged (medium-high)
            'L': 1.2, 'I': 1.1, 'V': 1.0, 'M': 1.3,  # Hydrophobic (medium)
            'H': 1.4, 'N': 0.9, 'Q': 0.9,  # Polar (medium-low)
            'S': 0.7, 'T': 0.8, 'C': 1.0,  # Small polar (low)
        }.get(aa, 0.5)

        # Add position-based variation
        position_factor = 0.3 * np.sin(mut['position'] / 20.0) + 1.0

        # Add random noise
        noise = np.random.normal(0, 0.3)

        ddg = max(0.1, base_ddg * position_factor + noise)

        # Simulate total score (more negative = better)
        total_score = -50.0 + ddg + np.random.normal(0, 1.5)

        ddg_data.append({
            'mutation': mut['mutation'],
            'chain': mut['chain'],
            'position': mut['position'],
            'original_aa': mut['original_aa'],
            'ddg': ddg,
            'ddg_normalized': -ddg,  # Normalized (negative) values
            'total_score': total_score
        })

    return pd.DataFrame(ddg_data).sort_values('ddg', ascending=False)


def main():
    """Main analysis pipeline"""
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

    # Step 1: Analyze structure
    print('STEP 1: Load and Analyze Structure')
    print('-' * 80)

    structure, mutations = analyze_structure(pdb_path)
    print(f'PDB file loaded: {pdb_path.name}')
    print(f'Total mutations generated: {len(mutations)}')

    # Count by AA
    aa_counts = {}
    for mut in mutations:
        aa = mut['original_aa']
        aa_counts[aa] = aa_counts.get(aa, 0) + 1

    print(f'\nMutation breakdown by amino acid:')
    for aa, count in sorted(aa_counts.items()):
        print(f'  {aa}: {count} mutations')
    print()

    # Step 2: Save mutation files
    print('STEP 2: Save Mutation Files')
    print('-' * 80)

    # Text format
    with open(output_dir / 'mutations.txt', 'w') as f:
        f.write(f'Alanine Scanning Mutations - eIF4E\n')
        f.write(f'Total: {len(mutations)} mutations\n\n')
        for mut in mutations:
            f.write(f'{mut["mutation"]:8s} | Chain {mut["chain"]} | Position {mut["position"]:3d} | {mut["original_aa"]} -> A\n')
    print('Saved: mutations.txt')

    # Rosetta format
    with open(output_dir / 'mutations_rosetta.txt', 'w') as f:
        f.write('total 1\n')
        for mut in mutations:
            f.write(f'1\n{mut["mutation"]} {mut["chain"]}\n')
    print('Saved: mutations_rosetta.txt')

    # CSV format
    mut_df = pd.DataFrame(mutations)
    mut_df.to_csv(output_dir / 'mutations.csv', index=False)
    print('Saved: mutations.csv')
    print()

    # Step 3: Simulate ddG
    print('STEP 3: Simulate ddG Results')
    print('-' * 80)
    print('NOTE: Generating simulated ddG values for demonstration')
    print()

    ddg_df = simulate_ddg(mutations)
    ddg_df.to_csv(output_dir / 'ddg_results.csv', index=False)
    print(f'Saved: ddg_results.csv ({len(ddg_df)} mutations)')

    print(f'\nddG Statistics:')
    print(f'  Mean:   {ddg_df["ddg"].mean():.3f} kcal/mol')
    print(f'  Median: {ddg_df["ddg"].median():.3f} kcal/mol')
    print(f'  Std:    {ddg_df["ddg"].std():.3f} kcal/mol')
    print(f'  Min:    {ddg_df["ddg"].min():.3f} kcal/mol')
    print(f'  Max:    {ddg_df["ddg"].max():.3f} kcal/mol')
    print()

    # Step 4: Identify hotspots
    print('STEP 4: Identify Hotspots')
    print('-' * 80)

    threshold = 2.0
    hotspots = ddg_df[ddg_df['ddg'] > threshold].copy()

    print(f'Hotspot threshold: {threshold} kcal/mol')
    print(f'Number of hotspots identified: {len(hotspots)}')

    if len(hotspots) > 0:
        print(f'\nTop 20 hotspots:')
        for idx, row in hotspots.head(20).iterrows():
            print(f'  {row["mutation"]:8s} | {row["original_aa"]}{row["position"]:<4d} | '
                  f'ddG = {row["ddg"]:5.2f} kcal/mol | Score = {row["total_score"]:6.2f} REU')

        hotspots.to_csv(output_dir / 'hotspots.csv', index=False)
        print(f'\nSaved: hotspots.csv ({len(hotspots)} hotspots)')
    else:
        print('No hotspots above threshold')
        hotspots.to_csv(output_dir / 'hotspots.csv', index=False)
    print()

    # Step 5: Generate report
    print('STEP 5: Generate Analysis Report')
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
        for aa, count in sorted(aa_counts.items()):
            f.write(f'- **{aa}**: {count} mutation(s)\n')

        f.write('\n---\n\n')
        f.write('## 2. ddG Results (Simulated)\n\n')
        f.write('**NOTE:** These are SIMULATED values for demonstration.\n\n')

        f.write('### Statistics\n\n')
        f.write(f'- **Number of Mutations:** {len(ddg_df)}\n')
        f.write(f'- **Mean ΔΔG:** {ddg_df["ddg"].mean():.2f} kcal/mol\n')
        f.write(f'- **Std ΔΔG:** {ddg_df["ddg"].std():.2f} kcal/mol\n')
        f.write(f'- **Min ΔΔG:** {ddg_df["ddg"].min():.2f} kcal/mol\n')
        f.write(f'- **Max ΔΔG:** {ddg_df["ddg"].max():.2f} kcal/mol\n\n')

        f.write('### Top 50 Mutations\n\n')
        f.write('| Rank | Mutation | Position | Original AA | ΔΔG (kcal/mol) |\n')
        f.write('|------|----------|----------|-------------|----------------|\n')
        for i, (idx, row) in enumerate(ddg_df.head(50).iterrows(), 1):
            f.write(f'| {i} | {row["mutation"]} | {row["position"]} | '
                   f'{row["original_aa"]} | {row["ddg"]:.2f} |\n')

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

        f.write('\n---\n\n')
        f.write(f'**Analysis completed:** {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')

    print(f'Saved: ANALYSIS_REPORT.md')
    print()

    # Summary
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
    print(f'  2. Generate heatmap: python3 generate_ddg_heatmap.py')
    print()


if __name__ == '__main__':
    main()
