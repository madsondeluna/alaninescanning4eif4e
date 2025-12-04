#!/usr/bin/env python3
"""
Alanine Scanning Analysis - gamma.pdb
High Accuracy Configuration

Protocol: Flex ddG (Barlow et al. 2018)
Score Function: REF2015 (Alford et al. 2017)

Configuration (optimized for small proteins):
- nstruct: 50 (robust ensemble statistics)
- repack_radius: 10.0 Å (extended shell for small proteins)
- max_minimization_iter: 100 (stability vs FastRelax)
- score_function: REF2015
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from rosetta_scan.protocols.flex_ddg_pyrosetta import (
    FlexDdGPyRosetta,
    FlexDdGConfig
)
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
import pandas as pd


def generate_mutations(pdb_path: str) -> list:
    """
    Generate alanine scanning mutation list.

    Excludes:
    - GLY (no Cβ)
    - ALA (already alanine)
    - PRO (conformational constraint)

    Args:
        pdb_path: Path to PDB file

    Returns:
        List of mutation dictionaries
    """
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)

    mutations = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True):
                    resname = residue.get_resname()
                    if resname in aa_dict:
                        aa = aa_dict[resname]
                        pos = residue.get_id()[1]

                        # Skip GLY, ALA, PRO
                        if aa not in ['G', 'A', 'P']:
                            mutations.append({
                                'mutation': f'{aa}{pos}A',
                                'position': pos,
                                'chain': chain.id,
                                'original_aa': aa
                            })

    return mutations


def main():
    """
    Main analysis pipeline.

    Executes Flex ddG protocol with high accuracy parameters
    optimized for publication-quality results.
    """
    print("="*80)
    print("ALANINE SCANNING ANALYSIS - HIGH ACCURACY CONFIGURATION")
    print("="*80)
    print()
    print("Protocol: Flex ddG (Barlow et al. 2018)")
    print("Score Function: REF2015 (Alford et al. 2017)")
    print()
    print("Parameters:")
    print("  nstruct:                  50")
    print("  repack_radius:            10.0 Angstrom")
    print("  max_minimization_iter:    100")
    print("  minimization_algorithm:   MinMover (LBFGS)")
    print("="*80)

    # File paths
    script_dir = Path(__file__).parent
    pdb_path = script_dir / 'gamma.pdb'
    output_dir = script_dir / 'analysis_results_gamma'
    output_dir.mkdir(exist_ok=True)

    # Verify PDB exists
    if not pdb_path.exists():
        print(f"\nERROR: {pdb_path} not found")
        sys.exit(1)

    # Step 1: Generate mutations
    print(f"\nStep 1: Generating mutation list")
    print("-" * 80)
    mutations = generate_mutations(str(pdb_path))
    print(f"Structure: {pdb_path.name}")
    print(f"Total mutations: {len(mutations)}")

    print("\nMutation list:")
    for i, mut in enumerate(mutations, 1):
        print(f"  {i:2d}. {mut['mutation']:6s} (position {mut['position']:2d}, chain {mut['chain']})")

    # Step 2: Configure protocol (HIGH ACCURACY)
    config = FlexDdGConfig(
        nstruct=50,                   # Robust ensemble
        repack_radius=10.0,           # Extended shell for small proteins
        max_minimization_iter=100     # Stability
    )

    # Estimate computational cost
    time_per_mutation_min = 10  # Conservative estimate
    total_time_min = len(mutations) * time_per_mutation_min
    total_time_h = total_time_min / 60

    print(f"\nStep 2: Flex ddG Protocol Execution")
    print("-" * 80)
    print(f"Estimated time: {total_time_h:.1f} hours ({total_time_min} minutes)")
    print()

    # Execute protocol
    print("Initializing Flex ddG protocol...")
    print("="*80)

    protocol = FlexDdGPyRosetta(config)

    try:
        results = protocol.run_flex_ddg(
            str(pdb_path),
            mutations,
            str(output_dir)
        )

        # Step 3: Statistical analysis
        print("\n" + "="*80)
        print("STATISTICAL ANALYSIS")
        print("="*80)

        ddg_mean = results['ddg'].mean()
        ddg_std = results['ddg'].std()
        ddg_max = results['ddg'].max()
        ddg_min = results['ddg'].min()
        ddg_sem = results['ddg_std'].mean()  # Mean standard error

        print(f"\nDescriptive Statistics:")
        print(f"  Mean ddG:        {ddg_mean:+7.2f} +/- {ddg_std:5.2f} kcal/mol")
        print(f"  Maximum ddG:     {ddg_max:+7.2f} kcal/mol")
        print(f"  Minimum ddG:     {ddg_min:+7.2f} kcal/mol")
        print(f"  Mean SEM:        {ddg_sem:7.2f} kcal/mol")

        print(f"\nConvergence Assessment:")
        if ddg_sem < 0.5:
            status = "Excellent (SEM < 0.5 kcal/mol)"
        elif ddg_sem < 1.0:
            status = "Good (SEM < 1.0 kcal/mol)"
        else:
            status = "Moderate (SEM > 1.0 kcal/mol)"
        print(f"  Status: {status}")

        # Identify hotspots (ddG > +2.0 kcal/mol)
        hotspot_threshold = 2.0
        hotspots = results[results['ddg'] > hotspot_threshold].copy()
        hotspots = hotspots.sort_values('ddg', ascending=False)

        print(f"\nHotspot Identification:")
        print(f"  Threshold: ddG > +{hotspot_threshold:.1f} kcal/mol")
        print(f"  Hotspots identified: {len(hotspots)}/{len(mutations)} ({len(hotspots)/len(mutations)*100:.1f}%)")

        if len(hotspots) > 0:
            print(f"\n{'='*80}")
            print("HOTSPOT RESIDUES")
            print("="*80)
            print(f"{'Rank':<6} {'Mutation':<10} {'Position':<10} {'ddG (kcal/mol)':<18} {'SEM':<12}")
            print("-"*80)

            for rank, (_, row) in enumerate(hotspots.iterrows(), 1):
                print(f"{rank:<6} {row['mutation']:<10} {row['position']:<10} "
                      f"{row['ddg']:+6.2f} +/- {row['ddg_std']:4.2f}     {row['ddg_std']:5.2f}")

            # Save hotspots
            hotspots_file = output_dir / 'hotspots.csv'
            hotspots.to_csv(hotspots_file, index=False)

        # Save complete results
        results_file = output_dir / 'ddg_results.csv'
        results.to_csv(results_file, index=False)

        # Generate summary report
        summary_file = output_dir / 'analysis_summary.txt'
        with open(summary_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("ALANINE SCANNING ANALYSIS SUMMARY\n")
            f.write("="*80 + "\n\n")

            f.write("STRUCTURE\n")
            f.write("-" * 80 + "\n")
            f.write(f"PDB file: {pdb_path.name}\n")
            f.write(f"Total mutations tested: {len(mutations)}\n\n")

            f.write("PROTOCOL\n")
            f.write("-" * 80 + "\n")
            f.write(f"Method: Flex ddG (Barlow et al. 2018)\n")
            f.write(f"Score function: REF2015 (Alford et al. 2017)\n")
            f.write(f"nstruct: {config.nstruct}\n")
            f.write(f"Repack radius: {config.repack_radius} A\n")
            f.write(f"Max minimization iterations: {config.max_minimization_iter}\n\n")

            f.write("STATISTICAL RESULTS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Mean ddG: {ddg_mean:+.2f} +/- {ddg_std:.2f} kcal/mol\n")
            f.write(f"Maximum ddG: {ddg_max:+.2f} kcal/mol\n")
            f.write(f"Minimum ddG: {ddg_min:+.2f} kcal/mol\n")
            f.write(f"Mean SEM: {ddg_sem:.2f} kcal/mol\n\n")

            f.write("HOTSPOT ANALYSIS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Threshold: ddG > +{hotspot_threshold:.1f} kcal/mol\n")
            f.write(f"Hotspots identified: {len(hotspots)}/{len(mutations)} ({len(hotspots)/len(mutations)*100:.1f}%)\n\n")

            if len(hotspots) > 0:
                f.write("HOTSPOT RESIDUES (ranked by ddG)\n")
                f.write("-" * 80 + "\n")
                for rank, (_, row) in enumerate(hotspots.iterrows(), 1):
                    f.write(f"{rank:2d}. {row['mutation']:8s} | ddG = {row['ddg']:+.2f} +/- {row['ddg_std']:.2f} kcal/mol\n")

            f.write("\n" + "="*80 + "\n")
            f.write("REFERENCES\n")
            f.write("="*80 + "\n")
            f.write("[1] Barlow KA, et al. (2018) Flex ddG. J Phys Chem B. 122(21):5389-5399\n")
            f.write("[2] Alford RF, et al. (2017) REF2015. J Chem Theory Comput. 13(6):3031-3048\n")

        print(f"\n{'='*80}")
        print(f"ANALYSIS COMPLETED SUCCESSFULLY")
        print(f"{'='*80}")
        print(f"\nOutput files:")
        print(f"  {results_file}")
        print(f"  {hotspots_file if len(hotspots) > 0 else '(no hotspots)'}")
        print(f"  {summary_file}")
        print("="*80)

    except Exception as e:
        print(f"\nERROR during analysis:")
        print(f"  {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
