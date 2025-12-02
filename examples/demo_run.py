#!/usr/bin/env python
"""
Complete Demo: Alanine Scanning with Example Structure

This script demonstrates the complete workflow using a small example protein.
It shows all steps from mutation generation to visualization.
"""

import sys
from pathlib import Path

# Add src to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from rosetta_scan.protocols.alanine_scanner import AlanineScan
from rosetta_scan.analysis.parser import ResultParser
from rosetta_scan.analysis.visualizer import ResultVisualizer
import pandas as pd


def print_header(text):
    """Print formatted header."""
    print("\n" + "=" * 80)
    print(f"  {text}")
    print("=" * 80 + "\n")


def print_section(text):
    """Print section header."""
    print(f"\n{'─' * 80}")
    print(f"▶ {text}")
    print(f"{'─' * 80}\n")


def main():
    """Run complete demonstration."""

    print_header(" ROSETTA ALANINE SCANNING - COMPLETE DEMO")

    # Setup
    pdb_file = Path(__file__).parent / "example_protein.pdb"
    output_dir = Path(__file__).parent / "demo_output"
    output_dir.mkdir(exist_ok=True)

    print(f" Working directory: {output_dir}")
    print(f" Input structure: {pdb_file.name}")

    # ========================================================================
    # STEP 1: ANALYZE STRUCTURE
    # ========================================================================

    print_section("STEP 1: Structure Analysis")

    print("Loading PDB structure...")
    scanner = AlanineScan(pdb_file)

    # Count residues per chain
    chains_info = {}
    for chain in scanner.structure[0]:
        n_residues = len([r for r in chain if r.id[0] == ' '])
        chains_info[chain.id] = n_residues

    print(f"\n✓ Structure loaded successfully!")
    print(f"\n  Chains found: {list(chains_info.keys())}")
    for chain_id, n_res in chains_info.items():
        print(f"    • Chain {chain_id}: {n_res} residues")

    # ========================================================================
    # STEP 2: GENERATE ALANINE MUTATIONS
    # ========================================================================

    print_section("STEP 2: Generate Alanine Mutations")

    # Generate mutations for all chains
    print("Generating alanine scanning mutations...")
    mutations = scanner.generate_mutations(
        chains=['A', 'B']
    )

    print(f"\n✓ Generated {len(mutations)} mutations")

    # Show mutation summary
    summary = scanner.get_mutation_summary()
    print(f"\nMutation breakdown by chain:")
    for chain, count in summary['chains'].items():
        print(f"  • Chain {chain}: {count} mutations")

    print(f"\nAmino acid distribution:")
    for aa, count in sorted(summary['amino_acid_distribution'].items()):
        print(f"  • {aa}: {count} mutations")

    # Show first few mutations
    print(f"\nFirst 10 mutations to be tested:")
    print(f"  {'#':<4} {'Mutation':<12} {'Chain':<8} {'Position':<10} {'Original→Ala':<15}")
    print(f"  {'-'*60}")
    for i, mut in enumerate(mutations[:10], 1):
        print(f"  {i:<4} {mut.to_rosetta_format():<12} {mut.chain:<8} "
              f"{mut.position:<10} {mut.original_aa}→A")

    if len(mutations) > 10:
        print(f"  ... and {len(mutations) - 10} more mutations")

    # Save mutations
    mutation_file = output_dir / "mutations.txt"
    scanner.save_mutation_report(mutation_file, format='txt')

    mutation_rosetta = output_dir / "mutations_rosetta.txt"
    scanner.save_mutation_report(mutation_rosetta, format='rosetta')

    mutation_csv = output_dir / "mutations.csv"
    scanner.save_mutation_report(mutation_csv, format='csv')

    print(f"\n Mutations saved to:")
    print(f"  • Text report: {mutation_file.name}")
    print(f"  • Rosetta format: {mutation_rosetta.name}")
    print(f"  • CSV format: {mutation_csv.name}")

    # ========================================================================
    # STEP 3: INTERFACE ANALYSIS
    # ========================================================================

    print_section("STEP 3: Interface Analysis")

    print("Identifying interface residues (cutoff = 8.0 Å)...")
    interface_mutations = scanner.generate_mutations(
        chains=['A', 'B'],
        interface_only=True,
        interface_cutoff=8.0
    )

    print(f"\n✓ Found {len(interface_mutations)} interface residues")

    if interface_mutations:
        print(f"\nInterface hotspot candidates:")
        print(f"  {'Chain':<8} {'Position':<10} {'Residue':<10} {'Mutation':<12}")
        print(f"  {'-'*50}")
        for mut in interface_mutations[:8]:
            print(f"  {mut.chain:<8} {mut.position:<10} {mut.original_aa:<10} "
                  f"{mut.to_rosetta_format():<12}")

        if len(interface_mutations) > 8:
            print(f"  ... and {len(interface_mutations) - 8} more")

    interface_file = output_dir / "interface_mutations.txt"
    scanner.save_mutation_report(interface_file, format='txt')
    print(f"\n Interface mutations saved to: {interface_file.name}")

    # ========================================================================
    # STEP 4: SIMULATE DDG RESULTS (for demo purposes)
    # ========================================================================

    print_section("STEP 4: Simulated ddG Results")

    print("NOTE: This is a demo with simulated ddG values.")
    print("In a real run, you would use:")
    print("  rosetta-scan run example_protein.pdb mutations_rosetta.txt --nstruct 35\n")

    # Create simulated results
    import numpy as np
    np.random.seed(42)

    simulated_results = []
    for mut in mutations:
        # Simulate ddG values (real values would come from Rosetta)
        # Some residues are "hotspots" with high ddG
        base_ddg = np.random.normal(0.5, 0.8)

        # Make some positions hotspots
        if mut.position in [2, 4] and mut.chain == 'A':
            base_ddg += np.random.uniform(1.5, 3.0)
        elif mut.position in [2, 3] and mut.chain == 'B':
            base_ddg += np.random.uniform(1.0, 2.5)

        simulated_results.append({
            'mutation': mut.to_rosetta_format(),
            'chain': mut.chain,
            'position': mut.position,
            'original_aa': mut.original_aa,
            'ddg': base_ddg,
            'total_score': base_ddg - 50.0,  # Simulated total score
        })

    results_df = pd.DataFrame(simulated_results)

    print("✓ Simulated results generated")
    print(f"\nStatistics:")
    print(f"  • Total mutations: {len(results_df)}")
    print(f"  • Mean ΔΔG: {results_df['ddg'].mean():.2f} kcal/mol")
    print(f"  • Std ΔΔG: {results_df['ddg'].std():.2f} kcal/mol")
    print(f"  • Min ΔΔG: {results_df['ddg'].min():.2f} kcal/mol")
    print(f"  • Max ΔΔG: {results_df['ddg'].max():.2f} kcal/mol")

    # ========================================================================
    # STEP 5: IDENTIFY HOTSPOTS
    # ========================================================================

    print_section("STEP 5: Hotspot Identification")

    threshold = 1.5
    hotspots = results_df[results_df['ddg'] >= threshold].sort_values('ddg', ascending=False)

    print(f"Identifying hotspots (threshold = {threshold} kcal/mol)...")
    print(f"\n✓ Found {len(hotspots)} hotspot residues")

    if len(hotspots) > 0:
        print(f"\n TOP HOTSPOTS (ranked by ΔΔG):")
        print(f"\n  {'Rank':<6} {'Mutation':<12} {'Chain':<8} {'Position':<10} {'ΔΔG (kcal/mol)':<15} {'Impact':<10}")
        print(f"  {'─'*70}")

        for i, (_, row) in enumerate(hotspots.head(10).iterrows(), 1):
            ddg = row['ddg']
            if ddg > 2.5:
                impact = " Critical"
            elif ddg > 2.0:
                impact = " High"
            else:
                impact = " Medium"

            print(f"  {i:<6} {row['mutation']:<12} {row['chain']:<8} "
                  f"{row['position']:<10} {ddg:<15.2f} {impact:<10}")

    # Save results
    results_file = output_dir / "ddg_results.csv"
    results_df.to_csv(results_file, index=False)

    hotspots_file = output_dir / "hotspots.csv"
    hotspots.to_csv(hotspots_file, index=False)

    print(f"\n Results saved to:")
    print(f"  • All results: {results_file.name}")
    print(f"  • Hotspots: {hotspots_file.name}")

    # ========================================================================
    # STEP 6: GENERATE VISUALIZATIONS
    # ========================================================================

    print_section("STEP 6: Generate Visualizations")

    print("Creating publication-quality plots...")

    visualizer = ResultVisualizer(results_df)
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(exist_ok=True)

    # Generate plots
    plots_created = []

    try:
        print("  • Creating ΔΔG distribution plot...")
        visualizer.plot_ddg_distribution(plot_dir / "ddg_distribution.png")
        plots_created.append("ddg_distribution.png")

        print("  • Creating top hotspots plot...")
        visualizer.plot_top_hotspots(n=10, output_path=plot_dir / "top_hotspots.png")
        plots_created.append("top_hotspots.png")

        print("  • Creating per-chain analysis...")
        visualizer.plot_per_chain_analysis(plot_dir / "chain_analysis.png")
        plots_created.append("chain_analysis.png")

        print("  • Creating hotspot heatmap...")
        visualizer.plot_hotspot_heatmap(plot_dir / "hotspot_heatmap.png", threshold=threshold)
        plots_created.append("hotspot_heatmap.png")

        print("  • Creating position scan plots...")
        for chain in results_df['chain'].unique():
            visualizer.plot_position_scan(
                chain,
                plot_dir / f"position_scan_chain_{chain}.png",
                threshold=threshold
            )
            plots_created.append(f"position_scan_chain_{chain}.png")

        print(f"\n✓ Created {len(plots_created)} visualization plots")
        print(f"\n Plots saved to: {plot_dir}/")
        for plot in plots_created:
            print(f"  • {plot}")

    except Exception as e:
        print(f"\n  Note: Visualization requires matplotlib backend.")
        print(f"    Plots can still be generated in non-interactive mode.")
        print(f"    Error: {e}")

    # ========================================================================
    # STEP 7: EXPORT PYMOL SCRIPT
    # ========================================================================

    print_section("STEP 7: PyMOL Visualization")

    pymol_script = output_dir / "visualize_hotspots.pml"

    print("Generating PyMOL script for 3D visualization...")

    # Create PyMOL script
    with open(pymol_script, 'w') as f:
        f.write("# PyMOL script to visualize alanine scanning hotspots\n")
        f.write(f"# Generated from: {pdb_file.name}\n\n")

        f.write(f"# Load structure\n")
        f.write(f"load {pdb_file.absolute()}\n\n")

        f.write("# Basic visualization\n")
        f.write("hide everything\n")
        f.write("show cartoon\n")
        f.write("color grey80, all\n")
        f.write("util.cbc\n\n")

        f.write("# Highlight hotspots\n")
        for i, (_, row) in enumerate(hotspots.iterrows(), 1):
            chain = row['chain']
            pos = row['position']
            ddg = row['ddg']

            # Color by intensity
            if ddg > 2.5:
                color = "red"
            elif ddg > 2.0:
                color = "orange"
            else:
                color = "yellow"

            f.write(f"# Hotspot {i}: {row['mutation']} (ΔΔG = {ddg:.2f})\n")
            f.write(f"select hot_{i}, chain {chain} and resi {pos}\n")
            f.write(f"show sticks, hot_{i}\n")
            f.write(f"color {color}, hot_{i}\n")
            f.write(f"label hot_{i}, '{row['mutation']}\\n{ddg:.1f}'\n\n")

        f.write("# Zoom to hotspots\n")
        if len(hotspots) > 0:
            f.write("zoom hot_*\n")

        f.write("\n# Set viewing angle\n")
        f.write("set_view (1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,20,0)\n")

    print(f"✓ PyMOL script created: {pymol_script.name}")
    print(f"\nTo visualize in PyMOL, run:")
    print(f"  pymol {pymol_script}")

    # ========================================================================
    # STEP 8: SUMMARY REPORT
    # ========================================================================

    print_section("STEP 8: Generate Summary Report")

    report_file = output_dir / "analysis_report.txt"

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("ROSETTA ALANINE SCANNING - ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Structure: {pdb_file.name}\n")
        f.write(f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("STRUCTURE INFORMATION\n")
        f.write("-" * 80 + "\n")
        for chain_id, n_res in chains_info.items():
            f.write(f"Chain {chain_id}: {n_res} residues\n")
        f.write(f"\n")

        f.write("MUTATION SUMMARY\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total mutations tested: {len(mutations)}\n")
        f.write(f"Interface mutations: {len(interface_mutations)}\n\n")

        f.write("DDG STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Mean ΔΔG: {results_df['ddg'].mean():.2f} kcal/mol\n")
        f.write(f"Std ΔΔG: {results_df['ddg'].std():.2f} kcal/mol\n")
        f.write(f"Min ΔΔG: {results_df['ddg'].min():.2f} kcal/mol\n")
        f.write(f"Max ΔΔG: {results_df['ddg'].max():.2f} kcal/mol\n\n")

        f.write(f"HOTSPOTS (threshold = {threshold} kcal/mol)\n")
        f.write("-" * 80 + "\n")
        f.write(f"Number of hotspots: {len(hotspots)}\n\n")

        if len(hotspots) > 0:
            f.write(f"{'Rank':<6} {'Mutation':<12} {'Chain':<8} {'Position':<10} {'ΔΔG':<12}\n")
            f.write("-" * 80 + "\n")
            for i, (_, row) in enumerate(hotspots.iterrows(), 1):
                f.write(f"{i:<6} {row['mutation']:<12} {row['chain']:<8} "
                       f"{row['position']:<10} {row['ddg']:<12.2f}\n")

        f.write("\n" + "=" * 80 + "\n")

    print(f"✓ Summary report saved: {report_file.name}")

    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================

    print_header(" DEMO COMPLETED SUCCESSFULLY")

    print(" All output files generated in: demo_output/\n")

    print(" Files created:")
    print("  Mutations:")
    print("    • mutations.txt (human-readable)")
    print("    • mutations_rosetta.txt (Rosetta format)")
    print("    • mutations.csv (spreadsheet)")
    print("    • interface_mutations.txt (interface only)")
    print("\n  Results:")
    print("    • ddg_results.csv (all ΔΔG values)")
    print("    • hotspots.csv (filtered hotspots)")
    print("    • analysis_report.txt (summary)")
    print("\n  Visualizations:")
    for plot in plots_created:
        print(f"    • plots/{plot}")
    print("\n  3D Visualization:")
    print("    • visualize_hotspots.pml (PyMOL script)")

    print("\n Key Findings:")
    print(f"  • Tested {len(mutations)} alanine mutations")
    print(f"  • Identified {len(hotspots)} hotspot residues")
    print(f"  • Top hotspot: {hotspots.iloc[0]['mutation']} "
          f"(ΔΔG = {hotspots.iloc[0]['ddg']:.2f} kcal/mol)")

    print("\n Next Steps:")
    print("  1. Review plots in the plots/ directory")
    print("  2. Examine hotspots.csv for detailed results")
    print("  3. Visualize in PyMOL: pymol visualize_hotspots.pml")
    print("  4. Run with real Rosetta for actual ddG values")

    print("\n" + "=" * 80 + "\n")


if __name__ == "__main__":
    main()
