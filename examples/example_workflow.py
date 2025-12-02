#!/usr/bin/env python
"""
Example workflow for Rosetta Alanine Scanning

This script demonstrates how to use the rosetta_scan package
programmatically for a complete alanine scanning analysis.
"""

from pathlib import Path
from rosetta_scan.protocols.flex_ddg import FlexDdGProtocol, FlexDdGConfig
from rosetta_scan.protocols.alanine_scanner import AlanineScan
from rosetta_scan.analysis.parser import ResultParser
from rosetta_scan.analysis.visualizer import ResultVisualizer


def main():
    """Run complete alanine scanning workflow."""

    # Configuration
    pdb_file = "protein.pdb"  # Replace with your PDB file
    output_dir = Path("workflow_results")
    output_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("Rosetta Alanine Scanning - Example Workflow")
    print("=" * 80)

    # Step 1: Generate alanine mutations
    print("\n[Step 1] Generating alanine mutations...")

    scanner = AlanineScan(pdb_file)

    # Option A: Scan all chains
    mutations = scanner.generate_mutations()

    # Option B: Scan specific chains
    # mutations = scanner.generate_mutations(chains=['A', 'B'])

    # Option C: Scan interface only
    # mutations = scanner.generate_mutations(
    #     chains=['A', 'B'],
    #     interface_only=True,
    #     interface_cutoff=8.0
    # )

    # Option D: Scan specific residue range
    # mutations = scanner.generate_mutations(
    #     residue_range={'A': (1, 100), 'B': (50, 150)}
    # )

    # Display summary
    summary = scanner.get_mutation_summary()
    print(f"  Total mutations: {summary['total_mutations']}")
    print(f"  Chains: {list(summary['chains'].keys())}")

    # Save mutations
    mutation_file = output_dir / "mutations.txt"
    scanner.save_mutation_report(mutation_file, format='rosetta')
    print(f"  Saved to: {mutation_file}")

    # Step 2: Configure and run Flex ddG
    print("\n[Step 2] Running Flex ddG protocol...")

    # Create configuration
    config = FlexDdGConfig(
        nstruct=5,  # Use 35 for production runs
        iterations=3,
        interface_ddg=False,  # Set True for protein-protein interfaces
        num_processors=1,
        # rosetta_path="/path/to/rosetta"  # Uncomment and set path
    )

    # Initialize protocol
    protocol = FlexDdGProtocol(config)

    # Run calculations
    ddg_dir = output_dir / "ddg_results"
    try:
        protocol.run_flex_ddg(
            pdb_path=pdb_file,
            mutations=scanner.get_rosetta_mutation_list(),
            output_dir=ddg_dir
        )
        print(f"  Results saved to: {ddg_dir}")
    except Exception as e:
        print(f"  Error running Flex ddG: {e}")
        print("  Make sure Rosetta is properly installed and configured.")
        return

    # Step 3: Parse and analyze results
    print("\n[Step 3] Analyzing results...")

    parser = ResultParser(ddg_dir)
    results_df = parser.parse_results()

    print(f"  Parsed {len(results_df)} results")
    print(f"  Mean ΔΔG: {results_df['ddg'].mean():.2f} kcal/mol")
    print(f"  Std ΔΔG: {results_df['ddg'].std():.2f} kcal/mol")

    # Identify hotspots
    hotspots = parser.identify_hotspots(threshold=1.0)
    print(f"  Identified {len(hotspots)} hotspots (ΔΔG > 1.0 kcal/mol)")

    # Display top hotspots
    if len(hotspots) > 0:
        print("\n  Top 5 hotspots:")
        for i, (_, row) in enumerate(hotspots.head(5).iterrows(), 1):
            print(f"    {i}. {row['mutation']}: {row['ddg']:.2f} kcal/mol")

    # Save results
    results_df.to_csv(output_dir / "results.csv", index=False)
    parser.export_summary_report(output_dir / "summary_report.txt")
    parser.export_pymol_script(output_dir / "visualize_hotspots.pml")

    # Step 4: Generate visualizations
    print("\n[Step 4] Generating visualizations...")

    visualizer = ResultVisualizer(results_df)
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(exist_ok=True)

    # Create all plots
    visualizer.create_dashboard(plot_dir, threshold=1.0)
    print(f"  Plots saved to: {plot_dir}")

    print("\n" + "=" * 80)
    print("Workflow complete!")
    print(f"All results saved to: {output_dir}")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
