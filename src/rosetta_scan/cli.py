"""
Command-Line Interface for Rosetta Alanine Scanning

Elegant CLI using Click for running Flex ddG protocols.
"""

import click
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich import print as rprint
import yaml
import json

from .protocols.flex_ddg import FlexDdGProtocol, FlexDdGConfig
from .protocols.alanine_scanner import AlanineScan
from .analysis.visualizer import ResultVisualizer
from .analysis.parser import ResultParser

console = Console()


@click.group()
@click.version_option(version="0.1.0")
def main():
    """
    🧬 Rosetta Alanine Scanning - Elegant Flex ddG Protocol

    A sophisticated tool for computational alanine scanning and
    hotspot identification using Rosetta's Flex ddG protocol.
    """
    pass


@main.command()
@click.argument('pdb_file', type=click.Path(exists=True))
@click.option(
    '--chains',
    '-c',
    multiple=True,
    help='Chain IDs to scan (e.g., -c A -c B). Scans all if not specified.'
)
@click.option(
    '--range',
    '-r',
    type=str,
    help='Residue range per chain (e.g., "A:1-100,B:50-150")'
)
@click.option(
    '--interface-only',
    is_flag=True,
    help='Only scan interface residues'
)
@click.option(
    '--interface-cutoff',
    type=float,
    default=8.0,
    help='Distance cutoff for interface definition (Å)'
)
@click.option(
    '--output',
    '-o',
    type=click.Path(),
    default='mutations',
    help='Output directory for mutation files'
)
@click.option(
    '--format',
    type=click.Choice(['txt', 'csv', 'rosetta']),
    default='txt',
    help='Output format'
)
def scan(pdb_file, chains, range, interface_only, interface_cutoff, output, format):
    """
    Generate alanine scanning mutations from PDB structure.

    Example:
        rosetta-scan scan protein.pdb -c A -c B --interface-only
    """
    console.print("\n[bold cyan]🔬 Alanine Scanning Mutation Generator[/bold cyan]\n")

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:

        # Load structure
        task = progress.add_task("Loading structure...", total=None)
        scanner = AlanineScan(pdb_file)
        progress.update(task, completed=True)

        # Parse residue ranges
        residue_range = None
        if range:
            residue_range = {}
            for part in range.split(','):
                chain, rng = part.split(':')
                start, end = map(int, rng.split('-'))
                residue_range[chain] = (start, end)

        # Generate mutations
        task = progress.add_task("Generating mutations...", total=None)
        chains_list = list(chains) if chains else None

        mutations = scanner.generate_mutations(
            chains=chains_list,
            residue_range=residue_range,
            interface_only=interface_only,
            interface_cutoff=interface_cutoff
        )
        progress.update(task, completed=True)

    # Display summary
    summary = scanner.get_mutation_summary()

    table = Table(title="Mutation Summary", show_header=True)
    table.add_column("Property", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Total Mutations", str(summary['total_mutations']))
    table.add_row("Structure File", str(summary['structure_file']))

    for chain, count in summary['chains'].items():
        table.add_row(f"Chain {chain}", str(count))

    console.print("\n")
    console.print(table)

    # Save mutations
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"mutations.{format}"
    scanner.save_mutation_report(output_file, format=format)

    console.print(f"\n✅ Mutations saved to: [bold]{output_file}[/bold]\n")


@main.command()
@click.argument('pdb_file', type=click.Path(exists=True))
@click.argument('mutation_file', type=click.Path(exists=True))
@click.option(
    '--config',
    type=click.Path(exists=True),
    help='YAML configuration file'
)
@click.option(
    '--nstruct',
    type=int,
    default=35,
    help='Number of structures per mutation'
)
@click.option(
    '--iterations',
    type=int,
    default=3,
    help='Number of backrub iterations'
)
@click.option(
    '--interface',
    is_flag=True,
    help='Run interface ddG mode'
)
@click.option(
    '--output',
    '-o',
    type=click.Path(),
    default='ddg_results',
    help='Output directory'
)
@click.option(
    '--rosetta-path',
    type=click.Path(exists=True),
    envvar='ROSETTA',
    help='Path to Rosetta installation'
)
def run(pdb_file, mutation_file, config, nstruct, iterations, interface, output, rosetta_path):
    """
    Run Flex ddG protocol for mutations.

    Example:
        rosetta-scan run protein.pdb mutations.txt --nstruct 35
    """
    console.print("\n[bold cyan]🧬 Running Flex ddG Protocol[/bold cyan]\n")

    # Load configuration
    if config:
        flex_config = FlexDdGConfig.from_yaml(config)
    else:
        flex_config = FlexDdGConfig(
            nstruct=nstruct,
            iterations=iterations,
            interface_ddg=interface,
            rosetta_path=rosetta_path
        )

    # Load mutations
    with open(mutation_file, 'r') as f:
        content = f.read()

    # Parse mutations (handle different formats)
    mutations = []
    for line in content.split('\n'):
        line = line.strip()
        if line and not line.startswith('total') and line != '1':
            mutations.append(line)

    # Display configuration
    table = Table(title="Protocol Configuration")
    table.add_column("Parameter", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Structure", str(pdb_file))
    table.add_row("Mutations", str(len(mutations)))
    table.add_row("Nstruct", str(flex_config.nstruct))
    table.add_row("Iterations", str(flex_config.iterations))
    table.add_row("Interface Mode", "Yes" if flex_config.interface_ddg else "No")

    console.print(table)
    console.print()

    # Initialize protocol
    with console.status("[bold green]Initializing Rosetta protocol..."):
        protocol = FlexDdGProtocol(flex_config)

    console.print("✅ Protocol initialized\n")

    # Run protocol
    try:
        with console.status("[bold green]Running Flex ddG calculations..."):
            result_dir = protocol.run_flex_ddg(
                pdb_path=pdb_file,
                mutations=mutations,
                output_dir=output
            )

        console.print(f"✅ Calculations complete!\n")
        console.print(f"📁 Results saved to: [bold]{result_dir}[/bold]\n")

    except Exception as e:
        console.print(f"[bold red]❌ Error:[/bold red] {str(e)}\n")
        raise


@main.command()
@click.argument('results_dir', type=click.Path(exists=True))
@click.option(
    '--output',
    '-o',
    type=click.Path(),
    help='Output file for analysis (CSV/JSON)'
)
@click.option(
    '--plot',
    is_flag=True,
    help='Generate visualization plots'
)
@click.option(
    '--threshold',
    type=float,
    default=1.0,
    help='ddG threshold for hotspot identification (kcal/mol)'
)
def analyze(results_dir, output, plot, threshold):
    """
    Analyze Flex ddG results and identify hotspots.

    Example:
        rosetta-scan analyze ddg_results/ --plot --threshold 1.5
    """
    console.print("\n[bold cyan]📊 Analyzing Flex ddG Results[/bold cyan]\n")

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console
    ) as progress:

        # Parse results
        task = progress.add_task("Parsing results...", total=None)
        parser = ResultParser(results_dir)
        results_df = parser.parse_results()
        progress.update(task, completed=True)

        # Identify hotspots
        task = progress.add_task("Identifying hotspots...", total=None)
        hotspots = parser.identify_hotspots(threshold=threshold)
        progress.update(task, completed=True)

    # Display summary
    table = Table(title="Analysis Summary")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Total Mutations", str(len(results_df)))
    table.add_row("Hotspots (ddG > threshold)", str(len(hotspots)))
    table.add_row("Mean ddG", f"{results_df['ddg'].mean():.2f} kcal/mol")
    table.add_row("Std ddG", f"{results_df['ddg'].std():.2f} kcal/mol")

    console.print(table)
    console.print()

    # Display hotspots
    if len(hotspots) > 0:
        hotspot_table = Table(title="Identified Hotspots")
        hotspot_table.add_column("Mutation", style="cyan")
        hotspot_table.add_column("ddG (kcal/mol)", style="red")
        hotspot_table.add_column("Chain", style="green")

        for _, row in hotspots.head(10).iterrows():
            hotspot_table.add_row(
                row['mutation'],
                f"{row['ddg']:.2f}",
                row.get('chain', 'N/A')
            )

        console.print(hotspot_table)
        console.print()

    # Save results
    if output:
        output_path = Path(output)
        if output_path.suffix == '.csv':
            results_df.to_csv(output_path, index=False)
        elif output_path.suffix == '.json':
            results_df.to_json(output_path, orient='records', indent=2)

        console.print(f"💾 Results saved to: [bold]{output_path}[/bold]\n")

    # Generate plots
    if plot:
        with console.status("[bold green]Generating visualizations..."):
            visualizer = ResultVisualizer(results_df)
            plot_dir = Path(results_dir) / "plots"
            plot_dir.mkdir(exist_ok=True)

            visualizer.plot_ddg_distribution(plot_dir / "ddg_distribution.png")
            visualizer.plot_hotspot_heatmap(plot_dir / "hotspot_heatmap.png")
            visualizer.plot_per_chain_analysis(plot_dir / "chain_analysis.png")

        console.print(f"📈 Plots saved to: [bold]{plot_dir}[/bold]\n")


@main.command()
@click.option(
    '--output',
    '-o',
    type=click.Path(),
    default='config.yaml',
    help='Output configuration file'
)
def init_config(output):
    """
    Generate example configuration file.

    Example:
        rosetta-scan init-config -o my_config.yaml
    """
    console.print("\n[bold cyan]📝 Generating Configuration File[/bold cyan]\n")

    config = FlexDdGConfig()
    config.save_yaml(output)

    console.print(f"✅ Configuration saved to: [bold]{output}[/bold]\n")
    console.print("Edit this file to customize protocol parameters.\n")


@main.command()
@click.argument('pdb_file', type=click.Path(exists=True))
@click.option(
    '--chains',
    '-c',
    multiple=True,
    help='Chain IDs to scan'
)
@click.option(
    '--interface-only',
    is_flag=True,
    help='Only scan interface residues'
)
@click.option(
    '--nstruct',
    type=int,
    default=35,
    help='Number of structures per mutation'
)
@click.option(
    '--output',
    '-o',
    type=click.Path(),
    default='alanine_scan',
    help='Output directory'
)
@click.option(
    '--config',
    type=click.Path(exists=True),
    help='Configuration file'
)
@click.option(
    '--rosetta-path',
    type=click.Path(exists=True),
    envvar='ROSETTA',
    help='Path to Rosetta installation'
)
def pipeline(pdb_file, chains, interface_only, nstruct, output, config, rosetta_path):
    """
    Run complete alanine scanning pipeline (scan + run + analyze).

    This command combines mutation generation, Flex ddG calculations,
    and result analysis into a single streamlined workflow.

    Example:
        rosetta-scan pipeline protein.pdb -c A -c B --interface-only
    """
    console.print("\n[bold cyan]🔬 Alanine Scanning Pipeline[/bold cyan]\n")

    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Generate mutations
    console.print("[bold]Step 1/3:[/bold] Generating mutations\n")

    scanner = AlanineScan(pdb_file)
    chains_list = list(chains) if chains else None

    mutations = scanner.generate_mutations(
        chains=chains_list,
        interface_only=interface_only
    )

    mutation_file = output_dir / "mutations.txt"
    scanner.save_mutation_report(mutation_file, format='rosetta')

    summary = scanner.get_mutation_summary()
    console.print(f"✅ Generated {summary['total_mutations']} mutations\n")

    # Step 2: Run Flex ddG
    console.print("[bold]Step 2/3:[/bold] Running Flex ddG protocol\n")

    if config:
        flex_config = FlexDdGConfig.from_yaml(config)
    else:
        flex_config = FlexDdGConfig(
            nstruct=nstruct,
            rosetta_path=rosetta_path
        )

    protocol = FlexDdGProtocol(flex_config)

    ddg_dir = output_dir / "ddg_results"

    with console.status("[bold green]Running calculations (this may take a while)..."):
        protocol.run_flex_ddg(
            pdb_path=pdb_file,
            mutations=scanner.get_rosetta_mutation_list(),
            output_dir=ddg_dir
        )

    console.print("✅ Flex ddG calculations complete\n")

    # Step 3: Analyze results
    console.print("[bold]Step 3/3:[/bold] Analyzing results\n")

    parser = ResultParser(ddg_dir)
    results_df = parser.parse_results()
    hotspots = parser.identify_hotspots(threshold=1.0)

    # Save analysis
    results_df.to_csv(output_dir / "results.csv", index=False)

    # Generate plots
    visualizer = ResultVisualizer(results_df)
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(exist_ok=True)

    visualizer.plot_ddg_distribution(plot_dir / "ddg_distribution.png")
    visualizer.plot_hotspot_heatmap(plot_dir / "hotspot_heatmap.png")

    console.print(f"✅ Analysis complete\n")
    console.print(f"📊 Found {len(hotspots)} hotspots\n")
    console.print(f"📁 All results saved to: [bold]{output_dir}[/bold]\n")


if __name__ == '__main__':
    main()
