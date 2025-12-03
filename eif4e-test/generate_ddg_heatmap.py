#!/usr/bin/env python3
"""
Standalone script to generate ΔΔG heatmap visualization
Uses twilight_shifted colormap with total_score values
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys


def load_ddg_results(csv_path):
    """
    Load ΔΔG results from CSV file

    Args:
        csv_path: Path to ddg_results.csv

    Returns:
        pandas.DataFrame with columns: mutation, chain, position, original_aa, ddg, total_score
    """
    try:
        df = pd.read_csv(csv_path)
        print(f"Loaded {len(df)} mutations from {csv_path}")
        return df
    except FileNotFoundError:
        print(f"Error: File not found: {csv_path}")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading CSV: {e}")
        sys.exit(1)




def create_heatmap_matrix(df, value_column='total_score'):
    """
    Create matrix for heatmap visualization

    Args:
        df: DataFrame with mutation data
        value_column: Column name to use for heatmap values ('ddg' or 'total_score')

    Returns:
        tuple: (matrix, row_labels, col_labels)
    """
    # Get unique chains and positions
    chains = sorted(df['chain'].unique())
    positions = sorted(df['position'].unique())

    # Create matrix
    matrix = np.full((len(chains), len(positions)), np.nan)

    # Fill matrix with selected values
    for idx, row in df.iterrows():
        chain_idx = chains.index(row['chain'])
        pos_idx = positions.index(row['position'])
        matrix[chain_idx, pos_idx] = row[value_column]

    # Row labels: Chain IDs
    row_labels = [f"Chain {c}" for c in chains]

    # Column labels: Position + Original AA
    col_labels = []
    for pos in positions:
        # Find original AA for this position (from first mutation at this position)
        aa_row = df[df['position'] == pos].iloc[0]
        col_labels.append(f"{aa_row['original_aa']}{pos}")

    return matrix, row_labels, col_labels


def generate_heatmap(ddg_df, output_path, title="Total Score Heatmap - Alanine Scanning",
                    value_column='total_score', colormap='twilight_shifted'):
    """
    Generate heatmap with total_score values using twilight_shifted colormap

    Args:
        ddg_df: DataFrame with ΔΔG results
        output_path: Path to save PNG file
        title: Plot title
        value_column: Column to visualize ('total_score' or 'ddg')
        colormap: Matplotlib colormap name
    """
    # Create matrix
    matrix, row_labels, col_labels = create_heatmap_matrix(ddg_df, value_column=value_column)

    # Get min and max values for colormap range
    valid_values = matrix[~np.isnan(matrix)]
    if len(valid_values) > 0:
        min_val = valid_values.min()
        max_val = valid_values.max()
    else:
        min_val, max_val = 0, 1

    # Create figure
    fig, ax = plt.subplots(figsize=(14, max(4, len(row_labels) * 1.5)))

    # Generate heatmap with twilight_shifted colormap
    im = ax.imshow(matrix, cmap=colormap, aspect='auto',
                   vmin=min_val, vmax=max_val, interpolation='nearest')

    # Set ticks and labels
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha='right')
    ax.set_yticklabels(row_labels)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    if value_column == 'total_score':
        cbar.set_label('Total Score (REU)', rotation=270, labelpad=20, fontsize=12)
    else:
        cbar.set_label('ΔΔG (kcal/mol)', rotation=270, labelpad=20, fontsize=12)

    # Annotate cells with values
    for i in range(len(row_labels)):
        for j in range(len(col_labels)):
            if not np.isnan(matrix[i, j]):
                # Determine text color based on value relative to range
                normalized_val = (matrix[i, j] - min_val) / (max_val - min_val) if max_val != min_val else 0.5
                text_color = 'white' if normalized_val < 0.5 else 'black'
                ax.text(j, i, f'{matrix[i, j]:.2f}',
                       ha='center', va='center', color=text_color,
                       fontsize=9, weight='bold')

    # Add title and labels
    ax.set_title(title, fontsize=14, weight='bold', pad=20)
    ax.set_xlabel('Residue Position', fontsize=12)
    ax.set_ylabel('Chain', fontsize=12)

    # Add grid
    ax.set_xticks(np.arange(len(col_labels)) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(row_labels)) - 0.5, minor=True)
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5)

    # Tight layout
    plt.tight_layout()

    # Save figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to: {output_path}")

    # Close figure to free memory
    plt.close()


def generate_statistics_table(ddg_df):
    """
    Print statistics table for ΔΔG values

    Args:
        ddg_df: DataFrame with ΔΔG results
    """
    print("\n" + "="*60)
    print("ΔΔG STATISTICS")
    print("="*60)

    print(f"\nTotal mutations: {len(ddg_df)}")
    print(f"Mean ΔΔG: {ddg_df['ddg'].mean():.3f} kcal/mol")
    print(f"Median ΔΔG: {ddg_df['ddg'].median():.3f} kcal/mol")
    print(f"Std ΔΔG: {ddg_df['ddg'].std():.3f} kcal/mol")
    print(f"Min ΔΔG: {ddg_df['ddg'].min():.3f} kcal/mol")
    print(f"Max ΔΔG: {ddg_df['ddg'].max():.3f} kcal/mol")

    # Top 5 mutations
    print("\nTop 5 mutations by ΔΔG:")
    print("-" * 60)
    top5 = ddg_df.nlargest(5, 'ddg')[['mutation', 'chain', 'position', 'original_aa', 'ddg']]
    for idx, row in top5.iterrows():
        print(f"  {row['mutation']:8s} | Chain {row['chain']} | {row['original_aa']}{row['position']:<3d} | ΔΔG = {row['ddg']:6.3f} kcal/mol")

    # Bottom 5 mutations
    print("\nBottom 5 mutations by ΔΔG:")
    print("-" * 60)
    bottom5 = ddg_df.nsmallest(5, 'ddg')[['mutation', 'chain', 'position', 'original_aa', 'ddg']]
    for idx, row in bottom5.iterrows():
        print(f"  {row['mutation']:8s} | Chain {row['chain']} | {row['original_aa']}{row['position']:<3d} | ΔΔG = {row['ddg']:6.3f} kcal/mol")

    print("\n" + "="*60 + "\n")


def main():
    """
    Main execution function
    """
    print("\n" + "="*60)
    print("TOTAL SCORE HEATMAP GENERATOR")
    print("Colormap: twilight_shifted")
    print("="*60 + "\n")

    # Define paths
    script_dir = Path(__file__).parent
    results_dir = script_dir / 'analysis_results'
    csv_path = results_dir / 'ddg_results.csv'
    plots_dir = results_dir / 'plots'
    plots_dir.mkdir(exist_ok=True)

    output_path = plots_dir / 'total_score_heatmap.png'

    # Check if input file exists
    if not csv_path.exists():
        print(f"Error: Input file not found: {csv_path}")
        print("Please run the analysis first to generate ddg_results.csv")
        sys.exit(1)

    # Load data
    print(f"Loading data from: {csv_path}")
    ddg_df = load_ddg_results(csv_path)

    # Print statistics
    generate_statistics_table(ddg_df)

    # Generate heatmap
    print(f"Generating heatmap...")
    generate_heatmap(ddg_df, output_path,
                    title="Total Score Heatmap - Alanine Scanning",
                    value_column='total_score',
                    colormap='twilight_shifted')

    print(f"\nVisualization complete!")
    print(f"Output saved to: {output_path}")
    print("\n" + "="*60 + "\n")


if __name__ == '__main__':
    main()
