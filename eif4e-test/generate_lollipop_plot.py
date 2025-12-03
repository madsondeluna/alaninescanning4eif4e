#!/usr/bin/env python3
"""
Generate lollipop plot for eIF4E Alanine Scanning results
Similar style to mgrB analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# -----------------------------------------------------------------------------
# Load and Prepare Data
# -----------------------------------------------------------------------------
def load_data():
    """Load ddG results from CSV file"""
    script_dir = Path(__file__).parent
    csv_path = script_dir / 'analysis_results' / 'ddg_results.csv'

    if not csv_path.exists():
        print(f"Error: File not found: {csv_path}")
        print("Please run the analysis first to generate ddg_results.csv")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    return df

# -----------------------------------------------------------------------------
# Chart Customization Parameters
# -----------------------------------------------------------------------------
# --- Figure Size & Style ---
figure_width = 28
figure_height = 8
marker_size = 50

# --- Font Sizes ---
title_fontsize = 24
axis_label_fontsize = 20
tick_label_fontsize = 16

# --- Y-Axis Limits ---
y_axis_min = -4
y_axis_max = 0

# --- Color Scheme ---
# Color based on ddG magnitude (higher = more critical)
color_low = '#1f77b4'      # Blue for low ddG
color_medium = '#ff7f0e'   # Orange for medium ddG
color_high = '#d62728'     # Red for high ddG

# Thresholds for color assignment
threshold_medium = 1.5
threshold_high = 2.5

# -----------------------------------------------------------------------------
# Generate the Lollipop Plot
# -----------------------------------------------------------------------------
def generate_lollipop_plot():
    """Generate and save lollipop plot"""

    print("\n" + "="*60)
    print("LOLLIPOP PLOT GENERATOR - eIF4E Alanine Scanning")
    print("="*60 + "\n")

    # Load data
    df = load_data()
    print(f"Loaded {len(df)} mutations")

    # Sort by position
    df_sorted = df.sort_values(by='position').reset_index(drop=True)

    # Use ddg_normalized column if available, otherwise create it
    if 'ddg_normalized' not in df_sorted.columns:
        df_sorted['ddg_normalized'] = -df_sorted['ddg']

    # Define colors based on ddG magnitude (using original positive values for thresholds)
    colors = []
    for ddg in df_sorted['ddg']:
        if ddg >= threshold_high:
            colors.append(color_high)
        elif ddg >= threshold_medium:
            colors.append(color_medium)
        else:
            colors.append(color_low)

    # Print statistics
    print(f"\nddG Statistics:")
    print(f"  Mean:   {df_sorted['ddg'].mean():.3f} kcal/mol")
    print(f"  Median: {df_sorted['ddg'].median():.3f} kcal/mol")
    print(f"  Min:    {df_sorted['ddg'].min():.3f} kcal/mol")
    print(f"  Max:    {df_sorted['ddg'].max():.3f} kcal/mol")

    # Count by color category
    n_high = sum(1 for c in colors if c == color_high)
    n_medium = sum(1 for c in colors if c == color_medium)
    n_low = sum(1 for c in colors if c == color_low)

    print(f"\nColor Distribution:")
    print(f"  High impact (red):    {n_high} mutations (ddG >= {threshold_high})")
    print(f"  Medium impact (orange): {n_medium} mutations ({threshold_medium} <= ddG < {threshold_high})")
    print(f"  Low impact (blue):    {n_low} mutations (ddG < {threshold_medium})")
    print()

    # Set the figure size
    plt.figure(figsize=(figure_width, figure_height))

    # Create the stems of the lollipop chart (using normalized negative values)
    (markerline, stemlines, baseline) = plt.stem(
        df_sorted['position'],
        df_sorted['ddg_normalized'],
        linefmt='grey',
        markerfmt='o',
        basefmt='black'
    )

    # Make default markers invisible
    plt.setp(markerline, 'markerfacecolor', 'none', 'markeredgecolor', 'none')

    # Draw custom-colored markers (using normalized negative values)
    plt.scatter(
        df_sorted['position'],
        df_sorted['ddg_normalized'],
        c=colors,
        s=marker_size,
        zorder=3,
        edgecolors='black',
        linewidths=0.5
    )

    # Apply Y-axis limits
    if y_axis_min is not None and y_axis_max is not None:
        plt.ylim(bottom=y_axis_min, top=y_axis_max)

    # Set titles and labels
    plt.title('eIF4E - Alanine Scan Energy Impact', fontsize=title_fontsize, weight='bold')
    plt.xlabel('Residue Position', fontsize=axis_label_fontsize)
    plt.ylabel('ΔΔG (kcal/mol)', fontsize=axis_label_fontsize)

    # Set tick properties
    max_pos = df_sorted['position'].max()
    plt.xticks(np.arange(0, max_pos + 1, 10), fontsize=tick_label_fontsize)
    plt.yticks(fontsize=tick_label_fontsize)

    # Add horizontal lines for thresholds (negative values)
    plt.axhline(y=-threshold_medium, color='orange', linestyle='--', alpha=0.3, linewidth=1)
    plt.axhline(y=-threshold_high, color='red', linestyle='--', alpha=0.3, linewidth=1)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color_high, edgecolor='black', label=f'High impact (ΔΔG ≤ -{threshold_high})'),
        Patch(facecolor=color_medium, edgecolor='black', label=f'Medium impact (-{threshold_high} < ΔΔG ≤ -{threshold_medium})'),
        Patch(facecolor=color_low, edgecolor='black', label=f'Low impact (ΔΔG > -{threshold_medium})')
    ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=14, framealpha=0.9)

    # Add grid
    plt.grid(axis='y', linestyle=':', color='gray', alpha=0.5)

    # Adjust layout
    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent / 'analysis_results' / 'plots'
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / 'eif4e_lollipop_plot.png'

    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")

    # Also save PDF version
    pdf_path = output_dir / 'eif4e_lollipop_plot.pdf'
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"PDF saved to: {pdf_path}")

    plt.close()

    print("\nVisualization complete!")
    print("="*60 + "\n")


# -----------------------------------------------------------------------------
# Variant: Plot with Total Score instead of ddG
# -----------------------------------------------------------------------------
def generate_lollipop_plot_score():
    """Generate lollipop plot using total_score instead of ddG"""

    print("\n" + "="*60)
    print("LOLLIPOP PLOT (Total Score) - eIF4E Alanine Scanning")
    print("="*60 + "\n")

    # Load data
    df = load_data()
    print(f"Loaded {len(df)} mutations")

    # Sort by position
    df_sorted = df.sort_values(by='position').reset_index(drop=True)

    # Define colors based on total_score (more negative = better)
    # Invert: worse scores (less negative) get red
    score_range = df_sorted['total_score'].max() - df_sorted['total_score'].min()
    threshold_score_high = df_sorted['total_score'].quantile(0.75)  # Top 25% (worst)
    threshold_score_medium = df_sorted['total_score'].quantile(0.50)  # Middle 50%

    colors = []
    for score in df_sorted['total_score']:
        if score >= threshold_score_high:
            colors.append(color_high)  # Worst structures
        elif score >= threshold_score_medium:
            colors.append(color_medium)
        else:
            colors.append(color_low)  # Best structures

    # Print statistics
    print(f"\nTotal Score Statistics:")
    print(f"  Mean:   {df_sorted['total_score'].mean():.3f} REU")
    print(f"  Median: {df_sorted['total_score'].median():.3f} REU")
    print(f"  Min:    {df_sorted['total_score'].min():.3f} REU (best)")
    print(f"  Max:    {df_sorted['total_score'].max():.3f} REU (worst)")
    print()

    # Set the figure size
    plt.figure(figsize=(figure_width, figure_height))

    # Create the stems of the lollipop chart
    (markerline, stemlines, baseline) = plt.stem(
        df_sorted['position'],
        df_sorted['total_score'],
        linefmt='grey',
        markerfmt='o',
        basefmt='grey'
    )

    # Make default markers invisible
    plt.setp(markerline, 'markerfacecolor', 'none', 'markeredgecolor', 'none')

    # Draw custom-colored markers
    plt.scatter(
        df_sorted['position'],
        df_sorted['total_score'],
        c=colors,
        s=marker_size,
        zorder=3,
        edgecolors='black',
        linewidths=0.5
    )

    # Set titles and labels
    plt.title('eIF4E - Alanine Scan Total Score', fontsize=title_fontsize, weight='bold')
    plt.xlabel('Residue Position', fontsize=axis_label_fontsize)
    plt.ylabel('Total Score (REU)', fontsize=axis_label_fontsize)

    # Set tick properties
    max_pos = df_sorted['position'].max()
    plt.xticks(np.arange(0, max_pos + 1, 10), fontsize=tick_label_fontsize)
    plt.yticks(fontsize=tick_label_fontsize)

    # Add horizontal lines for quartiles
    plt.axhline(y=threshold_score_medium, color='orange', linestyle='--', alpha=0.3, linewidth=1)
    plt.axhline(y=threshold_score_high, color='red', linestyle='--', alpha=0.3, linewidth=1)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color_high, edgecolor='black', label='Worst structures (top 25%)'),
        Patch(facecolor=color_medium, edgecolor='black', label='Medium structures (25-75%)'),
        Patch(facecolor=color_low, edgecolor='black', label='Best structures (bottom 25%)')
    ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=14, framealpha=0.9)

    # Add grid
    plt.grid(axis='y', linestyle=':', color='gray', alpha=0.5)

    # Adjust layout
    plt.tight_layout()

    # Save figure
    output_dir = Path(__file__).parent / 'analysis_results' / 'plots'
    output_path = output_dir / 'eif4e_lollipop_plot_score.png'

    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")

    pdf_path = output_dir / 'eif4e_lollipop_plot_score.pdf'
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"PDF saved to: {pdf_path}")

    plt.close()

    print("\nVisualization complete!")
    print("="*60 + "\n")


# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    # Generate both plots
    generate_lollipop_plot()        # ddG version
    generate_lollipop_plot_score()  # Total score version
