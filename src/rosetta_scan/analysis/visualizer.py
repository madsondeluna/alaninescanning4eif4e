"""
Visualization utilities for alanine scanning results.
"""

from pathlib import Path
from typing import Optional, Union, List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


class ResultVisualizer:
    """Create elegant visualizations for alanine scanning results."""

    def __init__(self, results_df: pd.DataFrame):
        """
        Initialize visualizer.

        Args:
            results_df: DataFrame with results from ResultParser
        """
        self.results_df = results_df
        self._setup_style()

    def _setup_style(self):
        """Setup matplotlib style."""
        sns.set_style("whitegrid")
        sns.set_context("paper", font_scale=1.5)

        # Color palette
        self.colors = {
            'primary': '#2E86AB',
            'secondary': '#A23B72',
            'accent': '#F18F01',
            'positive': '#C73E1D',
            'negative': '#6A4C93'
        }

    def plot_ddg_distribution(
        self,
        output_path: Optional[Union[str, Path]] = None,
        figsize: tuple = (10, 6)
    ):
        """
        Plot ddG distribution histogram.

        Args:
            output_path: Path to save figure
            figsize: Figure size (width, height)
        """
        fig, ax = plt.subplots(figsize=figsize)

        # Plot histogram
        ax.hist(
            self.results_df['ddg'],
            bins=30,
            color=self.colors['primary'],
            alpha=0.7,
            edgecolor='black'
        )

        # Add vertical line at mean
        mean_ddg = self.results_df['ddg'].mean()
        ax.axvline(
            mean_ddg,
            color=self.colors['accent'],
            linestyle='--',
            linewidth=2,
            label=f'Mean: {mean_ddg:.2f} kcal/mol'
        )

        # Styling
        ax.set_xlabel('ΔΔG (kcal/mol)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Frequency', fontsize=14, fontweight='bold')
        ax.set_title('Alanine Scanning ΔΔG Distribution', fontsize=16, fontweight='bold')
        ax.legend(fontsize=12)
        ax.grid(alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_hotspot_heatmap(
        self,
        output_path: Optional[Union[str, Path]] = None,
        figsize: tuple = (12, 8),
        threshold: float = 1.0
    ):
        """
        Plot heatmap of ddG values by position and chain.

        Args:
            output_path: Path to save figure
            figsize: Figure size
            threshold: Threshold for highlighting hotspots
        """
        if 'chain' not in self.results_df.columns or 'position' not in self.results_df.columns:
            raise ValueError("Chain and position information required")

        # Pivot data
        pivot_data = self.results_df.pivot_table(
            values='ddg',
            index='position',
            columns='chain',
            aggfunc='mean'
        )

        fig, ax = plt.subplots(figsize=figsize)

        # Create custom colormap (blue-white-red)
        cmap = LinearSegmentedColormap.from_list(
            'custom',
            ['#6A4C93', '#FFFFFF', '#C73E1D']
        )

        # Plot heatmap
        sns.heatmap(
            pivot_data,
            cmap=cmap,
            center=0,
            cbar_kws={'label': 'ΔΔG (kcal/mol)'},
            ax=ax,
            linewidths=0.5,
            linecolor='lightgray'
        )

        ax.set_xlabel('Chain', fontsize=14, fontweight='bold')
        ax.set_ylabel('Residue Position', fontsize=14, fontweight='bold')
        ax.set_title('Alanine Scanning Hotspot Map', fontsize=16, fontweight='bold')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_per_chain_analysis(
        self,
        output_path: Optional[Union[str, Path]] = None,
        figsize: tuple = (14, 6)
    ):
        """
        Plot per-chain ddG analysis.

        Args:
            output_path: Path to save figure
            figsize: Figure size
        """
        if 'chain' not in self.results_df.columns:
            raise ValueError("Chain information required")

        fig, axes = plt.subplots(1, 2, figsize=figsize)

        # Box plot
        sns.boxplot(
            data=self.results_df,
            x='chain',
            y='ddg',
            palette='Set2',
            ax=axes[0]
        )

        axes[0].set_xlabel('Chain', fontsize=12, fontweight='bold')
        axes[0].set_ylabel('ΔΔG (kcal/mol)', fontsize=12, fontweight='bold')
        axes[0].set_title('ΔΔG Distribution by Chain', fontsize=14, fontweight='bold')
        axes[0].grid(alpha=0.3, axis='y')

        # Violin plot
        sns.violinplot(
            data=self.results_df,
            x='chain',
            y='ddg',
            palette='Set2',
            ax=axes[1]
        )

        axes[1].set_xlabel('Chain', fontsize=12, fontweight='bold')
        axes[1].set_ylabel('ΔΔG (kcal/mol)', fontsize=12, fontweight='bold')
        axes[1].set_title('ΔΔG Density by Chain', fontsize=14, fontweight='bold')
        axes[1].grid(alpha=0.3, axis='y')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_top_hotspots(
        self,
        n: int = 20,
        output_path: Optional[Union[str, Path]] = None,
        figsize: tuple = (10, 8)
    ):
        """
        Plot top N hotspots as bar chart.

        Args:
            n: Number of top hotspots to show
            output_path: Path to save figure
            figsize: Figure size
        """
        # Get top hotspots by absolute ddG
        top_data = self.results_df.nlargest(n, 'ddg').copy()
        top_data = top_data.sort_values('ddg')

        fig, ax = plt.subplots(figsize=figsize)

        # Color bars by sign
        colors = [self.colors['positive'] if x > 0 else self.colors['negative']
                 for x in top_data['ddg']]

        # Plot horizontal bar chart
        ax.barh(
            range(len(top_data)),
            top_data['ddg'],
            color=colors,
            alpha=0.7,
            edgecolor='black'
        )

        # Set labels
        labels = []
        for _, row in top_data.iterrows():
            mutation = row.get('mutation', 'Unknown')
            chain = row.get('chain', '')
            labels.append(f"{mutation} ({chain})")

        ax.set_yticks(range(len(top_data)))
        ax.set_yticklabels(labels, fontsize=10)
        ax.set_xlabel('ΔΔG (kcal/mol)', fontsize=12, fontweight='bold')
        ax.set_title(f'Top {n} Hotspot Residues', fontsize=14, fontweight='bold')
        ax.axvline(0, color='black', linewidth=0.8)
        ax.grid(alpha=0.3, axis='x')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_position_scan(
        self,
        chain: str,
        output_path: Optional[Union[str, Path]] = None,
        figsize: tuple = (14, 6),
        threshold: float = 1.0
    ):
        """
        Plot ddG values along sequence for a specific chain.

        Args:
            chain: Chain ID to plot
            output_path: Path to save figure
            figsize: Figure size
            threshold: Threshold line for hotspots
        """
        if 'chain' not in self.results_df.columns or 'position' not in self.results_df.columns:
            raise ValueError("Chain and position information required")

        # Filter by chain
        chain_data = self.results_df[self.results_df['chain'] == chain].copy()
        chain_data = chain_data.sort_values('position')

        fig, ax = plt.subplots(figsize=figsize)

        # Plot line
        ax.plot(
            chain_data['position'],
            chain_data['ddg'],
            marker='o',
            linewidth=2,
            markersize=6,
            color=self.colors['primary'],
            alpha=0.7
        )

        # Highlight hotspots
        hotspots = chain_data[chain_data['ddg'].abs() >= threshold]
        ax.scatter(
            hotspots['position'],
            hotspots['ddg'],
            color=self.colors['accent'],
            s=100,
            zorder=5,
            label=f'Hotspots (|ΔΔG| ≥ {threshold})'
        )

        # Add threshold lines
        ax.axhline(threshold, color='red', linestyle='--', alpha=0.5, label=f'Threshold: ±{threshold}')
        ax.axhline(-threshold, color='red', linestyle='--', alpha=0.5)
        ax.axhline(0, color='black', linewidth=0.8)

        ax.set_xlabel('Residue Position', fontsize=12, fontweight='bold')
        ax.set_ylabel('ΔΔG (kcal/mol)', fontsize=12, fontweight='bold')
        ax.set_title(f'Alanine Scan Along Chain {chain}', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def create_dashboard(
        self,
        output_dir: Union[str, Path],
        threshold: float = 1.0
    ):
        """
        Create complete dashboard with all visualizations.

        Args:
            output_dir: Directory to save plots
            threshold: Hotspot threshold
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate all plots
        self.plot_ddg_distribution(output_dir / 'ddg_distribution.png')
        self.plot_top_hotspots(output_path=output_dir / 'top_hotspots.png')

        if 'chain' in self.results_df.columns:
            self.plot_per_chain_analysis(output_dir / 'chain_analysis.png')

            if 'position' in self.results_df.columns:
                self.plot_hotspot_heatmap(output_dir / 'hotspot_heatmap.png', threshold=threshold)

                # Plot per chain position scans
                for chain in self.results_df['chain'].unique():
                    self.plot_position_scan(
                        chain,
                        output_dir / f'position_scan_chain_{chain}.png',
                        threshold=threshold
                    )
