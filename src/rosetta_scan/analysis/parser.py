"""
Result parsing utilities for Flex ddG output.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Union
import pandas as pd
import numpy as np


class ResultParser:
    """Parse and process Rosetta Flex ddG results."""

    def __init__(self, results_dir: Union[str, Path]):
        """
        Initialize result parser.

        Args:
            results_dir: Directory containing Rosetta output files
        """
        self.results_dir = Path(results_dir)
        self.results_df: Optional[pd.DataFrame] = None

    def parse_results(self) -> pd.DataFrame:
        """
        Parse Flex ddG results from score files.

        Returns:
            DataFrame with parsed results
        """
        # Find score files
        score_files = list(self.results_dir.glob("*.sc"))

        if not score_files:
            raise FileNotFoundError(
                f"No score files (*.sc) found in {self.results_dir}"
            )

        all_results = []

        for score_file in score_files:
            results = self._parse_score_file(score_file)
            all_results.extend(results)

        # Create DataFrame
        self.results_df = pd.DataFrame(all_results)

        # Add derived columns
        if len(self.results_df) > 0:
            self._add_derived_columns()

        return self.results_df

    def _parse_score_file(self, score_file: Path) -> List[Dict]:
        """Parse a single score file."""
        results = []

        with open(score_file, 'r') as f:
            lines = f.readlines()

        # Find header line
        header_idx = None
        for i, line in enumerate(lines):
            if line.startswith('SCORE:') and 'description' in line:
                header_idx = i
                break

        if header_idx is None:
            return results

        # Parse header
        header = lines[header_idx].split()[1:]  # Skip 'SCORE:'

        # Parse data lines
        for line in lines[header_idx + 1:]:
            if not line.startswith('SCORE:'):
                continue

            fields = line.split()[1:]  # Skip 'SCORE:'

            if len(fields) != len(header):
                continue

            # Create result dict
            result = {}
            for key, value in zip(header, fields):
                # Try to convert to float
                try:
                    result[key] = float(value)
                except ValueError:
                    result[key] = value

            results.append(result)

        return results

    def _add_derived_columns(self):
        """Add derived columns to results DataFrame."""
        # Extract mutation info
        if 'description' in self.results_df.columns:
            self.results_df['mutation'] = self.results_df['description'].apply(
                self._extract_mutation
            )

            self.results_df['chain'] = self.results_df['mutation'].apply(
                lambda x: x[0] if x else None
            )

            self.results_df['position'] = self.results_df['mutation'].apply(
                lambda x: int(re.search(r'\d+', x).group()) if x and re.search(r'\d+', x) else None
            )

        # Rename ddG column if exists
        ddg_columns = ['ddg', 'dG', 'total_score']
        for col in ddg_columns:
            if col in self.results_df.columns:
                if col != 'ddg':
                    self.results_df.rename(columns={col: 'ddg'}, inplace=True)
                break

    def _extract_mutation(self, description: str) -> Optional[str]:
        """Extract mutation string from description."""
        # Look for pattern like A123A
        match = re.search(r'[A-Z]\d+[A-Z]', description)
        if match:
            return match.group()
        return None

    def identify_hotspots(
        self,
        threshold: float = 1.0,
        metric: str = 'ddg'
    ) -> pd.DataFrame:
        """
        Identify hotspot residues based on ddG threshold.

        Args:
            threshold: ddG threshold (kcal/mol)
            metric: Column to use for threshold

        Returns:
            DataFrame with hotspot residues
        """
        if self.results_df is None:
            raise ValueError("No results loaded. Call parse_results() first.")

        if metric not in self.results_df.columns:
            raise ValueError(f"Metric '{metric}' not found in results")

        # Filter by threshold
        hotspots = self.results_df[
            self.results_df[metric].abs() >= threshold
        ].copy()

        # Sort by ddG magnitude
        hotspots['abs_ddg'] = hotspots[metric].abs()
        hotspots = hotspots.sort_values('abs_ddg', ascending=False)

        return hotspots

    def get_per_chain_statistics(self) -> pd.DataFrame:
        """
        Calculate statistics per chain.

        Returns:
            DataFrame with per-chain statistics
        """
        if self.results_df is None:
            raise ValueError("No results loaded. Call parse_results() first.")

        if 'chain' not in self.results_df.columns:
            raise ValueError("Chain information not available")

        stats = self.results_df.groupby('chain')['ddg'].agg([
            'count',
            'mean',
            'std',
            'min',
            'max'
        ]).reset_index()

        stats.columns = ['chain', 'n_mutations', 'mean_ddg', 'std_ddg', 'min_ddg', 'max_ddg']

        return stats

    def get_top_mutations(
        self,
        n: int = 10,
        by: str = 'ddg',
        ascending: bool = False
    ) -> pd.DataFrame:
        """
        Get top N mutations by specified metric.

        Args:
            n: Number of top mutations to return
            by: Column to sort by
            ascending: Sort order

        Returns:
            DataFrame with top mutations
        """
        if self.results_df is None:
            raise ValueError("No results loaded. Call parse_results() first.")

        return self.results_df.nlargest(n, by) if not ascending else \
               self.results_df.nsmallest(n, by)

    def export_pymol_script(
        self,
        output_path: Union[str, Path],
        threshold: float = 1.0,
        color_scheme: str = 'red_blue'
    ):
        """
        Export PyMOL script to visualize hotspots.

        Args:
            output_path: Path to save PyMOL script
            threshold: ddG threshold for highlighting
            color_scheme: Color scheme ('red_blue', 'gradient')
        """
        if self.results_df is None:
            raise ValueError("No results loaded. Call parse_results() first.")

        hotspots = self.identify_hotspots(threshold=threshold)

        script_lines = [
            "# PyMOL script for visualizing alanine scanning hotspots",
            "# Generated by Rosetta Alanine Scanning",
            "",
            "# Hide everything",
            "hide everything",
            "show cartoon",
            "",
            "# Color by chain",
            "util.cbc",
            "",
            "# Highlight hotspots",
        ]

        for _, row in hotspots.iterrows():
            chain = row.get('chain', 'A')
            position = row.get('position', 0)
            ddg = row.get('ddg', 0)

            # Color based on ddG
            if ddg > 0:
                color = "red"
            else:
                color = "blue"

            script_lines.extend([
                f"select hotspot_{chain}{position}, chain {chain} and resi {position}",
                f"show sticks, hotspot_{chain}{position}",
                f"color {color}, hotspot_{chain}{position}",
            ])

        script_lines.extend([
            "",
            "# Label hotspots",
            "label hotspot_*, '%s%s' % (resn, resi)",
            "",
            "# Zoom to hotspots",
            "zoom hotspot_*",
        ])

        with open(output_path, 'w') as f:
            f.write('\n'.join(script_lines))

    def export_summary_report(
        self,
        output_path: Union[str, Path],
        threshold: float = 1.0
    ):
        """
        Export comprehensive summary report.

        Args:
            output_path: Path to save report
            threshold: Hotspot threshold
        """
        if self.results_df is None:
            raise ValueError("No results loaded. Call parse_results() first.")

        hotspots = self.identify_hotspots(threshold=threshold)

        with open(output_path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("Rosetta Alanine Scanning - Summary Report\n")
            f.write("=" * 80 + "\n\n")

            # Overall statistics
            f.write("OVERALL STATISTICS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Total mutations analyzed: {len(self.results_df)}\n")
            f.write(f"Mean ddG: {self.results_df['ddg'].mean():.2f} kcal/mol\n")
            f.write(f"Std ddG: {self.results_df['ddg'].std():.2f} kcal/mol\n")
            f.write(f"Min ddG: {self.results_df['ddg'].min():.2f} kcal/mol\n")
            f.write(f"Max ddG: {self.results_df['ddg'].max():.2f} kcal/mol\n")
            f.write("\n")

            # Hotspots
            f.write(f"IDENTIFIED HOTSPOTS (ddG >= {threshold} kcal/mol)\n")
            f.write("-" * 80 + "\n")
            f.write(f"Number of hotspots: {len(hotspots)}\n\n")

            if len(hotspots) > 0:
                f.write(f"{'Rank':<6}{'Mutation':<12}{'Chain':<8}{'Position':<10}{'ddG (kcal/mol)':<15}\n")
                f.write("-" * 80 + "\n")

                for i, (_, row) in enumerate(hotspots.head(20).iterrows(), 1):
                    f.write(
                        f"{i:<6}{row['mutation']:<12}{row.get('chain', 'N/A'):<8}"
                        f"{row.get('position', 'N/A'):<10}{row['ddg']:<15.2f}\n"
                    )

            f.write("\n")

            # Per-chain statistics
            if 'chain' in self.results_df.columns:
                f.write("PER-CHAIN STATISTICS\n")
                f.write("-" * 80 + "\n")

                chain_stats = self.get_per_chain_statistics()

                f.write(f"{'Chain':<8}{'N':<8}{'Mean ddG':<15}{'Std ddG':<15}{'Min ddG':<15}{'Max ddG':<15}\n")
                f.write("-" * 80 + "\n")

                for _, row in chain_stats.iterrows():
                    f.write(
                        f"{row['chain']:<8}{row['n_mutations']:<8}{row['mean_ddg']:<15.2f}"
                        f"{row['std_ddg']:<15.2f}{row['min_ddg']:<15.2f}{row['max_ddg']:<15.2f}\n"
                    )

            f.write("\n")
            f.write("=" * 80 + "\n")
