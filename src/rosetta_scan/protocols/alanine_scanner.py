"""
Alanine Scanning Mutagenesis Module

Generates systematic alanine mutations for hotspot identification.
"""

from pathlib import Path
from typing import List, Optional, Set, Dict, Union
from dataclasses import dataclass
from Bio import PDB
from Bio.PDB.Structure import Structure
from Bio.PDB.Polypeptide import is_aa, three_to_one
import re


@dataclass
class MutationSite:
    """Represents a mutation site."""

    chain: str
    position: int
    original_aa: str
    target_aa: str = "A"

    def to_rosetta_format(self) -> str:
        """
        Convert to Rosetta mutation format.

        Returns:
            String like "A123A" (chain, position, new_aa)
        """
        return f"{self.chain}{self.position}{self.target_aa}"

    def __str__(self) -> str:
        """Human-readable representation."""
        return f"{self.original_aa}{self.position}{self.target_aa} (chain {self.chain})"


class AlanineScan:
    """
    Generate alanine scanning mutations for hotspot identification.

    This class analyzes protein structures and generates systematic
    alanine mutations for computational hotspot mapping.
    """

    # Amino acids to exclude from scanning (already Ala, Gly, Pro)
    EXCLUDED_AA = {'ALA', 'GLY', 'PRO'}
    EXCLUDED_1LETTER = {'A', 'G', 'P'}

    def __init__(self, pdb_path: Union[str, Path]):
        """
        Initialize alanine scanner.

        Args:
            pdb_path: Path to PDB structure file
        """
        self.pdb_path = Path(pdb_path)
        self.structure = self._load_structure()
        self.mutation_sites: List[MutationSite] = []

    def _load_structure(self) -> Structure:
        """Load PDB structure using BioPython."""
        parser = PDB.PDBParser(QUIET=True)
        return parser.get_structure("protein", str(self.pdb_path))

    def generate_mutations(
        self,
        chains: Optional[List[str]] = None,
        residue_range: Optional[Dict[str, tuple]] = None,
        exclude_residues: Optional[Set[int]] = None,
        interface_only: bool = False,
        interface_cutoff: float = 8.0
    ) -> List[MutationSite]:
        """
        Generate alanine scanning mutations.

        Args:
            chains: List of chain IDs to scan. If None, scan all chains.
            residue_range: Dict mapping chain IDs to (start, end) tuples.
                          Example: {'A': (1, 100), 'B': (50, 150)}
            exclude_residues: Set of residue positions to exclude
            interface_only: Only scan interface residues if True
            interface_cutoff: Distance cutoff (Å) for interface definition

        Returns:
            List of MutationSite objects
        """
        self.mutation_sites = []
        exclude_residues = exclude_residues or set()

        # Get chains to scan
        if chains is None:
            chains = [chain.id for chain in self.structure[0]]

        # Build interface residues if needed
        interface_residues = None
        if interface_only:
            interface_residues = self._identify_interface_residues(
                chains, interface_cutoff
            )

        # Scan each chain
        for chain in self.structure[0]:
            if chain.id not in chains:
                continue

            # Get residue range for this chain
            if residue_range and chain.id in residue_range:
                start, end = residue_range[chain.id]
            else:
                start, end = None, None

            # Scan residues
            for residue in chain:
                # Skip heteroatoms
                if residue.id[0] != ' ':
                    continue

                res_num = residue.id[1]

                # Check range
                if start is not None and res_num < start:
                    continue
                if end is not None and res_num > end:
                    continue

                # Check exclusion list
                if res_num in exclude_residues:
                    continue

                # Check if amino acid
                if not is_aa(residue, standard=True):
                    continue

                # Get residue name
                res_name = residue.resname

                # Skip excluded amino acids
                if res_name in self.EXCLUDED_AA:
                    continue

                # Check interface if needed
                if interface_only and interface_residues:
                    if (chain.id, res_num) not in interface_residues:
                        continue

                # Convert to one-letter code
                try:
                    original_aa = three_to_one(res_name)
                except KeyError:
                    continue

                # Create mutation site
                mutation = MutationSite(
                    chain=chain.id,
                    position=res_num,
                    original_aa=original_aa,
                    target_aa='A'
                )

                self.mutation_sites.append(mutation)

        return self.mutation_sites

    def _identify_interface_residues(
        self,
        chains: List[str],
        cutoff: float
    ) -> Set[tuple]:
        """
        Identify interface residues between chains.

        Args:
            chains: List of chain IDs
            cutoff: Distance cutoff in Ångströms

        Returns:
            Set of (chain_id, residue_number) tuples
        """
        interface_residues = set()

        # Get atoms for each chain
        chain_atoms = {}
        for chain in self.structure[0]:
            if chain.id in chains:
                atoms = [atom for residue in chain for atom in residue
                        if residue.id[0] == ' ']
                chain_atoms[chain.id] = atoms

        # Find interface residues
        chain_ids = list(chain_atoms.keys())
        for i, chain_a in enumerate(chain_ids):
            for chain_b in chain_ids[i+1:]:
                for atom_a in chain_atoms[chain_a]:
                    for atom_b in chain_atoms[chain_b]:
                        distance = atom_a - atom_b
                        if distance <= cutoff:
                            # Add both residues
                            res_a = atom_a.get_parent()
                            res_b = atom_b.get_parent()
                            interface_residues.add((chain_a, res_a.id[1]))
                            interface_residues.add((chain_b, res_b.id[1]))

        return interface_residues

    def get_rosetta_mutation_list(self) -> List[str]:
        """
        Get mutations in Rosetta format.

        Returns:
            List of mutation strings (e.g., ["A123A", "B234A"])
        """
        return [mut.to_rosetta_format() for mut in self.mutation_sites]

    def save_mutation_report(
        self,
        output_path: Union[str, Path],
        format: str = 'txt'
    ):
        """
        Save mutation report to file.

        Args:
            output_path: Path to output file
            format: Output format ('txt', 'csv', or 'rosetta')
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format == 'txt':
            self._save_txt_report(output_path)
        elif format == 'csv':
            self._save_csv_report(output_path)
        elif format == 'rosetta':
            self._save_rosetta_format(output_path)
        else:
            raise ValueError(f"Unknown format: {format}")

    def _save_txt_report(self, output_path: Path):
        """Save human-readable text report."""
        with open(output_path, 'w') as f:
            f.write(f"Alanine Scanning Mutation Report\n")
            f.write(f"{'='*50}\n\n")
            f.write(f"Structure: {self.pdb_path.name}\n")
            f.write(f"Total mutations: {len(self.mutation_sites)}\n\n")
            f.write(f"{'Chain':<8}{'Position':<10}{'Original':<10}{'Target':<10}\n")
            f.write(f"{'-'*50}\n")

            for mut in self.mutation_sites:
                f.write(f"{mut.chain:<8}{mut.position:<10}"
                       f"{mut.original_aa:<10}{mut.target_aa:<10}\n")

    def _save_csv_report(self, output_path: Path):
        """Save CSV format report."""
        with open(output_path, 'w') as f:
            f.write("chain,position,original_aa,target_aa,mutation\n")
            for mut in self.mutation_sites:
                f.write(f"{mut.chain},{mut.position},{mut.original_aa},"
                       f"{mut.target_aa},{mut}\n")

    def _save_rosetta_format(self, output_path: Path):
        """Save Rosetta-compatible mutation file."""
        with open(output_path, 'w') as f:
            f.write("total 1\n")
            for mut in self.mutation_sites:
                f.write("1\n")
                f.write(f"{mut.to_rosetta_format()}\n")

    def filter_by_solvent_accessibility(
        self,
        sasa_threshold: float = 0.2,
        probe_radius: float = 1.4
    ) -> List[MutationSite]:
        """
        Filter mutations to include only surface-exposed residues.

        Args:
            sasa_threshold: Relative SASA threshold (0-1)
            probe_radius: Probe radius for SASA calculation (Å)

        Returns:
            Filtered list of MutationSite objects
        """
        # Calculate SASA
        calculator = PDB.SASA.ShrakeRupley()
        calculator.compute(self.structure, level="R")

        filtered_mutations = []

        for mut in self.mutation_sites:
            # Get residue
            chain = self.structure[0][mut.chain]
            residue = chain[(' ', mut.position, ' ')]

            # Check SASA
            if hasattr(residue, 'sasa'):
                # Get max SASA for this amino acid type
                max_sasa = self._get_max_sasa(mut.original_aa)
                relative_sasa = residue.sasa / max_sasa

                if relative_sasa >= sasa_threshold:
                    filtered_mutations.append(mut)
            else:
                # Include if SASA not calculated
                filtered_mutations.append(mut)

        self.mutation_sites = filtered_mutations
        return filtered_mutations

    def _get_max_sasa(self, aa: str) -> float:
        """
        Get maximum SASA for amino acid type.

        Values from Tien et al. (2013) PLoS ONE 8(11): e80635
        """
        max_sasa = {
            'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0,
            'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0,
            'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0,
            'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0,
            'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0
        }
        return max_sasa.get(aa, 200.0)

    def get_mutation_summary(self) -> Dict:
        """
        Get summary statistics for mutations.

        Returns:
            Dictionary with summary information
        """
        chains = {}
        for mut in self.mutation_sites:
            if mut.chain not in chains:
                chains[mut.chain] = 0
            chains[mut.chain] += 1

        aa_counts = {}
        for mut in self.mutation_sites:
            if mut.original_aa not in aa_counts:
                aa_counts[mut.original_aa] = 0
            aa_counts[mut.original_aa] += 1

        return {
            'total_mutations': len(self.mutation_sites),
            'chains': chains,
            'amino_acid_distribution': aa_counts,
            'structure_file': str(self.pdb_path)
        }
