"""
Rosetta Flex ddG Protocol Implementation

This module implements the Flex ddG protocol for calculating binding/folding
free energy changes upon mutation.

Reference:
Barlow KA, et al. (2018) Flex ddG: Rosetta Ensemble-Based Estimation of
Changes in Protein-Protein Binding Affinity upon Mutation.
J Phys Chem B. 122(21):5389-5399.
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Dict, List, Union
from dataclasses import dataclass
import yaml


@dataclass
class FlexDdGConfig:
    """Configuration for Flex ddG protocol."""

    # Basic settings
    nstruct: int = 35  # Number of structures per mutation (recommended: 35)
    iterations: int = 3  # Backrub iterations
    repack_radius: float = 8.0  # Repacking shell radius (Å)

    # Protocol settings
    use_backrub: bool = True
    use_talaris2014: bool = False  # Use ref2015 by default
    interface_ddg: bool = False  # Set True for protein-protein interfaces

    # Computational
    num_processors: int = 1
    memory_gb: int = 4

    # Paths
    rosetta_path: Optional[str] = None
    rosetta_database: Optional[str] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'nstruct': self.nstruct,
            'iterations': self.iterations,
            'repack_radius': self.repack_radius,
            'use_backrub': self.use_backrub,
            'use_talaris2014': self.use_talaris2014,
            'interface_ddg': self.interface_ddg,
            'num_processors': self.num_processors,
            'memory_gb': self.memory_gb,
            'rosetta_path': self.rosetta_path,
            'rosetta_database': self.rosetta_database,
        }

    @classmethod
    def from_yaml(cls, yaml_path: Union[str, Path]) -> 'FlexDdGConfig':
        """Load configuration from YAML file."""
        with open(yaml_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        return cls(**config_dict)

    def save_yaml(self, yaml_path: Union[str, Path]):
        """Save configuration to YAML file."""
        with open(yaml_path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)


class FlexDdGProtocol:
    """
    Elegant wrapper for Rosetta Flex ddG protocol.

    This class provides a clean interface for running Flex ddG calculations
    with proper error handling, logging, and result management.
    """

    def __init__(self, config: Optional[FlexDdGConfig] = None):
        """
        Initialize Flex ddG protocol.

        Args:
            config: Protocol configuration. Uses defaults if None.
        """
        self.config = config or FlexDdGConfig()
        self._validate_rosetta_installation()

    def _validate_rosetta_installation(self):
        """Validate that Rosetta is properly installed and accessible."""
        if self.config.rosetta_path is None:
            # Try to find Rosetta in common locations
            rosetta_path = os.environ.get('ROSETTA')
            if rosetta_path:
                self.config.rosetta_path = rosetta_path
            else:
                raise EnvironmentError(
                    "Rosetta path not found. Set ROSETTA environment variable "
                    "or provide rosetta_path in configuration."
                )

        if self.config.rosetta_database is None:
            db_path = Path(self.config.rosetta_path) / "main" / "database"
            if db_path.exists():
                self.config.rosetta_database = str(db_path)
            else:
                raise FileNotFoundError(f"Rosetta database not found at {db_path}")

    def generate_mutation_file(
        self,
        mutations: List[str],
        output_path: Union[str, Path] = "mutations.txt"
    ):
        """
        Generate mutation file for Rosetta.

        Args:
            mutations: List of mutations in format "A123A" (chain, position, new_aa)
            output_path: Path to save mutation file

        Example:
            mutations = ["A123A", "B234A", "A456A"]
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as f:
            f.write("total 1\n")  # Header for single mutation per run
            for mutation in mutations:
                f.write(f"1\n")
                f.write(f"{mutation}\n")

    def prepare_input_structure(
        self,
        pdb_path: Union[str, Path],
        output_dir: Union[str, Path] = "prepared_structures"
    ) -> Path:
        """
        Prepare and relax input structure using Rosetta's relax protocol.

        Args:
            pdb_path: Path to input PDB file
            output_dir: Directory for output

        Returns:
            Path to relaxed structure
        """
        pdb_path = Path(pdb_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        relax_script = output_dir / "relax.xml"
        self._write_relax_script(relax_script)

        cmd = [
            str(Path(self.config.rosetta_path) / "main" / "source" / "bin" / "relax.default.linuxgccrelease"),
            "-in:file:s", str(pdb_path),
            "-parser:protocol", str(relax_script),
            "-out:path:all", str(output_dir),
            "-out:suffix", "_relaxed",
            "-database", self.config.rosetta_database,
            "-nstruct", "1",
        ]

        if not self.config.use_talaris2014:
            cmd.extend(["-beta"])

        return self._run_rosetta_command(cmd, output_dir)

    def _write_relax_script(self, script_path: Path):
        """Write RosettaScripts XML for relaxation."""
        script = """<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="ref15" weights="ref2015.wts"/>
    </SCOREFXNS>
    <RESIDUE_SELECTORS>
    </RESIDUE_SELECTORS>
    <TASKOPERATIONS>
        <RestrictToRepacking name="repack_only"/>
    </TASKOPERATIONS>
    <FILTERS>
    </FILTERS>
    <MOVERS>
        <FastRelax name="relax" scorefxn="ref15" task_operations="repack_only"/>
    </MOVERS>
    <PROTOCOLS>
        <Add mover="relax"/>
    </PROTOCOLS>
</ROSETTASCRIPTS>
"""
        with open(script_path, 'w') as f:
            f.write(script)

    def run_flex_ddg(
        self,
        pdb_path: Union[str, Path],
        mutations: List[str],
        output_dir: Union[str, Path] = "ddg_results",
        cleanup: bool = True
    ) -> Path:
        """
        Run Flex ddG protocol for given mutations.

        Args:
            pdb_path: Path to input PDB structure
            mutations: List of mutations to evaluate
            output_dir: Directory for results
            cleanup: Remove intermediate files if True

        Returns:
            Path to results directory
        """
        pdb_path = Path(pdb_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate mutation file
        mut_file = output_dir / "mutations.txt"
        self.generate_mutation_file(mutations, mut_file)

        # Prepare command
        flex_ddg_exe = Path(self.config.rosetta_path) / "main" / "source" / "bin" / "flex_ddg.default.linuxgccrelease"

        cmd = [
            str(flex_ddg_exe),
            "-in:file:s", str(pdb_path),
            "-ddg:mut_file", str(mut_file),
            "-ddg:iterations", str(self.config.iterations),
            "-ddg:dump_pdbs", "false",
            "-ddg:bbnbrs", str(self.config.iterations),
            "-database", self.config.rosetta_database,
            "-out:path:all", str(output_dir),
            "-nstruct", str(self.config.nstruct),
        ]

        if self.config.interface_ddg:
            cmd.append("-ddg:interface_ddg")

        if not self.config.use_talaris2014:
            cmd.append("-beta")

        # Run protocol
        result_dir = self._run_rosetta_command(cmd, output_dir)

        if cleanup:
            self._cleanup_intermediate_files(output_dir)

        return result_dir

    def _run_rosetta_command(self, cmd: List[str], working_dir: Path) -> Path:
        """
        Execute Rosetta command with proper error handling.

        Args:
            cmd: Command and arguments
            working_dir: Working directory

        Returns:
            Path to output directory
        """
        log_file = working_dir / "rosetta.log"

        try:
            with open(log_file, 'w') as log:
                process = subprocess.run(
                    cmd,
                    cwd=working_dir,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                    text=True
                )
            return working_dir

        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Rosetta command failed. Check log file: {log_file}\n"
                f"Command: {' '.join(cmd)}"
            ) from e

        except FileNotFoundError as e:
            raise FileNotFoundError(
                f"Rosetta executable not found. Check installation path.\n"
                f"Looking for: {cmd[0]}"
            ) from e

    def _cleanup_intermediate_files(self, output_dir: Path):
        """Remove intermediate files to save space."""
        patterns = ["*.pdb", "*.fasc"]
        for pattern in patterns:
            for file in output_dir.glob(pattern):
                if not file.name.endswith("_relaxed.pdb"):
                    file.unlink()

    def parse_results(self, results_dir: Union[str, Path]) -> Dict:
        """
        Parse Flex ddG results from output directory.

        Args:
            results_dir: Directory containing Rosetta output

        Returns:
            Dictionary with parsed results
        """
        results_dir = Path(results_dir)
        score_file = None

        # Find score file
        for sf in results_dir.glob("*.sc"):
            score_file = sf
            break

        if score_file is None:
            raise FileNotFoundError(f"No score file found in {results_dir}")

        # Parse score file
        results = {
            'mutations': [],
            'ddg_values': [],
            'scores': []
        }

        with open(score_file, 'r') as f:
            lines = f.readlines()

            # Skip header lines
            data_lines = [l for l in lines if not l.startswith('SEQUENCE:') and
                         not l.startswith('SCORE:') and l.strip()]

            for line in data_lines:
                fields = line.split()
                if len(fields) >= 2:
                    results['mutations'].append(fields[0])
                    results['ddg_values'].append(float(fields[1]))
                    results['scores'].append(fields)

        return results
