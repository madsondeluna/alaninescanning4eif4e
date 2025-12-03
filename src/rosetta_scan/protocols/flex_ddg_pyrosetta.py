"""
Protocolo Flex ddG usando PyRosetta

Implementação do protocolo Flex ddG para cálculo de ΔΔG ao introduzir mutações
pontuais para alanina (alanine scanning).

Baseado em:
Barlow KA, et al. (2018) Flex ddG: Rosetta Ensemble-Based Estimation of
Changes in Protein-Protein Binding Affinity upon Mutation.
J Phys Chem B. 122(21):5389-5399.

Protocolo:
1. Relaxar estrutura wild-type (FastRelax)
2. Para cada mutação:
   a. Introduzir mutação
   b. Relaxar estrutura mutante (FastRelax)
   c. Calcular ΔΔG = E_mutante - E_wild-type
3. Repetir para nstruct estruturas (ensemble)
4. ΔΔG final = média do ensemble
"""

import pyrosetta
from pyrosetta import rosetta
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    RestrictToRepacking,
    IncludeCurrent,
    RestrictToRepackingRLT
)
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueIndexSelector,
    NeighborhoodResidueSelector
)
from pyrosetta.rosetta.core.scoring import get_score_function
from typing import List, Dict, Optional
from pathlib import Path
import pandas as pd
import numpy as np
from dataclasses import dataclass


@dataclass
class FlexDdGConfig:
    """
    Configuração do protocolo Flex ddG.

    Atributos:
        nstruct: Número de estruturas no ensemble por mutação
        repack_radius: Raio de repacking ao redor da mutação (Å)
        max_minimization_iter: Iterações máximas de minimização
    """
    nstruct: int = 35
    repack_radius: float = 8.0
    max_minimization_iter: int = 200


class FlexDdGPyRosetta:
    """
    Implementação do protocolo Flex ddG usando PyRosetta.

    Este protocolo calcula mudanças na energia livre (ΔΔG) ao introduzir
    mutações pontuais, considerando flexibilidade conformacional através
    de relaxamento completo da estrutura.

    Equação fundamental:
        ΔΔG = ΔG_mutante - ΔG_wild-type

    Onde ΔG é calculado usando a função de energia REF2015 do Rosetta:
        E_total = Σ w_i × E_i

    Componentes principais:
        - fa_atr: van der Waals atrativo
        - fa_rep: van der Waals repulsivo
        - fa_sol: Solvatação
        - fa_elec: Eletrostático
        - hbond: Ligações de hidrogênio
        - hbond_sc: Ligações H cadeia lateral

    Critério de hotspot (literatura):
        ΔΔG ≥ +2.0 kcal/mol (≈ redução de 30x na afinidade)
    """

    def __init__(self, config: Optional[FlexDdGConfig] = None):
        """
        Inicializa o protocolo Flex ddG.

        Args:
            config: Configuração do protocolo. Se None, usa padrões.
        """
        self.config = config or FlexDdGConfig()

        # Inicializar PyRosetta com flags apropriadas
        pyrosetta.init(
            "-ignore_unrecognized_res "
            "-ignore_zero_occupancy false "
            "-ex1 "
            "-ex2 "
            "-use_input_sc "
            "-mute all"
        )

        # Carregar função de energia REF2015
        self.scorefxn = get_score_function()

        print(f"PyRosetta inicializado")
        print(f"Função de energia: REF2015")
        print(f"Configuração:")
        print(f"  - nstruct: {self.config.nstruct}")
        print(f"  - repack_radius: {self.config.repack_radius} Å")
        print(f"  - max_minimization_iter: {self.config.max_minimization_iter}")

    def run_flex_ddg(
        self,
        pdb_path: str,
        mutations: List[Dict],
        output_dir: str = "results"
    ) -> pd.DataFrame:
        """
        Executa protocolo Flex ddG para lista de mutações.

        Args:
            pdb_path: Caminho para arquivo PDB
            mutations: Lista de dicts com informações de mutações
            output_dir: Diretório de saída

        Returns:
            DataFrame com resultados de ΔΔG
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)

        print(f"\nCarregando estrutura: {pdb_path}")
        pose_wt = pyrosetta.pose_from_pdb(pdb_path)

        # Relaxar estrutura wild-type
        print("Relaxando estrutura wild-type...")
        self._relax_structure(pose_wt)
        E_wt = self.scorefxn(pose_wt)

        print(f"Energia wild-type: {E_wt:.2f} REU")
        print(f"\nProcessando {len(mutations)} mutações...")

        results = []

        for i, mut_info in enumerate(mutations, 1):
            print(f"\n[{i}/{len(mutations)}] {mut_info['mutation']}", flush=True)

            ddg_result = self._calculate_single_mutation(
                pose_wt,
                mut_info,
                E_wt
            )

            results.append({
                'mutation': mut_info['mutation'],
                'position': mut_info['position'],
                'chain': mut_info.get('chain', 'A'),
                'original_aa': mut_info['original_aa'],
                'ddg': ddg_result['ddg'],
                'ddg_std': ddg_result['std'],
                'E_wt': E_wt,
                'E_mut': ddg_result['E_mut']
            })

            print(f"  ΔΔG = {ddg_result['ddg']:+.2f} ± {ddg_result['std']:.2f} kcal/mol")

        df_results = pd.DataFrame(results)

        # Salvar resultados
        output_file = output_dir / "ddg_results.csv"
        df_results.to_csv(output_file, index=False)
        print(f"\nResultados salvos em: {output_file}")

        # Identificar hotspots (ΔΔG > +2.0 kcal/mol)
        hotspots = df_results[df_results['ddg'] > 2.0].copy()
        hotspots = hotspots.sort_values('ddg', ascending=False)

        if len(hotspots) > 0:
            hotspots_file = output_dir / "hotspots.csv"
            hotspots.to_csv(hotspots_file, index=False)
            print(f"\nHotspots identificados: {len(hotspots)}")
            print(f"Hotspots salvos em: {hotspots_file}")

            print("\nTop 10 hotspots:")
            for idx, row in hotspots.head(10).iterrows():
                print(f"  {row['mutation']:8s} ΔΔG = {row['ddg']:+.2f} kcal/mol")

        return df_results

    def _calculate_single_mutation(
        self,
        pose_wt: pyrosetta.Pose,
        mut_info: Dict,
        E_wt: float
    ) -> Dict:
        """
        Calcula ΔΔG para uma única mutação usando ensemble.

        Protocolo Flex ddG:
        1. Para cada estrutura no ensemble (nstruct):
           a. Clonar wild-type
           b. Introduzir mutação
           c. Relaxar estrutura mutante
           d. Calcular ΔΔG_i = E_mutante - E_wild-type
        2. ΔΔG = média(ΔΔG_i)
        3. std = desvio_padrão(ΔΔG_i)

        Args:
            pose_wt: Estrutura wild-type relaxada
            mut_info: Dict com informação da mutação
            E_wt: Energia wild-type

        Returns:
            Dict com 'ddg', 'std', 'E_mut'
        """
        position = mut_info['position']
        original_aa = mut_info['original_aa']

        ddg_values = []
        E_mut_values = []

        # Gerar ensemble de estruturas mutantes
        for n in range(self.config.nstruct):
            # Clonar pose
            pose_mut = pose_wt.clone()

            # Introduzir mutação para alanina
            self._mutate_residue(pose_mut, position, 'ALA')

            # Relaxar estrutura mutante
            self._relax_structure(pose_mut)

            # Calcular energia mutante
            E_mut = self.scorefxn(pose_mut)

            # Calcular ΔΔG para esta estrutura
            ddg_i = E_mut - E_wt

            ddg_values.append(ddg_i)
            E_mut_values.append(E_mut)

            if (n + 1) % 5 == 0:
                print(f"  Estrutura {n+1}/{self.config.nstruct}: ΔΔG = {ddg_i:+.2f} kcal/mol", flush=True)

        # Calcular estatísticas do ensemble
        ddg_mean = np.mean(ddg_values)
        ddg_std = np.std(ddg_values)
        E_mut_mean = np.mean(E_mut_values)

        return {
            'ddg': ddg_mean,
            'std': ddg_std,
            'E_mut': E_mut_mean
        }

    def _relax_structure(self, pose: pyrosetta.Pose):
        """
        Relaxa estrutura usando FastRelax.

        FastRelax é o protocolo padrão do Rosetta para relaxamento de estruturas,
        alternando entre repacking de cadeias laterais e minimização.

        Args:
            pose: Estrutura a ser relaxada (modificada in-place)
        """
        # Configurar FastRelax
        fast_relax = FastRelax()
        fast_relax.set_scorefxn(self.scorefxn)
        fast_relax.max_iter(self.config.max_minimization_iter)

        # Aplicar relaxamento
        fast_relax.apply(pose)

    def _mutate_residue(
        self,
        pose: pyrosetta.Pose,
        position: int,
        target_aa: str
    ):
        """
        Introduz mutação pontual e reempacota ao redor.

        Args:
            pose: Estrutura (modificada in-place)
            position: Posição do resíduo a mutar
            target_aa: Código de 3 letras do aminoácido alvo (ex: 'ALA')
        """
        # Converter código de 3 letras para código Rosetta
        aa_name3_to_aa = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }

        target_aa_1letter = aa_name3_to_aa[target_aa]

        # Criar mutação
        mutate = rosetta.protocols.simple_moves.MutateResidue(
            target=position,
            new_res=target_aa_1letter
        )
        mutate.apply(pose)

        # Reempacotar ao redor da mutação
        self._repack_around_mutation(pose, position)

    def _repack_around_mutation(
        self,
        pose: pyrosetta.Pose,
        position: int
    ):
        """
        Reempacota cadeias laterais ao redor da mutação.

        Args:
            pose: Estrutura (modificada in-place)
            position: Posição central
        """
        # Selecionar resíduo mutado
        mutated_selector = ResidueIndexSelector(str(position))

        # Selecionar vizinhos dentro do raio
        neighborhood_selector = NeighborhoodResidueSelector(
            mutated_selector,
            self.config.repack_radius,
            include_focus_in_subset=True
        )

        # Configurar TaskFactory
        tf = TaskFactory()
        tf.push_back(RestrictToRepackingRLT())
        tf.push_back(IncludeCurrent())

        # Criar PackRotamersMover
        packer = PackRotamersMover(self.scorefxn)
        task = tf.create_task_and_apply_taskoperations(pose)

        # Restringir apenas aos resíduos selecionados
        subset = neighborhood_selector.apply(pose)
        for i in range(1, pose.total_residue() + 1):
            if not subset[i]:
                task.nonconst_residue_task(i).prevent_repacking()

        packer.task(task)
        packer.apply(pose)
