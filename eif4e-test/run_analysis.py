#!/usr/bin/env python3
"""
Análise de Alanine Scanning para eIF4E usando PyRosetta

Este script executa o protocolo Flex ddG para calcular ΔΔG de mutações
pontuais para alanina na proteína eIF4E.
"""

import sys
from pathlib import Path

# Adicionar diretório src ao path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from rosetta_scan.protocols.flex_ddg_pyrosetta import (
    FlexDdGPyRosetta,
    FlexDdGConfig
)
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
import pandas as pd


def gerar_mutacoes(pdb_path: str) -> list:
    """
    Gera lista de mutações para alanina a partir do PDB.

    Args:
        pdb_path: Caminho para arquivo PDB

    Returns:
        Lista de dicts com informações de mutações
    """
    # Dicionário de conversão 3-letras para 1-letra
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)

    mutations = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True):
                    resname = residue.get_resname()
                    if resname in aa_dict:
                        aa = aa_dict[resname]
                        pos = residue.get_id()[1]

                        # Skip GLY, ALA, PRO
                        if aa not in ['G', 'A', 'P']:
                            mutations.append({
                                'mutation': f'{aa}{pos}A',
                                'position': pos,
                                'chain': chain.id,
                                'original_aa': aa
                            })

    return mutations


def main():
    """Pipeline principal de análise."""
    print("="*80)
    print("ALANINE SCANNING - eIF4E")
    print("Protocolo Flex ddG com PyRosetta")
    print("="*80)

    # Caminhos
    script_dir = Path(__file__).parent
    pdb_path = script_dir / 'model.pdb'
    output_dir = script_dir / 'analysis_results'

    # Verificar PDB
    if not pdb_path.exists():
        print(f"ERRO: {pdb_path} não encontrado")
        sys.exit(1)

    # Gerar mutações
    print(f"\nEtapa 1: Gerando mutações")
    print("-" * 80)
    mutations = gerar_mutacoes(str(pdb_path))
    print(f"Total de mutações: {len(mutations)}")

    # Contar por tipo de aminoácido
    aa_counts = {}
    for mut in mutations:
        aa = mut['original_aa']
        aa_counts[aa] = aa_counts.get(aa, 0) + 1

    print("\nDistribuição por aminoácido:")
    for aa, count in sorted(aa_counts.items()):
        print(f"  {aa}: {count}")

    # Configurar protocolo
    # Para teste rápido: nstruct=5
    # Para produção: nstruct=35
    config = FlexDdGConfig(
        nstruct=5,                    # Alterar para 35 em produção
        repack_radius=8.0,
        max_minimization_iter=200
    )

    print(f"\nEtapa 2: Executando Flex ddG")
    print("-" * 80)
    print("ATENÇÃO: Este processo pode demorar várias horas")
    print(f"Tempo estimado: ~{len(mutations) * 2} minutos com nstruct={config.nstruct}")
    print()

    # Executar protocolo
    protocol = FlexDdGPyRosetta(config)
    results = protocol.run_flex_ddg(
        str(pdb_path),
        mutations,
        str(output_dir)
    )

    # Análise de resultados
    print("\n" + "="*80)
    print("RESULTADOS")
    print("="*80)

    # Estatísticas
    ddg_mean = results['ddg'].mean()
    ddg_std = results['ddg'].std()
    ddg_max = results['ddg'].max()
    ddg_min = results['ddg'].min()

    print(f"\nEstatísticas:")
    print(f"  ΔΔG médio: {ddg_mean:+.2f} ± {ddg_std:.2f} kcal/mol")
    print(f"  ΔΔG máximo: {ddg_max:+.2f} kcal/mol")
    print(f"  ΔΔG mínimo: {ddg_min:+.2f} kcal/mol")

    # Identificar hotspots (ΔΔG > +2.0 kcal/mol)
    hotspots = results[results['ddg'] > 2.0].copy()
    hotspots = hotspots.sort_values('ddg', ascending=False)

    print(f"\nHotspots (ΔΔG > +2.0 kcal/mol): {len(hotspots)}")

    if len(hotspots) > 0:
        print("\nTop 10 hotspots:")
        print("-" * 80)
        for _, row in hotspots.head(10).iterrows():
            print(f"  {row['mutation']:8s} | ΔΔG = {row['ddg']:+.2f} ± {row['ddg_std']:.2f} kcal/mol")

        # Salvar hotspots
        hotspots_file = output_dir / 'hotspots.csv'
        hotspots.to_csv(hotspots_file, index=False)
        print(f"\nHotspots salvos em: {hotspots_file}")

    # Salvar resumo
    with open(output_dir / 'resumo.txt', 'w') as f:
        f.write("ANÁLISE DE ALANINE SCANNING - eIF4E\n")
        f.write("="*80 + "\n\n")
        f.write(f"Total de mutações: {len(mutations)}\n")
        f.write(f"ΔΔG médio: {ddg_mean:+.2f} ± {ddg_std:.2f} kcal/mol\n")
        f.write(f"ΔΔG máximo: {ddg_max:+.2f} kcal/mol\n")
        f.write(f"ΔΔG mínimo: {ddg_min:+.2f} kcal/mol\n")
        f.write(f"\nHotspots identificados: {len(hotspots)}\n\n")

        if len(hotspots) > 0:
            f.write("Top 10 hotspots:\n")
            f.write("-" * 80 + "\n")
            for _, row in hotspots.head(10).iterrows():
                f.write(f"{row['mutation']:8s} | ΔΔG = {row['ddg']:+.2f} ± {row['ddg_std']:.2f} kcal/mol\n")

    print(f"\nResumo salvo em: {output_dir}/resumo.txt")
    print("\nAnálise completa!")


if __name__ == '__main__':
    main()
