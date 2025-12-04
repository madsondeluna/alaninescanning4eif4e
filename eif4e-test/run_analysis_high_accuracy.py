#!/usr/bin/env python3
"""
Análise de Alanine Scanning para eIF4E - ALTA ACURÁCIA

Este script usa configuração otimizada baseada em:
- Barlow KA, et al. (2018) Flex ddG. J Phys Chem B. 122(21):5389-5399
- Kellogg EH, et al. (2011) Conformational sampling. Proteins. 79(3):830-8
- Smith CA, Kortemme T (2008) Backrub simulation. Structure. 16(7):1126-33

Configuração para alta acurácia:
- nstruct=50 (ensemble robusto)
- backrub_iterations=5 (flexibilidade completa)
- repack_radius=10.0 Å (proteína pequena)
- max_minimization_iter=500 (convergência completa)

Tempo estimado: ~18-26 horas para 132 mutações
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

    Exclui:
    - GLY (sem Cβ)
    - ALA (já é alanina)
    - PRO (restrição conformacional)

    Args:
        pdb_path: Caminho para arquivo PDB

    Returns:
        Lista de dicts com informações de mutações
    """
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
    """Pipeline de análise com configuração de alta acurácia."""
    print("="*80)
    print("ALANINE SCANNING - eIF4E (ALTA ACURÁCIA)")
    print("Protocolo Flex ddG com PyRosetta")
    print("="*80)
    print()
    print("CONFIGURAÇÃO:")
    print("  • nstruct = 50 (ensemble robusto)")
    print("  • backrub iterations = 5 (flexibilidade completa)")
    print("  • repack radius = 10.0 Å (proteína pequena)")
    print("  • minimization iter = 500 (convergência total)")
    print("  • score function = REF2015 (mais acurada)")
    print()
    print("REFERÊNCIAS:")
    print("  [1] Barlow KA, et al. (2018) J Phys Chem B. 122(21):5389-5399")
    print("  [2] Kellogg EH, et al. (2011) Proteins. 79(3):830-8")
    print("  [3] Smith CA, Kortemme T (2008) Structure. 16(7):1126-33")
    print("="*80)

    # Caminhos
    script_dir = Path(__file__).parent
    pdb_path = script_dir / 'model.pdb'
    output_dir = script_dir / 'analysis_results_high_accuracy'
    output_dir.mkdir(exist_ok=True)

    # Verificar PDB
    if not pdb_path.exists():
        print(f"\nERRO: {pdb_path} não encontrado")
        sys.exit(1)

    # Gerar mutações
    print(f"\nEtapa 1: Gerando mutações")
    print("-" * 80)
    mutations = gerar_mutacoes(str(pdb_path))
    print(f"Total de mutações: {len(mutations)}")

    # Estatísticas de aminoácidos
    aa_counts = {}
    for mut in mutations:
        aa = mut['original_aa']
        aa_counts[aa] = aa_counts.get(aa, 0) + 1

    print("\nDistribuição por aminoácido:")
    for aa in sorted(aa_counts.keys()):
        count = aa_counts[aa]
        bar = '█' * (count // 2)
        print(f"  {aa}: {count:3d} {bar}")

    # Configuração de alta acurácia
    config = FlexDdGConfig(
        nstruct=50,                    # Ensemble robusto
        repack_radius=10.0,            # Proteína pequena
        max_minimization_iter=500      # Convergência completa
    )

    # Estimativa de tempo
    tempo_por_mut = 10  # minutos (estimativa conservadora)
    tempo_total_min = len(mutations) * tempo_por_mut
    tempo_total_h = tempo_total_min / 60

    print(f"\nEtapa 2: Executando Flex ddG")
    print("-" * 80)
    print(f"⚠️  ATENÇÃO: Este processo levará VÁRIAS HORAS")
    print(f"⏱️  Tempo estimado: ~{tempo_total_h:.1f} horas ({tempo_total_min} min)")
    print(f"💾  Resultados em: {output_dir}")
    print()
    print("Dicas:")
    print("  • Execute overnight ou durante fim de semana")
    print("  • Monitore progresso com: tail -f analysis_results_high_accuracy/resumo.txt")
    print("  • Processo pode ser pausado e retomado")
    print()

    # Confirmar execução
    resposta = input("Continuar com análise? [s/N]: ")
    if resposta.lower() not in ['s', 'sim', 'y', 'yes']:
        print("\nAnálise cancelada pelo usuário")
        sys.exit(0)

    # Executar protocolo
    print("\n🚀 Iniciando análise...")
    print("="*80)

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

    # Estatísticas gerais
    ddg_mean = results['ddg'].mean()
    ddg_std_overall = results['ddg'].std()
    ddg_max = results['ddg'].max()
    ddg_min = results['ddg'].min()

    # Média do desvio padrão (convergência)
    ddg_std_mean = results['ddg_std'].mean()

    print(f"\nEstatísticas Gerais:")
    print(f"  ΔΔG médio:  {ddg_mean:+.2f} ± {ddg_std_overall:.2f} kcal/mol")
    print(f"  ΔΔG máximo: {ddg_max:+.2f} kcal/mol")
    print(f"  ΔΔG mínimo: {ddg_min:+.2f} kcal/mol")
    print(f"\nConvergência:")
    print(f"  Desvio padrão médio: {ddg_std_mean:.2f} kcal/mol")
    if ddg_std_mean < 0.5:
        print(f"  ✓ Excelente convergência")
    elif ddg_std_mean < 1.0:
        print(f"  ✓ Boa convergência")
    else:
        print(f"  ⚠ Convergência moderada (considere aumentar nstruct)")

    # Identificar hotspots (ΔΔG > +2.0 kcal/mol)
    hotspots = results[results['ddg'] > 2.0].copy()
    hotspots = hotspots.sort_values('ddg', ascending=False)

    # Classificação por magnitude
    critical = results[results['ddg'] > 2.5]
    high = results[(results['ddg'] > 2.0) & (results['ddg'] <= 2.5)]
    moderate = results[(results['ddg'] > 1.5) & (results['ddg'] <= 2.0)]

    print(f"\nClassificação de Hotspots:")
    print(f"  🔴 CRÍTICO  (ΔΔG > +2.5): {len(critical):3d} ({len(critical)/len(results)*100:.1f}%)")
    print(f"  🟠 ALTO     (ΔΔG > +2.0): {len(high):3d} ({len(high)/len(results)*100:.1f}%)")
    print(f"  🟡 MODERADO (ΔΔG > +1.5): {len(moderate):3d} ({len(moderate)/len(results)*100:.1f}%)")
    print(f"  ⚪ FRACO    (ΔΔG < +1.5): {len(results) - len(critical) - len(high) - len(moderate):3d}")

    print(f"\nTotal de hotspots (ΔΔG > +2.0): {len(hotspots)}")

    if len(hotspots) > 0:
        print(f"\n{'='*80}")
        print("TOP 15 HOTSPOTS")
        print("="*80)
        print(f"{'Rank':<6} {'Mutação':<10} {'ΔΔG':<12} {'Desvio':<12} {'Classe':<10}")
        print("-"*80)

        for rank, (_, row) in enumerate(hotspots.head(15).iterrows(), 1):
            ddg = row['ddg']
            std = row['ddg_std']

            if ddg > 2.5:
                classe = "🔴 CRÍTICO"
            else:
                classe = "🟠 ALTO"

            print(f"{rank:<6} {row['mutation']:<10} {ddg:+6.2f} kcal/mol  ± {std:4.2f}      {classe}")

        # Salvar hotspots
        hotspots_file = output_dir / 'hotspots.csv'
        hotspots.to_csv(hotspots_file, index=False)
        print(f"\n✓ Hotspots salvos em: {hotspots_file}")

    # Análise por região (para eIF4E)
    print(f"\n{'='*80}")
    print("ANÁLISE POR REGIÃO (eIF4E)")
    print("="*80)

    # Definir regiões funcionais da eIF4E
    regions = {
        'Cap-binding (22-42)': (22, 42),
        'Core (43-70)': (43, 70),
        'Loop flexível (71-85)': (71, 85),
        'β-sheet (86-120)': (86, 120),
        'Interface eIF4G (145-157)': (145, 157)
    }

    for region_name, (start, end) in regions.items():
        region_muts = results[(results['position'] >= start) & (results['position'] <= end)]
        if len(region_muts) > 0:
            region_mean = region_muts['ddg'].mean()
            region_hotspots = len(region_muts[region_muts['ddg'] > 2.0])
            print(f"\n{region_name}:")
            print(f"  Mutações: {len(region_muts)}")
            print(f"  ΔΔG médio: {region_mean:+.2f} kcal/mol")
            print(f"  Hotspots: {region_hotspots} ({region_hotspots/len(region_muts)*100:.1f}%)")

            # Top 3 da região
            top3 = region_muts.nlargest(3, 'ddg')
            if len(top3) > 0:
                print(f"  Top 3: ", end="")
                print(", ".join([f"{row['mutation']} ({row['ddg']:+.1f})"
                                for _, row in top3.iterrows()]))

    # Salvar resultados completos
    results_file = output_dir / 'ddg_results.csv'
    results.to_csv(results_file, index=False)

    # Salvar resumo detalhado
    with open(output_dir / 'resumo.txt', 'w') as f:
        f.write("="*80 + "\n")
        f.write("ANÁLISE DE ALANINE SCANNING - eIF4E (ALTA ACURÁCIA)\n")
        f.write("="*80 + "\n\n")

        f.write("CONFIGURAÇÃO:\n")
        f.write(f"  nstruct = {config.nstruct}\n")
        f.write(f"  repack_radius = {config.repack_radius} Å\n")
        f.write(f"  max_minimization_iter = {config.max_minimization_iter}\n\n")

        f.write("ESTATÍSTICAS:\n")
        f.write(f"  Total de mutações: {len(mutations)}\n")
        f.write(f"  ΔΔG médio: {ddg_mean:+.2f} ± {ddg_std_overall:.2f} kcal/mol\n")
        f.write(f"  ΔΔG máximo: {ddg_max:+.2f} kcal/mol\n")
        f.write(f"  ΔΔG mínimo: {ddg_min:+.2f} kcal/mol\n")
        f.write(f"  Convergência (std médio): {ddg_std_mean:.2f} kcal/mol\n\n")

        f.write("CLASSIFICAÇÃO:\n")
        f.write(f"  CRÍTICO  (ΔΔG > +2.5): {len(critical)}\n")
        f.write(f"  ALTO     (ΔΔG > +2.0): {len(high)}\n")
        f.write(f"  MODERADO (ΔΔG > +1.5): {len(moderate)}\n\n")

        f.write("TOP 15 HOTSPOTS:\n")
        f.write("-" * 80 + "\n")
        for rank, (_, row) in enumerate(hotspots.head(15).iterrows(), 1):
            f.write(f"{rank:2d}. {row['mutation']:8s} | ΔΔG = {row['ddg']:+.2f} ± {row['ddg_std']:.2f} kcal/mol\n")

    print(f"\n{'='*80}")
    print(f"✓ Análise completa!")
    print(f"✓ Resultados em: {output_dir}")
    print(f"  • ddg_results.csv - Todos os resultados")
    print(f"  • hotspots.csv - Apenas hotspots")
    print(f"  • resumo.txt - Resumo executivo")
    print("="*80)


if __name__ == '__main__':
    main()
