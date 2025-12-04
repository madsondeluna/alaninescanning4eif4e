#!/usr/bin/env python3
"""
Teste simples para debug do crash.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

import pyrosetta
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap

# Inicializar PyRosetta
print("1. Inicializando PyRosetta...")
pyrosetta.init("-ignore_unrecognized_res -mute all")
print("   ✓ PyRosetta inicializado")

# Carregar estrutura
print("\n2. Carregando estrutura...")
pose = pyrosetta.pose_from_pdb("model.pdb")
print(f"   ✓ Estrutura carregada: {pose.total_residue()} resíduos")

# Scorefxn
print("\n3. Configurando score function...")
scorefxn = get_score_function()
energy = scorefxn(pose)
print(f"   ✓ Energia inicial: {energy:.2f} REU")

# Tentar mutar primeiro resíduo LEU2
print("\n4. Tentando mutar L2A...")
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

mutate = MutateResidue(target=2, new_res='A')
print("   - Aplicando mutação...")
mutate.apply(pose)
print("   ✓ Mutação aplicada")

energy_mut = scorefxn(pose)
print(f"   ✓ Energia após mutação: {energy_mut:.2f} REU")

# Tentar minimizar
print("\n5. Tentando minimizar...")
movemap = MoveMap()
movemap.set_bb(False)  # DESLIGAR backbone primeiro
movemap.set_chi(True)  # Apenas side chains

min_mover = MinMover()
min_mover.movemap(movemap)
min_mover.score_function(scorefxn)
min_mover.min_type('lbfgs_armijo_nonmonotone')
min_mover.tolerance(0.01)
min_mover.max_iter(50)  # Poucas iterações

print("   - Aplicando MinMover (50 iter, apenas chi)...")
min_mover.apply(pose)
print("   ✓ Minimização completa")

energy_final = scorefxn(pose)
print(f"   ✓ Energia final: {energy_final:.2f} REU")

print("\n✅ TESTE PASSOU! Nenhum crash detectado.")
print(f"\nΔΔG = {energy_final - energy:.2f} kcal/mol")
