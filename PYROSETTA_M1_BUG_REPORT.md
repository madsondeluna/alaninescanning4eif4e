# PyRosetta M1 (ARM) Bug Report - Segmentation Fault

## Summary

PyRosetta-4 2025.47 for M1 (ARM) crashes with segmentation fault (exit code 139) when performing basic operations required for Flex ddG protocol.

## Environment

```
Platform: macOS M1 (ARM64)
PyRosetta: PyRosetta-4 2025 [Rosetta PyRosetta4.Release.python310.m1 2025.47+release.8bb54e8a6dc3e1c027e4f028bdace6bd4691c823 2025-11-20T15:59:05]
Python: 3.10.19
Installation: conda (miniconda base environment)
```

## Bug Description

### Operations Tested (All Failed)

1. **FastRelax** - Segfault
   - Operation: `FastRelax().apply(pose)`
   - Location: Structure relaxation

2. **MinMover** - Segfault
   - Operation: `MinMover().apply(pose)`
   - Location: Energy minimization

3. **PackRotamersMover** - Segfault
   - Operation: `PackRotamersMover().apply(pose)`
   - Location: Side-chain repacking

4. **MutateResidue** - Segfault (CRITICAL)
   - Operation: `MutateResidue(target=pos, new_res='A').apply(pose)`
   - Location: Point mutation
   - **This is the root cause** - cannot even perform basic mutation

### Failure Point

```python
# This crashes consistently:
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

pose = pyrosetta.pose_from_pdb("gamma.pdb")
mutate = MutateResidue(target=3, new_res='A')
mutate.apply(pose)  # <-- SEGFAULT HERE (exit code 139)
```

## Test Cases

### Structures Tested

| Structure | Size | Result |
|-----------|------|--------|
| model.pdb | 157 residues | Segfault |
| gamma.pdb | 12 residues | Segfault |

### Mutations Attempted

| Mutation | Position | Result |
|----------|----------|--------|
| L2A | 2 | Segfault |
| R3A | 3 | Segfault |

**All mutations fail at the same point:** MutateResidue().apply()

## Workarounds Attempted

### 1. Replace FastRelax with MinMover
**Result:** Segfault in MinMover

### 2. Disable Minimization
**Result:** Segfault in PackRotamersMover

### 3. Disable Minimization + Repacking
**Result:** Segfault in MutateResidue

### 4. Minimalist Test (mutate + score only)
**Result:** Segfault in MutateResidue

**Conclusion:** No workaround possible. Core mutation functionality is broken.

## Stack Trace Pattern

```
PyRosetta inicializado
Função de energia: REF2015
Carregando estrutura: gamma.pdb
Relaxando estrutura wild-type...
Energia wild-type: 83.37 REU
Processando 9 mutações...
[1/9] R3A
[Exit code: 139]
```

Crash occurs immediately when attempting first mutation.

## Impact

### Unusable Features

- Alanine scanning
- Flex ddG protocol
- Point mutation analysis
- Any workflow requiring MutateResidue
- Any workflow requiring minimization or repacking

### Usable Features

- Structure loading
- Energy calculation (scorefxn)
- Basic pose operations (clone, etc.)

## Recommendations

### For Users

**Option 1: Use Rosetta C++ (Recommended for Production)**

```bash
# Install Rosetta C++ from RosettaCommons
# https://www.rosettacommons.org/software/license-and-download

# Use flex_ddg native executable
$ROSETTA/main/source/bin/flex_ddg.macosclangrelease \
  -s input.pdb \
  -ddg:mut_file mutations.txt \
  -nstruct 50

# More stable, faster, and fully featured
```

**Option 2: Use Intel x86 PyRosetta via Rosetta**

```bash
# Install x86 PyRosetta (runs under Rosetta 2)
# May have performance penalty but more stable
arch -x86_64 conda create -n pyrosetta_x86 python=3.10
arch -x86_64 conda activate pyrosetta_x86
arch -x86_64 conda install -c conda-forge pyrosetta
```

**Option 3: Use Cloud Computing (x86)**

- AWS EC2 (x86 instances)
- Google Cloud Compute
- Run analysis on Intel-based machines

### For PyRosetta Developers

**Bug Location:**

The segfault originates in the ARM-compiled binary during pose manipulation operations. Likely causes:

1. Memory alignment issues (ARM vs x86)
2. SIMD instruction incompatibility
3. Pointer arithmetic bugs in ARM build
4. Missing ARM-specific optimizations

**Suggested Fix:**

- Recompile PyRosetta M1 version with:
  - Address sanitizer (`-fsanitize=address`)
  - Debug symbols
  - ARM-specific optimizations review
  - Memory alignment checks

**Testing Priority:**

1. MutateResidue - CRITICAL (breaks all mutation workflows)
2. PackRotamersMover - HIGH (breaks repacking workflows)
3. MinMover - HIGH (breaks minimization workflows)
4. FastRelax - MEDIUM (alternative: MinMover, but also broken)

## Reproducible Test Case

```python
#!/usr/bin/env python3
"""
Minimal reproducible test for PyRosetta M1 bug.
"""

import pyrosetta
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

# Initialize
pyrosetta.init("-ignore_unrecognized_res -mute all")

# Load any small structure
pose = pyrosetta.pose_from_pdb("gamma.pdb")  # 12 residues

print(f"Loaded: {pose.total_residue()} residues")
print("Attempting mutation R3A...")

# This will segfault on PyRosetta M1:
mutate = MutateResidue(target=3, new_res='A')
mutate.apply(pose)  # <-- CRASH

print("SUCCESS")  # Never reached
```

**Expected:** Mutation applied, "SUCCESS" printed
**Actual:** Segmentation fault (exit code 139)

## Temporary Solution Implemented

For this project, we've documented the issue and recommend users:

1. Use Rosetta C++ flex_ddg for production analyses
2. Use PyRosetta on Intel hardware (not M1)
3. Wait for PyRosetta team to fix M1 build

The Python framework remains useful for:
- Mutation list generation
- Result parsing and visualization
- Workflow automation (calling C++ Rosetta)

## References

- PyRosetta Documentation: http://www.pyrosetta.org/
- Flex ddG Paper: Barlow et al. (2018) J Phys Chem B. 122(21):5389-5399
- Issue Reported: [Add GitHub issue link when created]

## Date

2025-12-04

## Reporter

Madson Luna - madsondeluna@gmail.com
