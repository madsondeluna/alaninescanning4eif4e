# Critical Bug: Segmentation Fault in MutateResidue on PyRosetta M1 (ARM) Build

## Summary

PyRosetta-4 2025.47 for M1 (ARM) crashes with segmentation fault (exit code 139) when performing core mutation operations required for any ddG-based protocol. The bug prevents all mutation-based workflows (Flex ddG, alanine scanning, saturation mutagenesis) from functioning on Apple Silicon hardware.

## Environment

```
Platform:        macOS (Apple Silicon M1/ARM64)
OS Version:      Darwin 25.1.0
PyRosetta:       PyRosetta-4 2025 [Rosetta PyRosetta4.Release.python310.m1 2025.47+release.8bb54e8a6dc3e1c027e4f028bdace6bd4691c823 2025-11-20T15:59:05]
Python:          3.10.19
Installation:    conda install -c https://conda.rosettacommons.org pyrosetta
Conda env:       miniconda base environment
```

## Bug Description

### Critical Failure Point

The segmentation fault occurs consistently in `MutateResidue().apply()`, which is the fundamental operation for all mutation-based protocols. The crash happens regardless of:

- Protein size (tested 12-157 residues)
- Mutation type (tested multiple residue types → Ala)
- Protocol complexity (fails even in minimal test case)
- Configuration parameters (fails with all settings)

### Operations Affected

All tested operations crash with exit code 139:

| Operation | Status | Location |
|-----------|--------|----------|
| `MutateResidue().apply()` | CRASH | Point mutation (CRITICAL) |
| `PackRotamersMover().apply()` | CRASH | Side-chain repacking |
| `MinMover().apply()` | CRASH | Energy minimization |
| `FastRelax().apply()` | CRASH | Structure relaxation |
| `pose_from_pdb()` | Works | Structure loading |
| `scorefxn(pose)` | Works | Energy calculation |

## Reproducible Test Case

### Minimal Example (12 residues)

```python
#!/usr/bin/env python3
"""
Minimal reproducible test for PyRosetta M1 bug.
Expected: Mutation applied successfully
Actual:   Segmentation fault (exit code 139)
"""

import pyrosetta
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

# Initialize PyRosetta
pyrosetta.init("-ignore_unrecognized_res -mute all")

# Load small test structure (12 residues)
pose = pyrosetta.pose_from_pdb("gamma.pdb")

print(f"Loaded: {pose.total_residue()} residues")
print("Attempting mutation R3A...")

# THIS CRASHES ON M1:
mutate = MutateResidue(target=3, new_res='A')
mutate.apply(pose)  # <-- SEGMENTATION FAULT

print("SUCCESS")  # Never reached
```

### Test Structures

| Structure | Size | PDB Source | Result |
|-----------|------|------------|--------|
| gamma.pdb | 12 residues | Simple test peptide | Segfault |
| model.pdb | 157 residues | eIF4E (AlphaFold) | Segfault |

Both structures are valid PDB files that load successfully with `pose_from_pdb()`.

### Stack Trace

```
PyRosetta inicializado
Função de energia: REF2015
Carregando estrutura: gamma.pdb
Estrutura carregada: 12 residues
Energia wild-type: 83.37 REU
Processando 9 mutações...
[1/9] R3A
[Exit code: 139]
```

No additional error messages or warnings are produced before the crash.

## Workarounds Attempted (All Failed)

### 1. Disable FastRelax → Use MinMover
**Result:** Segfault in MinMover

### 2. Disable All Minimization
**Result:** Segfault in PackRotamersMover

### 3. Disable Minimization + Repacking
**Result:** Segfault in MutateResidue

### 4. Minimalist Test (mutation + score only)
**Result:** Segfault in MutateResidue

### Conclusion
No workaround possible. The core `MutateResidue` operation is broken in the ARM build.

## Expected Behavior

The mutation should be applied to the pose without crashing:

```python
mutate = MutateResidue(target=3, new_res='A')
mutate.apply(pose)
print(f"Mutation successful: {pose.residue(3).name()}")  # Should print "ALA"
```

## Actual Behavior

```
[Exit code: 139]
Segmentation fault
```

The program crashes immediately during `mutate.apply(pose)` with no error message.

## Impact Assessment

### Unusable Features on M1

- ✗ Alanine scanning
- ✗ Flex ddG protocol
- ✗ Cartesian ddG protocol
- ✗ Point mutation analysis
- ✗ Saturation mutagenesis
- ✗ Protein design workflows
- ✗ Any protocol requiring `MutateResidue`
- ✗ Any protocol requiring minimization or repacking

### Usable Features on M1

- ✓ Structure loading (`pose_from_pdb`)
- ✓ Energy calculation (`scorefxn`)
- ✓ Basic pose operations (`clone`, `dump_pdb`)

## Root Cause Analysis

The segmentation fault likely originates from:

1. **Memory alignment issues** - ARM requires stricter memory alignment than x86
2. **SIMD instruction incompatibility** - x86 SSE/AVX vs ARM NEON
3. **Pointer arithmetic bugs** - Different pointer sizes or alignment assumptions
4. **Missing ARM-specific optimizations** - Code may assume x86 architecture

The bug appears to be in the compiled C++ binary, not in the Python interface layer, as:
- Python code executes correctly up to the C++ call
- No Python exceptions are raised
- The crash occurs during C++ pose manipulation

## Diagnostic Suggestions

For PyRosetta developers to debug this issue:

1. **Recompile with debugging tools:**
   ```bash
   # Add to compilation flags
   -fsanitize=address      # Detect memory errors
   -fsanitize=undefined    # Detect undefined behavior
   -g                      # Debug symbols
   -O0                     # Disable optimizations
   ```

2. **Check memory alignment:**
   ```cpp
   // Verify proper alignment for ARM
   static_assert(alignof(Pose) % 16 == 0);
   static_assert(alignof(Residue) % 16 == 0);
   ```

3. **Test SIMD code paths:**
   - Review all SIMD intrinsics for ARM compatibility
   - Consider using platform-agnostic SIMD libraries (e.g., xsimd, highway)

4. **Validate pointer operations:**
   - Check all pointer arithmetic in pose manipulation code
   - Verify no assumptions about pointer size or alignment

## Current Workarounds for Users

Since the M1 build is fundamentally broken for mutation workflows, users have three options:

### Option 1: Use Rosetta C++ (Recommended)

```bash
# Download from RosettaCommons
# https://www.rosettacommons.org/software/license-and-download

# Use native flex_ddg executable
$ROSETTA/main/source/bin/flex_ddg.macosclangrelease \
  -s input.pdb \
  -ddg:mut_file mutations.txt \
  -nstruct 50
```

**Advantages:** Stable, fast, full-featured, native ARM support

### Option 2: Use Intel x86 PyRosetta via Rosetta 2

```bash
# Install x86 PyRosetta (runs under Rosetta 2 emulation)
arch -x86_64 conda create -n pyrosetta_x86 python=3.10
arch -x86_64 conda activate pyrosetta_x86
arch -x86_64 conda install -c https://conda.rosettacommons.org pyrosetta
```

**Advantages:** Python interface maintained, more stable than ARM build
**Disadvantages:** Performance penalty from emulation

### Option 3: Use x86 Cloud Computing

- AWS EC2 (x86 instances)
- Google Cloud Compute
- Any Intel/AMD-based server

**Advantages:** Native x86 performance, scalable
**Disadvantages:** Requires cloud account, data transfer

## Testing Priority

For PyRosetta team testing:

1. **CRITICAL** - `MutateResidue` - Breaks all mutation workflows
2. **HIGH** - `PackRotamersMover` - Breaks repacking workflows
3. **HIGH** - `MinMover` - Breaks minimization workflows
4. **MEDIUM** - `FastRelax` - Alternative is MinMover (also broken)

## Additional Information

### PDB Test Files

I can provide the test PDB files used:
- `gamma.pdb` (12 residues) - minimal test case
- `model.pdb` (157 residues) - eIF4E structure

### Full Protocol Code

The complete Flex ddG protocol implementation is available if needed for testing. The code works correctly on Intel hardware and with Rosetta C++.

### Reproducibility

This bug is **100% reproducible** on:
- macOS M1 (ARM64)
- PyRosetta-4 2025.47
- Multiple protein structures
- Multiple mutation types

## System Information

```bash
# Architecture
uname -m
# arm64

# macOS Version
sw_vers
# ProductName:        macOS
# ProductVersion:     15.1
# BuildVersion:       24B83

# Python
python --version
# Python 3.10.19

# PyRosetta
python -c "import pyrosetta; print(pyrosetta.version())"
# PyRosetta-4 2025 [Rosetta PyRosetta4.Release.python310.m1 2025.47+release.8bb54e8a6dc3e1c027e4f028bdace6bd4691c823 2025-11-20T15:59:05]
```

## References

This bug was discovered while implementing protocols from:

1. Barlow KA, et al. (2018) **Flex ddG: Rosetta Ensemble-Based Estimation of Changes in Protein–Protein Binding Affinity upon Mutation.** *J Phys Chem B.* 122(21):5389-5399. DOI: 10.1021/acs.jpcb.8b01639

2. Kellogg EH, et al. (2011) **Role of conformational sampling in computing mutation-induced changes in protein structure and stability.** *Proteins.* 79(3):830-8. DOI: 10.1002/prot.22921

3. Smith CA, Kortemme T (2008) **Backrub-like backbone simulation recapitulates natural protein conformational variability and improves mutant side-chain prediction.** *Structure.* 16(7):1126-33. DOI: 10.1016/j.str.2008.05.006

## Contact

**Reporter:** Madson Luna
**Email:** madsondeluna@gmail.com
**Date:** 2024-12-04
**Issue Type:** Bug - Segmentation Fault
**Severity:** Critical (Core functionality broken on M1)

---

## Request for PyRosetta Team

Could you please:

1. Confirm if this is a known issue with the M1 build
2. Provide timeline for a fix if possible
3. Suggest any additional debugging steps I can perform
4. Confirm if the x86 build under Rosetta 2 is the recommended workaround until this is fixed

Thank you for your work on PyRosetta! Despite this M1 issue, the framework is excellent and works perfectly on Intel hardware.
