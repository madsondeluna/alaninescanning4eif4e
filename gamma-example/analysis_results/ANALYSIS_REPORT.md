# Alanine Scanning Analysis Report - Gamma Peptide

**Date:** 2025-12-02 19:27:44
**PDB File:** gamma.pdb
**Analysis Type:** Complete peptide alanine scanning

---

## 1. Structure Information

- **Chain:** A
- **Total Residues:** 12
- **Sequence:** GLY → GLY → ARG → CYS → ARG → GLY → PHE → ARG → ARG → ARG → CYS → PHE

### Peptide Composition

- **ARG**: 5 residue(s)
- **CYS**: 2 residue(s)
- **GLY**: 3 residue(s)
- **PHE**: 2 residue(s)


---

## 2. Alanine Scanning Parameters

- **Target:** Complete peptide (all scannable residues)
- **Chain Analyzed:** A
- **Total Mutations Generated:** 9
- **Excluded Residues:** GLY (already glycine)

### Mutation Breakdown

- **C**: 2 mutation(s)
- **F**: 2 mutation(s)
- **R**: 5 mutation(s)


### Complete Mutation List

| # | Mutation | Position | Original AA | Target AA |
|---|----------|----------|-------------|-----------|
| 1 | A3R | 3 | R | A |
| 2 | A4C | 4 | C | A |
| 3 | A5R | 5 | R | A |
| 4 | A7F | 7 | F | A |
| 5 | A8R | 8 | R | A |
| 6 | A9R | 9 | R | A |
| 7 | A10R | 10 | R | A |
| 8 | A11C | 11 | C | A |
| 9 | A12F | 12 | F | A |


---

## 3. ddG Results (Simulated)

**NOTE:** These are SIMULATED values for demonstration. For real results, run with Rosetta Flex ddG.

### Statistics

- **Number of Mutations:** 9
- **Mean ΔΔG:** 0.60 kcal/mol
- **Std ΔΔG:** 0.32 kcal/mol
- **Min ΔΔG:** 0.18 kcal/mol
- **Max ΔΔG:** 1.28 kcal/mol

### All Results

| Rank | Mutation | Position | Original AA | ΔΔG (kcal/mol) |
|------|----------|----------|-------------|----------------|
| 1 | A7F | 7 | F | 1.28 |
| 2 | A10R | 10 | R | 0.72 |
| 3 | A4C | 4 | C | 0.71 |
| 4 | A3R | 3 | R | 0.68 |
| 5 | A9R | 9 | R | 0.62 |
| 6 | A5R | 5 | R | 0.50 |
| 7 | A8R | 8 | R | 0.35 |
| 8 | A12F | 12 | F | 0.33 |
| 9 | A11C | 11 | C | 0.18 |


---

## 4. Hotspot Analysis

**Hotspot Threshold:** ΔΔG > 1.5 kcal/mol

**Number of Hotspots:** 0

No significant hotspots identified.


---

## 5. Biological Insights

### Peptide Characteristics

The gamma peptide is a **12-residue peptide** with the following notable features:

1. **High Arginine Content (5 residues):** Suggests strong positive charge, potentially important for:
   - Electrostatic interactions
   - DNA/RNA binding
   - Membrane interactions

2. **Cysteine Residues (2 residues):** May form disulfide bonds:
   - Positions: 
   - Potentially critical for structural stability

3. **Aromatic Residues (2 PHE):** Contribute to:
   - Hydrophobic core
   - π-π stacking interactions
   - Binding interfaces

### Functional Predictions

Based on the alanine scanning results, the most critical residues for peptide function appear to be:



---

## 6. Recommendations

### For Experimental Validation

1. **Top Priority Mutations:**


2. **Mutagenesis Strategy:**
   - Use site-directed mutagenesis
   - Test single mutations first
   - Measure binding affinity (if applicable)
   - Assess structural stability (CD, thermal shift)

3. **Controls:**
   - Wild-type peptide
   - Non-hotspot mutations (negative controls)

### For Computational Follow-up

1. **Molecular Dynamics Simulations:**
   - Simulate wild-type and top mutants
   - Assess conformational changes
   - Calculate binding free energies

2. **Rosetta Flex ddG:**
   - Run protocol with actual Rosetta software
   - Use nstruct=35 for statistical rigor
   - Compare with these simulated results

---

## 7. Output Files

All analysis files are saved in `gamma-example/analysis_results/`:

### Mutation Files
- `mutations.txt` - Human-readable mutation list
- `mutations_rosetta.txt` - Rosetta input format
- `mutations.csv` - Spreadsheet format

### Results Files
- `ddg_results.csv` - Complete ddG values
- `hotspots.csv` - Filtered hotspot residues
- `ANALYSIS_REPORT.md` - This report

### Visualizations (`plots/` subdirectory)
- `ddg_distribution.png` - ΔΔG histogram
- `top_hotspots.png` - Bar chart of top mutations
- `position_scan_chain_A.png` - ΔΔG along sequence
- `hotspot_heatmap.png` - 2D heatmap

---

## 8. Notes

- This analysis used **SIMULATED** ddG values for demonstration
- For publication-quality results, use Rosetta Flex ddG protocol
- Recommended Rosetta parameters:
  - `nstruct`: 35 (minimum for statistical significance)
  - `iterations`: 3 (backrub iterations)
  - `score function`: REF2015 (recommended)

---

**Analysis completed:** 2025-12-02 19:27:44

