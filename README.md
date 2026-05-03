# AlleleSelect

**Allele-Selective ASO Design Pipeline for CACNA1A Gain-of-Function Mutations**

The Xiu Lab | [thexiulab.org](https://www.thexiulab.org/current-projects/alleleselect) | [github.com/axshoe/alleleselect](https://github.com/axshoe/alleleselect)

---

AlleleSelect is a Python CLI tool for computational design of allele-selective antisense oligonucleotides (ASOs) targeting gain-of-function (GoF) CACNA1A variants, with a pre-run demo for R192Q (c.575G>A, FHM1).

For heterozygous GoF mutations like R192Q, the therapeutic constraint is non-negotiable: you cannot silence both copies of CACNA1A. The wildtype allele is keeping cerebellar and cortical synapses functional. AlleleSelect fills the gap left by general ASO design tools by computing allele selectivity ratios (ASR = ΔG_mutant − ΔG_wildtype) for every candidate window, ranking candidates by preferential mutant binding, and integrating accessibility, off-target, and splice-risk scoring into a single reproducible pipeline.

## What AlleleSelect does that general ASO tools do not

- Computes delta-delta-G between mutant and wildtype binding for every window using SantaLucia 1998 nearest-neighbor thermodynamics (implemented from scratch)
- Applies Peyret 1999 internal mismatch penalties for wildtype binding
- Scores mRNA accessibility via RNAfold partition function
- BLASTn off-target check against GENCODE v44 human transcriptome
- Flags candidates within 15 nt of CACNA1A splice sites
- Annotates recommended gapmer modification patterns (MOE/LNA, PS backbone)
- Pre-computed demo output for R192Q in `demo/R192Q_output/`

## Installation

```bash
# Clone
git clone https://github.com/axshoe/alleleselect.git
cd alleleselect

# Install Python dependencies
pip install -r requirements.txt

# Install package
pip install -e .
```

**Optional (recommended for full scoring):**
- ViennaRNA (RNAfold): https://www.tbi.univie.ac.at/RNA/
- BLAST+: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- GENCODE v44 transcriptome FASTA (for off-target): see SETUP.md

## Quick Start

```bash
# R192Q (pre-run demo output already in repo)
alleleselect --variant c.575G>A --output demo/R192Q_output/ --no-blast --no-rnafold

# Full run (requires ViennaRNA + BLAST + GENCODE)
alleleselect --variant c.575G>A --transcript ENST00000360228.10 --output R192Q_full/
```

## Run tests

```bash
python -m pytest tests/ -v
```

## Pre-computed R192Q output

`demo/R192Q_output/candidates.csv` and `demo/R192Q_output/report.html` contain pre-computed AlleleSelect results for CACNA1A c.575G>A (R192Q). These are ready for sharing with wet-lab collaborators.

## Output columns

| Column | Description |
|--------|-------------|
| priority_rank | Overall rank (1 = best) |
| ASO_sequence | 5'→3' DNA sequence |
| allele_selectivity_ratio_kcal_mol | ΔG_mutant − ΔG_wildtype; more negative = more selective |
| mRNA_accessibility_score | Mean unpaired probability from RNAfold (0–1) |
| off_target_count | BLASTn hits ≥80% identity over ≥14 nt (−1 = not checked) |
| splice_risk | Y if within 15 nt of annotated splice site |
| recommended_gapmer_pattern | MOE or LNA gapmer design |

## Contact / collaboration

Questions: angie.xiu27@gmail.com | The Xiu Lab: [thexiulab.org](https://thexiulab.org)

## Limitations

AlleleSelect is a research-grade computational design tool. All candidates require experimental synthesis and validation before any in vivo testing. Thermodynamic predictions use simplified nearest-neighbor models that do not account for all aspects of cellular ASO behavior (protein binding, delivery, off-target effects not captured by BLASTn, etc.).

## License

MIT License. See LICENSE.

---