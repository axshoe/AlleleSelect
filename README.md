# AlleleSelect

**Allele-Selective ASO Design Pipeline for Dominant Neurological Mutations**

The Xiu Lab | [thexiulab.org](https://www.thexiulab.org/current-projects/alleleselect) | [github.com/axshoe/alleleselect](https://github.com/axshoe/alleleselect)

---

AlleleSelect is a Python CLI tool for computational design of allele-selective gapmer antisense oligonucleotides (ASOs) targeting dominant gain-of-function mutations. It was built for CACNA1A FHM1 variants and extended to any coding SNP in any gene.

For heterozygous gain-of-function mutations, the therapeutic constraint is non-negotiable: you cannot silence both copies of the gene. AlleleSelect fills the gap left by general ASO design tools by computing allele selectivity ratios (ASR = ΔG_mutant − ΔG_wildtype) for every candidate window and integrating accessibility, off-target, splice-risk, and RNase H cleavage site scoring into a single reproducible pipeline.

Pre-computed demo output is available for six CACNA1A FHM1 variants (R192Q, S218L, G293R, R1349Q, A454T, T501M), ATXN1/SCA1, and COL6A3/UCMD (junction mode).

## What AlleleSelect does

- Allele selectivity ratio using Sugimoto 1995 RNA:DNA hybrid nearest-neighbor thermodynamics with Peyret 1999 mismatch corrections (v7 default; SantaLucia 1998 available via `--rna-params santalucia`)
- SNP position scoring: triangular weighting from gap center (Ostergaard 2013) or sequence-dependent RNase H cleavage site prediction (Kielpinski 2017, `--rnase-h-scoring`)
- Differential RNAfold accessibility (mutant vs. wildtype separately)
- BLASTn off-target screening against GENCODE v44 human transcriptome
- Splice site proximity flagging
- Gapmer modification pattern annotation (MOE/LNA, PS backbone)
- Engineered mismatch mode: deliberate extra mismatch against wildtype (`--extra-mismatch`)
- Junction mode for exon-skipping mutations (`--junction-mode`)

## Installation

```bash
git clone https://github.com/axshoe/alleleselect.git
cd alleleselect
pip install -e .
```

**After any `git pull`, reinstall to pick up updates:**
```bash
git pull
pip install -e .
```

**Optional dependencies (recommended for full scoring):**
- ViennaRNA (RNAfold): https://www.tbi.univie.ac.at/RNA/
- BLAST+: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
- GENCODE v44 transcriptome FASTA: see SETUP.md

## Quick start

```bash
# R192Q — no external dependencies required
alleleselect --variant c.575G>A --transcript ENST00000360228.10 --gene CACNA1A \
  --output demo/R192Q/ --no-blast

# Full run with BLAST, engineered mismatch candidates, and v8 RNase H scoring
alleleselect --variant c.575G>A --transcript ENST00000360228.10 --gene CACNA1A \
  --output demo/R192Q_v8/ --extra-mismatch --rnase-h-scoring

# Fixed-length 20-mer 5-10-5 architecture (e.g. for wet-lab screen comparison)
alleleselect --variant c.XXXN>Y --transcript ENST00000436367.6 --gene ATXN1 \
  --no-splice-check --fixed-length 20 --gapmer-architecture 5-10-5 \
  --output atxn1_20mer/

# Junction mode (exon-skipping mutations)
alleleselect --junction-mode --transcript ENST00000359218.9 --gene COL6A3 \
  --skip-exon 16 --output col6a3_junction/
```

## Run tests

```bash
python -m pytest tests/ -v
```

## Key flags

| Flag | Description |
|------|-------------|
| `--variant` | HGVS coding variant (e.g. c.575G>A) |
| `--transcript` | Ensembl transcript ID |
| `--gene` | Gene name for output labeling |
| `--output` | Output directory |
| `--no-blast` | Skip BLASTn (faster, for testing) |
| `--extra-mismatch` | Generate engineered mismatch candidates |
| `--rnase-h-scoring` | v8: sequence-dependent RNase H cleavage site scoring |
| `--fixed-length INT` | Generate only candidates of specified length |
| `--gapmer-architecture W-G-W` | Filter to specific wing-gap-wing architecture |
| `--rna-params` | `sugimoto` (default) or `santalucia` |
| `--top-n-blast INT` | Number of candidates to BLASTn (default 50) |
| `--verbose` | Detailed progress output |

## Output columns

| Column | Description |
|--------|-------------|
| priority_rank | Overall rank (1 = best) |
| ASO_sequence | 5' to 3' DNA sequence |
| allele_selectivity_ratio_kcal_mol | ΔG_mutant − ΔG_wildtype; more negative = more selective |
| mRNA_accessibility_score | Mean unpaired probability from RNAfold (0 to 1) |
| off_target_count | BLASTn hits at >=80% identity over >=14 nt (-1 = not screened) |
| min_off_target_mismatches | Minimum edit distance to any off-target |
| splice_risk | Y if within 15 nt of annotated splice site |
| recommended_gapmer_pattern | MOE or LNA gapmer design string |
| composite_score | Weighted rank score (ASR 40%, position 35%, accessibility 25%) |

## Validation

External experimental validation on ATXN1/SCA1 (Scholten, LUMC, 2026): Pearson r = 0.91 between AlleleSelect composite score and experimental allele discrimination across 20 gapmers (p < 0.0001, n = 20). Top 4 of 5 AlleleSelect predictions matched the high-discrimination cluster experimentally.

## Contact

Questions: angie.xiu27@gmail.com | The Xiu Lab: [thexiulab.org](https://thexiulab.org)

## Limitations

AlleleSelect is a research-grade computational design tool. All candidates require experimental synthesis and validation. Thermodynamic predictions use nearest-neighbor models and do not account for all aspects of cellular ASO behavior.

## License

MIT License. See LICENSE.