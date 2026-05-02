# AlleleSelect — SETUP.md

**Allele-Selective ASO Design Pipeline for CACNA1A Gain-of-Function Mutations**
The Xiu Lab | thexiulab.org | github.com/axshoe/alleleselect

---

## What You Need Before Starting

| Requirement | Required? | Notes |
|---|---|---|
| Python 3.9+ | Yes | Check: `python --version` |
| Git | Yes | For cloning and pushing to GitHub |
| pip | Yes | Comes with Python |
| ViennaRNA (RNAfold) | Recommended | For mRNA accessibility scoring |
| BLAST+ | Recommended | For off-target screening |
| GENCODE v44 FASTA | Recommended | ~4 GB, needed for BLAST |
| Internet connection | Yes (first run) | Fetches CDS from Ensembl API |

If ViennaRNA or BLAST+ are not installed, the pipeline runs with `--no-rnafold` and `--no-blast` flags. All thermodynamic scoring still works; accessibility scores default to 0.5 (neutral) and off-target counts to -1 (unknown).

---

## Step 1: Install Python

If you do not have Python 3.9+:
- Windows: Download from https://python.org/downloads — check "Add Python to PATH" during install
- Verify: open PowerShell and run `python --version`

---

## Step 2: Clone the Repository

Open PowerShell:

```powershell
cd C:\Users\<YourName>\Documents    # or wherever you keep projects
git clone https://github.com/axshoe/alleleselect.git
cd alleleselect
```

---

## Step 3: Install Python Dependencies

```powershell
pip install -r requirements.txt
```

This installs: requests, biopython, numpy, scipy, pytest.

If you get a permissions error, try:
```powershell
pip install -r requirements.txt --user
```

---

## Step 4: Install the Package (Editable Mode)

```powershell
pip install -e .
```

This makes the `alleleselect` command available in your terminal. Verify:
```powershell
alleleselect --help
```

---

## Step 5: Install ViennaRNA (Recommended)

ViennaRNA provides RNAfold, used for mRNA secondary structure and accessibility scoring.

1. Download the Windows installer from: https://www.tbi.univie.ac.at/RNA/#download
2. Run the installer, accept defaults.
3. After installation, open a new PowerShell window and verify:
   ```powershell
   RNAfold --version
   ```
4. If `RNAfold` is not found, add the ViennaRNA bin directory to your PATH:
   - Search Windows: "Environment Variables"
   - Edit PATH → Add the ViennaRNA install directory (usually `C:\Program Files\ViennaRNA\bin`)

---

## Step 6: Install BLAST+ (Recommended)

BLAST+ is needed for transcriptome-wide off-target screening.

1. Download from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
   - Choose the Windows installer: `ncbi-blast-...-win64.exe`
2. Run the installer, accept defaults.
3. Verify in a new PowerShell window:
   ```powershell
   blastn -version
   makeblastdb -version
   ```
4. If not found, add BLAST+ to PATH (usually `C:\Program Files\NCBI\blast-...\bin`).

---

## Step 7: Download GENCODE v44 Transcriptome (For BLAST Off-Target)

This is a large download (~4 GB compressed). Only needed if using BLAST off-target screening.

1. Download: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
2. Create a data directory:
   ```powershell
   mkdir $HOME\alleleselect_data
   ```
3. Move the downloaded file to `C:\Users\<YourName>\alleleselect_data\gencode.v44.transcripts.fa.gz`
4. Decompress (use 7-Zip or PowerShell):
   ```powershell
   # Using 7-Zip (if installed):
   & "C:\Program Files\7-Zip\7z.exe" e "$HOME\alleleselect_data\gencode.v44.transcripts.fa.gz" -o"$HOME\alleleselect_data\"
   ```
5. Set the environment variable so AlleleSelect finds the file:
   ```powershell
   $env:ALLELESELECT_GENCODE_FASTA = "$HOME\alleleselect_data\gencode.v44.transcripts.fa"
   $env:ALLELESELECT_BLAST_DB = "$HOME\alleleselect_data\gencode_v44_db"
   ```
   To make these permanent, add them via System → Environment Variables in Windows Settings.

6. Build the BLAST database (one time, takes 5-10 min):
   ```powershell
   alleleselect --variant c.575G>A --output test_db_build/ --no-rnafold
   ```
   The first run with BLAST enabled will automatically build the database.

---

## Step 8: Run the Tests

```powershell
python -m pytest tests/ -v
```

All 27 tests should pass. If any fail, check that the package was installed correctly (`pip install -e .`) and that your Python version is 3.9+.

---

## Step 9: Run the Pipeline

### Minimal run (no external tools required):
```powershell
alleleselect --variant c.575G>A --output my_R192Q_run/ --no-blast --no-rnafold
```

### Full run (requires ViennaRNA + BLAST + GENCODE):
```powershell
alleleselect --variant c.575G>A --transcript ENST00000360228.10 --output my_R192Q_full/ --verbose
```

### View pre-computed demo output:
Open `demo/R192Q_output/report.html` in any browser. No additional setup required.

### Other variants:
```powershell
# Example for a different CACNA1A GoF variant
alleleselect --variant c.1748C>T --output my_variant_run/ --no-blast --no-rnafold
```

### All options:
```powershell
alleleselect --help
```

---

## Step 10: Open in PyCharm

1. Open PyCharm → File → Open → navigate to the `alleleselect/` folder → Open
2. PyCharm may prompt to create a virtual environment. Accept and point it to your Python 3.9+ interpreter.
3. In the terminal at the bottom of PyCharm, run:
   ```
   pip install -e .
   pip install -r requirements.txt
   python -m pytest tests/ -v
   ```
4. The project uses standard Python package structure with `alleleselect/` as the main package. PyCharm will recognize this automatically.

---

## Step 11: Push to GitHub

```powershell
cd alleleselect

# Initialize (if not already a git repo from clone):
git init
git remote add origin https://github.com/axshoe/alleleselect.git

# First push:
git add .
git commit -m "Initial AlleleSelect release — R192Q pre-computed output included"
git branch -M main
git push -u origin main

# Subsequent pushes:
git add .
git commit -m "Your message here"
git push
```

If you get an authentication error on push, use a GitHub Personal Access Token (Settings → Developer settings → Personal access tokens) instead of your password.

---

## Output Files Reference

After running the pipeline, the output directory contains:

| File | Description |
|---|---|
| `candidates.csv` | All ranked candidates with thermodynamic scores |
| `report.html` | Interactive HTML report with scatter plot and sortable table |

The CSV columns are:
- `priority_rank`: Overall rank (1 = best)
- `ASO_sequence`: The candidate sequence (5'→3', DNA)
- `allele_selectivity_ratio_kcal_mol`: Core metric; more negative = more selective for mutant
- `mRNA_accessibility_score`: From RNAfold; >0.65 preferred
- `off_target_count`: BLASTn hits; 0 is ideal (-1 = not checked)
- `splice_risk`: Y if within 15 nt of a CACNA1A splice site
- `recommended_gapmer_pattern`: e.g. "5MOE-10DNA-5MOE (all-PS backbone)"

---

## Environment Variables

| Variable | Purpose | Default |
|---|---|---|
| `ALLELESELECT_GENCODE_FASTA` | Path to GENCODE transcriptome FASTA | `~/alleleselect_data/gencode.v44.transcripts.fa` |
| `ALLELESELECT_BLAST_DB` | Path prefix for BLAST database | `~/alleleselect_data/gencode_v44_db` |
| `ALLELESELECT_GENCODE_GTF` | Path to GENCODE GTF (optional, for splice sites) | `~/alleleselect_data/gencode.v44.annotation.gtf.gz` |

---

## Troubleshooting

**`alleleselect: command not found`**
Run `pip install -e .` again from inside the `alleleselect/` directory.

**`EnsemblFetchError: Network error`**
Check your internet connection. Ensembl REST API is required for CDS fetching.

**`RNAfold not found`**
Install ViennaRNA and ensure the bin directory is on your PATH. Use `--no-rnafold` to skip.

**`makeblastdb not found`**
Install BLAST+ and ensure the bin directory is on your PATH. Use `--no-blast` to skip.

**`GENCODE FASTA not found`**
Download and decompress the file. Set `ALLELESELECT_GENCODE_FASTA` env variable to its path.

**Tests fail**
Run `pip install -e .` first. Check Python is 3.9+ with `python --version`.

---

*AlleleSelect v1.0.0 | The Xiu Lab | thexiulab.org | angie.xiu27@gmail.com*
