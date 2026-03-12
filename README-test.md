###############################################################################################
  ██████╗ ███████╗███╗   ██╗███████╗██╗   ██╗██╗███████╗
 ██╔════╝ ██╔════╝████╗  ██║██╔════╝██║   ██║██║╚══███╔╝
 ██║  ███╗█████╗  ██╔██╗ ██║█████╗  ██║   ██║██║  ███╔╝
 ██║   ██║██╔══╝  ██║╚██╗██║██╔══╝  ╚██╗ ██╔╝██║ ███╔╝
 ╚██████╔╝███████╗██║ ╚████║███████╗  ╚████╔╝ ██║███████╗
  ╚═════╝ ╚══════╝╚═╝  ╚═══╝╚══════╝   ╚═══╝  ╚═╝╚══════╝
###############################################################################################

# GeneViz: A User-Friendly Toolkit for Visualizing Pairwise Genomic Microsynteny across Any Species

[![Python](https://img.shields.io/badge/Python-3.6+-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

> **GitHub**: https://github.com/yourusername/GeneViz  
> **License**: MIT | **Release Date**: 2026-03-12

---

## 📖 Table of Contents

- [What does GeneViz do?](#-what-does-geneviz-do)
- [Key Features](#-key-features)
- [Dependencies](#-dependencies)
- [Installation](#-installation)
  - [Option 1: Direct script usage (simplest)](#option-1-direct-script-usage-simplest)
  - [Option 2: Install via pip (from GitHub)](#option-2-install-via-pip-from-github)
  - [Option 3: Install via pip (from PyPI)](#option-3-install-via-pip-from-pypi)
  - [Option 4: Conda environment (recommended for reproducibility)](#option-4-conda-environment-recommended-for-reproducibility)
- [Input File Formats](#-input-file-formats)
- [Command Line Arguments](#-command-line-arguments)
- [Usage Examples](#-usage-examples)
  - [Single‑species comparison](#single‑species-comparison)
  - [Cross‑species comparison](#cross‑species-comparison)
  - [Including TE tracks](#including-te-tracks)
  - [Automatic reverse complement](#automatic-reverse-complement)
  - [Highlighting a core region](#highlighting-a-core-region)
  - [Adjusting plot appearance](#adjusting-plot-appearance)
- [Output Files](#-output-files)
- [Tips & Troubleshooting](#-tips--troubleshooting)
- [Packaging & Publishing (for developers)](#-packaging--publishing-for-developers)
- [Citation](#-citation)
- [License](#-license)

---

## 🔬 What does GeneViz do?

GeneViz is a lightweight yet powerful Python script that generates **publication‑quality pairwise microsynteny plots** for any two genomic regions, within the same species or across different species. It combines:

- **BLAST‑based homology** ribbons (colored by bitscore)
- **Gene structure** (exons/introns with direction arrows)
- **Transposable element (TE)** tracks with motif labels
- **Automatic orientation detection** (reverse‑complement when needed)
- **Core region highlighting** (e.g., conserved domains)

All you need are GFF annotation files, a genome FASTA, and BLAST output (optional – the script can run BLAST for you). GeneViz produces SVG and PDF vector figures ready for publication.

---

## ✨ Key Features

- **Single‑species or cross‑species** – compare genes on the same genome or between two different genomes.
- **Automatic reverse complement** – detects if the second region should be flipped based on BLAST hit directions.
- **TE track** – draws transposable elements above/below the chromosome, extracts motif names from GFF attributes.
- **Gene structure** – shows exons (colored boxes) and introns (lines), with an arrow indicating transcription direction.
- **Core range highlighting** – mark a specific interval (e.g., protein domain) in red within genes.
- **BLAST ribbon coloring** – polygons connecting homologous regions, colored by bitscore (gray → purple → red).
- **Vector output** – SVG and PDF (via `cairosvg`) for easy editing.
- **Highly customizable** – dozens of command‑line options to control sizes, colors, gaps, etc.
- **No complex installation** – just a single Python script; all Python dependencies are listed in `requirements.txt`.

---

## 📦 Dependencies

### Python packages
- `biopython`
- `pandas`
- `cairosvg` (for PDF conversion)

### External command‑line tools (must be installed separately and in `PATH`)
- `samtools` (≥1.9)
- `blastn` (NCBI BLAST+ ≥2.10)

Install Python dependencies with:
```bash
pip install biopython pandas cairosvg
```
Install external tools via conda (recommended):
```bash
conda install -c bioconda samtools blast
```
or via your system package manager (e.g., apt-get install samtools ncbi-blast+ on Ubuntu).
## 🛠 Installation
### Option 1: Direct script usage (simplest)

1.Download the genevis.py script from the GitHub repository.

2.Install Python dependencies: pip install biopython pandas cairosvg.

3.Ensure samtools and blastn are installed and accessible.

4.Run the script with:
```bash
python Geneviz.py --help
```

### Option 2: Install via pip (from GitHub)
Clone the repository and install in editable mode (this makes the genevis command available):
```bash
git clone https://github.com/yourusername/GeneViz.git
cd GeneViz
pip install -e .
```
Now you can run genevis from anywhere:
```bash
Geneviz --help
```
### Option 3: Install via pip (from PyPI)
If the package is uploaded to PyPI:
```bash
pip install GeneViz
genevis --help
```
## 📁 Input File Formats
GeneViz requires several input files. All are plain text; please follow the exact formats.
### 1. Gene pair & Region specification
#### - Gene pair specification:
The two gene IDs are supplied directly via `--gene1` and `--gene2`. No separate pair file is needed. Use:
```bash
--gene1 GeneID-1
--gene2 GeneID-2
```
#### - Region specification：
To mark a specific sub-region on each chromosome (e.g., a conserved domain), use:
```bash
--region1 Chr:start-end #for the first region
--region2 Chr:start-end #for the second region
```
Genes overlapping these intervals will be drawn in red (target color).

### 2. GFF3 annotation file(s)
- For single‑species mode: `--gff`
- For cross‑species mode: `--gff1` and `--gff2`
The GFF must contain `gene` features with an `ID` attribute (or `Name`). Exon/intron structures are extracted from `mRNA` and `exon` features. Example:
```text
Chr1    .       gene    1000    2000    .       +       .       ID=GeneA
Chr1    .       mRNA    1000    2000    .       +       .       ID=GeneA.1;Parent=GeneA
Chr1    .       exon    1000    1200    .       +       .       Parent=GeneA.1
Chr1    .       exon    1500    2000    .       +       .       Parent=GeneA.1
```
### 3. Genome FASTA file(s)
- For single‑species mode: `--fasta`
- For cross‑species mode: `--fasta1` and `--fasta2`
The FASTA must be indexed with `samtools faidx` (run `samtools faidx genome.fasta` beforehand). GeneViz will extract the region around each gene using `samtools`.
### 4. TE GFF file(s) (optional)
- For single‑species mode: `--te_gff`
- For cross‑species mode: `--te_gff1` and `--te_gff2`

If the attribute contains `Target "Motif:xxx"`, the motif name is displayed; otherwise "unknown". We recommend using the RepeatMasker annotation output file (.out.gff format) as input. Example:
```text
##gff-version 3
##sequence-region Chr1 1 43270923
Chr1	RepeatMasker	dispersed_repeat	2762	2924	 9.2	-	.	ID=1;Target "Motif:Harbinger-N11_OS" 1 171
Chr1	RepeatMasker	dispersed_repeat	18790	18926	 2.9	+	.	ID=2;Target "Motif:LINE1-N1F_OS" 1480 1616
Chr1	RepeatMasker	dispersed_repeat	20731	20772	 7.3	-	.	ID=3;Target "Motif:MuDR-N18H_OS" 2627 2669
Chr1	RepeatMasker	dispersed_repeat	22455	22538	 9.9	-	.	ID=4;Target "Motif:hAT-N43B_OS" 2135 2228
Chr1	RepeatMasker	dispersed_repeat	27884	27974	10.0	+	.	ID=5;Target "Motif:EnSpm-11_OS" 4064 4156
Chr1	RepeatMasker	dispersed_repeat	27895	28006	 7.2	+	.	ID=6;Target "Motif:HELITRON7_OS" 1599 1713
Chr1	RepeatMasker	dispersed_repeat	28732	28974	 8.7	-	.	ID=7;Target "Motif:Harbinger-N9_OS" 1 231
Chr1	RepeatMasker	dispersed_repeat	30292	30684	 1.0	-	.	ID=8;Target "Motif:ECR" 1 395
```
### 5. BLAST output (optional)
If you already have BLAST results, provide them with --blast_result. Otherwise, GeneViz will run BLAST for you using the extracted sequences. The format must be outfmt 6 (tabular) with columns:
```text
