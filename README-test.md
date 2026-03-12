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
### 1. Gene pair specification
The two gene IDs are supplied directly via --gene1 and --gene2. No separate pair file is needed.
### 2. GFF3 annotation file(s)
- For single‑species mode: `--gff`
- For cross‑species mode: `--gff1` and `--gff2`
