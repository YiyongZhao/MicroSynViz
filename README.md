```
###############################################################################################
 ███╗   ███╗██╗ ██████╗██████╗  ██████╗ ███████╗██╗   ██╗███╗   ██╗██╗   ██╗██╗███████╗
 ████╗ ████║██║██╔════╝██╔══██╗██╔═══██╗██╔════╝╚██╗ ██╔╝████╗  ██║██║   ██║██║╚══███╔╝
 ██╔████╔██║██║██║     ██████╔╝██║   ██║███████╗ ╚████╔╝ ██╔██╗ ██║██║   ██║██║  ███╔╝
 ██║╚██╔╝██║██║██║     ██╔══██╗██║   ██║╚════██║  ╚██╔╝  ██║╚██╗██║╚██╗ ██╔╝██║ ███╔╝
 ██║ ╚═╝ ██║██║╚██████╗██║  ██║╚██████╔╝███████║   ██║   ██║ ╚████║ ╚████╔╝ ██║███████╗
 ╚═╝     ╚═╝╚═╝ ╚═════╝╚═╝  ╚═╝ ╚═════╝╚══════╝   ╚═╝   ╚═╝  ╚═══╝  ╚═══╝  ╚═╝╚══════╝
###############################################################################################
```

# MicroSynViz: Visualizing Pairwise Genomic Microsynteny Within and Between Species

[![Python](https://img.shields.io/badge/Python-3.8+-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **GitHub**: https://github.com/YiyongZhao/MicroSynViz
> **License**: MIT

---

## Table of Contents

- [What does MicroSynViz do?](#what-does-microsynviz-do)
- [Key Features](#key-features)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Input Files](#input-files)
- [Quick Start](#quick-start)
- [Command-line Options](#command-line-options)
- [Output Files](#output-files)
- [Examples](#examples)
- [License](#license)

---

## What does MicroSynViz do?

MicroSynViz provides a reproducible workflow for visualizing **pairwise genomic microsynteny** between any two genomic regions within or across species. It integrates BLASTN-based sequence homology with genome annotation to generate publication-quality linkview plots, enabling intuitive classification of gene duplication types (transposed, segmental, tandem) and mechanistic inference of transposable element involvement.

All you need are GFF annotation files, a genome FASTA, and optionally a TE annotation. MicroSynViz produces SVG and PDF vector figures ready for publication.

---

## Key Features

- **Single-species and cross-species** pairwise microsynteny comparison
- **Two input modes**: gene IDs (`--gene1/--gene2`) or genomic regions (`--region1/--region2`)
- **Automatic reverse complement** detection based on BLAST alignment orientation
- **TE track** with motif labels for transposable element visualization
- **Gene structure** (exons/introns) visualization from GFF3 annotation
- **Core range highlighting** for focused synteny regions
- **SVG and PDF output** (vector graphics via CairoSVG)
- **Bezier curve** rendering for smooth homology ribbons

---

## Dependencies

### Python packages (installed automatically via pip)

- Python >= 3.8
- pandas
- BioPython
- CairoSVG

### External tools (must be in PATH)

- **samtools** >= 1.10 — for FASTA region extraction
- **blastn** (NCBI BLAST+) — for sequence homology search

Install external tools via conda:

```bash
conda install -c bioconda samtools blast
```

---

## Installation

### Option 1: Install via pip (from GitHub)

```bash
pip install git+https://github.com/YiyongZhao/MicroSynViz.git
MicroSynViz --help
```

### Option 2: Install via conda environment

```bash
git clone https://github.com/YiyongZhao/MicroSynViz.git
cd MicroSynViz
conda env create -f environment.yml
conda activate MicroSynViz
pip install -e .
MicroSynViz --help
```

### Option 3: Direct script usage

```bash
git clone https://github.com/YiyongZhao/MicroSynViz.git
cd MicroSynViz
pip install pandas biopython cairosvg
python -m microsynviz --help
```

---

## Input Files

MicroSynViz requires several input files. All are plain text; please follow the exact formats.

### 1. GFF3 Annotation

Standard GFF3 format with `gene`, `mRNA`, and `exon` features. The `ID=` attribute in column 9 is used to identify genes.

```
Chr1    MSU_osa1r7    gene    1000    2000    .    +    .    ID=LOC_Os01g01010;Name=...
Chr1    MSU_osa1r7    mRNA    1000    2000    .    +    .    ID=LOC_Os01g01010.1;Parent=LOC_Os01g01010
Chr1    MSU_osa1r7    exon    1000    1200    .    +    .    Parent=LOC_Os01g01010.1
```

### 2. Genome FASTA

Must be indexed with `samtools faidx`:

```bash
samtools faidx genome.fasta
```

### 3. TE Annotation GFF (optional)

GFF3 with transposable element annotations. The `Name=` attribute is displayed as the TE label.

### 4. BLAST Result (optional)

If not provided, MicroSynViz runs BLASTN automatically. To provide pre-computed results, use BLAST tabular format (`-outfmt 6`):

```bash
blastn -query seq1.fa -subject seq2.fa -outfmt 6 -evalue 1e-5 > blast.txt
```

---

## Quick Start

### Single-species: Compare two genes

```bash
MicroSynViz \
    --gene1 LOC_Os06g50440 \
    --gene2 LOC_Os06g50789 \
    --gff annotation.gff \
    --fasta genome.fa \
    --te_gff TEs.gff \
    --auto_complementary \
    --bezier \
    --output my_synteny
```

### Single-species: Compare two regions

```bash
MicroSynViz \
    --region1 Chr1:1000-5000 \
    --region2 Chr2:3000-8000 \
    --gff annotation.gff \
    --fasta genome.fa \
    --te_gff TEs.gff \
    --extend 5000 \
    --output region_synteny
```

### Cross-species comparison

```bash
MicroSynViz \
    --gene1 GeneA \
    --gene2 GeneB \
    --gff1 speciesA.gff \
    --fasta1 speciesA.fa \
    --gff2 speciesB.gff \
    --fasta2 speciesB.fa \
    --te_gff1 speciesA_te.gff \
    --te_gff2 speciesB_te.gff \
    --auto_complementary \
    --bezier \
    --output cross_species_synteny
```

---

## Command-line Options

Run `MicroSynViz --help` to see all options.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--gene1` | str | — | Gene ID for region 1 |
| `--gene2` | str | — | Gene ID for region 2 |
| `--region1` | str | — | Genomic region (Chr:start-end) for region 1 |
| `--region2` | str | — | Genomic region (Chr:start-end) for region 2 |
| `--gff` | str | — | GFF3 annotation (single species) |
| `--gff1` / `--gff2` | str | — | GFF3 annotations (cross-species) |
| `--fasta` | str | — | Genome FASTA (single species) |
| `--fasta1` / `--fasta2` | str | — | Genome FASTAs (cross-species) |
| `--te_gff` | str | — | TE annotation GFF (single species) |
| `--te_gff1` / `--te_gff2` | str | — | TE annotation GFFs (cross-species) |
| `--blast_result` | str | — | Pre-computed BLAST result (skip auto-BLAST) |
| `--extend` | int | 50000 | Flanking extension in bp |
| `--evalue` | float | 1e-5 | BLAST e-value threshold |
| `--threads` | int | 4 | BLAST threads |
| `--min_identity` | float | 0 | Minimum BLAST identity % |
| `--min_length` | int | 0 | Minimum alignment length |
| `--auto_complementary` | flag | — | Auto-detect reverse complement |
| `--bezier` | flag | — | Use Bezier curves for ribbons |
| `--output` | str | MicroSynViz_result | Output file prefix |

---

## Output Files

All output files are named with the prefix given by `--output` (default: `MicroSynViz_result`).

| File | Description |
|------|-------------|
| `*_linkview.svg` | Vector SVG figure |
| `*_linkview.pdf` | Vector PDF figure (converted from SVG) |
| `*_genes.fasta` | Extracted gene sequences |
| `*_blast.txt` | BLAST tabular results |

---

## Examples

Example data and expected outputs are provided in the `example-test/` directory:

- **Single-species**: `example-test/Single-species_comparison/`
- **Cross-species**: `example-test/Cross-species_comparison/`

To run an example:

```bash
cd example-test/Single-species_comparison
bash run_example.sh
```

---


## License

MIT License. See [LICENSE](LICENSE) for details.
