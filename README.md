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

### How is MicroSynViz different?

Unlike genome-wide synteny tools (MCScanX, SynVisio, JupiterPlot), MicroSynViz focuses on **base-pair resolution visualization of individual gene pairs**, integrating gene structure (exons/introns), transposable element annotations, and BLASTN homology ribbons in a single publication-ready figure. This makes it particularly suited for studying the mechanisms of gene duplication (transposed, segmental, tandem) and the role of TEs in genomic rearrangements.

| Feature | MCScanX | SynVisio | MicroSynViz |
|---------|---------|----------|-------------|
| Scale | Genome-wide | Genome-wide | Gene-pair level |
| Gene structure (exon/intron) | No | No | **Yes** |
| TE track overlay | No | No | **Yes** |
| Auto BLAST + visualize | No | No | **Yes** |
| Publication-ready vector output | No | Screenshot | **SVG + PDF** |
| One-command workflow | No | No | **Yes** |

---

## Key Features

- **Per-region input**: each region gets its own genome + annotations — no mode switching
- **Unified GFF parsing**: gene and TE annotations are auto-detected from any GFF file
- **Two input modes**: gene IDs (`--gene1/--gene2`) or genomic regions (`--region1/--region2`)
- **Automatic reverse complement** detection based on BLAST alignment orientation
- **Flexible ribbon coloring**: color by bitscore, identity, or e-value (`--color_by`)
- **BLAST filtering**: minimum identity and alignment length thresholds
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

---

## Input Files

### GFF3 Annotation

Standard GFF3 format with `gene`, `mRNA`, `exon`, and/or TE features. MicroSynViz automatically detects gene structure and TE features from the same file — no need to separate them.

```
Chr1    MSU_osa1r7    gene    1000    2000    .    +    .    ID=LOC_Os01g01010
Chr1    MSU_osa1r7    exon    1000    1200    .    +    .    Parent=LOC_Os01g01010.1
Chr1    RepeatMasker  transposable_element  3000  3500  .  +  .  ID=TE001;Target="Motif:hAT"
```

### Genome FASTA

Must be indexed with `samtools faidx`:

```bash
samtools faidx genome.fasta
```

---

## Quick Start

### Single-species: Compare two genes

```bash
MicroSynViz \
    --gene1 LOC_Os06g50440 \
    --gene2 LOC_Os06g50789 \
    --g1 genome.fa --gffs1 annotation.gff TEs.gff \
    --g2 genome.fa --gffs2 annotation.gff TEs.gff \
    --auto_complementary --bezier \
    --output my_synteny
```

### Cross-species comparison

```bash
MicroSynViz \
    --gene1 GeneA \
    --gene2 GeneB \
    --g1 rice.fa --gffs1 rice.gff rice_te.gff \
    --g2 maize.fa --gffs2 maize.gff maize_te.gff \
    --auto_complementary --bezier \
    --output cross_species
```

### Region mode

```bash
MicroSynViz \
    --region1 Chr1:1000-5000 \
    --region2 Chr2:3000-8000 \
    --g1 genome.fa --gffs1 annotation.gff \
    --g2 genome.fa --gffs2 annotation.gff \
    --extend 5000 \
    --output region_synteny
```

> **Design principle**: each region is self-contained with `--g1/--gffs1` and `--g2/--gffs2`. For single-species, simply provide the same files for both. GFF files are auto-parsed for both gene structure and TE features.

---

## Command-line Options

Run `MicroSynViz --help` to see all options.

### Input

| Parameter | Type | Description |
|-----------|------|-------------|
| `--gene1` | str | Gene ID for region 1 |
| `--gene2` | str | Gene ID for region 2 |
| `--region1` | str | Genomic region 1 (Chr:start-end) |
| `--region2` | str | Genomic region 2 (Chr:start-end) |
| `--g1` | FASTA | Genome FASTA for region 1 |
| `--g2` | FASTA | Genome FASTA for region 2 |
| `--gffs1` | GFF+ | GFF file(s) for region 1 (gene, TE, etc.) |
| `--gffs2` | GFF+ | GFF file(s) for region 2 (gene, TE, etc.) |
| `--extend` | int | Flanking extension in bp (default: 3000) |

### BLAST

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--evalue` | float | 1e-5 | E-value threshold |
| `--identity` | float | 50 | Minimum identity % |
| `--alignment_length` | int | 5 | Minimum alignment length (bp) |
| `--color_by` | choice | bitscore | Ribbon color: `bitscore`, `identity`, or `evalue` |
| `--threads` | int | 8 | BLAST threads |
| `--blast_result` | str | — | Pre-computed BLAST (skip auto-BLAST) |
| `--auto_complementary` | flag | — | Auto-detect reverse complement |

### Appearance

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--bezier` | flag | — | Bezier curves for ribbons |
| `--svg_width` | int | 2000 | SVG canvas width |
| `--svg_height` | int | 800 | SVG canvas height |
| `--output` | str | MicroSynViz_result | Output file prefix |

---

## Output Files

| File | Description |
|------|-------------|
| `*_linkview.svg` | Vector SVG figure |
| `*_linkview.pdf` | Vector PDF (via CairoSVG) |
| `*_genes.fasta` | Extracted sequences |
| `*_blast.txt` | BLAST tabular results |

---

## Examples

Example data is provided in `example-test/`:

- **Single-species**: `example-test/Single-species_comparison/`
- **Cross-species**: `example-test/Cross-species_comparison/`

```bash
cd example-test/Single-species_comparison
bash run_example.sh
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
