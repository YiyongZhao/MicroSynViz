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
| Publication-ready vector output | No | Screenshot | **PDF / SVG** |
| One-command workflow | No | No | **Yes** |

---

## Key Features

- **Per-region input**: each region gets its own genome + annotations — no mode switching
- **Format-agnostic**: auto-detects GFF3, GTF, and BED — mix freely in annotation inputs
- **Two input modes**: gene IDs (`--gene1/--gene2`) or regions (`--region1/--region2` in `SeqID:start-end` format)
- **Flexible FASTA input**: genome FASTA (chromosome-level) or CDS/transcript FASTA (gene-level) — annotations optional for CDS mode
- **Automatic reverse complement** detection based on BLAST alignment orientation
- **Flexible ribbon coloring**: color by bitscore, identity, or e-value (`--color_by`)
- **BLAST filtering**: minimum identity and alignment length thresholds
- **PDF or SVG output** (vector graphics via CairoSVG)
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

### Annotation Files (GFF3 / GTF / BED)

MicroSynViz **auto-detects** the file format. You can mix formats freely in `--annos1` / `--annos2`.

**GFF3** — standard format with `gene`, `mRNA`, `exon`, TE features:

```
Chr1    MSU_osa1r7    gene    1000    2000    .    +    .    ID=LOC_Os01g01010
Chr1    MSU_osa1r7    exon    1000    1200    .    +    .    Parent=LOC_Os01g01010.1
Chr1    RepeatMasker  transposable_element  3000  3500  .  +  .  ID=TE001;Target="Motif:hAT"
```

**GTF** — Ensembl/GENCODE style:

```
Chr1    ensembl    gene    1000    2000    .    +    .    gene_id "ENSG001"; gene_name "ABC";
Chr1    ensembl    transcript    1000    2000    .    +    .    gene_id "ENSG001"; transcript_id "ENST001";
Chr1    ensembl    exon    1000    1200    .    +    .    gene_id "ENSG001"; transcript_id "ENST001";
```

**BED** — BED3/BED6/BED12 (0-based coordinates, auto-converted):

```
Chr1    999    2000    LOC_Os01g01010    0    +    999    2000    0    3    200,150,300    0,400,800
```

BED12 block columns are parsed as exons; BED3/BED6 are treated as single-exon genes.

### Genome FASTA

Must be indexed with `samtools faidx`:

```bash
samtools faidx genome.fasta
```

### Pre-computed BLAST Result (optional)

By default, MicroSynViz runs BLASTN automatically. If you want to provide your own BLAST result (e.g., with custom parameters), use `--blast_result`:

The file must be in **BLAST tabular format** (`-outfmt 6`), with 12 tab-separated columns:

```
qseqid  sseqid  pident  length  mismatch  gapopen  qstart  qend  sstart  send  evalue  bitscore
```

Generate it with:

```bash
blastn -query region1.fa -subject region2.fa \
    -outfmt 6 -evalue 1e-5 -num_threads 8 \
    -out my_blast.txt
```

Then pass it to MicroSynViz:

```bash
MicroSynViz --gene1 A --gene2 B \
    --fa1 genome.fa --annos1 anno.gff \
    --fa2 genome.fa --annos2 anno.gff \
    --blast_result my_blast.txt
```

> **Note**: When `--blast_result` is provided, `--evalue` and `--threads` are ignored (BLAST is skipped). Post-filters `--identity` and `--alignment_length` still apply.

---

## Quick Start

> **CDS/Transcript FASTA mode**: When using CDS or transcript FASTA files (where each sequence header is a gene ID), annotation files (`--annos1`/`--annos2`) are **optional**. Gene IDs are matched directly against FASTA headers. Use `--extend 0` to avoid extracting beyond sequence boundaries.

### Single-species: Compare two genes

```bash
MicroSynViz \
    --gene1 LOC_Os06g50440 \
    --gene2 LOC_Os06g50789 \
    --fa1 genome.fa --annos1 annotation.gff TEs.gff \
    --fa2 genome.fa --annos2 annotation.gff TEs.gff \
    --auto_complementary --bezier \
    --output my_synteny
```

### Cross-species comparison

```bash
MicroSynViz \
    --gene1 GeneA \
    --gene2 GeneB \
    --fa1 rice.fa --annos1 rice.gff rice_te.gff \
    --fa2 maize.fa --annos2 maize.gff maize_te.gff \
    --auto_complementary --bezier \
    --output cross_species
```

### Region mode (genome FASTA)

```bash
MicroSynViz \
    --region1 Chr1:1000-5000 \
    --region2 Chr2:3000-8000 \
    --fa1 genome.fa --annos1 annotation.gff \
    --fa2 genome.fa --annos2 annotation.gff \
    --extend 5000 \
    --output region_synteny
```

### CDS/Transcript mode (no annotations needed)

```bash
MicroSynViz \
    --gene1 GENE_A --gene2 GENE_B \
    --fa1 cds.fa --fa2 cds.fa \
    --extend 0 --bezier \
    --output cds_comparison
```

### Region mode (CDS/transcript FASTA)

```bash
MicroSynViz \
    --region1 LOC_Os06g50440:1-1000 \
    --region2 LOC_Os06g50789:1-1500 \
    --fa1 cds.fa --fa2 cds.fa \
    --extend 0 \
    --output cds_comparison
```

> **Design principle**: each region is self-contained with `--fa1/--annos1` and `--fa2/--annos2`. For single-species, simply provide the same files for both. GFF files are auto-parsed for both gene structure and TE features.

---

## Command-line Options

Run `MicroSynViz --help` to see all options.

### Input

| Parameter | Type | Description |
|-----------|------|-------------|
| `--gene1` | str | Gene ID for region 1. Searched in `--annos1` files; falls back to FASTA header if not found |
| `--gene2` | str | Gene ID for region 2. Searched in `--annos2` files; falls back to FASTA header if not found |
| `--region1` | str | Region 1 in `SeqID:start-end` format. Genome: `Chr1:1000-5000`; CDS: `GeneID:1-1000` |
| `--region2` | str | Region 2 in `SeqID:start-end` format. Genome: `Chr2:3000-8000`; CDS: `GeneID:1-2000` |
| `--fa1` | FASTA | FASTA file for region 1 — genome (chromosome-level) or CDS/transcript (gene-level). Auto-indexed if `.fai` missing |
| `--fa2` | FASTA | FASTA file for region 2 — genome (chromosome-level) or CDS/transcript (gene-level). Auto-indexed if `.fai` missing |
| `--annos1` | GFF3/GTF/BED | Annotation file(s) for region 1. Multiple files accepted (e.g., `genes.gff TEs.gff`); rendered from center outward by order. **Optional** for CDS FASTA |
| `--annos2` | GFF3/GTF/BED | Annotation file(s) for region 2. Multiple files accepted; rendered from center outward by order. **Optional** for CDS FASTA |
| `--extend` | int | Flanking extension in bp (default: 3000). Use `0` for CDS FASTA |

### BLAST

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--evalue` | float | 1e-5 | E-value threshold |
| `--identity` | float | 50 | Minimum identity % |
| `--alignment_length` | int | 5 | Minimum alignment length (bp) |
| `--color_by` | choice | bitscore | Ribbon color: `bitscore`, `identity`, or `evalue` |
| `--threads` | int | 8 | BLAST threads |
| `--blast_result` | FILE | — | Pre-computed BLAST tabular result (`-outfmt 6`); skips auto-BLAST |
| `--auto_complementary` | flag | — | Auto-detect reverse complement |

### Appearance

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--bezier` | flag | — | Use Bezier curves for smooth homology ribbons |
| `--ribbon_opacity` | float | 0.2 | Opacity of BLAST ribbons (0=invisible, 1=opaque) |
| `--chro_axis` | flag | — | Draw tick marks on chromosome bars |
| `--no_scale` | flag | — | Omit the scale bar |
| `--output` | str | MicroSynViz_result | Output file prefix |
| `--quiet` / `-q` | flag | — | Suppress informational output (warnings/errors only) |

---

## Output Files

| File | Description |
|------|-------------|
| `*_linkview.pdf` | Publication-ready vector PDF figure |
| `*_genes.fasta` | Extracted sequences |
| `*_blast.txt` | BLAST tabular results (outfmt 6) |

---

## Examples

Example data is provided in `example-test/`:

- **Single-species**: `example-test/Single-species_comparison/`
```bash
cd example-test/Single-species_comparison
bash run_example.sh
```
- **Cross-species**: `example-test/Cross-species_comparison/`
```bash
cd example-test/Cross-species_comparison
bash run_example.sh
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
