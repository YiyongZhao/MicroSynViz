```
###############################################################################################
  ██████╗ ███████╗███╗   ██╗███████╗██╗   ██╗██╗███████╗
 ██╔════╝ ██╔════╝████╗  ██║██╔════╝██║   ██║██║╚══███╔╝
 ██║  ███╗█████╗  ██╔██╗ ██║█████╗  ██║   ██║██║  ███╔╝
 ██║   ██║██╔══╝  ██║╚██╗██║██╔══╝  ╚██╗ ██╔╝██║ ███╔╝
 ╚██████╔╝███████╗██║ ╚████║███████╗  ╚████╔╝ ██║███████╗
  ╚═════╝ ╚══════╝╚═╝  ╚═══╝╚══════╝   ╚═══╝  ╚═╝╚══════╝
###############################################################################################
```

**GeneViz: A User-Friendly Toolkit for Visualizing Pairwise Genomic Microsynteny across Any Species**

[![PyPI](https://img.shields.io/pypi/v/GeneViz)](https://pypi.org/project/GeneViz/)
[![Python](https://img.shields.io/badge/Python-3.8--3.12-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> Pypi: https://pypi.org/project/GeneViz
> Github: https://github.com/YiyongZhao/GeneViz
> License: MIT license | Release Date: 2026-3

---

## What does GeneViz do?

GeneViz provides a reproducible workflow for visualizing and interpreting **pairwise genomic microsynteny** between any two genomic regions across any species. It integrates BLASTN-based sequence homology with genome annotation to generate publication-quality linkview plots, enabling intuitive classification of gene duplication types (transposed, segmental, tandem) and mechanistic inference of transposable element involvement.

---

## Table of Contents

- [Module Features](#module-features)
- [Getting Started](#getting-started)
- [Installation](#installation)
- [Example Input Files](#example-input-files)
- [GeneViz Result Files](#geneviz-result-files)
- [Command Line Options](#command-line-options)
- [Bug Reports](#bug-reports)
- [License](#license)

---

## Module Features

GeneViz integrates modular tools for microsynteny visualization.

- **LinkView_Visualizer**: Generates publication-quality pairwise linkview plots displaying two genomic regions, BLASTN homology ribbons, gene annotations, and flanking TE tracks. Supports any species and any chromosome.
- **MicroSynteny_Plotter**: Batch-processes multiple gene pairs and produces standardized linkview PDFs for systematic comparison.

---

## Getting Started

**Option A (recommended): clone + conda + editable install**

```bash
git clone https://github.com/YiyongZhao/GeneViz.git
cd GeneViz
conda env create -f environment.yml
conda activate GeneViz
python -m pip install -e .
GeneViz -h
```

**Option B: install from PyPI**

```bash
pip install GeneViz
GeneViz -h
```

**Quick start**

```bash
GeneViz LinkView_Visualizer \
    --input_pairs example_data/LinkView_Visualizer/gene_pairs.tsv \
    --blast_result example_data/LinkView_Visualizer/blastn.out \
    --gff_a example_data/LinkView_Visualizer/speciesA.gff \
    --gff_b example_data/LinkView_Visualizer/speciesB.gff \
    --lens_a example_data/LinkView_Visualizer/speciesA.lens \
    --lens_b example_data/LinkView_Visualizer/speciesB.lens \
    --output_dir results/
```

---

## Installation

> Important: Python 3.13 is NOT supported. Please use Python 3.8-3.12.

```bash
conda create -n GeneViz python=3.12 -y
conda activate GeneViz
pip install GeneViz
```

Required dependencies: Python 3.8-3.12, matplotlib >= 3.5, numpy, pandas, biopython, scipy

---

## Example Input Files

**gene_pairs.tsv** - Two-column TSV: one gene pair per line

```
geneA_id    geneB_id
LOC_Os07g14000    LOC_Os01g23940
LOC_Os11g03280    LOC_Os11g25130
```

**blastn.out** - Standard BLASTN tabular output (-outfmt 6)

**speciesA.gff / speciesB.gff** - Standard GFF3 genome annotation

**speciesA.lens / speciesB.lens** - Chromosome length file (chr TAB length)

```
Chr1    43270923
Chr2    35937250
```

---

## GeneViz Result Files

| Output | Description |
|--------|-------------|
| `linkview_*.pdf` | Publication-quality pairwise linkview plot (vector PDF) |
| `classification_summary.tsv` | Gene pair classification: Transposed / Segmental / Tandem |
| `te_annotation_summary.tsv` | Flanking TE types for mechanism inference |

---

## Command Line Options

### LinkView_Visualizer

Description: Generates pairwise linkview plots showing BLASTN homology ribbons, gene annotations, and TE tracks for any two genomic regions.

Required parameters:

| Parameter | Description |
|-----------|-------------|
| --input_pairs | TSV file with gene pairs (geneA TAB geneB per line) |
| --blast_result | BLASTN tabular output (-outfmt 6) |
| --gff_a | GFF3 annotation for genome A |
| --gff_b | GFF3 annotation for genome B |
| --lens_a | Chromosome lengths for genome A (chr TAB length) |
| --lens_b | Chromosome lengths for genome B (chr TAB length) |

Optional parameters:

| Parameter | Description |
|-----------|-------------|
| --window | Flanking window size in bp, default = 50000 |
| --highlight_color | Color for focal gene pair ribbon, default = #E74C3C |
| --output_dir | Output directory, default = current working directory |

Usage:
```bash
GeneViz LinkView_Visualizer \
    --input_pairs gene_pairs.tsv \
    --blast_result blastn.out \
    --gff_a A.gff --gff_b B.gff \
    --lens_a A.lens --lens_b B.lens \
    [--window 50000] [--output_dir DIR]
```

---

### MicroSynteny_Plotter

Description: Batch-processes multiple gene pairs and generates standardized linkview PDFs in bulk.

Usage:
```bash
GeneViz MicroSynteny_Plotter \
    --input_list batch_pairs.txt \
    --blast_result blastn.out \
    --gff_a A.gff --gff_b B.gff \
    --lens_a A.lens --lens_b B.lens \
    [--output_dir DIR]
```

---

### GD_Classifier

Description: Classifies gene duplication type based on BLASTN ribbon pattern. Transposed = single focal ribbon with no flanking connections. Segmental = multiple parallel ribbons across the region.

Usage:
```bash
GeneViz GD_Classifier \
    --input_pairs gene_pairs.tsv \
    --blast_result blastn.out \
    --gff_a A.gff --gff_b B.gff \
    [--window 50000] [--output_dir DIR]
```

---

### TE_Annotator

Description: Overlays TE annotations onto linkview plots and infers likely transposition mechanism (Pack-MULE, Pack-Helitron, LTR retrotransposition) based on flanking TE family composition.

Usage:
```bash
GeneViz TE_Annotator \
    --input_pairs gene_pairs.tsv \
    --te_gff te_annotation.gff \
    --gff_a A.gff --gff_b B.gff \
    [--window 50000] [--output_dir DIR]
```

---

## Bug Reports

Report bugs or request features via [GitHub Issues](https://github.com/YiyongZhao/GeneViz/issues).
Questions: yiyong.zhao@yale.edu

---

## Contributing

Contributions welcome! Please check [Contribution Guidelines](CONTRIBUTING.md).

---

## License

GeneViz is licensed under the [MIT LICENSE](LICENSE).
