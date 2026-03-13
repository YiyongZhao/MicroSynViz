#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MicroSynViz: Visualizing Pairwise Genomic Microsynteny Within and Between Species
================================================================================================
A professional tool for visualizing gene synteny with TE tracks,
gene structures, and BLAST-based homology links. Designed for
scientific publication figures.

Features:
- Single-species and cross-species comparisons
- Automatic detection of reverse complement based on BLAST alignment
- TE track with motif labels
- Gene structure (exons/introns) visualization
- Core range highlighting
- SVG and PDF output (vector graphics)
- Two input modes: gene IDs or genomic regions (with automatic extension)

Dependencies:
    Python 3.8+
    BioPython
    pandas
    cairosvg
    samtools (command-line)
    blastn (command-line)

Usage example (gene mode, single species):
    MicroSynViz --gene1 GeneA --gene2 GeneB \\
        --fa1 genome.fa --annos1 genes.gff TEs.gff \\
        --fa2 genome.fa --annos2 genes.gff TEs.gff \\
        --auto_complementary --output my_synteny

Usage example (gene mode, cross species):
    MicroSynViz --gene1 GeneA --gene2 GeneB \\
        --fa1 species1.fa --annos1 sp1_genes.gff sp1_TEs.gff \\
        --fa2 species2.fa --annos2 sp2_genes.gff sp2_TEs.gff \\
        --auto_complementary --output cross_synteny

Usage example (region mode, genome FASTA):
    MicroSynViz --region1 "Chr1:1000-2000" --region2 "Chr2:3000-4000" \\
        --fa1 genome.fa --annos1 genes.gff \\
        --fa2 genome.fa --annos2 genes.gff \\
        --extend 5000 --output region_synteny

Usage example (CDS FASTA, no annotations needed):
    MicroSynViz --gene1 GENE_A --gene2 GENE_B \\
        --fa1 cds.fa --fa2 cds.fa --extend 0 --output cds_compare

    MicroSynViz --region1 "LOC01:1-1000" --region2 "LOC02:1-1500" \\
        --fa1 cds.fa --fa2 cds.fa --extend 0 --output cds_region

All parameters use long option style (--xxxx).
"""

import os
import sys
import argparse
import subprocess
import re
import math
import cairosvg
from Bio import SeqIO
from io import StringIO
import pandas as pd
import shutil
import tempfile
import logging

logger = logging.getLogger("MicroSynViz")

__version__ = "1.0.0"  # Updated version with region mode

# -------------------------- XML escaping for SVG text -------------------------
def xml_escape(s):
    """Escape special characters for XML/SVG text content."""
    if s is None:
        return ""
    return (str(s)
            .replace('&', '&amp;')
            .replace('<', '&lt;')
            .replace('>', '&gt;')
            .replace('"', '&quot;')
            .replace("'", '&apos;'))

# -------------------------- SVG style definitions -------------------------
style_css = {
    'simple': '''
        <defs><style>
            .chro{
                stroke: #000000;
                fill: none;
            }
            .line {
                stroke: #000000;
                stroke-width: 3.5px;
            }
            .leader_line {
                stroke: #000000;
                stroke-width: 1px;
            }
            * {
                font-family: Arial, Helvetica, sans-serif;
            }
            text {
                paint-order: stroke;
                stroke: white;
                stroke-width: 5;
                stroke-linecap: round;
                stroke-linejoin: round;
            }
            .feature_label {
                font-size: 20;
                text-anchor: middle;
                fill: black;
            }
            .target_gene_label {
                font-size: 23;
                text-anchor: middle;
                fill: black; /* Carrot Orange */
                font-weight: bold;
            }
            .pos_label {
                font-size: 20;
                fill: black;
            }
            .scale-text {
                font-size: 18;
                fill: black;
            }
        </style></defs>
    ''',
}

# -------------------------- Custom Exceptions -------------------------
class MicroSynVizError(Exception):
    """Base exception for MicroSynViz."""
    pass

class FormatError(MicroSynVizError):
    def __init__(self, fmt_type, filename, line=""):
        super().__init__(f"Invalid {fmt_type} format in {filename}: {line[:100]}")

class InputError(MicroSynVizError):
    """Missing or invalid input files/parameters."""
    pass

class ExternalToolError(MicroSynVizError):
    """External tool (samtools, blastn) not found or failed."""
    pass

class GeneNotFoundError(MicroSynVizError):
    """Gene ID not found in any annotation file."""
    pass

# -------------------------- Chromosome/Drawing Helper Class -------------------------
class Chro():
    """Represents a chromosome track in the SVG drawing."""
    def __init__(self, name, real_name, start, end, is_full):
        self.name = real_name if real_name else name
        self.start = start
        self.end = end
        self.length = end - start + 1
        self.is_full = is_full
        self.left = None
        self.right = None
        self.top = None
        self.level = None

    def coordinate(self, pos, is_up, is_start=False):
        """Convert genomic position to SVG coordinates."""
        if self.left and self.right and self.top:
            if is_up:
                y = self.top
            else:
                y = self.top + args.chro_thickness
            if is_start:
                x = self.left + (pos - self.start) * scale
            else:
                x = self.left + (pos - self.start + 1) * scale
            return (x, y)
        else:
            return None

# -------------------------- GFF Parsing (Gene Features) -------------------------
gene_info = {}      # gene_id -> {'chro':, 'strand':, 'gene':(s,e), 'mRNA':{mrna_id:...}}
mRNA_info = {}      # mrna_id -> {'pos':(s,e), 'exon':[(s,e),...]}

def parse_attributes(attr):
    """Parse GFF3/GTF attribute column into a dictionary.

    Supports both formats:
      GFF3: ID=gene01;Name=ABC      (key=value separated by ;)
      GTF:  gene_id "gene01"; transcript_id "t01";  (key "value" separated by ;)
    """
    res = {}
    if not attr.endswith(';'):
        attr += ';'
    contentFlag = False
    lastPos = 0
    for i, c in enumerate(attr):
        if c == '"':
            contentFlag = not contentFlag
        elif c == ';' and not contentFlag:
            item = attr[lastPos:i].strip()
            lastPos = i + 1
            if not item:
                continue
            # GFF3 format: key=value
            if '=' in item and '"' not in item.split('=', 1)[0]:
                tmpArr = item.split('=', 1)
                key = tmpArr[0].strip()
                val = tmpArr[1].strip(' "') if len(tmpArr) >= 2 else None
            else:
                # GTF format: key "value"
                parts = item.split(None, 1)
                if len(parts) == 2:
                    key = parts[0].strip()
                    val = parts[1].strip(' "')
                elif len(parts) == 1:
                    key = parts[0].strip()
                    val = None
                else:
                    continue
            res[key] = val
            # GTF -> GFF3 key mapping
            if key == 'gene_id' and 'ID' not in res:
                res['ID'] = val
            elif key == 'transcript_id':
                if 'ID' not in res or res.get('gene_id') == res.get('ID'):
                    res['ID'] = val
                res.setdefault('Parent', res.get('gene_id', val))
    return res

def parse_gene_gff(gff_file):
    """
    Parse gene GFF file and populate gene_info and mRNA_info.
    Supports both single and multi-exon genes.
    """
    global gene_info, mRNA_info
    gene_type_name_set = {'gene', 'GENE'}
    mRNA_type_name_set = {'mRNA', 'mrna', 'MRNA', 'transcript'}
    exon_type_name_set = {'exon', 'EXON', 'CDS', 'cds'}

    logger.info(f" Parsing gene GFF: {gff_file}")
    sample_count = 0
    with open(gff_file, 'r') as GFF:
        for line in GFF:
            if line.startswith('#'):
                continue
            line = line.strip()
            if line == '':
                continue
            items = line.split('\t')
            if len(items) != 9:
                continue  # skip non-standard lines silently
            seqid, source, feat, start, end, score, strand, phase, attr = items
            start, end = int(start), int(end)
            att = parse_attributes(attr)

            if feat in gene_type_name_set:
                gene_id = att.get('ID')
                if gene_id is None:
                    continue
                gene_info.setdefault(gene_id, {})
                gene_info[gene_id]['chro'] = seqid
                gene_info[gene_id]['strand'] = strand
                gene_info[gene_id]['gene'] = (start, end)
                if sample_count < 5:
                    logger.debug(f" Gene: {gene_id} at {seqid}:{start}-{end} strand={strand}")
                    sample_count += 1
            elif feat in mRNA_type_name_set:
                if 'Parent' not in att:
                    # If mRNA lacks Parent, treat as gene itself
                    gene_id = f"{att['ID']}_gene"
                    gene_info.setdefault(gene_id, {})
                    gene_info[gene_id]['chro'] = seqid
                    gene_info[gene_id]['strand'] = strand
                    gene_info[gene_id]['gene'] = (start, end)
                    att['Parent'] = gene_id
                parent = att['Parent']
                mrna_id = att['ID']
                gene_info.setdefault(parent, {})
                gene_info[parent].setdefault('mRNA', {})
                mRNA_info.setdefault(mrna_id, {})
                mRNA_info[mrna_id]['pos'] = (start, end)
                gene_info[parent]['mRNA'][mrna_id] = mRNA_info[mrna_id]
            elif feat in exon_type_name_set:
                if 'Parent' not in att:
                    raise FormatError('gff_no_parent', gff_file, line)
                parent = att['Parent']
                if parent not in mRNA_info:
                    continue
                mRNA_info[parent].setdefault('exon', [])
                mRNA_info[parent]['exon'].append((start, end))

# -------------------------- TE GFF Parsing -------------------------
te_info = []   # list of (file_id, chro, start, end, motif_name, strand)

def parse_te_gff(te_gff_file, source_order=0):
    """
    Parse TE GFF file and populate te_info.
    Recognizes common TE feature types and extracts motif name if available.
    """
    global te_info
    te_type_name_set = {'dispersed_repeat', 'transposable_element', 'repeat_region',
                        'LTR_retrotransposon', 'DNA_transposon', 'TE', 'LINE', 'SINE', 'helitron'}
    id_pattern = re.compile(r'ID=([^;]+)')
    motif_pattern = re.compile(r'Target\s+"Motif:([^"\s]+)')

    logger.info(f" Parsing TE file: {te_gff_file}")
    sample_count = 0
    te_count = 0
    with open(te_gff_file, 'r') as GFF:
        for line in GFF:
            if line.startswith('#'):
                continue
            line = line.strip()
            if line == '':
                continue
            items = line.split('\t')
            if len(items) != 9:
                continue  # skip non-standard lines silently
            seqid, source, feat, start, end, score, strand, phase, attr = items
            start, end = int(start), int(end)
            if feat not in te_type_name_set:
                continue
            id_match = id_pattern.search(attr)
            if id_match:
                file_id = id_match.group(1)
            else:
                file_id = f"{seqid}_{start}_{end}"
            motif_match = motif_pattern.search(attr)
            if motif_match:
                motif_name = motif_match.group(1)
            else:
                motif_name = "unknown"
            te_info.append((file_id, seqid, start, end, motif_name, strand, source_order))
            te_count += 1
            if sample_count < 5:
                logger.debug(f" TE {file_id}: {seqid}:{start}-{end} motif={motif_name} strand={strand}")
                sample_count += 1
    logger.info(f" TE parsing complete. TE entries stored: {len(te_info)}")



# -------------------------- BED Parsing -------------------------
def detect_format(filepath):
    """Auto-detect annotation file format: gff/gtf or bed.
    
    Heuristic: read first non-comment, non-empty line.
    - 9 tab-separated columns → GFF/GTF
    - 3-12 tab-separated columns with int cols[1], int cols[2] → BED
    """
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) == 9:
                return 'gff'
            if 3 <= len(cols) <= 12:
                try:
                    int(cols[1])
                    int(cols[2])
                    return 'bed'
                except ValueError:
                    pass
            return 'gff'  # fallback
    return 'gff'


def parse_bed_genes(bed_file):
    """Parse BED12 file for gene structure (exons from blocks).
    
    BED12 columns: chrom, start, end, name, score, strand,
                   thickStart, thickEnd, rgb, blockCount, blockSizes, blockStarts
    BED6 or BED3: treated as single-exon genes.
    """
    global gene_info, mRNA_info
    logger.info(f" Parsing BED (gene mode): {bed_file}")
    count = 0
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 3:
                continue
            chrom = cols[0]
            start = int(cols[1]) + 1  # BED is 0-based, convert to 1-based
            end = int(cols[2])
            name = cols[3] if len(cols) > 3 else f"{chrom}_{start}_{end}"
            strand = cols[5] if len(cols) > 5 else '+'

            gene_info.setdefault(name, {})
            gene_info[name]['chro'] = chrom
            gene_info[name]['strand'] = strand
            gene_info[name]['gene'] = (start, end)

            # BED12: parse blocks as exons
            if len(cols) >= 12:
                block_count = int(cols[9])
                block_sizes = [int(x) for x in cols[10].rstrip(',').split(',') if x]
                block_starts = [int(x) for x in cols[11].rstrip(',').split(',') if x]
                mrna_id = f"{name}.1"
                mRNA_info.setdefault(mrna_id, {})
                mRNA_info[mrna_id]['pos'] = (start, end)
                mRNA_info[mrna_id]['exon'] = []
                for i in range(block_count):
                    exon_start = int(cols[1]) + block_starts[i] + 1  # 0-based to 1-based
                    exon_end = exon_start + block_sizes[i] - 1
                    mRNA_info[mrna_id]['exon'].append((exon_start, exon_end))
                gene_info[name].setdefault('mRNA', {})
                gene_info[name]['mRNA'][mrna_id] = mRNA_info[mrna_id]
            else:
                # BED3/6: single exon = whole gene
                mrna_id = f"{name}.1"
                mRNA_info.setdefault(mrna_id, {})
                mRNA_info[mrna_id]['pos'] = (start, end)
                mRNA_info[mrna_id]['exon'] = [(start, end)]
                gene_info[name].setdefault('mRNA', {})
                gene_info[name]['mRNA'][mrna_id] = mRNA_info[mrna_id]
            count += 1
    logger.info(f" BED gene parsing complete: {count} entries")


def parse_bed_te(bed_file, source_order=0):
    """Parse BED file for TE annotations."""
    global te_info
    logger.info(f" Parsing BED (TE mode): {bed_file}")
    count = 0
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 3:
                continue
            chrom = cols[0]
            start = int(cols[1]) + 1  # 0-based to 1-based
            end = int(cols[2])
            name = cols[3] if len(cols) > 3 else "unknown"
            strand = cols[5] if len(cols) > 5 else '+'
            file_id = f"{chrom}_{start}_{end}"
            te_info.append((file_id, chrom, start, end, name, strand, source_order))
            count += 1
    logger.info(f" BED TE parsing complete: {count} entries")


def parse_annotation(filepath, source_order=0):
    """Universal annotation parser: auto-detect format and parse for both gene and TE features.
    
    source_order: controls track layering (0 = closest to chromosome bar, 1 = next layer out, etc.)
    """
    fmt = detect_format(filepath)
    logger.info(f"  {filepath} [detected: {fmt.upper()}, layer={source_order}]\n")
    if fmt == 'bed':
        parse_bed_genes(filepath)
        parse_bed_te(filepath, source_order=source_order)
    else:
        # GFF/GTF: try both gene and TE parsing
        try:
            parse_gene_gff(filepath)
        except Exception as e:
            sys.stderr.write(f"  [WARNING] Gene parsing failed for {filepath}: {e}\n")
        try:
            parse_te_gff(filepath, source_order=source_order)
        except Exception:
            pass  # Not all GFFs have TE features


def extract_gene_from_bed(gene_id, bed_file):
    """Extract gene coordinates from a BED file."""
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 4:
                continue
            if cols[3] == gene_id:
                chrom = cols[0]
                start = int(cols[1]) + 1  # 0-based to 1-based
                end = int(cols[2])
                strand = cols[5] if len(cols) > 5 else '+'
                return chrom, start, end, strand
    return None

# -------------------------- Utility Functions -------------------------
def check_file_exists(file_path):
    """Raise InputError if file does not exist."""
    if not os.path.isfile(file_path):
        raise InputError(f"File not found: {file_path}")

def _find_gene_in_annos(gene_id, anno_files):
    """Search multiple annotation files for a gene. Returns (chr, start, end, strand) or None."""
    for anno_file in anno_files:
        try:
            fmt = detect_format(anno_file)
            if fmt == 'bed':
                result = extract_gene_from_bed(gene_id, anno_file)
                if result:
                    return result
            else:
                with open(anno_file, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) != 9:
                            continue
                        if parts[2] != 'gene':
                            continue
                        att = parse_attributes(parts[8])
                        if (att.get('ID') == gene_id or att.get('Name') == gene_id or
                                att.get('gene_id') == gene_id or att.get('gene_name') == gene_id):
                            return parts[0], int(parts[3]), int(parts[4]), parts[6]
        except Exception:
            continue
    return None


def extract_gene_coordinates(gene_id, anno_file):
    """Extract gene coordinates from GFF3, GTF, or BED file."""
    check_file_exists(anno_file)
    fmt = detect_format(anno_file)
    if fmt == 'bed':
        result = extract_gene_from_bed(gene_id, anno_file)
        if result:
            return result
        raise GeneNotFoundError(f"Gene '{gene_id}' not found in BED file {anno_file}")
    # GFF/GTF
    with open(anno_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue
            if parts[2] != 'gene':
                continue
            att = parse_attributes(parts[8])
            match = (att.get('ID') == gene_id or
                     att.get('Name') == gene_id or
                     att.get('gene_id') == gene_id or
                     att.get('gene_name') == gene_id)
            if match:
                return parts[0], int(parts[3]), int(parts[4]), parts[6]
    raise GeneNotFoundError(f"Gene '{gene_id}' not found in {anno_file}")



def _find_gene_in_fasta_index(gene_id, fasta_file):
    """Fallback: search FASTA index (.fai) for a sequence matching gene_id.
    
    This enables CDS/transcript FASTA input where sequence headers ARE the gene IDs.
    Auto-creates .fai index if missing.
    Returns (seq_name, 1, length, '+') or None.
    """
    fai_file = fasta_file + ".fai"
    
    # Auto-index if .fai doesn't exist
    if not os.path.isfile(fai_file):
        try:
            subprocess.check_call(["samtools", "faidx", fasta_file])
            logger.info(f" Auto-indexed: {fasta_file}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None
    
    if not os.path.isfile(fai_file):
        return None
    
    with open(fai_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 2 and cols[0] == gene_id:
                seq_len = int(cols[1])
                return (gene_id, 1, seq_len, '+')
    return None

def extract_sequence_samtools(chr_id, start, end, gene_id, genome_fa, extend=3000):
    """Extract genomic region using samtools faidx (gene mode)."""
    check_file_exists(genome_fa)
    region_start = max(1, start - extend)
    region_end = end + extend
    region = f"{chr_id}:{region_start}-{region_end}"
    cmd = ["samtools", "faidx", genome_fa, region]
    try:
        output = subprocess.check_output(cmd, encoding="utf-8")
    except subprocess.CalledProcessError:
        raise ExternalToolError(f"samtools failed to extract region: {region}")
    except FileNotFoundError:
        raise ExternalToolError("samtools not found in PATH. Install: conda install -c bioconda samtools")
    lines = output.strip().split("\n")
    if len(lines) < 2 or len(lines[1].strip()) == 0:
        raise InputError(f"Empty sequence extracted for {gene_id}")
    fixed_fasta = [f">{gene_id}"] + lines[1:]
    seq_str = "\n".join(fixed_fasta)
    # Adjust region_end based on actual extracted sequence length
    actual_len = sum(len(l.strip()) for l in lines[1:])
    actual_end = region_start + actual_len - 1
    if actual_end < region_end:
        region_end = actual_end
    return seq_str, region_start, region_end

def extract_region_samtools(chr_id, start, end, genome_fa, extend=0):
    """Extract genomic region using samtools faidx (region mode)."""
    check_file_exists(genome_fa)
    region_start = max(1, start - extend)
    region_end = end + extend
    region = f"{chr_id}:{region_start}-{region_end}"
    cmd = ["samtools", "faidx", genome_fa, region]
    try:
        output = subprocess.check_output(cmd, encoding="utf-8")
    except subprocess.CalledProcessError:
        raise ExternalToolError(f"samtools failed to extract region: {region}")
    except FileNotFoundError:
        raise ExternalToolError("samtools not found in PATH. Install: conda install -c bioconda samtools")
    lines = output.strip().split("\n")
    if len(lines) < 2 or len(lines[1].strip()) == 0:
        raise InputError(f"Empty sequence extracted for region {region}")
    # Adjust region_end based on actual extracted sequence length
    actual_len = sum(len(l.strip()) for l in lines[1:])
    actual_end = region_start + actual_len - 1
    if actual_end < region_end:
        region_end = actual_end
    # Use region identifier as FASTA header
    header = f">{chr_id}:{region_start}-{region_end}"
    fixed_fasta = [header] + lines[1:]
    seq_str = "\n".join(fixed_fasta)
    return seq_str, region_start, region_end

def run_blastn(seq_fasta, blast_out, evalue=1e-5, threads=8):
    """Run BLASTn with tabular output format 6."""
    check_file_exists(seq_fasta)
    tmp_dir = tempfile.mkdtemp(prefix="microsynviz_blast_")
    tmp_db = os.path.join(tmp_dir, "blast_db")
    makeblastdb_cmd = ["makeblastdb", "-in", seq_fasta, "-dbtype", "nucl", "-out", tmp_db]
    try:
        subprocess.run(makeblastdb_cmd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        raise ExternalToolError(f"Failed to build BLAST database: {e}")
    blast_cmd = ["blastn", "-query", seq_fasta, "-db", tmp_db,
                 "-evalue", str(evalue), "-out", blast_out,
                 "-outfmt", "6", "-num_threads", str(threads)]
    try:
        subprocess.run(blast_cmd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        raise ExternalToolError(f"BLASTn alignment failed: {e}")
    finally:
        # Cleanup temporary database directory
        shutil.rmtree(tmp_dir, ignore_errors=True)

def read_sequence_lengths(fasta_file):
    """Return dictionary of sequence lengths from FASTA."""
    check_file_exists(fasta_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return {gene_id: len(seq) for gene_id, seq in seq_dict.items()}

def parse_blast_results(blast_file, min_identity=0, min_length=0):
    """Parse BLAST tabular output (outfmt 6) into list of dicts.
    
    Args:
        blast_file: Path to BLAST outfmt 6 file
        min_identity: Minimum percent identity to keep (default: 0)
        min_length: Minimum alignment length in bp to keep (default: 0)
    """
    blast_hits = []
    if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
        blast_cols = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        df = pd.read_csv(blast_file, sep="\t", names=blast_cols)
        n_total = len(df)
        if min_identity > 0:
            df = df[df["pident"] >= min_identity]
        if min_length > 0:
            df = df[df["length"] >= min_length]
        n_kept = len(df)
        if n_total > 0:
            logger.info(f"BLAST hits: {n_total} total, {n_kept} passed filters "
                           f"(identity >= {min_identity}%, length >= {min_length} bp)\n")
        blast_hits = df.to_dict("records")
    return blast_hits

def parse_region(region_str):
    """Parse region string like 'Chr1:1000-2000'."""
    pattern = re.compile(r'^([^:]+):(\d+)-(\d+)$')
    m = pattern.match(region_str)
    if not m:
        raise ValueError(f"Invalid region format: {region_str}. Expected e.g., Chr1:1-5000")
    chr_name = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))
    if start > end:
        raise ValueError(f"Start must be <= end in region: {region_str}")
    return chr_name, start, end

# -------------------------- SVG Generation (Core) -------------------------
def generate_svg(gene1, gene2, len1, len2, blast_hits,
                 seq1_start, seq1_end, seq2_start, seq2_end,
                 chr1, chr2, out_prefix, bezier=False, target_set=None,
                 core_ranges=None, revcomp2=False):
    """
    Generate SVG visualization of synteny.
    Returns path to the generated SVG file.
    """
    global scale, args

    chro_lst = {}
    chro1 = Chro(gene1, gene1, seq1_start, seq1_end, [True, True])
    chro2 = Chro(gene2, gene2, seq2_start, seq2_end, [True, True])
    chro_lst[gene1] = chro1
    chro_lst[gene2] = chro2

    max_len = max(len1, len2)
    scale = args.svg_width * (1 - args.svg_space) / max_len

    space_vertical = args.svg_height / (2 + 1)
    top1 = space_vertical
    top2 = space_vertical * 2

    left1 = (args.svg_width - len1 * scale) / 2
    left2 = (args.svg_width - len2 * scale) / 2

    chro1.level = 0
    chro1.top = top1
    chro1.left = left1
    chro1.right = left1 + len1 * scale

    chro2.level = 1
    chro2.top = top2
    chro2.left = left2
    chro2.right = left2 + len2 * scale

    svg_content_parts = []
    svg_content_parts.append(re.sub(r'\s+', '', style_css['simple']))

    # Gradient definition for color legend
    svg_content_parts.append('''
        <defs>
            <linearGradient id="grayPurpleRedGradient" x1="0%" y1="100%" x2="0%" y2="0%">
                <stop offset="0%" stop-color="rgb(200,200,200)"/>
                <stop offset="50%" stop-color="rgb(128,0,128)"/>
                <stop offset="100%" stop-color="rgb(255,0,0)"/>
            </linearGradient>
        </defs>
    ''')

    # Draw chromosome rectangles
    svg_content_parts.append(
        f'<rect x="{chro1.left}" y="{chro1.top}" width="{chro1.right - chro1.left}" '
        f'height="{args.chro_thickness}" class="chro"/>'
    )
    svg_content_parts.append(
        f'<rect x="{chro2.left}" y="{chro2.top}" width="{chro2.right - chro2.left}" '
        f'height="{args.chro_thickness}" class="chro"/>'
    )

    # ========== Plot BLAST hits ==========
    relevant_hits = []
    for hit in blast_hits:
        qid = hit['qseqid']
        sid = hit['sseqid']
        if (qid == gene1 and sid == gene2) or (qid == gene2 and sid == gene1):
            relevant_hits.append(hit)

    # Remove duplicate hits (keep highest bitscore)
    unique_hits = {}
    for hit in relevant_hits:
        qid = hit['qseqid']
        sid = hit['sseqid']
        if qid == gene1 and sid == gene2:
            g1_start_abs = seq1_start + hit['qstart'] - 1
            g1_end_abs   = seq1_start + hit['qend'] - 1
            g2_start_abs = seq2_start + hit['sstart'] - 1
            g2_end_abs   = seq2_start + hit['send'] - 1
        else:
            g1_start_abs = seq1_start + hit['sstart'] - 1
            g1_end_abs   = seq1_start + hit['send'] - 1
            g2_start_abs = seq2_start + hit['qstart'] - 1
            g2_end_abs   = seq2_start + hit['qend'] - 1
        key = (g1_start_abs, g1_end_abs, g2_start_abs, g2_end_abs)
        if key not in unique_hits or hit['bitscore'] > unique_hits[key]['bitscore']:  # always dedup by bitscore
            unique_hits[key] = hit

    relevant_hits = list(unique_hits.values())

    # Score range for legend (based on --color_by)
    color_metric = args.color_by
    legend_min = legend_max = None
    if relevant_hits:
        scores = [hit[color_metric if color_metric != 'identity' else 'pident'] for hit in relevant_hits]
        if color_metric == 'evalue':
            # Use -log10(evalue) for color mapping (higher = better match)
            scores = [-math.log10(max(s, 1e-300)) for s in scores]
        legend_min = min(scores)
        legend_max = max(scores)

        # Color endpoints
        low_rgb = (200, 200, 200)   # gray
        mid_rgb = (128, 0, 128)     # purple
        high_rgb = (255, 0, 0)       # red

        for hit in relevant_hits:
            qid = hit['qseqid']
            sid = hit['sseqid']
            if qid == gene1 and sid == gene2:
                g1_start_abs = seq1_start + hit['qstart'] - 1
                g1_end_abs   = seq1_start + hit['qend'] - 1
                g2_start_abs = seq2_start + hit['sstart'] - 1
                g2_end_abs   = seq2_start + hit['send'] - 1
            else:
                g1_start_abs = seq1_start + hit['sstart'] - 1
                g1_end_abs   = seq1_start + hit['send'] - 1
                g2_start_abs = seq2_start + hit['qstart'] - 1
                g2_end_abs   = seq2_start + hit['qend'] - 1

            # Normalize score for color mapping
            if color_metric == 'evalue':
                score_val = -math.log10(max(hit['evalue'], 1e-300))
            elif color_metric == 'identity':
                score_val = hit['pident']
            else:
                score_val = hit['bitscore']
            if legend_max > legend_min:
                norm = (score_val - legend_min) / (legend_max - legend_min)
            else:
                norm = 0.5

            # Interpolate color
            if norm <= 0.5:
                t = norm / 0.5
                r = int(low_rgb[0] + (mid_rgb[0] - low_rgb[0]) * t)
                g = int(low_rgb[1] + (mid_rgb[1] - low_rgb[1]) * t)
                b = int(low_rgb[2] + (mid_rgb[2] - low_rgb[2]) * t)
            else:
                t = (norm - 0.5) / 0.5
                r = int(mid_rgb[0] + (high_rgb[0] - mid_rgb[0]) * t)
                g = int(mid_rgb[1] + (high_rgb[1] - mid_rgb[1]) * t)
                b = int(mid_rgb[2] + (high_rgb[2] - mid_rgb[2]) * t)

            fill_color = f'rgb({r},{g},{b})'

            # Get polygon vertices
            v1 = chro_lst[gene1].coordinate(g1_start_abs, is_up=False, is_start=True)
            v2 = chro_lst[gene2].coordinate(g2_start_abs, is_up=True,  is_start=True)
            v3 = chro_lst[gene2].coordinate(g2_end_abs,   is_up=True,  is_start=False)
            v4 = chro_lst[gene1].coordinate(g1_end_abs,   is_up=False, is_start=False)

            if not (v1 and v2 and v3 and v4):
                continue

            # Apply gap
            v1 = (v1[0], v1[1] + args.gap)
            v2 = (v2[0], v2[1] - args.gap)
            v3 = (v3[0], v3[1] - args.gap)
            v4 = (v4[0], v4[1] + args.gap)

            if bezier:
                if v1[1] > v2[1]:
                    multipliers = [1, 2, 2, 1]
                else:
                    multipliers = [2, 1, 1, 2]
                bezier_coor1 = [v1[0], min(v1[1], v2[1]) + abs(v1[1] - v2[1]) / 3 * multipliers[0]]
                bezier_coor2 = [v2[0], min(v1[1], v2[1]) + abs(v1[1] - v2[1]) / 3 * multipliers[1]]
                bezier_coor3 = [v3[0], min(v3[1], v4[1]) + abs(v3[1] - v4[1]) / 3 * multipliers[2]]
                bezier_coor4 = [v4[0], min(v3[1], v4[1]) + abs(v3[1] - v4[1]) / 3 * multipliers[3]]
                svg_content_parts.append(
                    f'<path d="M{v1[0]},{v1[1]} C{bezier_coor1[0]},{bezier_coor1[1]} {bezier_coor2[0]},{bezier_coor2[1]} {v2[0]},{v2[1]} L{v3[0]},{v3[1]} C{bezier_coor3[0]},{bezier_coor3[1]} {bezier_coor4[0]},{bezier_coor4[1]} {v4[0]},{v4[1]} Z" fill="{fill_color}" opacity="{args.ribbon_opacity}"/>'
                )
            else:
                svg_content_parts.append(
                    f'<path d="M{v1[0]},{v1[1]} L{v2[0]},{v2[1]} L{v3[0]},{v3[1]} L{v4[0]},{v4[1]} Z" fill="{fill_color}" opacity="{args.ribbon_opacity}"/>'
                )
    else:
        logger.info(" No relevant BLAST hits found between the two genes/regions.")

    # -------------------------- Draw TE track --------------------------
    if te_info:
        logger.info(f" Starting TE track drawing for {gene1} and {gene2}")
        debug_level = int(os.environ.get('DEBUG_TE', '0'))

        def draw_te_track(gene_id, chromosome, chro_obj, position, revcomp):
            te_svg = []
            chro_start = chro_obj.start
            chro_end = chro_obj.end
            track_h = args.te_track_height
            track_offset = args.te_track_offset

            if position == 'top':
                base_y = chro_obj.top - track_offset - track_h
            else:
                base_y = chro_obj.top + args.chro_thickness + track_offset

            te_candidates = []
            matched_chr_count = 0
            for entry in te_info:
                file_id, te_chr, te_start, te_end, motif_name, strand = entry[:6]
                source_order = entry[6] if len(entry) > 6 else 0
                if te_chr.lower() != chromosome.lower():
                    continue
                matched_chr_count += 1

                # Apply reverse complement if needed
                if revcomp:
                    mapped_start = chro_start + chro_end - te_end
                    mapped_end   = chro_start + chro_end - te_start
                    if mapped_start > mapped_end:
                        mapped_start, mapped_end = mapped_end, mapped_start
                else:
                    mapped_start = te_start
                    mapped_end   = te_end

                # Trim to displayed region
                if mapped_end < chro_start or mapped_start > chro_end:
                    continue
                draw_start = max(mapped_start, chro_start)
                draw_end   = min(mapped_end,   chro_end)
                if draw_start > draw_end:
                    continue

                # Convert to SVG x coordinates
                x1, _ = chro_obj.coordinate(draw_start, is_up=True, is_start=True)
                x2, _ = chro_obj.coordinate(draw_end,   is_up=True, is_start=False)
                if x1 is None or x2 is None or x2 - x1 <= 0:
                    continue

                # Adjust strand if reversed
                if revcomp:
                    actual_strand = '+' if strand == '-' else '-'
                else:
                    actual_strand = strand

                te_candidates.append((x1, x2, motif_name, actual_strand, source_order))
                if debug_level >= 1:
                    logger.debug(f" TE {file_id}: original ({te_start},{te_end}) mapped to ({mapped_start},{mapped_end}) -> draw ({draw_start},{draw_end}) strand={actual_strand}")

            # Sort and layout to avoid overlaps
            te_candidates.sort(key=lambda t: t[0])
            used_layers = []
            used_id_offsets = []
            base_id_offset = args.te_offset_base
            id_step = args.te_offset_step

            for (x1, x2, motif_name, strand, src_order) in te_candidates:
                y_layer_offset = 0
                while True:
                    conflict = False
                    for used_entry in used_layers:
                        used_offset, ux1, ux2 = used_entry[0], used_entry[1], used_entry[2]
                        if max(x1, ux1) <= min(x2, ux2) and used_offset == y_layer_offset:
                            conflict = True
                            break
                    if not conflict:
                        break
                    y_layer_offset += track_h + 2
                    if y_layer_offset > 10 * (track_h + 2):
                        y_layer_offset = 0
                        break
                used_layers.append((y_layer_offset, x1, x2, src_order))

                if position == 'top':
                    rect_y = base_y - y_layer_offset
                else:
                    rect_y = base_y + y_layer_offset

                # Draw TE rectangle
                te_svg.append(
                    f'<rect x="{x1}" y="{rect_y}" width="{x2-x1}" height="{track_h}" '
                    f'fill="#1E90FF" fill-opacity="0.5" stroke="#00008B" stroke-width="0.5"/>'
                )

                # Draw strand arrow
                arrow_size = track_h * 0.6
                arrow_y = rect_y + (track_h - arrow_size) / 2
                if strand == '+':
                    arrow_x = x2 - arrow_size
                    points = f"{arrow_x},{arrow_y} {x2},{arrow_y+arrow_size/2} {arrow_x},{arrow_y+arrow_size}"
                else:
                    arrow_x = x1 + arrow_size
                    points = f"{x1},{arrow_y} {arrow_x},{arrow_y+arrow_size/2} {x1},{arrow_y+arrow_size}"
                te_svg.append(f'<polygon points="{points}" fill="white" />')

                # Draw motif label with leader line
                mid_x = (x1 + x2) / 2
                if position == 'top':
                    rect_edge = rect_y
                else:
                    rect_edge = rect_y + track_h

                id_offset = base_id_offset
                while True:
                    conflict = False
                    for (used_offset, ux1, ux2) in used_id_offsets:
                        if max(x1, ux1) <= min(x2, ux2) and abs(used_offset - id_offset) < id_step:
                            conflict = True
                            break
                    if not conflict:
                        break
                    id_offset += id_step
                    if id_offset > base_id_offset + 10 * id_step:
                        id_offset = base_id_offset
                        break
                used_id_offsets.append((id_offset, x1, x2))

                if position == 'top':
                    line_end_y = rect_edge - id_offset
                else:
                    line_end_y = rect_edge + id_offset

                te_svg.append(
                    f'<line x1="{mid_x}" y1="{rect_edge}" x2="{mid_x}" y2="{line_end_y}" '
                    f'stroke="black" stroke-width="1"/>'
                )
                escaped_motif = xml_escape(motif_name)
                label_y = line_end_y - 2 if position == 'top' else line_end_y + 12
                te_svg.append(
                    f'<text x="{mid_x}" y="{label_y}" fill="black" font-size="8" '
                    f'text-anchor="middle">{escaped_motif}</text>'
                )

            logger.info(f" TE track for {gene_id}: matched chromosome count = {matched_chr_count}, {len(te_candidates)} TE(s) drawn.")
            return ''.join(te_svg)

        svg_content_parts.append(draw_te_track(gene1, chr1, chro1, 'top', revcomp=False))
        svg_content_parts.append(draw_te_track(gene2, chr2, chro2, 'bottom', revcomp2))
    else:
        logger.info(" te_info is empty, skipping TE drawing.")

    # -------------------------- Draw all gene structures --------------------------
    if gene_info:
        def draw_all_gene_structures(chromosome, chro_obj, target_set, core_range, revcomp):
            gene_svg = []
            chro_start = chro_obj.start
            chro_end = chro_obj.end
            core_start, core_end = core_range if core_range else (None, None)
            chromosome_lower = chromosome.lower()

            # Debug info
            if not hasattr(draw_all_gene_structures, "chromosomes_printed"):
                all_chromosomes = set(ginfo['chro'].lower() for ginfo in gene_info.values() if 'chro' in ginfo)
                logger.debug(f" All chromosome names in gene_info (lowercase): {sorted(all_chromosomes)}")
                draw_all_gene_structures.chromosomes_printed = True

            genes_on_chr = [gid for gid, ginfo in gene_info.items() if ginfo['chro'].lower() == chromosome_lower]
            logger.debug(f" Genes on {chromosome} (lower={chromosome_lower}): {genes_on_chr}")

            # Adjust core range if reversed
            if revcomp and core_start is not None and core_end is not None:
                core_start, core_end = chro_start + chro_end - core_end, chro_start + chro_end - core_start
                if core_start > core_end:
                    core_start, core_end = core_end, core_start
                logger.debug(f" Mapped core range for revcomp: ({core_start}, {core_end})")

            COLOR_TARGET = "#FF0000"
            COLOR_OTHER  = "#228B22"

            for gid, ginfo in gene_info.items():
                if ginfo['chro'].lower() != chromosome_lower:
                    continue
                if 'gene' not in ginfo:
                    continue
                g_start, g_end = ginfo['gene']
                logger.debug(f" candidate gene {gid} original ({g_start},{g_end})")

                if revcomp:
                    mapped_g_start = chro_start + chro_end - g_end
                    mapped_g_end   = chro_start + chro_end - g_start
                    if mapped_g_start > mapped_g_end:
                        mapped_g_start, mapped_g_end = mapped_g_end, mapped_g_start
                    strand = '+' if ginfo['strand'] == '-' else '-'
                    logger.debug(f" revcomp gene {gid}: mapped to ({mapped_g_start},{mapped_g_end}) strand={strand}")
                else:
                    mapped_g_start = g_start
                    mapped_g_end   = g_end
                    strand = ginfo['strand']

                # Check if gene overlaps display region
                if mapped_g_end < chro_start or mapped_g_start > chro_end:
                    logger.debug(f" gene {gid} out of display region after mapping, skipping")
                    continue
                gene_region_start = max(mapped_g_start, chro_start)
                gene_region_end = min(mapped_g_end, chro_end)

                # Get longest mRNA and its exons
                mrna_ids = list(ginfo.get('mRNA', {}).keys())
                if not mrna_ids:
                    logger.debug(f" gene {gid} has no mRNA info, skipping")
                    continue
                mrna_id = max(mrna_ids, key=lambda mid: mRNA_info[mid]['pos'][1] - mRNA_info[mid]['pos'][0])
                exons = mRNA_info[mrna_id].get('exon', [])
                if not exons:
                    logger.debug(f" gene {gid} has no exon info, skipping")
                    continue

                # Map exons to display coordinates
                mapped_exons = []
                for (ex_start, ex_end) in exons:
                    if revcomp:
                        mex_start = chro_start + chro_end - ex_end
                        mex_end   = chro_start + chro_end - ex_start
                        if mex_start > mex_end:
                            mex_start, mex_end = mex_end, mex_start
                    else:
                        mex_start, mex_end = ex_start, ex_end
                    if mex_end < chro_start or mex_start > chro_end:
                        continue
                    mex_start = max(mex_start, chro_start)
                    mex_end   = min(mex_end,   chro_end)
                    if mex_start <= mex_end:
                        mapped_exons.append((mex_start, mex_end))

                if not mapped_exons:
                    logger.debug(f" gene {gid} has no exon in display region after mapping")
                    continue
                mapped_exons.sort()
                logger.debug(f" gene {gid} mapped exons: {mapped_exons}")

                # Introns (gaps between exons)
                introns = []
                for i in range(len(mapped_exons) - 1):
                    intron_start = mapped_exons[i][1] + 1
                    intron_end = mapped_exons[i + 1][0] - 1
                    if intron_start <= intron_end:
                        introns.append((intron_start, intron_end))

                # Determine arrow position (direction of transcription)
                arrow_index = 0 if strand == '-' else len(mapped_exons) - 1
                arrow_pos = gene_region_end if strand == '+' else gene_region_start

                # Draw exons (possibly split by core range)
                for i, (ex_start, ex_end) in enumerate(mapped_exons):
                    if core_start is not None and core_end is not None:
                        segments = []
                        if ex_start < core_start:
                            left_end = min(ex_end, core_start - 1)
                            if left_end >= ex_start:
                                segments.append((ex_start, left_end, False))
                        inter_start = max(ex_start, core_start)
                        inter_end = min(ex_end, core_end)
                        if inter_start <= inter_end:
                            segments.append((inter_start, inter_end, True))
                        if ex_end > core_end:
                            right_start = max(ex_start, core_end + 1)
                            if right_start <= ex_end:
                                segments.append((right_start, ex_end, False))
                    else:
                        segments = [(ex_start, ex_end, False)]

                    for seg_start, seg_end, in_core in segments:
                        x1, y = chro_obj.coordinate(seg_start, is_up=True, is_start=True)
                        x2, y = chro_obj.coordinate(seg_end, is_up=True, is_start=False)
                        if x1 is None or x2 is None or x2 - x1 <= 0:
                            logger.debug(f" segment ({seg_start},{seg_end}) gave invalid coordinates")
                            continue

                        # Determine color
                        if gid in target_set:
                            if core_start is not None and core_end is not None:
                                fill = COLOR_TARGET if in_core else COLOR_OTHER
                            else:
                                fill = COLOR_TARGET
                        else:
                            fill = COLOR_OTHER

                        # Check if this segment contains the arrow (for transcription direction)
                        contains_arrow = (i == arrow_index) and (seg_start <= arrow_pos <= seg_end)

                        if contains_arrow:
                            y_top = y
                            y_bottom = y + args.chro_thickness
                            y_mid = y + args.chro_thickness / 2
                            length = x2 - x1
                            if length <= 0:
                                continue
                            arrow_len = length * 0.2
                            if strand == '+':
                                x_arrow_start = x2 - arrow_len
                                vertices = [
                                    (x1, y_top),
                                    (x_arrow_start, y_top),
                                    (x2, y_mid),
                                    (x_arrow_start, y_bottom),
                                    (x1, y_bottom)
                                ]
                            else:
                                x_arrow_end = x1 + arrow_len
                                vertices = [
                                    (x2, y_top),
                                    (x2, y_bottom),
                                    (x_arrow_end, y_bottom),
                                    (x1, y_mid),
                                    (x_arrow_end, y_top)
                                ]
                            gene_svg.append(
                                f'<path d="M{vertices[0][0]},{vertices[0][1]} L{vertices[1][0]},{vertices[1][1]} L{vertices[2][0]},{vertices[2][1]} L{vertices[3][0]},{vertices[3][1]} L{vertices[4][0]},{vertices[4][1]} Z" fill="{fill}" stroke="white" stroke-width="0.5"/>'
                            )
                        else:
                            gene_svg.append(
                                f'<rect x="{x1}" y="{y}" width="{x2-x1}" height="{args.chro_thickness}" fill="{fill}" stroke="white" stroke-width="0.5"/>'
                            )

                # Draw introns as thin lines
                for int_start, int_end in introns:
                    if int_end < chro_start or int_start > chro_end:
                        continue
                    int_start = max(int_start, chro_start)
                    int_end = min(int_end, chro_end)
                    if int_start > int_end:
                        continue

                    if core_start is not None and core_end is not None:
                        segments = []
                        if int_start < core_start:
                            left_end = min(int_end, core_start - 1)
                            if left_end >= int_start:
                                segments.append((int_start, left_end, False))
                        inter_start = max(int_start, core_start)
                        inter_end = min(int_end, core_end)
                        if inter_start <= inter_end:
                            segments.append((inter_start, inter_end, True))
                        if int_end > core_end:
                            right_start = max(int_start, core_end + 1)
                            if right_start <= int_end:
                                segments.append((right_start, int_end, False))
                    else:
                        segments = [(int_start, int_end, False)]

                    for seg_start, seg_end, in_core in segments:
                        x1, y = chro_obj.coordinate(seg_start, is_up=True, is_start=True)
                        x2, y = chro_obj.coordinate(seg_end, is_up=True, is_start=False)
                        y += args.chro_thickness * 0.5
                        if gid in target_set:
                            if core_start is not None and core_end is not None:
                                line_color = COLOR_TARGET if in_core else COLOR_OTHER
                            else:
                                line_color = COLOR_TARGET
                        else:
                            line_color = COLOR_OTHER
                        gene_svg.append(
                            f'<line x1="{x1}" y1="{y}" x2="{x2}" y2="{y}" stroke="{line_color}" stroke-width="1.5"/>'
                        )

            return ''.join(gene_svg)

        if target_set is None:
            target_set = {gene1, gene2}

        core_range1 = core_ranges.get(chr1, None) if core_ranges else None
        core_range2 = core_ranges.get(chr2, None) if core_ranges else None

        svg_content_parts.append(draw_all_gene_structures(chr1, chro1, target_set, core_range1, revcomp=False))
        svg_content_parts.append(draw_all_gene_structures(chr2, chro2, target_set, core_range2, revcomp2))

    # -------------------------- Draw gene ID labels --------------------------
    if gene_info:
        def draw_gene_labels_below(chromosome, chro_obj, target_set, revcomp):
            gene_svg = []
            chro_start = chro_obj.start
            chro_end = chro_obj.end
            chromosome_lower = chromosome.lower()
            for gid, ginfo in gene_info.items():
                if ginfo['chro'].lower() != chromosome_lower:
                    continue
                if 'gene' not in ginfo:
                    continue
                g_start, g_end = ginfo['gene']
                if revcomp:
                    mapped_g_start = chro_start + chro_end - g_end
                    mapped_g_end   = chro_start + chro_end - g_start
                    if mapped_g_start > mapped_g_end:
                        mapped_g_start, mapped_g_end = mapped_g_end, mapped_g_start
                else:
                    mapped_g_start = g_start
                    mapped_g_end   = g_end
                if mapped_g_end < chro_start or mapped_g_start > chro_end:
                    continue
                draw_start = max(mapped_g_start, chro_start)
                draw_end = min(mapped_g_end, chro_end)
                if draw_start > draw_end:
                    continue
                x1, y = chro_obj.coordinate(draw_start, is_up=True, is_start=True)
                x2, y = chro_obj.coordinate(draw_end, is_up=True, is_start=False)
                mid_x = (x1 + x2) / 2
                line_bottom_y = y + args.chro_thickness + 5
                label_y = line_bottom_y + 15
                gene_svg.append(
                    f'<line x1="{mid_x}" y1="{y + args.chro_thickness}" x2="{mid_x}" y2="{label_y-5}" class="leader_line"/>'
                )
                if gid in target_set:
                    escaped_gid = xml_escape(gid)
                    gene_svg.append(
                        f'<text x="{mid_x}" y="{label_y}" class="target_gene_label">{escaped_gid}</text>'
                    )
                else:
                    escaped_gid = xml_escape(gid)
                    gene_svg.append(
                        f'<text x="{mid_x}" y="{label_y}" class="feature_label">{escaped_gid}</text>'
                    )
            return ''.join(gene_svg)

        if target_set is None:
            target_set = {gene1, gene2}
        svg_content_parts.append(draw_gene_labels_below(chr1, chro1, target_set, revcomp=False))
        svg_content_parts.append(draw_gene_labels_below(chr2, chro2, target_set, revcomp2))

    # -------------------------- Scale bar --------------------------
    if not args.no_scale:
        L = max_len / 5
        if L < 1:
            L = 1
        else:
            exp = int(math.log10(L))
            mant = L / (10**exp)
            if mant < 2:
                mant = 1
            elif mant < 5:
                mant = 2
            else:
                mant = 5
            L = mant * (10**exp)
        scale_x = args.svg_width * 0.8
        scale_y = args.svg_height * 0.9
        svg_content_parts.append(
            f'<polyline points="{scale_x},{scale_y-5} {scale_x},{scale_y} {scale_x+L*scale},{scale_y} {scale_x+L*scale},{scale_y-5}" fill="none" stroke="black" class="scale"/>'
        )
        svg_content_parts.append(
            f'<text x="{scale_x + L*scale/3}" y="{scale_y-10}" fill="black" font-size="{args.label_font_size}" class="scale-text">{int(L)} bp</text>'
        )

    # -------------------------- Chromosome position labels --------------------------
    pos_label_y1 = chro1.top - args.chro_thickness - 0
    svg_content_parts.append(
        f'<text x="{chro1.left - args.pos_label_x_offset}" y="{pos_label_y1}" fill="black" font-size="12" class="pos_label">{chr1}: {seq1_start:,} - {seq1_end:,}</text>'
    )
    pos_label_y2 = chro2.top - args.chro_thickness - 0
    svg_content_parts.append(
        f'<text x="{chro2.left - args.pos_label_x_offset}" y="{pos_label_y2}" fill="black" font-size="12" class="pos_label">{chr2}: {seq2_start:,} - {seq2_end:,}</text>'
    )

    # -------------------------- Optional axis --------------------------
    if args.chro_axis:
        def add_axis(chro, position):
            axis_left = chro.left
            axis_right = chro.right
            axis_top = chro.top
            axis_start = chro.start
            axis_end = chro.end
            axis_svg = []
            axis_svg.append(
                f'<path d="M {axis_left} {axis_top} L {axis_right} {axis_top}" stroke="black" fill="none"/>'
            )
            step = (axis_end - axis_start) / 10
            for i in range(11):
                pos = axis_start + i * step
                x = axis_left + (pos - axis_start) * scale
                axis_svg.append(
                    f'<path d="M {x} {axis_top} L {x} {axis_top+5}" stroke="black" fill="none"/>'
                )
            return ''.join(axis_svg)

        svg_content_parts.append(add_axis(chro1, 'top'))
        svg_content_parts.append(add_axis(chro2, 'bottom'))

    # -------------------------- Color legend --------------------------
    if legend_min is not None and legend_max is not None:
        legend_width = 25
        legend_height = (top2 + args.chro_thickness - top1)
        legend_x = args.svg_width - 100
        legend_y = top1

        svg_content_parts.append(
            f'<rect x="{legend_x}" y="{legend_y}" width="{legend_width}" height="{legend_height}" fill="url(#grayPurpleRedGradient)" stroke="black" stroke-width="1"/>'
        )
        svg_content_parts.append(
            f'<text x="{legend_x + legend_width + 5}" y="{legend_y + 12}" fill="black" font-size="10" text-anchor="start">{legend_max:.1f}{"%" if color_metric == "identity" else ""}</text>'
        )
        svg_content_parts.append(
            f'<text x="{legend_x + legend_width + 5}" y="{legend_y + legend_height - 5}" fill="black" font-size="10" text-anchor="start">{legend_min:.1f}{"%" if color_metric == "identity" else ""}</text>'
        )
        svg_content_parts.append(
            f'<text x="{legend_x}" y="{legend_y - 10}" fill="black" font-size="11" text-anchor="start">'
            f'{"Bitscore" if color_metric == "bitscore" else "Identity (%)" if color_metric == "identity" else "-log10(E-value)"}</text>'
        )

    # Assemble SVG — physical size = A4 width (595.28pt), viewBox = virtual coords
    phys_w = 595.28
    phys_h = phys_w * args.svg_height / args.svg_width
    svg_content = f'<svg width="{phys_w}pt" height="{phys_h:.1f}pt\" viewBox=\"0 0 {args.svg_width} {args.svg_height}" xmlns="http://www.w3.org/2000/svg" version="1.1">'
    svg_content += ''.join(svg_content_parts)
    svg_content += '</svg>'

    svg_file = f"{out_prefix}_linkview.svg"
    with open(svg_file, 'w') as f:
        f.write(svg_content)

    return svg_file

# -------------------------- Main Function --------------------------
def main():
    global args, scale
    parser = argparse.ArgumentParser(
        description="MicroSynViz: Visualize Pairwise Genomic Microsynteny",
        epilog="For detailed documentation, see https://github.com/YiyongZhao/MicroSynViz"
    )
    # Target: gene IDs or genomic regions
    parser.add_argument("--gene1", help="Gene ID for region 1. Searched in --annos1 files (GFF/GTF/BED); falls back to FASTA header match if not found or no annotations provided")
    parser.add_argument("--gene2", help="Gene ID for region 2. Searched in --annos2 files; falls back to FASTA header match")
    parser.add_argument("--region1", help="Region 1 in SeqID:start-end format. "
                        "For genome FASTA: Chr1:1000-5000; for CDS FASTA: LOC_Os06g50440:1-1000")
    parser.add_argument("--region2", help="Region 2 in SeqID:start-end format. "
                        "For genome FASTA: Chr2:3000-8000; for CDS FASTA: AT1G01010:1-2000")
    parser.add_argument("--extend", type=int, default=3000, help="Bases to extend around genes/regions (default: 3000). Use 0 for CDS/transcript FASTA")

    # Per-region input files
    parser.add_argument("--fa1", default=None, metavar='FASTA',
                        help="FASTA file for region 1. Accepts genome (chromosome-level) or CDS/transcript (gene-level) FASTA. Auto-indexed if .fai missing. Takes priority over legacy flags.")
    parser.add_argument("--g1", default=None, metavar='FASTA',
                        help=argparse.SUPPRESS)  # Legacy alias for --fa1
    parser.add_argument("--g2", default=None, metavar='FASTA',
                        help=argparse.SUPPRESS)  # Legacy alias for --fa2
    parser.add_argument("--fa2", default=None, metavar='FASTA',
                        help="FASTA file for region 2. Accepts genome or CDS/transcript FASTA. Auto-indexed if .fai missing.")
    parser.add_argument("--annos1", nargs='+', default=None, metavar='FILE',
                        help="Annotation file(s) for region 1 (GFF3, GTF, or BED format). Multiple files accepted: first file renders closest to chromosome bar, subsequent files layer outward. Optional for CDS/transcript FASTA.")
    parser.add_argument("--annos2", nargs='+', default=None, metavar='FILE',
                        help="Annotation file(s) for region 2 (GFF3, GTF, or BED format). Multiple files accepted, layered outward by order. Optional for CDS/transcript FASTA.")

    # Legacy aliases (hidden, for backward compatibility)
    parser.add_argument("--gffs1", nargs='+', help=argparse.SUPPRESS)
    parser.add_argument("--gffs2", nargs='+', help=argparse.SUPPRESS)
    parser.add_argument("-g", "--genome", nargs='+', help=argparse.SUPPRESS)
    parser.add_argument("--gff", nargs='+', help=argparse.SUPPRESS)
    parser.add_argument("--te", nargs='*', default=None, help=argparse.SUPPRESS)
    parser.add_argument("--fasta", help=argparse.SUPPRESS)
    parser.add_argument("--te_gff", help=argparse.SUPPRESS)
    parser.add_argument("--gff1", help=argparse.SUPPRESS)
    parser.add_argument("--gff2", help=argparse.SUPPRESS)
    parser.add_argument("--fasta1", help=argparse.SUPPRESS)
    parser.add_argument("--fasta2", help=argparse.SUPPRESS)
    parser.add_argument("--te_gff1", help=argparse.SUPPRESS)
    parser.add_argument("--te_gff2", help=argparse.SUPPRESS)

    # BLAST options
    parser.add_argument("--evalue", type=float, default=1e-5, help="BLAST e-value threshold (default: 1e-5)")
    parser.add_argument("--identity", type=float, default=50, help="Minimum BLAST identity %% to display (default: 50)")
    parser.add_argument("--alignment_length", type=int, default=5, help="Minimum alignment length in bp to display (default: 5)")
    parser.add_argument("--color_by", choices=["bitscore", "identity", "evalue"], default="bitscore",
                        help="Metric for ribbon color gradient: bitscore (default), identity, or evalue")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for BLAST (default: 8)")
    parser.add_argument("--blast_result", metavar="FILE", help="Pre-computed BLAST result in tabular format (outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore). If provided, skips automatic BLAST.")
    parser.add_argument("--auto_complementary", action='store_true',
                        help='Enable automatic reverse complement based on BLAST alignment direction. '
                             'When enabled, if most BLAST hits are reverse complementary, the second sequence '
                             'will be reversed and BLAST will be re-run. Default: off (keep original orientation).')

    # SVG appearance
    parser.add_argument('--svg_height', default=800, type=float, help="SVG canvas height in virtual units (default: 800)")
    parser.add_argument('--svg_width', default=2000, type=float, help="SVG canvas width in virtual units (default: 2000)")
    parser.add_argument('--svg_space', default=0.2, type=float, help="Fraction of width used as left/right margins (default: 0.2)")
    parser.add_argument('--chro_thickness', default=15, type=int, help="Thickness of chromosome rectangles (default: 15)")
    parser.add_argument('--label_font_size', default=18, type=int, help="Font size for scale bar (default: 18)")
    parser.add_argument('--chro_axis', action="store_true", help="Draw tick marks on chromosomes")
    parser.add_argument('--no_scale', action="store_true", help="Omit scale bar")
    parser.add_argument('--gap', type=float, default=0, help="Gap between chromosome and hit polygons (default: 0)")
    parser.add_argument('--pos_label_x_offset', type=float, default=170, help="X offset for position labels (default: 170)")
    parser.add_argument('--ribbon_opacity', default=0.1, type=float, help="Opacity of BLAST homology ribbons. 0=invisible, 1=opaque (default: 0.1)")
    parser.add_argument('--bezier', action="store_true", help="Use Bezier curves for hit polygons")

    # TE track options
    parser.add_argument('--te_offset_base', type=float, default=15, help="Base offset for TE motif labels (default: 15)")
    parser.add_argument('--te_offset_step', type=float, default=15, help="Step increment to avoid overlapping TE labels (default: 15)")
    parser.add_argument('--te_track_height', type=float, default=12, help="Height of TE rectangles (default: 12)")
    parser.add_argument('--te_track_offset', type=float, default=30, help="Distance from chromosome to TE track (default: 30)")

    # Output
    parser.add_argument("--output", default="MicroSynViz_result", help="Output file prefix. Generates {prefix}_linkview.svg, {prefix}_linkview.pdf, {prefix}_genes.fasta, {prefix}_blast.txt (default: MicroSynViz_result)")

    # Verbosity
    parser.add_argument('--quiet', '-q', action='store_true',
                        help="Suppress informational output (only show warnings and errors)")

    # Version
    parser.add_argument('--version', action='version', version=f'MicroSynViz {__version__}')

    args = parser.parse_args()
    scale = 1.0

    # Reset global state (safety for repeated calls / testing)
    global gene_info, mRNA_info, te_info
    gene_info = {}
    mRNA_info = {}
    te_info = []

    # Setup logging
    if args.quiet:
        logging.basicConfig(level=logging.WARNING, format="%(message)s")
    else:
        logging.basicConfig(level=logging.DEBUG, format="%(message)s")

    def log(msg):
        """Log informational messages (respects --quiet)."""
        logger.info(msg.rstrip('\n'))

    # Check external dependencies
    for tool in ['samtools', 'blastn', 'makeblastdb']:
        if not shutil.which(tool):
            raise ExternalToolError(f"'{tool}' not found in PATH. Install via: conda install -c bioconda samtools blast")

    # Validate input mode
    gene_mode = args.gene1 is not None and args.gene2 is not None
    region_mode = args.region1 is not None and args.region2 is not None

    if not gene_mode and not region_mode:
        raise InputError("You must provide either --gene1 and --gene2, or --region1 and --region2.")

    if gene_mode and region_mode:
        sys.stderr.write("[WARNING] Both gene IDs and regions provided. Using gene IDs.\n")
        # Override: use gene mode
        region_mode = False

    # ===================== Resolve input files =====================
    # Primary: --fa1/--fa2 + --annos1/--annos2 (--g1/--g2 are legacy aliases)
    # Legacy fallback: --genome, --gff, --te, --fasta*, --gff*, --te_gff*, --gffs*
    fasta1 = getattr(args, 'fa1', None) or getattr(args, 'g1', None)
    fasta2 = getattr(args, 'fa2', None) or getattr(args, 'g2', None)
    annos1 = list(args.annos1) if getattr(args, 'annos1', None) else []
    annos2 = list(args.annos2) if getattr(args, 'annos2', None) else []

    # Legacy --gffs1/--gffs2
    if not annos1 and getattr(args, 'gffs1', None):
        annos1 = list(args.gffs1)
    if not annos2 and getattr(args, 'gffs2', None):
        annos2 = list(args.gffs2)

    # Legacy --genome/-g, --gff, --te, --fasta*, --gff*, --te_gff*
    if not fasta1 or not fasta2:
        genome_files = getattr(args, 'genome', None) or []
        if not genome_files:
            f1 = getattr(args, 'fasta1', None)
            f2 = getattr(args, 'fasta2', None)
            f0 = getattr(args, 'fasta', None)
            if f1 and f2:
                genome_files = [f1, f2]
            elif f0:
                genome_files = [f0]
        if genome_files:
            fasta1 = fasta1 or genome_files[0]
            fasta2 = fasta2 or (genome_files[1] if len(genome_files) >= 2 else genome_files[0])

    if not annos1 or not annos2:
        gff_files = list(args.gff) if getattr(args, 'gff', None) else []
        te_files = list(args.te) if getattr(args, 'te', None) else []
        if not gff_files:
            g1 = getattr(args, 'gff1', None)
            g2 = getattr(args, 'gff2', None)
            if g1 and g2:
                gff_files = [g1, g2]
            elif g1:
                gff_files = [g1]
        legacy_te = []
        if not te_files:
            for attr_name in ['te_gff', 'te_gff1', 'te_gff2']:
                v = getattr(args, attr_name, None)
                if v and v not in legacy_te:
                    legacy_te.append(v)
        all_te = te_files or legacy_te
        if not annos1 and gff_files:
            annos1 = [gff_files[0]] + (all_te[:1] if all_te else [])
        if not annos2 and gff_files:
            idx = min(1, len(gff_files) - 1)
            te_idx = min(1, len(all_te) - 1) if all_te else -1
            annos2 = [gff_files[idx]] + ([all_te[te_idx]] if te_idx >= 0 else [])

    # Validate
    if not fasta1 or not fasta2:
        raise InputError("Genome FASTA required. Use --fa1/--fa2 (or legacy --g1/--g2/--genome).")
    # annos are optional — CDS FASTA mode works without annotations
    if not annos1:
        annos1 = []
    if not annos2:
        annos2 = []

    check_file_exists(fasta1)
    check_file_exists(fasta2)
    for f in annos1 + annos2:
        check_file_exists(f)
    
    # Auto-index FASTA if .fai missing
    for fa in [fasta1, fasta2]:
        if not os.path.isfile(fa + ".fai"):
            try:
                subprocess.check_call(["samtools", "faidx", fa])
                logger.info(f" Auto-indexed: {fa}")
            except (subprocess.CalledProcessError, FileNotFoundError):
                pass  # Will fail later with a clear error

    dual_species = (fasta1 != fasta2)
    mode_label = "cross-species" if dual_species else "single-species"

    log(f"[Step 1/4] Extracting coordinates ({mode_label})...\n")
    log(f"  Region 1: genome={fasta1}, annotations={annos1}\n")
    log(f"  Region 2: genome={fasta2}, annotations={annos2}\n")

    if gene_mode:
        # gene1 searches annos1 only; gene2 searches annos2 only
        result1 = _find_gene_in_annos(args.gene1, annos1) if annos1 else None
        if not result1:
            result1 = _find_gene_in_fasta_index(args.gene1, fasta1)
            if result1:
                logger.info(f" gene1 '{args.gene1}' found in FASTA index (CDS/sequence mode)")
        if not result1:
            raise GeneNotFoundError(
                f"Gene '{args.gene1}' not found in annos1 {annos1} or FASTA index {fasta1}")
        chr1, start1, end1, strand1 = result1

        result2 = _find_gene_in_annos(args.gene2, annos2) if annos2 else None
        if not result2:
            result2 = _find_gene_in_fasta_index(args.gene2, fasta2)
            if result2:
                logger.info(f" gene2 '{args.gene2}' found in FASTA index (CDS/sequence mode)")
        if not result2:
            raise GeneNotFoundError(
                f"Gene '{args.gene2}' not found in annos2 {annos2} or FASTA index {fasta2}")
        chr2, start2, end2, strand2 = result2

        seq1, seq1_start, seq1_end = extract_sequence_samtools(chr1, start1, end1, args.gene1, fasta1, args.extend)
        seq2, seq2_start, seq2_end = extract_sequence_samtools(chr2, start2, end2, args.gene2, fasta2, args.extend)
    else:  # region mode
        chr1, start1, end1 = parse_region(args.region1)
        chr2, start2, end2 = parse_region(args.region2)
        strand1 = '+'
        strand2 = '+'
        seq1, seq1_start, seq1_end = extract_region_samtools(chr1, start1, end1, fasta1, args.extend)
        seq2, seq2_start, seq2_end = extract_region_samtools(chr2, start2, end2, fasta2, args.extend)

    # Parse annotation files (auto-detect GFF3/GTF/BED) — optional for CDS mode
    if annos1 or annos2:
        log("[Step 2/4] Parsing annotations...\n")
        # Parse annos1 and annos2 separately so source_order is per-region (0 = closest to chro bar)
        parsed_set = set()
        for idx, anno_path in enumerate(annos1):
            if anno_path not in parsed_set:
                parse_annotation(anno_path, source_order=idx)
                parsed_set.add(anno_path)
        for idx, anno_path in enumerate(annos2):
            if anno_path not in parsed_set:
                parse_annotation(anno_path, source_order=idx)
                parsed_set.add(anno_path)
        log(f"  Found {len(gene_info)} genes, {len(te_info)} TE entries.\n")
    else:
        log("[Step 2/4] No annotation files provided — CDS/sequence mode (no gene structures)\n")

    # Common steps
    if gene_mode:
        logger.debug(f" gene1: {args.gene1} on {chr1}:{start1}-{end1} strand={strand1}")
        logger.debug(f" gene2: {args.gene2} on {chr2}:{start2}-{end2} strand={strand2}")
        gene1_display = args.gene1
        gene2_display = args.gene2
    else:
        logger.debug(f" region1: {chr1}:{start1}-{end1}")
        logger.debug(f" region2: {chr2}:{start2}-{end2}")
        # Use region strings as identifiers
        gene1_display = f"{chr1}:{seq1_start}-{seq1_end}"
        gene2_display = f"{chr2}:{seq2_start}-{seq2_end}"

    revcomp2 = False
    seq_fasta = f"{args.output}_genes.fasta"
    blast_out = f"{args.output}_blast.txt"

    # Write initial FASTA (original orientation)
    with open(seq_fasta, "w", encoding="utf-8") as f:
        f.write(seq1 + "\n")
        f.write(seq2 + "\n")

    if args.blast_result:
        # Use pre-computed BLAST result (must be outfmt 6)
        check_file_exists(args.blast_result)
        blast_out = args.blast_result
        log(f"[Step 3a/4] Using pre-computed BLAST result: {blast_out}\n")
    else:
        log("[Step 3a/4] Running initial BLASTn alignment...\n")
        run_blastn(seq_fasta, blast_out, args.evalue, args.threads)
    blast_hits = parse_blast_results(blast_out, min_identity=args.identity, min_length=args.alignment_length)

    # If auto_complementary enabled, decide whether to reverse second sequence
    if args.auto_complementary:
        need_revcomp = False
        if blast_hits:
            # Only consider hits between the two sequences
            relevant_hits = [h for h in blast_hits 
                             if (h['qseqid'] == gene1_display and h['sseqid'] == gene2_display) 
                                or (h['qseqid'] == gene2_display and h['sseqid'] == gene1_display)]
            if relevant_hits:
                forward = 0
                reverse = 0
                for hit in relevant_hits:
                    qdir = hit['qstart'] < hit['qend']
                    sdir = hit['sstart'] < hit['send']
                    if qdir == sdir:
                        forward += 1
                    else:
                        reverse += 1
                if reverse > forward:
                    need_revcomp = True
                    log("[INFO] BLAST hits indicate reverse complement. Reversing second sequence and re-running BLAST.\n")
                else:
                    log("[INFO] BLAST hits indicate same direction. No reversal needed.\n")
            else:
                log("[INFO] No inter-sequence BLAST hits found. Skipping direction adjustment.\n")
        else:
            log("[INFO] No BLAST hits found. Skipping direction adjustment.\n")

        if need_revcomp:
            # Reverse complement second sequence
            seq_record = SeqIO.read(StringIO(seq2), "fasta")
            seq_record.seq = seq_record.seq.reverse_complement()
            seq2 = seq_record.format("fasta")
            # Swap coordinates to maintain left<right
            seq2_start, seq2_end = min(seq2_start, seq2_end), max(seq2_start, seq2_end)  # ensure left < right
            revcomp2 = True
            # Write updated FASTA
            with open(seq_fasta, "w", encoding="utf-8") as f:
                f.write(seq1 + "\n")
                f.write(seq2 + "\n")
            # Re-run BLAST
            log("[Step 3b/4] Re-running BLASTn with reversed sequence...\n")
            run_blastn(seq_fasta, blast_out, args.evalue, args.threads)
            blast_hits = parse_blast_results(blast_out, min_identity=args.identity, min_length=args.alignment_length)   # update
    else:
        log("[INFO] --auto_complementary not specified: keeping original orientation, skipping direction adjustment.\n")

    # Determine target genes and core ranges if regions are provided (only relevant in region mode)
    if region_mode and args.region1 and args.region2:
        target_genes = set()
        chr1_lower = chr1.lower()
        chr2_lower = chr2.lower()
        for gid, ginfo in gene_info.items():
            g_chr = ginfo.get('chro', '').lower()
            if g_chr == chr1_lower:
                g_start, g_end = ginfo.get('gene', (0,0))
                if g_start <= end1 and g_end >= start1:
                    target_genes.add(gid)
            if g_chr == chr2_lower:
                g_start, g_end = ginfo.get('gene', (0,0))
                if g_start <= end2 and g_end >= start2:
                    target_genes.add(gid)
        log(f"[Info] {len(target_genes)} genes found overlapping the specified regions and will be highlighted in red.\n")
        core_ranges = {chr1: (start1, end1), chr2: (start2, end2)}
    else:
        target_genes = None
        core_ranges = None

    seq_lengths = read_sequence_lengths(seq_fasta)
    len1 = seq_lengths[gene1_display]
    len2 = seq_lengths[gene2_display]

    log("[Step 4/4] Generating SVG plot...\n")
    svg_file = generate_svg(
        gene1=gene1_display, gene2=gene2_display,
        len1=len1, len2=len2,
        blast_hits=blast_hits,
        seq1_start=seq1_start, seq1_end=seq1_end,
        seq2_start=seq2_start, seq2_end=seq2_end,
        chr1=chr1, chr2=chr2,
        out_prefix=args.output,
        bezier=args.bezier,
        target_set=target_genes,
        core_ranges=core_ranges,
        revcomp2=revcomp2
    )

    # Convert SVG to PDF using cairosvg
    pdf_file = svg_file.replace('.svg', '.pdf')
    try:
        cairosvg.svg2pdf(url=svg_file, write_to=pdf_file)
        log(f"[INFO] PDF generated: {pdf_file}\n")
    except Exception as e:
        sys.stderr.write(f"[WARNING] Failed to convert SVG to PDF: {e}\n")

    log("\n[✅ Analysis Completed Successfully]\n")
    log(f"1. Gene sequences: {seq_fasta}\n")
    log(f"2. BLAST results: {blast_out}\n")
    log(f"3. Output SVG: {svg_file}\n")
    log(f"4. Output PDF: {pdf_file}\n")

if __name__ == "__main__":
    try:
        main()
    except MicroSynVizError as e:
        sys.stderr.write(f"[ERROR] {e}\n")
        sys.exit(1)
    except KeyboardInterrupt:
        sys.stderr.write("\n[Interrupted]\n")
        sys.exit(130)
    except Exception as e:
        sys.stderr.write(f"[FATAL] Unexpected error: {e}\n")
        sys.exit(2)
