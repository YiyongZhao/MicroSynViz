#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GeneViz: A User-Friendly Toolkit for Visualizing Pairwise Genomic Microsynteny across Any Species
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
    Python 3.6+
    BioPython
    pandas
    cairosvg
    samtools (command-line)
    blastn (command-line)

Usage example (gene mode):
    ./genevis.py --gene1 GeneA --gene2 GeneB \\
        --gff species.gff --fasta genome.fasta --te_gff TEs.gff \\
        --auto_complementary --output my_synteny

Usage example (region mode):
    ./genevis.py --region1 "Chr1:1000-2000" --region2 "Chr2:3000-4000" \\
        --gff species.gff --fasta genome.fasta --te_gff TEs.gff \\
        --extend 5000 --output region_synteny

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

__version__ = "1.1.0"  # Updated version with region mode

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
            .feature_label {
                font-size: 9px;
                text-anchor: middle;
                fill: #228B22;
            }
            .target_gene_label {
                font-size: 9px;
                text-anchor: middle;
                fill: #FF0000;
                font-weight: bold;
            }
        </style></defs>
    ''',
}

# -------------------------- Custom Exceptions -------------------------
class ArgumentError(Exception): pass
class FormatError(Exception): pass
class FatalError(Exception): pass

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
    """Parse GFF attribute column into a dictionary."""
    res = {}
    if not attr.endswith(';'):
        attr += ';'
    contentFlag = False
    lastPos = 0
    for i, c in enumerate(attr):
        if c == '"':
            contentFlag = not contentFlag
        elif c == ';' and not contentFlag:
            item = attr[lastPos:i]
            tmpArr = item.split('=')
            key = tmpArr[0]
            val = tmpArr[1].strip(' "') if len(tmpArr) >= 2 else None
            res[key] = val
            lastPos = i + 1
    return res

def parse_gene_gff(gff_file):
    """
    Parse gene GFF file and populate gene_info and mRNA_info.
    Supports both single and multi-exon genes.
    """
    global gene_info, mRNA_info
    gene_type_name_set = {'gene', 'GENE'}
    mRNA_type_name_set = {'mRNA', 'mrna', 'MRNA'}
    exon_type_name_set = {'exon', 'EXON'}

    print(f"[INFO] Parsing gene GFF: {gff_file}")
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
                raise FormatError('gff', gff_file, line)
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
                    print(f"[SAMPLE] Gene: {gene_id} at {seqid}:{start}-{end} strand={strand}")
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

def parse_te_gff(te_gff_file):
    """
    Parse TE GFF file and populate te_info.
    Recognizes common TE feature types and extracts motif name if available.
    """
    global te_info
    te_type_name_set = {'dispersed_repeat', 'transposable_element', 'repeat_region',
                        'LTR_retrotransposon', 'DNA_transposon', 'TE', 'LINE', 'SINE', 'helitron'}
    id_pattern = re.compile(r'ID=([^;]+)')
    motif_pattern = re.compile(r'Target\s+"Motif:([^"\s]+)')

    print(f"[INFO] Parsing TE file: {te_gff_file}")
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
                raise FormatError('te_gff', te_gff_file, line)
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
            te_info.append((file_id, seqid, start, end, motif_name, strand))
            te_count += 1
            if sample_count < 5:
                print(f"[SAMPLE] TE {file_id}: {seqid}:{start}-{end} motif={motif_name} strand={strand}")
                sample_count += 1
    print(f"[INFO] TE parsing complete. TE entries stored: {len(te_info)}")

# -------------------------- Utility Functions -------------------------
def check_file_exists(file_path):
    """Exit with error if file does not exist."""
    if not os.path.isfile(file_path):
        sys.stderr.write(f"[ERROR] File not found: {file_path}\n")
        sys.exit(1)

def extract_gene_coordinates(gene_id, gff_file):
    """Use awk to extract gene coordinates from GFF."""
    check_file_exists(gff_file)
    awk_cmd = (
        f"awk -v id=\"{gene_id}\" '$3==\"gene\" && "
        f"($9 ~ \"ID=\"id\";\" || $9 ~ \"Name=\"id\";\" || "
        f"$9 ~ \"ID=\"id\" \" || $9 ~ \"Name=\"id\" \") "
        f"{{print $1, $4, $5, $7}}' {gff_file}"
    )
    try:
        result = subprocess.check_output(awk_cmd, shell=True, encoding="utf-8").strip()
    except subprocess.CalledProcessError:
        sys.stderr.write(f"[ERROR] Failed to parse GFF file: {gff_file}\n")
        sys.exit(1)
    if not result:
        sys.stderr.write(f"[ERROR] Gene {gene_id} not found in GFF file\n")
        sys.exit(1)
    fields = result.split()
    if len(fields) != 4:
        sys.stderr.write(f"[ERROR] Unexpected output from awk for gene {gene_id}: {result}\n")
        sys.exit(1)
    chr_id, start, end, strand = fields
    return chr_id, int(start), int(end), strand

def extract_sequence_samtools(chr_id, start, end, gene_id, genome_fa, extend=3000):
    """Extract genomic region using samtools faidx (gene mode)."""
    check_file_exists(genome_fa)
    region_start = max(1, start - extend)
    region_end = end + extend
    region = f"{chr_id}:{region_start}-{region_end}"
    cmd = f"samtools faidx {genome_fa} {region}"
    try:
        output = subprocess.check_output(cmd, shell=True, encoding="utf-8")
    except subprocess.CalledProcessError:
        sys.stderr.write(f"[ERROR] samtools failed to extract region: {region}\n")
        sys.exit(1)
    lines = output.strip().split("\n")
    if len(lines) < 2 or len(lines[1].strip()) == 0:
        sys.stderr.write(f"[ERROR] Empty sequence for {gene_id}\n")
        sys.exit(1)
    fixed_fasta = [f">{gene_id}"] + lines[1:]
    seq_str = "\n".join(fixed_fasta)
    return seq_str, region_start, region_end

def extract_region_samtools(chr_id, start, end, genome_fa, extend=0):
    """Extract genomic region using samtools faidx (region mode)."""
    check_file_exists(genome_fa)
    region_start = max(1, start - extend)
    region_end = end + extend
    region = f"{chr_id}:{region_start}-{region_end}"
    cmd = f"samtools faidx {genome_fa} {region}"
    try:
        output = subprocess.check_output(cmd, shell=True, encoding="utf-8")
    except subprocess.CalledProcessError:
        sys.stderr.write(f"[ERROR] samtools failed to extract region: {region}\n")
        sys.exit(1)
    lines = output.strip().split("\n")
    if len(lines) < 2 or len(lines[1].strip()) == 0:
        sys.stderr.write(f"[ERROR] Empty sequence for region {region}\n")
        sys.exit(1)
    # Use region identifier as FASTA header
    header = f">{chr_id}:{region_start}-{region_end}"
    fixed_fasta = [header] + lines[1:]
    seq_str = "\n".join(fixed_fasta)
    return seq_str, region_start, region_end

def run_blastn(seq_fasta, blast_out, evalue=1e-5, threads=8):
    """Run BLASTn with tabular output format 6."""
    check_file_exists(seq_fasta)
    tmp_db = "tmp_blast_db"
    makeblastdb_cmd = f"makeblastdb -in {seq_fasta} -dbtype nucl -out {tmp_db}"
    try:
        subprocess.run(makeblastdb_cmd, shell=True, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        sys.stderr.write("[ERROR] Failed to build BLAST database\n")
        sys.exit(1)
    blast_cmd = (f"blastn -query {seq_fasta} -db {tmp_db} -evalue {evalue} "
                 f"-out {blast_out} -outfmt 6 -num_threads {threads}")
    try:
        subprocess.run(blast_cmd, shell=True, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        sys.stderr.write("[ERROR] BLASTn alignment failed\n")
        sys.exit(1)
    # Cleanup temporary database files
    for f in os.listdir("."):
        if f.startswith(tmp_db):
            os.remove(f)

def read_sequence_lengths(fasta_file):
    """Return dictionary of sequence lengths from FASTA."""
    check_file_exists(fasta_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return {gene_id: len(seq) for gene_id, seq in seq_dict.items()}

def parse_blast_results(blast_file):
    """Parse BLAST tabular output (outfmt 6) into list of dicts."""
    blast_hits = []
    if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
        blast_cols = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        df = pd.read_csv(blast_file, sep="\t", names=blast_cols)
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

    # Gradient definition for bitscore legend
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
        if key not in unique_hits or hit['bitscore'] > unique_hits[key]['bitscore']:
            unique_hits[key] = hit

    relevant_hits = list(unique_hits.values())

    # Bitscore range for legend
    legend_min = legend_max = None
    if relevant_hits:
        scores = [hit['bitscore'] for hit in relevant_hits]
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

            # Normalize bitscore
            if legend_max > legend_min:
                norm = (hit['bitscore'] - legend_min) / (legend_max - legend_min)
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
                    f'<path d="M{v1[0]},{v1[1]} C{bezier_coor1[0]},{bezier_coor1[1]} {bezier_coor2[0]},{bezier_coor2[1]} {v2[0]},{v2[1]} L{v3[0]},{v3[1]} C{bezier_coor3[0]},{bezier_coor3[1]} {bezier_coor4[0]},{bezier_coor4[1]} {v4[0]},{v4[1]} Z" fill="{fill_color}" opacity="0.6"/>'
                )
            else:
                svg_content_parts.append(
                    f'<path d="M{v1[0]},{v1[1]} L{v2[0]},{v2[1]} L{v3[0]},{v3[1]} L{v4[0]},{v4[1]} Z" fill="{fill_color}" opacity="0.6"/>'
                )
    else:
        print("[INFO] No relevant BLAST hits found between the two genes/regions.")

    # -------------------------- Draw TE track --------------------------
    if te_info:
        print(f"[INFO] Starting TE track drawing for {gene1} and {gene2}")
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
            for (file_id, te_chr, te_start, te_end, motif_name, strand) in te_info:
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

                te_candidates.append((x1, x2, motif_name, actual_strand))
                if debug_level >= 1:
                    print(f"[DEBUG] TE {file_id}: original ({te_start},{te_end}) mapped to ({mapped_start},{mapped_end}) -> draw ({draw_start},{draw_end}) strand={actual_strand}")

            # Sort and layout to avoid overlaps
            te_candidates.sort(key=lambda t: t[0])
            used_layers = []
            used_id_offsets = []
            base_id_offset = args.te_offset_base
            id_step = args.te_offset_step

            for (x1, x2, motif_name, strand) in te_candidates:
                y_layer_offset = 0
                while True:
                    conflict = False
                    for (used_offset, ux1, ux2) in used_layers:
                        if max(x1, ux1) <= min(x2, ux2) and used_offset == y_layer_offset:
                            conflict = True
                            break
                    if not conflict:
                        break
                    y_layer_offset += track_h + 2
                    if y_layer_offset > 10 * (track_h + 2):
                        y_layer_offset = 0
                        break
                used_layers.append((y_layer_offset, x1, x2))

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

            print(f"[INFO] TE track for {gene_id}: matched chromosome count = {matched_chr_count}, {len(te_candidates)} TE(s) drawn.")
            return ''.join(te_svg)

        svg_content_parts.append(draw_te_track(gene1, chr1, chro1, 'top', revcomp=False))
        svg_content_parts.append(draw_te_track(gene2, chr2, chro2, 'bottom', revcomp2))
    else:
        print("[INFO] te_info is empty, skipping TE drawing.")

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
                print(f"[DEBUG] All chromosome names in gene_info (lowercase): {sorted(all_chromosomes)}")
                draw_all_gene_structures.chromosomes_printed = True

            genes_on_chr = [gid for gid, ginfo in gene_info.items() if ginfo['chro'].lower() == chromosome_lower]
            print(f"[DEBUG] Genes on {chromosome} (lower={chromosome_lower}): {genes_on_chr}")

            # Adjust core range if reversed
            if revcomp and core_start is not None and core_end is not None:
                core_start, core_end = chro_start + chro_end - core_end, chro_start + chro_end - core_start
                if core_start > core_end:
                    core_start, core_end = core_end, core_start
                print(f"[DEBUG] Mapped core range for revcomp: ({core_start}, {core_end})")

            COLOR_TARGET = "#FF0000"
            COLOR_OTHER  = "#228B22"

            for gid, ginfo in gene_info.items():
                if ginfo['chro'].lower() != chromosome_lower:
                    continue
                if 'gene' not in ginfo:
                    continue
                g_start, g_end = ginfo['gene']
                print(f"[DEBUG] candidate gene {gid} original ({g_start},{g_end})")

                if revcomp:
                    mapped_g_start = chro_start + chro_end - g_end
                    mapped_g_end   = chro_start + chro_end - g_start
                    if mapped_g_start > mapped_g_end:
                        mapped_g_start, mapped_g_end = mapped_g_end, mapped_g_start
                    strand = '+' if ginfo['strand'] == '-' else '-'
                    print(f"[DEBUG] revcomp gene {gid}: mapped to ({mapped_g_start},{mapped_g_end}) strand={strand}")
                else:
                    mapped_g_start = g_start
                    mapped_g_end   = g_end
                    strand = ginfo['strand']

                # Check if gene overlaps display region
                if mapped_g_end < chro_start or mapped_g_start > chro_end:
                    print(f"[DEBUG] gene {gid} out of display region after mapping, skipping")
                    continue
                gene_region_start = max(mapped_g_start, chro_start)
                gene_region_end = min(mapped_g_end, chro_end)

                # Get longest mRNA and its exons
                mrna_ids = list(ginfo.get('mRNA', {}).keys())
                if not mrna_ids:
                    print(f"[DEBUG] gene {gid} has no mRNA info, skipping")
                    continue
                mrna_id = max(mrna_ids, key=lambda mid: mRNA_info[mid]['pos'][1] - mRNA_info[mid]['pos'][0])
                exons = mRNA_info[mrna_id].get('exon', [])
                if not exons:
                    print(f"[DEBUG] gene {gid} has no exon info, skipping")
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
                    print(f"[DEBUG] gene {gid} has no exon in display region after mapping")
                    continue
                mapped_exons.sort()
                print(f"[DEBUG] gene {gid} mapped exons: {mapped_exons}")

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
                            print(f"[DEBUG] segment ({seg_start},{seg_end}) gave invalid coordinates")
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

    # -------------------------- Bitscore legend --------------------------
    if legend_min is not None and legend_max is not None:
        legend_width = 25
        legend_height = (top2 + args.chro_thickness - top1)
        legend_x = args.svg_width - 100
        legend_y = top1

        svg_content_parts.append(
            f'<rect x="{legend_x}" y="{legend_y}" width="{legend_width}" height="{legend_height}" fill="url(#grayPurpleRedGradient)" stroke="black" stroke-width="1"/>'
        )
        svg_content_parts.append(
            f'<text x="{legend_x + legend_width + 5}" y="{legend_y + 12}" fill="black" font-size="10" text-anchor="start">{legend_max:.1f}</text>'
        )
        svg_content_parts.append(
            f'<text x="{legend_x + legend_width + 5}" y="{legend_y + legend_height - 5}" fill="black" font-size="10" text-anchor="start">{legend_min:.1f}</text>'
        )
        svg_content_parts.append(
            f'<text x="{legend_x}" y="{legend_y - 10}" fill="black" font-size="11" text-anchor="start">Bitscore</text>'
        )

    # Assemble SVG
    svg_content = f'<svg width="{args.svg_width}" height="{args.svg_height}" xmlns="http://www.w3.org/2000/svg" version="1.1">'
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
        description="GeneViz: Visualize Pairwise Genomic Microsynteny",
        epilog="For detailed documentation, see https://github.com/yourusername/GeneViz"
    )
    # Two mutually exclusive input modes: gene IDs or genomic regions
    parser.add_argument("--gene1", help="Gene ID for the first gene (required if --region1 not used)")
    parser.add_argument("--gene2", help="Gene ID for the second gene (required if --region2 not used)")
    parser.add_argument("--region1", help="First genomic region (e.g., Chr1:1000-2000); required if --gene1 not used")
    parser.add_argument("--region2", help="Second genomic region; required if --gene2 not used")

    # Region mode (optional, if not provided, uses full extracted region)
    parser.add_argument("--extend", type=int, default=3000, help="Bases to extend around genes/regions (default: 3000)")

    # Single species mode
    parser.add_argument("--gff", help="GFF file for single species mode")
    parser.add_argument("--fasta", help="Genome FASTA for single species mode")
    parser.add_argument("--te_gff", help="TE GFF for single species mode")

    # Dual species mode
    parser.add_argument("--gff1", help="GFF file for species 1")
    parser.add_argument("--gff2", help="GFF file for species 2")
    parser.add_argument("--fasta1", help="Genome FASTA for species 1")
    parser.add_argument("--fasta2", help="Genome FASTA for species 2")
    parser.add_argument("--te_gff1", help="TE GFF for species 1")
    parser.add_argument("--te_gff2", help="TE GFF for species 2")

    # BLAST options
    parser.add_argument("--evalue", type=float, default=1e-5, help="BLAST e-value threshold (default: 1e-5)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for BLAST (default: 8)")
    parser.add_argument("--blast_result", help="Pre-computed BLAST output (outfmt 6); if not given, BLAST is run automatically")
    parser.add_argument("--auto_complementary", action='store_true',
                        help='Enable automatic reverse complement based on BLAST alignment direction. '
                             'When enabled, if most BLAST hits are reverse complementary, the second sequence '
                             'will be reversed and BLAST will be re-run. Default: off (keep original orientation).')

    # SVG appearance
    parser.add_argument('--svg_height', default=800, type=int, help="SVG canvas height (default: 800)")
    parser.add_argument('--svg_width', default=2000, type=int, help="SVG canvas width (default: 2000)")
    parser.add_argument('--svg_space', default=0.2, type=float, help="Fraction of width used as left/right margins (default: 0.2)")
    parser.add_argument('--chro_thickness', default=15, type=int, help="Thickness of chromosome rectangles (default: 15)")
    parser.add_argument('--label_font_size', default=18, type=int, help="Font size for scale bar (default: 18)")
    parser.add_argument('--chro_axis', action="store_true", help="Draw tick marks on chromosomes")
    parser.add_argument('--no_scale', action="store_true", help="Omit scale bar")
    parser.add_argument('--gap', type=float, default=0, help="Gap between chromosome and hit polygons (default: 0)")
    parser.add_argument('--pos_label_x_offset', type=float, default=170, help="X offset for position labels (default: 170)")
    parser.add_argument('--bezier', action="store_true", help="Use Bezier curves for hit polygons")

    # TE track options
    parser.add_argument('--te_offset_base', type=float, default=15, help="Base offset for TE motif labels (default: 15)")
    parser.add_argument('--te_offset_step', type=float, default=15, help="Step increment to avoid overlapping TE labels (default: 15)")
    parser.add_argument('--te_track_height', type=float, default=12, help="Height of TE rectangles (default: 12)")
    parser.add_argument('--te_track_offset', type=float, default=30, help="Distance from chromosome to TE track (default: 30)")

    # Output
    parser.add_argument("--output", default="GeneViz_result", help="Output file prefix (default: GeneViz_result)")

    # Version
    parser.add_argument('--version', action='version', version=f'GeneViz {__version__}')

    args = parser.parse_args()
    scale = 1.0

    # Validate input mode
    gene_mode = args.gene1 is not None and args.gene2 is not None
    region_mode = args.region1 is not None and args.region2 is not None

    if not gene_mode and not region_mode:
        sys.stderr.write("[ERROR] You must provide either --gene1 and --gene2, or --region1 and --region2.\n")
        sys.exit(1)

    if gene_mode and region_mode:
        sys.stderr.write("[WARNING] Both gene IDs and regions provided. Using gene IDs.\n")
        # Override: use gene mode
        region_mode = False

    # Determine single vs dual species mode
    dual_species = (args.fasta1 is not None and args.fasta2 is not None)

    if dual_species:
        if args.fasta is not None or args.gff is not None:
            sys.stderr.write("[WARNING] Both single-species and dual-species arguments provided. Using dual-species mode.\n")
        fasta1 = args.fasta1
        fasta2 = args.fasta2
        gff1 = args.gff1
        gff2 = args.gff2
        te_gff1 = args.te_gff1
        te_gff2 = args.te_gff2

        if gff1 is None or gff2 is None:
            sys.stderr.write("[ERROR] In dual-species mode, --gff1 and --gff2 are required.\n")
            sys.exit(1)

        for f in [fasta1, fasta2, gff1, gff2]:
            check_file_exists(f)

        sys.stdout.write("[Step 1/4] Extracting coordinates (dual species)...\n")
        if gene_mode:
            chr1, start1, end1, strand1 = extract_gene_coordinates(args.gene1, gff1)
            chr2, start2, end2, strand2 = extract_gene_coordinates(args.gene2, gff2)
            seq1, seq1_start, seq1_end = extract_sequence_samtools(chr1, start1, end1, args.gene1, fasta1, args.extend)
            seq2, seq2_start, seq2_end = extract_sequence_samtools(chr2, start2, end2, args.gene2, fasta2, args.extend)
        else:  # region mode
            chr1, start1, end1 = parse_region(args.region1)
            chr2, start2, end2 = parse_region(args.region2)
            strand1 = '+'
            strand2 = '+'
            seq1, seq1_start, seq1_end = extract_region_samtools(chr1, start1, end1, fasta1, args.extend)
            seq2, seq2_start, seq2_end = extract_region_samtools(chr2, start2, end2, fasta2, args.extend)

        # Parse gene GFFs
        sys.stdout.write("[Step 2a/4] Parsing gene GFFs...\n")
        try:
            parse_gene_gff(gff1)
            parse_gene_gff(gff2)
            sys.stdout.write(f"[Info] Found {len(gene_info)} genes in total.\n")
        except Exception as e:
            sys.stderr.write(f"[WARNING] Failed to parse gene GFF: {e}\n")
            sys.stderr.write("Gene structures will not be drawn.\n")

        # Parse TE GFFs
        if te_gff1:
            sys.stdout.write("[Step 2b/4] Parsing TE GFF for species 1...\n")
            try:
                parse_te_gff(te_gff1)
            except Exception as e:
                sys.stderr.write(f"[WARNING] Failed to parse TE GFF1: {e}\n")
        if te_gff2:
            sys.stdout.write("[Step 2b/4] Parsing TE GFF for species 2...\n")
            try:
                parse_te_gff(te_gff2)
            except Exception as e:
                sys.stderr.write(f"[WARNING] Failed to parse TE GFF2: {e}\n")
        if te_gff1 or te_gff2:
            sys.stdout.write(f"[Info] Total TE entries: {len(te_info)}\n")
        else:
            sys.stdout.write("[Step 2b/4] No TE GFF provided.\n")

    else:  # Single species mode
        if args.fasta is None or args.gff is None:
            sys.stderr.write("[ERROR] In single-species mode, --fasta and --gff are required.\n")
            sys.exit(1)
        fasta = args.fasta
        gff = args.gff
        te_gff = args.te_gff

        check_file_exists(fasta)
        check_file_exists(gff)

        sys.stdout.write("[Step 1/4] Extracting coordinates (single species)...\n")
        if gene_mode:
            chr1, start1, end1, strand1 = extract_gene_coordinates(args.gene1, gff)
            chr2, start2, end2, strand2 = extract_gene_coordinates(args.gene2, gff)
            seq1, seq1_start, seq1_end = extract_sequence_samtools(chr1, start1, end1, args.gene1, fasta, args.extend)
            seq2, seq2_start, seq2_end = extract_sequence_samtools(chr2, start2, end2, args.gene2, fasta, args.extend)
        else:  # region mode
            chr1, start1, end1 = parse_region(args.region1)
            chr2, start2, end2 = parse_region(args.region2)
            strand1 = '+'
            strand2 = '+'
            seq1, seq1_start, seq1_end = extract_region_samtools(chr1, start1, end1, fasta, args.extend)
            seq2, seq2_start, seq2_end = extract_region_samtools(chr2, start2, end2, fasta, args.extend)

        sys.stdout.write("[Step 2a/4] Parsing gene GFF...\n")
        try:
            parse_gene_gff(gff)
            sys.stdout.write(f"[Info] Found {len(gene_info)} genes in gene GFF.\n")
        except Exception as e:
            sys.stderr.write(f"[WARNING] Failed to parse gene GFF: {e}\n")
            sys.stderr.write("Gene structures will not be drawn.\n")

        if te_gff:
            sys.stdout.write("[Step 2b/4] Parsing TE GFF...\n")
            try:
                parse_te_gff(te_gff)
                sys.stdout.write(f"[Info] Found {len(te_info)} TE elements in TE GFF.\n")
            except Exception as e:
                sys.stderr.write(f"[WARNING] Failed to parse TE GFF: {e}\n")
                sys.stderr.write("TE annotations will not be drawn.\n")
        else:
            sys.stdout.write("[Step 2b/4] No TE GFF provided.\n")

    # Common steps
    if gene_mode:
        print(f"[DEBUG] gene1: {args.gene1} on {chr1}:{start1}-{end1} strand={strand1}")
        print(f"[DEBUG] gene2: {args.gene2} on {chr2}:{start2}-{end2} strand={strand2}")
        gene1_display = args.gene1
        gene2_display = args.gene2
    else:
        print(f"[DEBUG] region1: {chr1}:{start1}-{end1}")
        print(f"[DEBUG] region2: {chr2}:{start2}-{end2}")
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

    sys.stdout.write("[Step 3a/4] Running initial BLASTn alignment...\n")
    run_blastn(seq_fasta, blast_out, args.evalue, args.threads)
    blast_hits = parse_blast_results(blast_out)

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
                    sys.stdout.write("[INFO] BLAST hits indicate reverse complement. Reversing second sequence and re-running BLAST.\n")
                else:
                    sys.stdout.write("[INFO] BLAST hits indicate same direction. No reversal needed.\n")
            else:
                sys.stdout.write("[INFO] No inter-sequence BLAST hits found. Skipping direction adjustment.\n")
        else:
            sys.stdout.write("[INFO] No BLAST hits found. Skipping direction adjustment.\n")

        if need_revcomp:
            # Reverse complement second sequence
            seq_record = SeqIO.read(StringIO(seq2), "fasta")
            seq_record.seq = seq_record.seq.reverse_complement()
            seq2 = seq_record.format("fasta")
            # Swap coordinates to maintain left<right
            seq2_start, seq2_end = seq2_end, seq2_start
            seq2_start, seq2_end = min(seq2_start, seq2_end), max(seq2_start, seq2_end)
            revcomp2 = True
            # Write updated FASTA
            with open(seq_fasta, "w", encoding="utf-8") as f:
                f.write(seq1 + "\n")
                f.write(seq2 + "\n")
            # Re-run BLAST
            sys.stdout.write("[Step 3b/4] Re-running BLASTn with reversed sequence...\n")
            run_blastn(seq_fasta, blast_out, args.evalue, args.threads)
            blast_hits = parse_blast_results(blast_out)   # update
    else:
        sys.stdout.write("[INFO] --auto_complementary not specified: keeping original orientation, skipping direction adjustment.\n")

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
        sys.stdout.write(f"[Info] {len(target_genes)} genes found overlapping the specified regions and will be highlighted in red.\n")
        core_ranges = {chr1: (start1, end1), chr2: (start2, end2)}
    else:
        target_genes = None
        core_ranges = None

    seq_lengths = read_sequence_lengths(seq_fasta)
    len1 = seq_lengths[gene1_display]
    len2 = seq_lengths[gene2_display]

    sys.stdout.write("[Step 4/4] Generating SVG plot...\n")
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
        sys.stdout.write(f"[INFO] PDF generated: {pdf_file}\n")
    except Exception as e:
        sys.stderr.write(f"[WARNING] Failed to convert SVG to PDF: {e}\n")

    sys.stdout.write("\n[✅ Analysis Completed Successfully]\n")
    sys.stdout.write(f"1. Gene sequences: {seq_fasta}\n")
    sys.stdout.write(f"2. BLAST results: {blast_out}\n")
    sys.stdout.write(f"3. Output SVG: {svg_file}\n")
    sys.stdout.write(f"4. Output PDF: {pdf_file}\n")

if __name__ == "__main__":
    main()
