"""Basic tests for MicroSynViz core module."""
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def test_import():
    import microsynviz
    assert hasattr(microsynviz, '__version__')
    assert microsynviz.__version__ == "1.0.0"


def test_parse_attributes_gff3():
    from microsynviz.core import parse_attributes
    result = parse_attributes("ID=gene01;Name=ABC;Parent=chr1")
    assert result['ID'] == 'gene01'
    assert result['Name'] == 'ABC'
    assert result['Parent'] == 'chr1'


def test_parse_attributes_gtf():
    from microsynviz.core import parse_attributes
    result = parse_attributes('gene_id "ENSG001"; transcript_id "ENST001";')
    assert result['gene_id'] == 'ENSG001'
    assert result.get('ID') is not None


def test_parse_region_valid():
    from microsynviz.core import parse_region
    chrom, start, end = parse_region("Chr1:1000-5000")
    assert chrom == "Chr1"
    assert start == 1000
    assert end == 5000


def test_parse_region_invalid():
    from microsynviz.core import parse_region
    with pytest.raises(ValueError):
        parse_region("bad_format")


def test_parse_region_start_gt_end():
    from microsynviz.core import parse_region
    with pytest.raises(ValueError):
        parse_region("Chr1:5000-1000")


def test_xml_escape():
    from microsynviz.core import xml_escape
    assert xml_escape('A<B>C&D"E') == "A&lt;B&gt;C&amp;D&quot;E"
    assert xml_escape(None) == ""


def test_detect_format_gff(tmp_path):
    from microsynviz.core import detect_format
    gff = tmp_path / "test.gff"
    gff.write_text("chr1\tmaker\tgene\t100\t200\t.\t+\t.\tID=gene01\n")
    assert detect_format(str(gff)) == 'gff'


def test_detect_format_bed(tmp_path):
    from microsynviz.core import detect_format
    bed = tmp_path / "test.bed"
    bed.write_text("chr1\t100\t200\tgene01\t0\t+\n")
    assert detect_format(str(bed)) == 'bed'


def test_format_error():
    from microsynviz.core import FormatError
    err = FormatError('gff', 'test.gff', 'bad_line_here')
    assert 'gff' in str(err)
    assert 'test.gff' in str(err)


def test_check_file_exists_raises():
    """check_file_exists raises InputError for missing file."""
    from microsynviz.core import check_file_exists, InputError
    with pytest.raises(InputError, match="File not found"):
        check_file_exists("/nonexistent/path/file.txt")


def test_exception_hierarchy():
    """All custom exceptions inherit from MicroSynVizError."""
    from microsynviz.core import (MicroSynVizError, FormatError,
                                   InputError, ExternalToolError, GeneNotFoundError)
    assert issubclass(FormatError, MicroSynVizError)
    assert issubclass(InputError, MicroSynVizError)
    assert issubclass(ExternalToolError, MicroSynVizError)
    assert issubclass(GeneNotFoundError, MicroSynVizError)


def test_find_gene_in_annos_not_found(tmp_path):
    """_find_gene_in_annos returns None when gene not found."""
    from microsynviz.core import _find_gene_in_annos
    gff = tmp_path / "empty.gff"
    gff.write_text("##gff-version 3\n")
    result = _find_gene_in_annos("NONEXISTENT", [str(gff)])
    assert result is None


def test_find_gene_in_annos_found(tmp_path):
    """_find_gene_in_annos finds gene across multiple files."""
    from microsynviz.core import _find_gene_in_annos
    gff1 = tmp_path / "a.gff"
    gff1.write_text("chr1\tmaker\tgene\t100\t500\t.\t+\t.\tID=GeneA;Name=GeneA\n")
    gff2 = tmp_path / "b.gff"
    gff2.write_text("chr2\tmaker\tgene\t200\t800\t.\t-\t.\tID=GeneB;Name=GeneB\n")
    
    r1 = _find_gene_in_annos("GeneA", [str(gff1), str(gff2)])
    assert r1 == ("chr1", 100, 500, "+")
    
    r2 = _find_gene_in_annos("GeneB", [str(gff1), str(gff2)])
    assert r2 == ("chr2", 200, 800, "-")


def test_parse_gene_gff_with_cds(tmp_path):
    """CDS features are recognized as exon-like."""
    from microsynviz.core import parse_gene_gff, gene_info, mRNA_info
    # Reset globals
    gene_info.clear()
    mRNA_info.clear()
    
    gff = tmp_path / "cds.gff"
    gff.write_text(
        "chr1\tmaker\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "chr1\tmaker\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tmaker\tCDS\t100\t200\t.\t+\t.\tID=cds1;Parent=t1\n"
        "chr1\tmaker\tCDS\t300\t500\t.\t+\t.\tID=cds2;Parent=t1\n"
    )
    parse_gene_gff(str(gff))
    assert 'g1' in gene_info
    assert 't1' in mRNA_info
    assert len(mRNA_info['t1'].get('exon', [])) == 2


def test_parse_te_gff(tmp_path):
    """TE GFF parsing recognizes transposable_element features."""
    from microsynviz.core import parse_te_gff, te_info
    te_info.clear()
    gff = tmp_path / "te.gff"
    gff.write_text(
        "chr1\tRepeatMasker\ttransposable_element\t500\t800\t.\t+\t.\tID=TE001;Target \"Motif:hAT\" 1 300\n"
        "chr1\tRepeatMasker\ttransposable_element\t1000\t1200\t.\t-\t.\tID=TE002\n"
    )
    parse_te_gff(str(gff))
    assert len(te_info) == 2
    assert te_info[0][4] == "hAT"  # motif name extracted
    assert te_info[1][4] == "unknown"  # no motif → unknown


def test_parse_bed_genes(tmp_path):
    """BED12 gene parsing extracts exon blocks."""
    from microsynviz.core import parse_bed_genes, gene_info, mRNA_info
    gene_info.clear()
    mRNA_info.clear()
    bed = tmp_path / "genes.bed"
    # BED12: 2 blocks (exons)
    bed.write_text("chr1\t99\t500\tMyGene\t0\t+\t99\t500\t0\t2\t100,150\t0,251\n")
    parse_bed_genes(str(bed))
    assert 'MyGene' in gene_info
    assert gene_info['MyGene']['gene'] == (100, 500)  # 0-based → 1-based
    mrna = mRNA_info.get('MyGene.1')
    assert mrna is not None
    assert len(mrna['exon']) == 2


def test_parse_annotation_gff(tmp_path):
    """parse_annotation auto-detects GFF and parses genes + TEs."""
    from microsynviz.core import parse_annotation, gene_info, mRNA_info, te_info
    gene_info.clear()
    mRNA_info.clear()
    te_info.clear()
    gff = tmp_path / "mixed.gff"
    gff.write_text(
        "chr1\tmaker\tgene\t100\t500\t.\t+\t.\tID=g1\n"
        "chr1\tmaker\tmRNA\t100\t500\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\tmaker\texon\t100\t300\t.\t+\t.\tID=e1;Parent=t1\n"
        "chr1\tRepMask\ttransposable_element\t600\t900\t.\t-\t.\tID=te1\n"
    )
    parse_annotation(str(gff))
    assert 'g1' in gene_info
    assert len(te_info) >= 1


def test_quiet_logging(tmp_path):
    """Verify logger exists and is named correctly."""
    from microsynviz.core import logger
    assert logger.name == "MicroSynViz"
