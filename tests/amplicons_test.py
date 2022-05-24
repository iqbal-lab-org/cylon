import collections
import filecmp
import os
import pysam
import pytest
from unittest import mock

from cylon import amplicons, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "amplicons")


def test_masked_seq_centre_coord():
    amplicon = amplicons.Amplicon("name", 0, 10, 1, 1)
    assert amplicon.masked_seq_centre_coord() is None
    amplicon.masked_seq = "ACT"
    assert amplicon.masked_seq_centre_coord() == 1
    amplicon.masked_seq = "ACTG"
    assert amplicon.masked_seq_centre_coord() == 2
    amplicon.masked_seq = "ACTGA"
    assert amplicon.masked_seq_centre_coord() == 2
    amplicon.masked_seq = "ACTGAT"
    assert amplicon.masked_seq_centre_coord() == 3


def test_ref_centre_coord():
    amplicon = amplicons.Amplicon("name", 0, 9, 1, 1)
    assert amplicon.ref_centre_coord() == 5
    amplicon = amplicons.Amplicon("name", 0, 10, 1, 1)
    assert amplicon.ref_centre_coord() == 5
    amplicon = amplicons.Amplicon("name", 0, 11, 1, 1)
    assert amplicon.ref_centre_coord() == 6
    amplicon = amplicons.Amplicon("name", 0, 12, 1, 1)
    assert amplicon.ref_centre_coord() == 6
    amplicon = amplicons.Amplicon("name", 0, 13, 1, 1)
    assert amplicon.ref_centre_coord() == 7
    amplicon = amplicons.Amplicon("name", 10, 19, 1, 1)
    assert amplicon.ref_centre_coord() == 15
    amplicon = amplicons.Amplicon("name", 10, 20, 1, 1)
    assert amplicon.ref_centre_coord() == 15
    amplicon = amplicons.Amplicon("name", 10, 21, 1, 1)
    assert amplicon.ref_centre_coord() == 16
    amplicon = amplicons.Amplicon("name", 10, 22, 1, 1)
    assert amplicon.ref_centre_coord() == 16
    amplicon = amplicons.Amplicon("name", 10, 23, 1, 1)
    assert amplicon.ref_centre_coord() == 17


def test_expected_overlap_length():
    amplicon1 = amplicons.Amplicon("name", 0, 10, 1, 1)
    amplicon2 = amplicons.Amplicon("name", 8, 20, 1, 1)
    amplicon3 = amplicons.Amplicon("name", 11, 20, 1, 1)
    assert amplicon1.expected_overlap_length(amplicon2) == 3
    assert amplicon1.expected_overlap_length(amplicon3) is None


def test_use_read_for_polishing():
    amplicon = amplicons.Amplicon("name", 50, 100, 1, 1)
    read = mock.Mock()
    # Test start of read is within X bp of start of amplicon
    read.reference_start = 48
    read.reference_end = 75
    assert amplicon.use_read_for_polishing(read, 2, None, wgs=False)
    assert not amplicon.use_read_for_polishing(read, 1, None, wgs=False)

    # Test end of read is within X bp of end of amplicon
    read.reference_start = 75
    read.reference_end = 103
    assert amplicon.use_read_for_polishing(read, 3, None, wgs=False)
    assert not amplicon.use_read_for_polishing(read, 2, None, wgs=False)

    # Test overlapping read when wgs is True
    read.reference_start = 1
    read.reference_end = 49
    assert not amplicon.use_read_for_polishing(read, None, 10, wgs=True)
    read.reference_end = 60
    assert amplicon.use_read_for_polishing(read, None, 11, wgs=True)
    assert not amplicon.use_read_for_polishing(read, None, 12, wgs=True)
    read.reference_start = 60
    read.reference_end = 69
    assert amplicon.use_read_for_polishing(read, None, 10, wgs=True)
    assert not amplicon.use_read_for_polishing(read, None, 11, wgs=True)
    read.reference_start = 20
    read.reference_end = 110
    assert amplicon.use_read_for_polishing(read, None, 51, wgs=True)
    assert not amplicon.use_read_for_polishing(read, None, 52, wgs=True)


def test_get_reads_for_polishing():
    reads_bam = os.path.join(data_dir, "get_reads_for_polishing.bam")
    bam = pysam.AlignmentFile(reads_bam, "rb")
    reads_out = "tmp.get_reads_for_polishing.reads.fa"
    utils.rm_rf(reads_out)
    amplicon = amplicons.Amplicon("amp1", 59, 419, 1, 1)

    got_reads, got_used, got_cov = amplicon.get_reads_for_polishing(
        "ref1",
        bam,
        reads_out,
        min_coverage=1,
        trim_ends=5,
        tolerance=1,
        min_output_length=300,
        target_depth=3,
    )
    assert got_reads == 6
    assert got_used == 4
    assert got_cov == pytest.approx(4.49, 4.50)
    expect_reads = os.path.join(data_dir, "get_reads_for_polishing.expect.fa")
    assert filecmp.cmp(reads_out, expect_reads, shallow=False)
    os.unlink(reads_out)

    amplicon = amplicons.Amplicon("amp1", 50, 100, 1, 1)
    got_reads, got_used, got_cov = amplicon.get_reads_for_polishing(
        "ref2",
        bam,
        reads_out,
        min_coverage=1,
        trim_ends=5,
        tolerance=1,
        min_output_length=30,
        target_depth=1,
    )
    assert got_reads == 0
    assert got_used == 0
    assert got_cov == 0
    assert not os.path.exists(reads_out)


def test_polish():
    ref_fasta = os.path.join(data_dir, "polish.ref.fa")
    ref_genome = utils.load_single_seq_fasta(ref_fasta)
    amplicon = amplicons.Amplicon("amplicon1", 60, 259, 1, 1)
    reads_bam = os.path.join(data_dir, "polish.bam")
    bam = pysam.AlignmentFile(reads_bam, "rb")
    outdir = "tmp.polish.out"
    utils.rm_rf(outdir)
    amplicon.polish(
        ref_genome,
        outdir,
        bam_to_slice_reads=bam,
        min_mean_coverage=3,
        racon_iterations=3,
        min_depth_for_not_N=3,
        min_read_length=100,
        max_polished_N_prop=0.5,
        debug=True,
    )
    assert (
        amplicon.masked_seq
        == "NNNNNNNNNNNNNNNNNNNNAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCATAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGNNNNNNNNNNNNNNNNNNNNN"
    )
    assert amplicon.assemble_success
    assert amplicon.polish_data["Polish success"]
    utils.rm_rf(outdir)

    # Same again, but this time use the fasta of reads instead of the BAM file.
    # Plus, this is giving untrimmed reads, so we get less masking. In the
    # previous run 20bp trimmed off all the reads
    reads_file = os.path.join(data_dir, "polish.reads.fa")
    amplicon.polish(
        ref_genome,
        outdir,
        reads_file=reads_file,
        min_mean_coverage=3,
        racon_iterations=3,
        min_depth_for_not_N=3,
        min_read_length=100,
        max_polished_N_prop=0.5,
        debug=True,
    )
    assert (
        amplicon.masked_seq
        == "CGTTAATCCTAGGGCAGTTAAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCATAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGAAATCTGCAACAAGCCCGGT"
    )
    assert amplicon.assemble_success
    assert amplicon.polish_data["Polish success"]
    utils.rm_rf(outdir)

    # The reads are such that there's a dip in coverage in the middle of the
    # amplicon. Setting min_depth_for_not_N higher makes this region get
    # masked, and then the amplicon should get failed
    amplicon = amplicons.Amplicon("amplicon1", 60, 259, 1, 1)
    amplicon.polish(
        ref_genome,
        outdir,
        bam_to_slice_reads=bam,
        min_mean_coverage=3,
        racon_iterations=3,
        min_depth_for_not_N=18,
        min_read_length=50,
        max_polished_N_prop=0.1,
        debug=True,
    )
    assert (
        amplicon.masked_seq
        == "NNNNNNNNNNNNNNNNNNNNAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGNNNNNNNNNNNNNNNNNNNNN"
    )
    assert not amplicon.assemble_success
    assert not amplicon.polish_data["Polish success"]
    utils.rm_rf(outdir)


def test_masked_overlap():
    Match = collections.namedtuple("Match", ("a", "b", "size"))
    amp1 = amplicons.Amplicon("a1", 50, 100, 1, 1)
    amp2 = amplicons.Amplicon("a2", 90, 150, 1, 1)
    assert amp1.masked_overlap(amp2, 10) is None
    amp2.masked_seq = "AAAAAAAAAAAAAAAAAATGCTGAACAGTCCCCCCC"
    assert amp1.masked_overlap(amp2, 10) is None
    amp1.masked_seq = "AAAAAAAAAAAAAAAAAATGCTGAACAGTCCCCCCC"
    amp2.masked_seq = None
    assert amp1.masked_overlap(amp2, 10) is None
    amp2.masked_seq = "CCTGCTGAACGGTTGATGCATCTCATGCTGACNNAGGTGTGGCCAAAAA"
    assert amp1.masked_overlap(amp2, 7) == Match(18, 2, 8)
    assert amp1.masked_overlap(amp2, 8) == Match(18, 2, 8)
    assert amp1.masked_overlap(amp2, 9) == None

    amp1.masked_seq = "CCTGCTGAACGGTTGATGCATCTCATGCTGACNNAGGTGTGGCCAAAAA"
    amp2.masked_seq = "NNAGGTGTGGCCTTTTTTTTTTTTTTTTTTTTTTTTTT"
    assert amp1.masked_overlap(amp2, 10) == Match(34, 2, 10)
    assert amp1.masked_overlap(amp2, 11) == None

    amp1.masked_seq = "NNAGGTGTGGCCTTTTTTTTTTTTTTTTTTTTTTTTTT"
    amp2.masked_seq = "AAAAAAAAAAAAAAAAANNNNNNNNNNN"
    assert amp1.masked_overlap(amp2, 0) == Match(2, 0, 1)
    assert amp1.masked_overlap(amp2, 2) == None


def test_load_amplicons_json_file():
    expect = [
        amplicons.Amplicon("name1", 42, 99, 5, 3),
        amplicons.Amplicon("name2", 90, 150, 2, 9),
    ]
    infile = os.path.join(data_dir, "load_amplicons_json_file.json")
    assert expect == amplicons.load_amplicons_json_file(infile)
