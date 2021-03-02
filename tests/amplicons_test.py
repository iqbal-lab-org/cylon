import collections
import filecmp
import os
import pysam
import pytest

from viridian import amplicons, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "amplicons")


def test_masked_seq_centre_coord():
    amplicon = amplicons.Amplicon("name", 0, 10)
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
    amplicon = amplicons.Amplicon("name", 0, 9)
    assert amplicon.ref_centre_coord() == 5
    amplicon = amplicons.Amplicon("name", 0, 10)
    assert amplicon.ref_centre_coord() == 5
    amplicon = amplicons.Amplicon("name", 0, 11)
    assert amplicon.ref_centre_coord() == 6
    amplicon = amplicons.Amplicon("name", 0, 12)
    assert amplicon.ref_centre_coord() == 6
    amplicon = amplicons.Amplicon("name", 0, 13)
    assert amplicon.ref_centre_coord() == 7
    amplicon = amplicons.Amplicon("name", 10, 19)
    assert amplicon.ref_centre_coord() == 15
    amplicon = amplicons.Amplicon("name", 10, 20)
    assert amplicon.ref_centre_coord() == 15
    amplicon = amplicons.Amplicon("name", 10, 21)
    assert amplicon.ref_centre_coord() == 16
    amplicon = amplicons.Amplicon("name", 10, 22)
    assert amplicon.ref_centre_coord() == 16
    amplicon = amplicons.Amplicon("name", 10, 23)
    assert amplicon.ref_centre_coord() == 17


def test_expected_overlap_length():
    amplicon1 = amplicons.Amplicon("name", 0, 10)
    amplicon2 = amplicons.Amplicon("name", 8, 20)
    amplicon3 = amplicons.Amplicon("name", 11, 20)
    assert amplicon1.expected_overlap_length(amplicon2) == 3
    assert amplicon1.expected_overlap_length(amplicon3) is None


def test_get_reads_for_polishing():
    reads_bam = os.path.join(data_dir, "get_reads_for_polishing.bam")
    bam = pysam.AlignmentFile(reads_bam, "rb")
    reads_out = "tmp.get_reads_for_polishing.reads.fa"
    utils.rm_rf(reads_out)
    amplicon = amplicons.Amplicon("amp1", 59, 419)

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

    amplicon = amplicons.Amplicon("amp1", 50, 100)
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
    amplicon = amplicons.Amplicon("amplicon1", 60, 259)
    reads_bam = os.path.join(data_dir, "polish.bam")
    bam = pysam.AlignmentFile(reads_bam, "rb")
    outdir = "tmp.polish.out"
    utils.rm_rf(outdir)
    amplicon.polish(
        ref_genome,
        bam,
        outdir,
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

    amplicon = amplicons.Amplicon("amplicon1", 60, 259)
    amplicon.polish(
        ref_genome,
        bam,
        outdir,
        min_mean_coverage=3,
        racon_iterations=3,
        min_depth_for_not_N=3,
        min_read_length=100,
        max_polished_N_prop=0.1,
        debug=True,
    )
    amplicon.masked_seq is None
    assert not amplicon.assemble_success
    assert not amplicon.polish_data["Polish success"]
    utils.rm_rf(outdir)


def test_masked_overlap():
    Match = collections.namedtuple("Match", ("a", "b", "size"))
    amp1 = amplicons.Amplicon("a1", 50, 100)
    amp2 = amplicons.Amplicon("a2", 90, 150)
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


def test_load_amplicons_bed_file():
    expect = [
        amplicons.Amplicon("name1", 42, 99),
        amplicons.Amplicon("name2", 85, 150),
    ]
    infile = os.path.join(data_dir, "load_amplicons_bed_file.bed")
    assert expect == amplicons.load_amplicons_bed_file(infile)


def test_amplicons_from_fasta_and_bed():
    fasta = os.path.join(data_dir, "amplicons_from_fasta_and_bed.fa")
    bed = os.path.join(data_dir, "amplicons_from_fasta_and_bed.bed")
    got = amplicons.load_amplicons_from_fasta_and_bed(fasta, bed)
    expect = [
        amplicons.Amplicon("name1", 42, 99),
        amplicons.Amplicon("name2", 85, 150),
        amplicons.Amplicon("name3", 140, 249),
    ]
    expect[0].masked_seq = "ACGT"
    expect[2].masked_seq = "GGG"
    expect[0].assemble_success = True
    expect[2].assemble_success = True
    assert got == expect
