import os
import pytest

import pyfastaq

from viridian import utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "utils")


def test_load_single_seq_fasta():
    expect = pyfastaq.sequences.Fasta("seq", "ACGT")
    infile = os.path.join(data_dir, "load_single_seq_fasta.ok.fa")
    assert expect == utils.load_single_seq_fasta(infile)
    infile = os.path.join(data_dir, "load_single_seq_fasta.bad1.fa")
    with pytest.raises(Exception):
        utils.load_single_seq_fasta(infile)
    infile = os.path.join(data_dir, "load_single_seq_fasta.bad2.fa")
    with pytest.raises(Exception):
        utils.load_single_seq_fasta(infile)


def test_mask_low_coverage():
    outprefix = "tmp.mask_low_coverage"
    expect_debug_files = [f"{outprefix}.{x}" for x in ["fa", "sam", "bam"]]
    for filename in expect_debug_files:
        utils.rm_rf(filename)
    ref_seq = "CGTTAATCCTAGGGCAGTTAAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCATAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGAAATCTGCAACAAGCCCGGT"
    reads_file = os.path.join(data_dir, "mask_low_coverage.reads.fa")
    got_masked = utils.mask_low_coverage(
        ref_seq, reads_file, outprefix, min_depth=4, debug=True
    )
    assert (
        got_masked
        == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCATAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGAAATCTGCAACAAGCCCGGT"
    )
    for filename in expect_debug_files:
        assert os.path.exists(filename)
        os.unlink(filename)
    got_masked = utils.mask_low_coverage(
        ref_seq, reads_file, outprefix, min_depth=1, debug=False
    )
    for filename in expect_debug_files:
        assert not os.path.exists(filename)
    assert (
        got_masked
        == "NNNNNNNNNNNNNNNNGTTAAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCATAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGAAATCTGCAACAAGCCCGGT"
    )
