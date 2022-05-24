import os
import pytest

from cylon import racon, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "racon")


def test_run_racon():
    # This has a SNP and two indels to fix. Also one position where
    # about 2/3 of the reads say A and the rest say T. Expect this to get
    # corrected to T.
    fa_to_polish = os.path.join(data_dir, "run_racon.to_polish.fa")
    seq_to_polish = utils.load_single_seq_fasta(fa_to_polish)
    reads = os.path.join(data_dir, "run_racon.reads.fa")
    pre_out = "tmp.run_racon"
    utils.rm_rf(f"{pre_out}.sam")
    utils.rm_rf(f"{pre_out}.to_polish.fa")
    polished1 = racon.run_racon(seq_to_polish, reads, pre_out, debug=True)
    assert polished1 != fa_to_polish
    assert (
        polished1
        == "CGTTAATCCTAGGGCAGTTAAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCATAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGAAATCTGCAACAAGCCCGGT"
    )
    # we used debug mode, so intermediate files should be left on disk
    assert os.path.exists(f"{pre_out}.sam")
    assert os.path.exists(f"{pre_out}.to_polish.fa")
    os.unlink(f"{pre_out}.sam")
    os.unlink(f"{pre_out}.to_polish.fa")
    # Another round of polishing shouldn't do anything
    polished2 = racon.run_racon(polished1, reads, pre_out, debug=False)
    assert polished1 == polished2
    # we didn't use debug mode so intermediate files should be deleted
    assert not os.path.exists(f"{pre_out}.sam")
    assert not os.path.exists(f"{pre_out}.to_polish.fa")


def test_run_racon_bad_data():
    fa_to_polish = os.path.join(data_dir, "run_racon.to_polish.fa")
    seq_to_polish = utils.load_single_seq_fasta(fa_to_polish)
    reads = os.path.join(data_dir, "run_racon_bad_reads.fa")
    pre_out = "tmp.run_racon"
    utils.rm_rf(f"{pre_out}.sam")
    utils.rm_rf(f"{pre_out}.to_polish.fa")
    polished = racon.run_racon(seq_to_polish, reads, pre_out, debug=True)
    utils.rm_rf(f"{pre_out}.sam")
    utils.rm_rf(f"{pre_out}.to_polish.fa")
    assert polished is None


def test_run_racon_iterations():
    # A bit hard to come with small artificial test data for this one.
    # We'll just use the same data as for test_run_racon. Should stop after
    # 2 iterations because only the first run corrects anything
    fa_to_polish = os.path.join(data_dir, "run_racon.to_polish.fa")
    seq_to_polish = utils.load_single_seq_fasta(fa_to_polish)
    reads = os.path.join(data_dir, "run_racon.reads.fa")
    outdir = "tmp.run_racon_iterations"
    utils.rm_rf(outdir)
    got_polished = racon.run_racon_iterations(
        seq_to_polish, reads, outdir, max_iterations=3, debug=True
    )
    for i in range(2):
        outprefix = os.path.join(outdir, f"racon.{i}")
        assert os.path.exists(f"{outprefix}.sam")
        assert os.path.exists(f"{outprefix}.polished.fa")
        assert os.path.exists(f"{outprefix}.to_polish.fa")
    assert len(os.listdir(outdir)) == 6
    assert (
        got_polished
        == "CGTTAATCCTAGGGCAGTTAAAAGCCCCATTTTGTACAGCTTTTTCTAGAACAGTCAGGGCGCGCTCCCAGGAGTTGCTTCGCTTCCAGCTAGAAATGATCATCGAACCTGGGTAAGGGCATAATACGAGAATGCTGCCCTATTGCCAGTGCTTAGAAATGGACTGGTGTTACGTCCACGAAATCTGCAACAAGCCCGGT"
    )
    utils.rm_rf(outdir)


def test_run_racon_iterations_bad_data():
    fa_to_polish = os.path.join(data_dir, "run_racon.to_polish.fa")
    seq_to_polish = utils.load_single_seq_fasta(fa_to_polish)
    reads = os.path.join(data_dir, "run_racon_bad_reads.fa")
    outdir = "tmp.run_racon_iterations"
    utils.rm_rf(outdir)
    got_polished = racon.run_racon_iterations(
        seq_to_polish, reads, outdir, max_iterations=3, debug=True
    )
    assert got_polished is None
    utils.rm_rf(outdir)
