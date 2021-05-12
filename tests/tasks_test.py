import os
import pytest
from unittest import mock

from viridian import utils, tasks

this_dir = os.path.dirname(os.path.abspath(__file__))
data_root = os.path.join(this_dir, "data")


def test_assemble():
    data_dir = os.path.join(data_root, "assemble")
    options = mock.Mock()
    outdir = "tmp.test_task_assemble"
    utils.rm_rf(outdir)
    options.bam = None
    options.ref_fasta = os.path.join(data_dir, "run_assembly_pipeline.ref.fa")
    options.amplicons_bed = os.path.join(
        data_dir, "run_assembly_pipeline.amplicons.bed"
    )
    options.outdir = outdir
    options.reads_to_map = os.path.join(data_dir, "run_assembly_pipeline.reads.fa")
    options.mates_to_map = None
    options.minimap_opts = "-t 1 -x map-ont"
    options.min_mean_coverage = 5
    options.target_coverage = 500
    options.read_end_trim = 1
    options.read_map_tolerance = 20
    options.min_read_length = 200
    options.racon_iterations = 3
    options.min_depth_for_not_N = 1
    options.min_amp_overlap_len = 20
    options.contig_map_end_allowance = 20
    options.amplicons_to_fail_file = None
    options.wgs = False
    options.debug = True
    got = tasks.assemble.run(options)
    expect_fa = os.path.join(data_dir, "run_assembly_pipeline.expect.fa")
    expect_seq = utils.load_single_seq_fasta(expect_fa)
    # expected fasta is the fasta used to generate the reads. But the amplicons
    # don't cover the whole genome, so we expect to miss the ends
    assert got == expect_seq[11:-10]
    consensus_from_file = utils.load_single_seq_fasta(
        os.path.join(outdir, "consensus.final_assembly.fa")
    )
    assert got == consensus_from_file.seq
    assert os.path.exists(os.path.join(options.outdir, "run_info.json"))
    utils.rm_rf(outdir)

    # Test the option amplicons_to_fail_file
    options.amplicons_to_fail_file = "tmp.amplicons_to_fail.txt"
    with open(options.amplicons_to_fail_file, "w") as f:
        print("a1", file=f)
    got = tasks.assemble.run(options)
    expect_fa = os.path.join(data_dir, "run_assembly_pipeline.expect.fa")
    expect_seq = utils.load_single_seq_fasta(expect_fa)
    assert got == expect_seq[351:-10]
    consensus_from_file = utils.load_single_seq_fasta(
        os.path.join(outdir, "consensus.final_assembly.fa")
    )
    assert got == consensus_from_file.seq
    assert os.path.exists(os.path.join(options.outdir, "run_info.json"))
    utils.rm_rf(outdir)
    os.unlink(options.amplicons_to_fail_file)


def test_amplicon_overlap():
    options = mock.Mock()
    data_dir = os.path.join(data_root, "tasks")
    options.ref_fasta = os.path.join(data_dir, "amplicon_overlap.ref.fa")
    options.amplicons_fasta = os.path.join(data_dir, "amplicon_overlap.amplicons.fa")
    options.amplicons_bed = os.path.join(data_dir, "amplicon_overlap.bed")
    options.outdir = "tmp.task_amplicon_overlap"
    utils.rm_rf(options.outdir)
    options.min_amp_overlap_len = 20
    options.contig_map_end_allowance = 20
    got = tasks.amplicon_overlap.run(options)
    expect = utils.load_single_seq_fasta(options.ref_fasta)
    assert got == expect.seq
    assert os.path.exists(os.path.join(options.outdir, "out.final_assembly.fa"))
    utils.rm_rf(options.outdir)
