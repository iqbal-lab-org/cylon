import os
import pytest
from unittest import mock

from cylon import utils, tasks

this_dir = os.path.dirname(os.path.abspath(__file__))
data_root = os.path.join(this_dir, "data")


def test_assemble():
    data_dir = os.path.join(data_root, "assemble")
    options = mock.Mock()
    outdir = "tmp.test_task_assemble"
    utils.rm_rf(outdir)
    options.bam = None
    options.ref_fasta = os.path.join(data_dir, "run_assembly_pipeline.ref.fa")
    options.amplicons_json = os.path.join(
        data_dir, "run_assembly_pipeline.amplicons.json"
    )
    options.outdir = outdir
    options.reads_to_map = os.path.join(data_dir, "run_assembly_pipeline.reads.fa")
    options.reads_per_amp_dir = None
    options.mates_to_map = None
    options.minimap_opts = "-t 1 -x map-ont"
    options.min_mean_coverage = 5
    options.target_coverage = 500
    options.read_end_trim = 1
    options.read_map_tolerance = 20
    options.min_read_length = 200
    options.racon_iterations = 3
    options.min_depth_for_not_N = 1
    options.max_amp_N_proportion = 0.5
    options.contig_map_end_allowance = 20
    options.amplicons_to_fail_file = None
    options.wgs = False
    options.debug = True
    got = tasks.assemble.run(options)
    expect_fa = os.path.join(data_dir, "run_assembly_pipeline.expect.fa")
    expect_seq = utils.load_single_seq_fasta(expect_fa)
    # expected fasta is the fasta used to generate the reads. But the amplicons
    # don't cover the whole genome, so we expect to miss the ends
    assert got == expect_seq[11:989]
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
    assert got == expect_seq[351:989]
    consensus_from_file = utils.load_single_seq_fasta(
        os.path.join(outdir, "consensus.final_assembly.fa")
    )
    assert got == consensus_from_file.seq
    assert os.path.exists(os.path.join(options.outdir, "run_info.json"))
    utils.rm_rf(outdir)
    os.unlink(options.amplicons_to_fail_file)
