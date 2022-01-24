import json
import os
import pytest

from viridian import assemble, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "assemble")


def test_run_assembly_pipeline():
    ref_fa = os.path.join(data_dir, "run_assembly_pipeline.ref.fa")
    reads_fa = os.path.join(data_dir, "run_assembly_pipeline.reads.fa")
    amplicon_json = os.path.join(data_dir, "run_assembly_pipeline.amplicons.json")
    outdir = "tmp.run_assembly_pipeline"
    utils.rm_rf(outdir)
    got = assemble.run_assembly_pipeline(
        ref_fa,
        amplicon_json,
        outdir,
        reads_fastaq=reads_fa,
        debug=True,
        min_mean_coverage=5,
        min_depth_for_not_N=1,
        read_end_trim=1,
    )

    expect_fa = os.path.join(data_dir, "run_assembly_pipeline.expect.fa")
    expect_seq = utils.load_single_seq_fasta(expect_fa)
    # expected fasta is the fasta used to generate the reads. But the amplicons
    # don't cover the whole genome, so we expect to miss the ends
    assert got == expect_seq[11:979]
    consensus_from_file = utils.load_single_seq_fasta(
        os.path.join(outdir, "consensus.final_assembly.fa")
    )
    assert got == consensus_from_file.seq
    utils.rm_rf(outdir)

    # Rerun, but test force failing the first amplicon
    got = assemble.run_assembly_pipeline(
        ref_fa,
        amplicon_json,
        outdir,
        reads_fastaq=reads_fa,
        debug=True,
        min_mean_coverage=5,
        min_depth_for_not_N=1,
        read_end_trim=1,
        amplicons_to_fail={"a1"},
    )

    expect_fa = os.path.join(data_dir, "run_assembly_pipeline.expect.fa")
    expect_seq = utils.load_single_seq_fasta(expect_fa)
    # This time, we should not have the first amplicon, and the returned
    # sequence should start with the second amplicon
    assert got == expect_seq[356:979]
    consensus_from_file = utils.load_single_seq_fasta(
        os.path.join(outdir, "consensus.final_assembly.fa")
    )
    assert got == consensus_from_file.seq

    # some checks of the contents of the json summary
    with open(os.path.join(outdir, "run_info.json")) as f:
        run_info = json.load(f)
    assert run_info["run_summary"]["made_consensus"] is True
    assert run_info["run_summary"]["amplicon_success"] == {
        "a1": False, "a2": True, "a3": True,
    }
    assert run_info["run_summary"]["successful_amplicons"] == 2
    assert run_info["run_summary"]["total_amplicons"] == 3
    assert run_info["run_summary"]["consensus_length"] == 623
    assert run_info["run_summary"]["consensus_N_count"] == 0
    utils.rm_rf(outdir)
