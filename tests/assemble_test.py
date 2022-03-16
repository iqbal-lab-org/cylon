import json
import os
import pytest
from unittest import mock

from viridian import assemble, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "assemble")


def test_load_and_check_reads_amp_dir():
    tmp_dir = "tmp.load_and_check_reads_amp_dir"
    utils.rm_rf(tmp_dir)
    os.mkdir(tmp_dir)
    json_file = os.path.join(tmp_dir, "manifest.json")
    json_data = {
        "a1": "a1.fasta",
        "a2": "a2.fasta",
    }
    with open(json_file, "w") as f:
        json.dump(json_data, f)

    amp1 = mock.Mock()
    amp1.name = "a1"
    amp2 = mock.Mock()
    amp2.name = "a2"
    amp3 = mock.Mock()
    amp3.name = "a3"

    with pytest.raises(Exception):
        assemble.load_and_check_reads_amp_dir(tmp_dir, set())

    with pytest.raises(Exception):
        assemble.load_and_check_reads_amp_dir(tmp_dir, [amp1])

    with pytest.raises(Exception):
        assemble.load_and_check_reads_amp_dir(tmp_dir, [amp1, amp2])

    for filename in json_data.values():
        with open(os.path.join(tmp_dir, filename), "w"):
            pass

    got = assemble.load_and_check_reads_amp_dir(tmp_dir, [amp1, amp2])
    assert got == {k: os.path.join(tmp_dir, v) for k, v in json_data.items()}

    with pytest.raises(Exception):
        assemble.load_and_check_reads_amp_dir(tmp_dir, [amp1, amp2, amp3])

    utils.rm_rf(tmp_dir)


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
    assert got == expect_seq[11:988]
    consensus_from_file = utils.load_single_seq_fasta(
        os.path.join(outdir, "consensus.final_assembly.fa")
    )
    assert got == consensus_from_file.seq
    utils.rm_rf(outdir)

    # rerun, but using a directory of one file of reads per amplicon. This is
    # what viridian workflow will be making as input
    reads_per_amp_dir = os.path.join(data_dir, "run_assembly_pipeline.reads_per_amp")
    got = assemble.run_assembly_pipeline(
        ref_fa,
        amplicon_json,
        outdir,
        reads_per_amp_dir=reads_per_amp_dir,
        debug=True,
        min_mean_coverage=5,
        min_depth_for_not_N=1,
        read_end_trim=1,
    )
    assert got == expect_seq[10:988]
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
    assert got == expect_seq[356:988]
    consensus_from_file = utils.load_single_seq_fasta(
        os.path.join(outdir, "consensus.final_assembly.fa")
    )
    assert got == consensus_from_file.seq

    # some checks of the contents of the json summary
    with open(os.path.join(outdir, "run_info.json")) as f:
        run_info = json.load(f)
    assert run_info["run_summary"]["made_consensus"] is True
    assert run_info["run_summary"]["amplicon_success"] == {
        "a1": False,
        "a2": True,
        "a3": True,
    }
    assert run_info["run_summary"]["successful_amplicons"] == 2
    assert run_info["run_summary"]["total_amplicons"] == 3
    assert run_info["run_summary"]["consensus_length"] == 632
    assert run_info["run_summary"]["consensus_N_count"] == 0
    utils.rm_rf(outdir)
