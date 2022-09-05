import argparse
import datetime
import json
import logging
import os
import socket
import sys
import tempfile

import pysam

from cylon import amplicon_overlapper, utils
from cylon import amplicons as amps
from cylon import __version__ as cylon_version


def map_reads(
    ref_fasta, reads, outfile, minimap_opts="-t 1 -x map-ont", mates_file=None
):
    if minimap_opts is None:
        minimap_opts = ""
    if mates_file is None:
        mates_file = ""
    sam = outfile + ".tmp.sam"
    utils.syscall(
        f"minimap2 -a {minimap_opts} {ref_fasta} {reads} {mates_file} > {sam}"
    )
    pysam.sort("-o", outfile, "--output-fmt", "BAM", sam)
    os.unlink(sam)
    pysam.index(outfile)


def polish_each_amplicon(
    ref_genome,
    amplicons,
    outdir,
    bam_to_slice_reads=None,
    amplicon_to_reads_file=None,
    min_mean_coverage=25,
    target_coverage=500,
    read_end_trim=20,
    read_map_tolerance=20,
    min_read_length=200,
    racon_iterations=3,
    min_depth_for_not_N=5,
    amplicons_to_fail=None,
    wgs=False,
    debug=False,
    minimap_opts=None,
):
    if amplicon_to_reads_file is None:
        amplicon_to_reads_file = {}
    if amplicons_to_fail is None:
        amplicons_to_fail = set()

    for i, amplicon in enumerate(amplicons):
        logging.debug(
            f"Start processing amplicon {amplicon.name} ({i+1}/{len(amplicons)})"
        )
        if amplicon.name in amplicons_to_fail:
            logging.debug(f"User chose to fail amplicon {amplicon.name}. Moving on")
            amplicon.force_polish_fail()
            continue

        logging.debug(f"Extracting reads and polishing amplicon {amplicon.name}")
        amplicon_dir = os.path.join(outdir, str(i + 1))
        amplicon.polish(
            ref_genome,
            amplicon_dir,
            bam_to_slice_reads=bam_to_slice_reads,
            reads_file=amplicon_to_reads_file.get(amplicon.name, None),
            min_mean_coverage=min_mean_coverage,
            target_coverage=target_coverage,
            read_end_trim=read_end_trim,
            read_map_tolerance=read_map_tolerance,
            min_read_length=min_read_length,
            racon_iterations=racon_iterations,
            min_depth_for_not_N=min_depth_for_not_N,
            wgs=wgs,
            debug=debug,
            minimap_opts=minimap_opts,
        )
        ok = "yes" if amplicon.assemble_success else "no"
        logging.debug(f"Finish polishing amplicon {amplicon.name}. Success: {ok}")
        if not debug:
            utils.rm_rf(amplicon_dir)

        logging.debug(
            f"Finish processing amplicon {amplicon.name} ({i+1}/{len(amplicons)})"
        )
        if i % 10 == 0:
            logging.info(f"Processed {i+1} of {len(amplicons)} amplicons")


def add_successful_amplicons_to_json_data(data, amplicons):
    data["run_summary"]["successful_amplicons"] = len(
        [a for a in amplicons if a.assemble_success]
    )
    data["run_summary"]["amplicon_success"] = {
        a.name: a.assemble_success for a in amplicons
    }


def add_consensus_length_N_count_to_json_data(data):
    if data["run_summary"]["made_consensus"]:
        data["run_summary"]["consensus_length"] = len(data["run_summary"]["consensus"])
        data["run_summary"]["consensus_N_count"] = data["run_summary"][
            "consensus"
        ].count("N")
    else:
        data["run_summary"]["consensus_length"] = None
        data["run_summary"]["consensus_N_count"] = None


def load_and_check_reads_amp_dir(reads_per_amp_dir, amplicons):
    json_file = os.path.join(reads_per_amp_dir, "manifest.json")
    with open(json_file) as f:
        manifest = json.load(f)

    if len(amplicons) < len(manifest):
        raise Exception(
            f"Amplicon scheme has {len(amplicons)} amplicons, which is less than number of reads files in {json_file}"
        )
    all_amp_names = set([x.name for x in amplicons])

    for amplicon_name in manifest:
        if amplicon_name not in all_amp_names:
            raise Exception(
                f"Amplicon '{amplicon_name}' in json {json_file} but not in amplicon scheme"
            )

        reads_file_full_path = os.path.join(reads_per_amp_dir, manifest[amplicon_name])
        if not os.path.exists(reads_file_full_path):
            raise Exception(
                f"Reads file for amplicon '{amplicon_name}' not found: {reads_file_full_path}"
            )
        manifest[amplicon_name] = reads_file_full_path

    return manifest


def run_assembly_pipeline(
    ref_fasta,
    amplicons_json,
    outdir,
    sorted_bam=None,
    reads_per_amp_dir=None,
    reads_fastaq=None,
    mates_fastaq=None,
    minimap_opts=None,
    min_mean_coverage=25,
    target_coverage=500,
    read_end_trim=20,
    read_map_tolerance=20,
    min_read_length=200,
    racon_iterations=3,
    min_depth_for_not_N=5,
    contig_map_end_allowance=20,
    amplicons_to_fail=None,
    wgs=False,
    debug=False,
    command_line_args=None,
):
    # Make a dict of the command line options to go in the JSON output file.
    # The tests don't use argparse (they use Mock), which means convert to dict
    # doesn't work. Don't care about that case anyway in the final output, so
    # just set to None
    if isinstance(command_line_args, argparse.Namespace):
        options_dict = {k: v for k, v in vars(command_line_args).items() if k != "func"}
    else:
        options_dict = None

    start_time = datetime.datetime.now()
    os.mkdir(outdir)
    json_out = os.path.join(outdir, "run_info.json")

    json_data = {
        "run_summary": {
            "total_amplicons": None,
            "successful_amplicons": None,
            "command": " ".join(sys.argv),
            "options": options_dict,
            "cwd": os.getcwd(),
            "version": cylon_version,
            "finished_running": False,
            "made_consensus": False,
            "consensus": None,
            "start_time": start_time.replace(microsecond=0).isoformat(),
            "end_time": None,
            "hostname": socket.gethostname(),
        },
        "amplicons": None,
    }
    with open(json_out, "w") as f:
        json.dump(json_data, f, indent=2, sort_keys=True)

    ref_genome = utils.load_single_seq_fasta(ref_fasta)
    logging.info(f"Loaded ref genome {ref_genome.id}")
    amplicons = amps.load_amplicons_json_file(amplicons_json)
    json_data["run_summary"]["total_amplicons"] = len(amplicons)
    logging.info(f"Loaded amplicons file {amplicons_json}")
    amplicon_to_reads_file = None

    if reads_per_amp_dir is not None:
        assert reads_fastaq is None
        assert sorted_bam is None
        amplicon_to_reads_file = load_and_check_reads_amp_dir(
            reads_per_amp_dir, amplicons
        )
        bam = None
    elif reads_fastaq is not None:
        assert sorted_bam is None
        assert reads_per_amp_dir is None
        logging.info("Reads in FASTA/FASTQ format provided. Mapping reads")
        sorted_bam = os.path.join(outdir, "map_reads.bam")
        map_reads(
            ref_fasta,
            reads_fastaq,
            sorted_bam,
            minimap_opts=minimap_opts,
            mates_file=mates_fastaq,
        )
        logging.info("Finished mapping reads")
        bam = pysam.AlignmentFile(sorted_bam, "rb")
    else:
        assert sorted_bam is not None
        bam = pysam.AlignmentFile(sorted_bam, "rb")

    if debug:
        polish_root_dir = os.path.join(outdir, "Amplicon_polish")
        os.mkdir(polish_root_dir)
    else:
        polish_root_dir = tempfile.mkdtemp(prefix="cylon_polish_")

    logging.info(f"Start polishing each amplicon. Directory: {polish_root_dir}")
    try:
        polish_each_amplicon(
            ref_genome,
            amplicons,
            polish_root_dir,
            bam_to_slice_reads=bam,
            amplicon_to_reads_file=amplicon_to_reads_file,
            min_mean_coverage=min_mean_coverage,
            target_coverage=target_coverage,
            read_end_trim=read_end_trim,
            read_map_tolerance=read_map_tolerance,
            min_read_length=min_read_length,
            racon_iterations=racon_iterations,
            min_depth_for_not_N=min_depth_for_not_N,
            amplicons_to_fail=amplicons_to_fail,
            wgs=wgs,
            debug=debug,
            minimap_opts=minimap_opts,
        )
    finally:
        if not debug:
            utils.rm_rf(polish_root_dir)

    logging.info("Finished polishing each amplicon")
    add_successful_amplicons_to_json_data(json_data, amplicons)
    if json_data["run_summary"]["successful_amplicons"] == 0:
        logging.warning("No amplicons successfully polished!")
        consensus = None
    else:
        logging.info("Start making consensus from polished amplicons")
        overlap_out = os.path.join(outdir, "consensus")
        consensus = amplicon_overlapper.assemble_amplicons(
            amplicons,
            ref_fasta,
            overlap_out,
            ref_map_end_allowance=contig_map_end_allowance,
            debug=debug,
        )
    json_data["run_summary"]["consensus"] = consensus

    # Need to recalculate successful amplicons because they can get failed
    # during overlapping. If two adjacent amplicons have no overlap, then they
    # both get failed.
    add_successful_amplicons_to_json_data(json_data, amplicons)

    if consensus is None:
        logging.warning("Did not make consensus sequence. Please see previous warnings")
    else:
        logging.info("Finished making consensus sequence.")
        json_data["run_summary"]["made_consensus"] = True

    add_consensus_length_N_count_to_json_data(json_data)
    json_data["amplicons"] = amps.amplicons_to_list_of_dicts(amplicons)
    json_data["run_summary"]["finished_running"] = True
    end_time = datetime.datetime.now()
    json_data["run_summary"]["end_time"] = end_time.replace(microsecond=0).isoformat()
    json_data["run_summary"]["run_time"] = str(end_time - start_time)
    with open(json_out, "w") as f:
        json.dump(json_data, f, indent=2, sort_keys=True)
    return consensus
