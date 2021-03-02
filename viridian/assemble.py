import logging
import os

import pysam

from viridian import amplicon_overlapper, utils
from viridian import amplicons as amps


def map_reads(ref_fasta, reads, outfile, minimap_opts="-t 1 -x map-ont", mates_file=None):
    if minimap_opts is None:
        minimap_opts = ""
    if mates_file is None:
        mates_file = ""
    sam = outfile + ".tmp.sam"
    utils.syscall(f"minimap2 -a {minimap_opts} {ref_fasta} {reads} {mates_file} > {sam}")
    pysam.sort("-o", outfile, "--output-fmt", "BAM", sam)
    os.unlink(sam)
    pysam.index(outfile)


def polish_each_amplicon(
    ref_genome,
    amplicons,
    bam,
    outdir,
    min_mean_coverage=25,
    target_coverage=500,
    read_end_trim=20,
    read_map_tolerance=20,
    min_read_length=200,
    racon_iterations=3,
    min_depth_for_not_N=5,
    amplicons_to_fail=None,
    debug=False,
):
    if amplicons_to_fail is None:
        amplicons_to_fail = set()

    for i, amplicon in enumerate(amplicons):
        logging.debug(
            f"Start processing amplicon {amplicon.name} ({i+1}/{len(amplicons)})"
        )
        if amplicon.name in amplicons_to_fail:
            logging.debug("User chose to fail amplicon {amplcon.name}. Moving on")
            amplicon.force_polish_fail()
            continue

        logging.debug(f"Extracting reads and polishing amplicon {amplicon.name}")
        amplicon_dir = os.path.join(outdir, str(i + 1))
        amplicon.polish(
            ref_genome,
            bam,
            amplicon_dir,
            min_mean_coverage=min_mean_coverage,
            target_coverage=target_coverage,
            read_end_trim=read_end_trim,
            read_map_tolerance=read_map_tolerance,
            min_read_length=min_read_length,
            racon_iterations=racon_iterations,
            min_depth_for_not_N=min_depth_for_not_N,
            debug=debug,
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


def run_assembly_pipeline(
    ref_fasta,
    amplicons_bed,
    outdir,
    sorted_bam=None,
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
    min_amp_overlap_len=20,
    contig_map_end_allowance=20,
    amplicons_to_fail=None,
    debug=False,
):
    ref_genome = utils.load_single_seq_fasta(ref_fasta)
    logging.info(f"Loaded ref genome {ref_genome.id}")
    amplicons = amps.load_amplicons_bed_file(amplicons_bed)
    logging.info(f"Loaded amplicons file {amplicons_bed}")
    os.mkdir(outdir)

    if sorted_bam is None:
        assert reads_fastaq is not None
        logging.info("Reads in FASTA/FASTQ format provided. Mapping reads")
        sorted_bam = os.path.join(outdir, "map_reads.bam")
        map_reads(ref_fasta, reads_fastaq, sorted_bam, minimap_opts=minimap_opts, mates_file=mates_fastaq)
        logging.info("Finished mapping reads")
    else:
        assert reads_fastaq is None

    bam = pysam.AlignmentFile(sorted_bam, "rb")
    polish_root_dir = os.path.join(outdir, "Amplicon_polish")
    os.mkdir(polish_root_dir)
    logging.info("Start polishing each amplicon")
    polish_each_amplicon(
        ref_genome,
        amplicons,
        bam,
        polish_root_dir,
        min_mean_coverage=min_mean_coverage,
        target_coverage=target_coverage,
        read_end_trim=read_end_trim,
        read_map_tolerance=read_map_tolerance,
        min_read_length=min_read_length,
        racon_iterations=racon_iterations,
        min_depth_for_not_N=min_depth_for_not_N,
        amplicons_to_fail=amplicons_to_fail,
        debug=debug,
    )
    if not debug:
        utils.rm_rf(polish_root_dir)

    logging.info("Finished polishing each amplicon")
    any_amplicon_ok = False
    for amplicon in amplicons:
        if amplicon.assemble_success:
            any_amplicon_ok = True
            break

    if not any_amplicon_ok:
        logging.warning("No amplicons successfully polished!")
        consensus = None
    else:
        logging.info("Start making consensus from polished amplicons")
        overlap_out = os.path.join(outdir, "consensus")
        consensus = amplicon_overlapper.assemble_amplicons(
            amplicons,
            ref_fasta,
            overlap_out,
            min_match_length=min_amp_overlap_len,
            ref_map_end_allowance=contig_map_end_allowance,
            debug=debug,
        )
    if consensus is None:
        logging.warning("Did not make consensus sequence. Please see previous warnings")
    else:
        logging.info("Finished making consensus sequence.")
    json_out = os.path.join(outdir, "amplicon_data.json")
    amps.amplicons_to_json(amplicons, json_out)
    return consensus
