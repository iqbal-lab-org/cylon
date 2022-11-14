#!/usr/bin/env python3

import argparse
import logging
import cylon


TECH_DEPENDENT_DEFAULTS = {
    "minimap_opts": {"illumina": "-t 1 -x sr", "ont": "-t 1 -x map-ont"},
    "target_coverage": {"illumina": 150, "ont": 250},
    "min_read_length": {"illumina": 50, "ont": 200},
}


def tech_dependent_usage_default_string(option):
    return "; ".join([f"{k}:{v}" for k, v in TECH_DEPENDENT_DEFAULTS[option].items()])


def set_tech_dependent_args(args):
    for option in TECH_DEPENDENT_DEFAULTS:
        try:
            if getattr(args, option) is None:
                setattr(args, option, TECH_DEPENDENT_DEFAULTS[option][args.tech])
        except:
            pass


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="cylon",
        usage="cylon <command> <options>",
        description="cylon: virus amplicon assembler",
    )

    parser.add_argument("--version", action="version", version=cylon.__version__)
    parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ------------------------ assemble ---------------------------------------
    subparser_assemble = subparsers.add_parser(
        "assemble",
        help="Reference-guided assembly from reads",
        usage="cylon assemble [options] <tech> <ref_fasta> <amplicons_json> <outdir>",
        description="Reference-guided assembly from reads",
        epilog="Required: --reads_per_amp_dir, or --bam, or --reads_to_map, or both --reads_to_map and --mates_to_map",
    )
    subparser_assemble.add_argument(
        "tech",
        help="Sequencing technology. Must be 'ont' or 'illumina'",
        choices=["ont", "illumina"],
        metavar="tech",
    )
    subparser_assemble.add_argument(
        "ref_fasta",
        help="FASTA file of reference genome",
    )
    subparser_assemble.add_argument(
        "amplicons_json", help="BED file of amplicon names and positions"
    )
    subparser_assemble.add_argument("outdir", help="Output directory")

    reads_group = subparser_assemble.add_argument_group(
        "Reads options. Must use: --reads_per_amp_dir; or --bam; or --reads_to_map; or --reads_to_map and --mates_to_map"
    )
    reads_group.add_argument(
        "--bam",
        help="Input reads in a sorted indexed BAM file",
        metavar="FILENAME",
    )
    reads_group.add_argument(
        "--reads_to_map",
        help="Input reads to be mapped, in FASTA or FASTQ format",
        metavar="FILENAME",
    )
    reads_group.add_argument(
        "--mates_to_map",
        help="Input mate reads to be mapped, in FASTA or FASTQ format. If you have paired reads, use this for second file of reads",
        metavar="FILENAME",
    )
    reads_group.add_argument(
        "--reads_per_amp_dir",
        help="Directory of reads, one reads file per amplicon (this option is for viridian_workflow)",
        metavar="DIRNAME",
    )
    subparser_assemble.add_argument(
        "--minimap_opts",
        help=f"Options string to pass to minimap2. Is used for initial mapping if reads provided with --reads_to_map, otherwise is ignored. Do not use -a or -o! This string is not sanity checked - it is up to you to provide valid options [{tech_dependent_usage_default_string('minimap_opts')}]",
        metavar="STRING",
    )
    subparser_assemble.add_argument(
        "--min_mean_coverage",
        help="Minimum mean read depth needed to polish an amplicon. Any amplicon less than this depth is considered failed [%(default)s]",
        default=25,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--target_coverage",
        help=f"Aim for this much coverage when extracting reads for polishing [{tech_dependent_usage_default_string('target_coverage')}]",
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--read_end_trim",
        help="Trim this many bases off the end of all reads [%(default)s]",
        default=0,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--read_map_tolerance",
        help="Max allowed distance read start or end can be outside amplicon coords [%(default)s]",
        default=5,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--min_read_length",
        help=f"Only use reads at least this long (after trimming with --read_end_trim). If the --wgs is used, then the reads must also overlap an 'amplicons' by at least this length to be used [{tech_dependent_usage_default_string('min_read_length')}]",
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--racon_iterations",
        help="Run racon up to this many times (stops if no more corrections made) [%(default)s]",
        default=10,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--min_depth_for_not_N",
        help="After polishing, each position with read depth less than this value will be masked with an N [%(default)s]",
        default=5,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--max_amp_N_proportion",
        help="After polishing, fail the amplicon if more than this proportion of the sequence is Ns [%(default)s]",
        default=0.5,
        type=float,
        metavar="FLOAT",
    )
    subparser_assemble.add_argument(
        "--contig_map_end_allowance",
        help="When mapping contigs ends to the reference to fill in failed amplicons with Ns, allow the mapping to start up to this distance away from the contig end [%(default)s]",
        default=20,
        type=int,
        metavar="INT",
    )
    subparser_assemble.add_argument(
        "--amplicons_to_fail_file",
        help="File of amplicon names to force to count as failed. One name in each line of the file. Names must exactly match those in amplicons_json file",
        metavar="FILENAME",
    )
    subparser_assemble.add_argument(
        "--wgs",
        help="Use this flag if the reads are not from amplicons (eg are WGS or SISPA). It drops the assumption that reads are sequenced from amplicons. This flag is essential if you do not have amplicon data.",
        action="store_true",
    )
    subparser_assemble.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output directory if it already exists",
    )
    subparser_assemble.set_defaults(func=cylon.tasks.assemble.run)

    args = parser.parse_args()
    set_tech_dependent_args(args)

    logging.basicConfig(
        format="[%(asctime)s cylon %(levelname)s] %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    log = logging.getLogger()
    if args.debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
