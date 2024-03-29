from cylon import assemble, utils


def run(options):
    if not utils.look_for_required_binaries_in_path():
        raise Exception(
            "At least one required program was not found in $PATH. Cannot continue"
        )

    if (
        len(
            [
                x
                for x in (options.bam, options.reads_to_map, options.reads_per_amp_dir)
                if x is not None
            ]
        )
        != 1
    ):
        raise Exception(
            "Must provide exactly one of: --bam, --reads_to_map, --reads_per_amp_dir"
        )

    if options.mates_to_map is not None and options.reads_to_map is None:
        raise Exception(
            "--mates_to_map was used, but --reads_to_map was not. --reads_to_map is required by --mates_to_map"
        )

    if options.force:
        utils.rm_rf(options.outdir)

    if options.amplicons_to_fail_file is None:
        amplicons_to_fail = None
    else:
        with open(options.amplicons_to_fail_file) as f:
            amplicons_to_fail = set([x.rstrip() for x in f])

    return assemble.run_assembly_pipeline(
        options.ref_fasta,
        options.amplicons_json,
        options.outdir,
        sorted_bam=options.bam,
        reads_per_amp_dir=options.reads_per_amp_dir,
        reads_fastaq=options.reads_to_map,
        mates_fastaq=options.mates_to_map,
        minimap_opts=options.minimap_opts,
        min_mean_coverage=options.min_mean_coverage,
        target_coverage=options.target_coverage,
        read_end_trim=options.read_end_trim,
        read_map_tolerance=options.read_map_tolerance,
        min_read_length=options.min_read_length,
        racon_iterations=options.racon_iterations,
        min_depth_for_not_N=options.min_depth_for_not_N,
        max_polished_N_prop=options.max_amp_N_proportion,
        contig_map_end_allowance=options.contig_map_end_allowance,
        amplicons_to_fail=amplicons_to_fail,
        wgs=options.wgs,
        debug=options.debug,
        command_line_args=options,
    )
