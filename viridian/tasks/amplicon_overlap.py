import os

from viridian import amplicon_overlapper
from viridian import amplicons as amps


def run(options):
    amplicons = amps.load_amplicons_from_fasta_and_bed(
        options.amplicons_fasta, options.amplicons_bed
    )
    os.mkdir(options.outdir)
    outprefix = os.path.join(options.outdir, "out")
    return amplicon_overlapper.assemble_amplicons(
        amplicons,
        options.ref_fasta,
        outprefix,
        min_match_length=options.min_amp_overlap_len,
        ref_map_end_allowance=options.contig_map_end_allowance,
    )
