import logging
import os

import pyfastaq

from cylon import utils


def run_racon(
    seq_to_polish,
    reads_filename,
    outprefix,
    minimap_opts="-t 1 -x map-ont",
    debug=False,
):
    if minimap_opts is None:
        minimap_opts = "-t 1 -x map-ont"
    fasta_to_polish = f"{outprefix}.to_polish.fa"
    with open(fasta_to_polish, "w") as f:
        if not isinstance(seq_to_polish, pyfastaq.sequences.Fasta):
            assert isinstance(seq_to_polish, str)
            print(">to_polish", file=f)
        print(seq_to_polish, file=f)
    sam = f"{outprefix}.sam"
    utils.syscall(
        f"minimap2 -a {minimap_opts} {fasta_to_polish} {reads_filename} > {sam}"
    )
    # Appears that racon can have errors at N * (window length). Default is
    # 500. So set it to be a bit more that the sequence we are correcting.
    window_length = len(seq_to_polish) + 100
    completed_process = utils.syscall(
        f"racon --window-length {window_length} --no-trimming {reads_filename} {sam} {fasta_to_polish}",
        allow_fail=True,
    )
    if completed_process.returncode != 0:
        return None

    try:
        polished_name, polished_seq = completed_process.stdout.rstrip().split("\n")
    except:
        return None

    if not debug:
        os.unlink(fasta_to_polish)
        os.unlink(sam)

    return polished_seq


def run_racon_iterations(
    seq_to_polish,
    reads_filename,
    outdir,
    max_iterations=3,
    debug=False,
    minimap_opts=None,
):
    reads_filename = os.path.abspath(reads_filename)
    os.mkdir(outdir)

    for iteration in range(max_iterations):
        logging.debug(f"Start racon iteration {iteration}")
        outprefix = os.path.join(outdir, f"racon.{iteration}")
        polished_seq = run_racon(
            seq_to_polish,
            reads_filename,
            outprefix,
            debug=debug,
            minimap_opts=minimap_opts,
        )
        if polished_seq is None:
            return None
        if debug:
            with open(f"{outprefix}.polished.fa", "w") as f:
                print(">" + "polished", file=f)
                print(polished_seq, file=f)
        if polished_seq != seq_to_polish:
            logging.debug(f"  racon iteration {iteration} finished. Made changes")
            seq_to_polish = polished_seq
        else:
            logging.debug(f"  racon iteration {iteration} finished. No changes")
            break

    return polished_seq
