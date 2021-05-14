import logging
import os

import pyfastaq

from viridian import utils


def run_racon(seq_to_polish, reads_filename, outprefix, debug=False):
    fasta_to_polish = f"{outprefix}.to_polish.fa"
    with open(fasta_to_polish, "w") as f:
        if not isinstance(seq_to_polish, pyfastaq.sequences.Fasta):
            assert isinstance(seq_to_polish, str)
            print(">to_polish", file=f)
        print(seq_to_polish, file=f)
    sam = f"{outprefix}.sam"
    utils.syscall(
        f"minimap2 -a -x map-ont -t 1 {fasta_to_polish} {reads_filename} > {sam}"
    )
    completed_process = utils.syscall(
        f"racon --no-trimming {reads_filename} {sam} {fasta_to_polish}", allow_fail=True
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
    seq_to_polish, reads_filename, outdir, max_iterations=3, debug=False
):
    reads_filename = os.path.abspath(reads_filename)
    os.mkdir(outdir)

    for iteration in range(max_iterations):
        logging.debug(f"Start racon iteration {iteration}")
        outprefix = os.path.join(outdir, f"racon.{iteration}")
        polished_seq = run_racon(seq_to_polish, reads_filename, outprefix, debug=debug)
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
