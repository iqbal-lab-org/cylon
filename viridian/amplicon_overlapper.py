import logging
import os

import pyfastaq
from viridian import utils


def get_amplicon_overlaps(amplicons, min_match_length):
    overlaps = []
    for i in range(0, len(amplicons) - 1, 1):
        overlaps.append(amplicons[i].masked_overlap(amplicons[i + 1], min_match_length))
    return overlaps


def amplicons_to_consensus_contigs(amplicons, min_match_length=20):
    if all(not x.assemble_success for x in amplicons):
        return None

    overlaps = get_amplicon_overlaps(amplicons, min_match_length)
    contigs = [[]]

    # force adjacent contigs that have no overlap to be fails because
    # we don't really know what happened
    indexes_to_clear = set()

    for i in range(0, len(amplicons) - 1):
        if (
            amplicons[i].assemble_success
            and amplicons[i + 1].assemble_success
            and overlaps[i] is None
        ):
            indexes_to_clear.add(i)
            indexes_to_clear.add(i + 1)

    for i in indexes_to_clear:
        amplicons[i].clear_seqs_because_overlap_fail()
        if i > 0:
            overlaps[i - 1] = None
        if i < len(overlaps):
            overlaps[i] = None

    if all(not x.assemble_success for x in amplicons):
        return None

    for i in range(0, len(amplicons)):
        if amplicons[i].assemble_success:
            if i == 0 or overlaps[i - 1] is None:
                start = 0
            else:
                start = max(overlaps[i - 1].b, 0)

            if i > len(overlaps) - 1 or overlaps[i] is None:
                end = len(amplicons[i].masked_seq)
            else:
                end = overlaps[i].a

            if len(contigs) == 0:
                contigs.append([])
            contigs[-1].append(amplicons[i].masked_seq[start:end])
        else:
            if len(contigs) > 0:
                if len(contigs[-1]) == 0:
                    contigs.pop()
                else:
                    contigs.append([])

    contigs = ["".join(x).strip("N") for x in contigs if len(x) > 0]
    return contigs


def _make_split_contigs_fasta(contigs, outfile):
    with open(outfile, "w") as f:
        for i, contig in enumerate(contigs):
            print(f">{i}.left", file=f)
            print(contig[: int(len(contig) / 2)], file=f)
            print(f">{i}.right", file=f)
            print(contig[int(len(contig) / 2) :], file=f)


def _map_split_contigs(to_map_fasta, ref_fasta, end_allowance=20):
    command = f"minimap2 -t 1 -c {ref_fasta} {to_map_fasta}"
    minimap2_out = utils.syscall(command)
    paf_lines = minimap2_out.stdout.strip().split("\n")
    mappings = {}
    for line in paf_lines:
        fields = line.rstrip().split()
        if fields[4] != "+":
            continue

        index, left_or_right = fields[0].split(".")
        if (left_or_right == "left" and int(fields[2]) > end_allowance) or (
            left_or_right == "right" and int(fields[3]) + end_allowance < int(fields[3])
        ):
            continue

        index = int(index)
        key = (index, left_or_right)

        if key not in mappings:
            mappings[key] = fields
        elif (left_or_right == "left" and int(fields[7]) < int(mappings[key][7])) or (
            left_or_right == "right" and int(fields[8]) > int(mappings[key][8])
        ):
            mappings[key] = fields

        logging.debug(
            "PAF mapping contigs to ref: " + " ".join([str(x) for x in fields])
        )
    return mappings


def _check_mappings(contigs, mappings):
    for i, contig in enumerate(contigs):
        left_key = (i, "left")
        right_key = (i, "right")
        if left_key not in mappings:
            logging.warning(
                f"Error mapping left half of contig to reference genome: {contig}"
            )
            return False
        elif right_key not in mappings:
            logging.warning(
                f"Error mapping right half of contig to reference genome: {contig}"
            )
            return False
    return True


def consensus_contigs_to_consensus(
    contigs, ref_fasta, outprefix, map_end_allowance=20, debug=False
):
    if contigs is None or len(contigs) == 0:
        logging.warning("No contigs were made. Aborting assembly")
        return None

    fa_to_map = f"{outprefix}.to_map.fa"
    _make_split_contigs_fasta(contigs, fa_to_map)
    mappings = _map_split_contigs(fa_to_map, ref_fasta, end_allowance=map_end_allowance)
    if not debug:
        os.unlink(fa_to_map)
    if not _check_mappings(contigs, mappings):
        logging.warning(
            "Errors aligning contigs to reference to make final sequence. Aborting assembly"
        )
        return None

    consensus = []

    for i, contig in enumerate(contigs):
        right_fields = mappings[(i, "right")]

        consensus.append(contig)
        if i < len(contigs) - 1:
            next_left_fields = mappings[(i + 1, "left")]
            contig_end_in_ref = int(right_fields[8])
            next_contig_start_in_ref = int(next_left_fields[7])
            if contig_end_in_ref < next_contig_start_in_ref:
                consensus.append("N" * (next_contig_start_in_ref - contig_end_in_ref))
            else:
                logging.warning(
                    f"Errors aligning contigs to reference to make final sequence. Order of contigs vs reference not as epxected. Cannot get gap between contigs: {contig[i]} and {contig[i+1]}. Aborting assembly"
                )
                return None

    return "".join(consensus)


def assemble_amplicons(
    amplicons,
    ref_fasta,
    outprefix,
    min_match_length=20,
    ref_map_end_allowance=20,
    debug=False,
):
    consensus_contigs = amplicons_to_consensus_contigs(
        amplicons, min_match_length=min_match_length
    )
    if consensus_contigs is None:
        logging.warning("Errors trying to overlap contigs. Aborting assembly")
        return

    with open(f"{outprefix}.contigs.fa", "w") as f:
        for i, contig in enumerate(consensus_contigs):
            print(f">{i}", file=f)
            print(contig, file=f)

    consensus_prefix = f"{outprefix}.consensus.tmp"
    consensus = consensus_contigs_to_consensus(
        consensus_contigs,
        ref_fasta,
        consensus_prefix,
        map_end_allowance=ref_map_end_allowance,
        debug=debug,
    )
    with open(f"{outprefix}.final_assembly.fa", "w") as f:
        seq = pyfastaq.sequences.Fasta("assembly", consensus)
        print(seq, file=f)
    return consensus
