import logging
import os
import re

import pyfastaq
import pymummer
from Bio import pairwise2

from cylon import utils



def global_align(seq1, seq2):
    """Returns global alignment strings from N-W alignment of the
    two sequences. Dashes for gaps"""
    match=1
    mismatch=-1
    gap_open=-5
    gap_extend=-3

    alignments = pairwise2.align.globalms(
        seq1,
        seq2,
        1, # match
        -1, # mismatch
        -5, #gap_open
        -3, #gap_extend
    )
    # Alignments is a list of tuples. Each tuple has length 5. Entries:
    # 0: seq1 alignment (ie with dashes for indels)
    # 1: seq2 alignemnt
    # 2: alignment score
    # 4, 5: don't know (not using them)
    if len(alignments) == 0:
        return None
    else:
        return alignments[0][0], alignments[0][1]


def make_split_contigs_fasta(contigs, outfile):
    with open(outfile, "w") as f:
        for d in contigs:
            length = min(250, int(len(d["seq"]) / 2))
            print(f">{d['name']}.left", file=f)
            print(d["seq"][:length], file=f)
            print(f">{d['name']}.right", file=f)
            print(d["seq"][-length:], file=f)


def map_split_contigs(to_map_fasta, ref_fasta, end_allowance=20, debug=False):
    tmp_nucmer_out = f"{to_map_fasta}.tmp.nucmer"
    runner = pymummer.nucmer.Runner(
        ref_fasta, to_map_fasta, tmp_nucmer_out, mincluster=20
    )
    runner.run()
    mappings = {}
    for hit in pymummer.coords_file.reader(tmp_nucmer_out):
        if not hit.on_same_strand():
            continue

        name, left_or_right = hit.qry_name.split(".")
        if (left_or_right == "left" and hit.qry_start > end_allowance) or (
            left_or_right == "right" and hit.qry_end + end_allowance < hit.qry_length
        ):
            continue

        key = (name, left_or_right)

        if key not in mappings or mappings[key].hit_length_qry < hit.hit_length_qry:
            mappings[key] = hit

        logging.debug(f"nucmer mapping contigs to ref {hit}")

    if not debug:
        os.unlink(tmp_nucmer_out)

    return mappings


def remove_contigs_without_both_ends_well_mapped(contigs, mappings):
    bad_contig_indexes = set()
    for i, contig in enumerate(contigs):
        left_key = (contig["name"], "left")
        right_key = (contig["name"], "right")
        if left_key not in mappings:
            logging.warning(
                f"Error mapping left half of contig to reference genome. Contig excluded from consensus: {contig['seq']}"
            )
            bad_contig_indexes.add(i)
        elif right_key not in mappings:
            logging.warning(
                f"Error mapping right half of contig to reference genome. Contig excluded from consensus: {contig['seq']}"
            )
            bad_contig_indexes.add(i)
        elif not mappings[left_key].ref_start < mappings[left_key].ref_end <= mappings[right_key].ref_start < mappings[right_key].ref_end:
            bad_contig_indexes.add(i)

    return [contigs[i] for i in range(len(contigs)) if i not in bad_contig_indexes]


def remove_out_of_order_or_contained_contigs(contigs, mappings):
    i = 0
    while i < len(contigs) - 1:
        this_left = mappings[(contigs[i]["name"], "left")]
        this_right = mappings[(contigs[i]["name"], "right")]
        next_left = mappings[(contigs[i+1]["name"], "left")]
        next_right = mappings[(contigs[i+1]["name"], "right")]

        if this_left.ref_start <= next_left.ref_start <= next_right.ref_end <= this_right.ref_end:
            contigs.pop(i+1)
            continue
        elif next_left.ref_start <= this_left.ref_start <= this_right.ref_end <= next_right.ref_end:
            contigs.pop(i)
            continue

        if next_left.ref_start < this_left.ref_start or this_right.ref_end > next_right.ref_end:
            if this_right.ref_end -  this_left.ref_start > next_right.ref_end - next_right.ref_start:
                contigs.pop(i+1)
            else:
                contigs.pop(i)
            continue

        i += 1

    return contigs


def split_contigs_on_gaps(contigs_in, gap_length=10):
    contigs_out = []
    regex = re.compile("N{" + str(gap_length) + ",}")
    for contig in contigs_in:
        contigs_out.extend([x for x in re.split(regex, contig) if len(x) > 0])
    return contigs_out


def make_trimmed_contigs(contigs_in, window=20):
    contigs_out = []
    contigs_after_gaps = split_contigs_on_gaps(contigs_in)
    for i, seq in enumerate(contigs_after_gaps):
        start = 0
        while start < len(seq) - window and seq[start : start + window].count("N") > 0:
            start += 1
        end = len(seq)
        while end > window and seq[end - window : end].count("N") > 0:
            end -= 1
        if start < end:
            logging.debug(f"contig trim {i}, {len(seq)}, {start}, {end}, {seq}")
            contigs_out.append({"name": str(i), "seq": seq[start:end]})
            logging.debug(
                f"contig trimmed {i}, {len(seq)}, {start}, {end}, {seq[start:end]}"
            )

    return contigs_out


def self_map_contigs(contigs, outprefix, debug=False):
    to_map_fasta = f"{outprefix}.fa"
    with open(to_map_fasta, "w") as f:
        for c in contigs:
            print(">" + c["name"], c["seq"], sep="\n", file=f)
    nucmer_out = f"{outprefix}.nucmer"
    runner = pymummer.nucmer.Runner(
        to_map_fasta, to_map_fasta, nucmer_out, mincluster=5, maxmatch=True,
    )
    runner.run()
    mappings = {}
    for hit in pymummer.coords_file.reader(nucmer_out):
        if not hit.on_same_strand() or int(hit.ref_name) + 1 != int(hit.qry_name):
            continue

        key = (int(hit.ref_name), int(hit.qry_name))
        if key not in mappings or hit.hit_length_qry > mappings[key].hit_length_qry:
            mappings[key] = hit

        logging.debug(f"nucmer mapping contigs to each other {hit}")

    if not debug:
        os.unlink(to_map_fasta)
        os.unlink(nucmer_out)

    return mappings


def consensus_contigs_to_consensus(
    contigs, ref_fasta, outprefix, map_end_allowance=20, debug=False, trim_end_window=20
):
    if contigs is None or len(contigs) == 0:
        logging.warning("No contigs were made. Aborting assembly")
        return None

    fa_to_map = f"{outprefix}.to_map.fa"
    contigs = make_trimmed_contigs(contigs, window=trim_end_window)
    make_split_contigs_fasta(contigs, fa_to_map)
    split_mappings = map_split_contigs(fa_to_map, ref_fasta, end_allowance=map_end_allowance)
    if not debug:
        os.unlink(fa_to_map)
    contigs = remove_contigs_without_both_ends_well_mapped(contigs, split_mappings)
    contigs = remove_out_of_order_or_contained_contigs(contigs, split_mappings)
    if len(contigs) == 0:
        logging.warning(
            "Errors aligning contigs to reference to make final sequence. Aborting assembly"
        )
        return None

    self_map_out = f"{outprefix}.self_map_contigs"
    self_mappings = self_map_contigs(contigs, self_map_out, debug=debug)
    consensus = []
    trim_next = 0

    for i, contig in enumerate(contigs):
        right_hit = split_mappings[(contig["name"], "right")]
        if i < len(contigs) - 1:
            next_left_hit = split_mappings[(contigs[i + 1]["name"], "left")]

            # If the contigs do not overlap each other based on ref coords
            if next_left_hit.ref_start > right_hit.ref_end:
                consensus.append(contig["seq"][trim_next:])
                consensus.append("N" * (next_left_hit.ref_start - right_hit.ref_end - 1))
                trim_next = 0
            else: # split mapping to ref coords say the contigs overlap
                self_hit = self_mappings.get((i, i+1), None)

                if self_hit is None:
                    # No overlapping match found by nucmer, but from reference
                    # mapping contig ends, they should overlap. If we're here
                    # then things may be bit dodgy. Likely a tiny overlap
                    # and/or enough disagreement between the contigs to
                    # mean we didn't get a match. Get the expected overlapping
                    # sequence from this and the next contig. If they are the
                    # same then all good, otherwise align them and use the
                    # first part from this contig, and the second part from the
                    # next contig

                    # Pull out the sequences that should match
                    overlap_len = right_hit.ref_end - next_left_hit.ref_start + 1
                    this_end_in_contig = len(contig["seq"]) - right_hit.qry_length + right_hit.qry_end
                    this_ol_start_in_contig = 1 + this_end_in_contig - overlap_len
                    this_ol_seq = contig["seq"][this_ol_start_in_contig:this_end_in_contig + 1]
                    next_ol_seq = contigs[i+1]["seq"][next_left_hit.qry_start:next_left_hit.qry_start+overlap_len]

                    # Nice case: the two sequences are same
                    if this_ol_seq == next_ol_seq:
                        consensus.append(contig["seq"][trim_next:right_hit.ref_end + 1])
                    else: # different overlapping sequences
                        match = global_align(this_ol_seq, next_ol_seq)
                        if match is None:
                            # not a lot can do here. Bail out, trim the contig
                            # ends and stick some Ns in there
                            consensus.append(contig["seq"][trim_next:right_hit.this_ol_start_in_contig])
                            consensus.append("N" * overlap_len)
                        else:
                            first_good = None
                            for i in range(len(match[0])):
                                if match[0][i] == match[1][i] and match[0] != "-":
                                    first_good = i
                                    break
                            if first_good is None:
                                # not a lot can do here. Bail out, trim the contig
                                # ends and stick some Ns in there
                                consensus.append(contig["seq"][trim_next:right_hit.this_ol_start_in_contig])
                                consensus.append("N" * overlap_len)
                            else:
                                new_seq = (match[0][:i] + match[1][i:]).replace("-", "")
                                consensus.append(contig["seq"][trim_next:this_ol_start_in_contig] + new_seq)
                    trim_next = next_left_hit.qry_start + overlap_len
                else: # overlap bewteen this/next contig found by nucmer
                    assert self_hit.qry_name == contigs[i+1]["name"]
                    consensus.append(contig["seq"][trim_next:self_hit.ref_start])
                    right_ns = contig["seq"][self_hit.ref_start:self_hit.ref_end].count("N")
                    next_left_ns = contigs[i+1]["seq"][self_hit.qry_start:self_hit.qry_end].count("N")
                    if right_ns < next_left_ns:
                        consensus.append(contig["seq"][self_hit.ref_start:self_hit.ref_end])
                        trim_next = self_hit.qry_end
                    else:
                        trim_next = self_hit.qry_start
        else: # add the last piece of the final contig
            consensus.append(contig["seq"][trim_next:])

    consensus = "".join([x for x in consensus if len(x) > 0])
    if any([x != "N" for x in consensus]):
        return consensus
    else:
        logging.warning(
            "Error making consensus, did not get any non-N bases"
        )
        return None


def _get_amplicon_ref_matches(amplicons, ref_fasta, outfile, debug=False):
    amp_fa = f"{outfile}.tmp.amps.fa"
    any_ok = False
    with open(amp_fa, "w") as f:
        for i, amplicon in enumerate(amplicons):
            if amplicon.assemble_success:
                any_ok = True
                print(f">{i}", amplicon.name, file=f)
                print(amplicon.final_seq, file=f)

    if not any_ok:
        logging.warning("No good amplicons to overlap. Aborting assembly")
        if not debug:
            os.unlink(amp_fa)
        return {}

    utils.syscall(f"minimap2 -x sr -t 1 {ref_fasta} {amp_fa} > {outfile}")
    if not debug:
        os.unlink(amp_fa)
    matches = {}

    with open(outfile) as f:
        for line in f:
            (
                amp_index,
                _,
                amp_start,
                amp_end,
                strand,
                _,
                _,
                ref_start,
                ref_end,
                match_len,
                *_,
            ) = line.rstrip().split("\t")
            if strand != "+":
                continue
            new_match = {
                "amp_start": int(amp_start),
                "amp_end": int(amp_end),
                "ref_start": int(ref_start),
                "ref_end": int(ref_end),
                "len": int(match_len),
            }
            i = int(amp_index)
            if i in matches and matches[i]["len"] > new_match["len"]:
                continue
            matches[i] = new_match
            amplicons[i].ref_match = new_match

    return matches


def _clear_no_ref_matches(amplicons, ref_matches):
    indexes_to_clear = set()
    for i in range(0, len(amplicons) - 1):
        if amplicons[i].assemble_success and i not in ref_matches:
            indexes_to_clear.add(i)

    for i in indexes_to_clear:
        logging.debug(
            f"Failing amplicon because no match to reference: {amplicons[i].name}"
        )
        amplicons[i].clear_seqs_because_no_ref_match()


def _get_starts_and_ends(amplicons, ref_matches):
    starts_ends = [[None, None] for _ in range(len(amplicons))]
    add_breaks = set()
    fails = set()
    min_overlap_len = 10
    for i in range(len(amplicons)):
        if i in ref_matches and (i + 1) in ref_matches:
            expect_ol_len = ref_matches[i]["ref_end"] - ref_matches[i + 1]["ref_start"]
            if expect_ol_len > 0:
                if expect_ol_len >= min_overlap_len:
                    perfect_ol = amplicons[i].final_overlap(
                        amplicons[i + 1],
                        min_overlap_len,
                        self_start=ref_matches[i]["amp_end"] - expect_ol_len - 5,
                        other_end=ref_matches[i + 1]["amp_end"] + 5,
                    )
                else:
                    perfect_ol = None

                if perfect_ol is None:
                    end = ref_matches[i]["amp_end"] - expect_ol_len
                    start = ref_matches[i + 1]["amp_start"]
                    if end < 0 or start > len(amplicons[i + 1]):
                        fails.add(i)
                    else:
                        starts_ends[i][1] = end
                        starts_ends[i + 1][0] = start
                else:
                    starts_ends[i][1] = perfect_ol.a
                    starts_ends[i + 1][0] = max(perfect_ol.b, 0)
            else:  # no overlap from the minimap mappings
                starts_ends[i][1] = len(amplicons[i].final_seq)
                starts_ends[i + 1][0] = 0
                add_breaks.add(i)
        elif i in ref_matches:  # this one matches, next one doesn't
            starts_ends[i][1] = len(amplicons[i].final_seq)
        elif (i + 1) in ref_matches:  # this doesn't match, next one does
            starts_ends[i + 1][0] = 0
            fails.add(i)
        else:  # neither matched the ref
            fails.add(i)

        if i in ref_matches and starts_ends[i][0] is None:
            starts_ends[i][0] = 0

    last_i = len(amplicons) - 1
    if (
        starts_ends[last_i][0] is not None
        and last_i in ref_matches
        and last_i not in fails
    ):
        starts_ends[last_i][1] = len(amplicons[last_i].final_seq)

    return starts_ends, add_breaks


def _contigs_from_starts_ends(amplicons, starts_ends, add_breaks):
    contigs = [[]]
    for i, amplicon in enumerate(amplicons):
        logging.debug(
            f"Overlapping amplicons. Processing {amplicon.name}. Assemble succes: {amplicon.assemble_success}"
        )
        if None not in starts_ends[i]:
            contigs[-1].append(
                amplicon.final_seq[starts_ends[i][0] : starts_ends[i][1]]
            )
        elif len(contigs[-1]) > 0:
            contigs.append([])

        if i in add_breaks:
            contigs.append([])

    return ["".join(x).strip("N") for x in contigs if len(x) > 0]


def amplicons_to_consensus_contigs(amplicons, ref_fasta, minimap_out, debug=False):
    ref_matches = _get_amplicon_ref_matches(
        amplicons, ref_fasta, minimap_out, debug=debug
    )
    if len(ref_matches) == 0:
        logging.warning(
            "No amplicon consensus sequences matched the reference. Aborting assembly"
        )
    if not debug:
        utils.rm_rf(minimap_out)

    _clear_no_ref_matches(amplicons, ref_matches)
    starts_ends, add_breaks = _get_starts_and_ends(amplicons, ref_matches)
    if all(not x.assemble_success for x in amplicons):
        return None
    return _contigs_from_starts_ends(amplicons, starts_ends, add_breaks)


def assemble_amplicons(
    amplicons,
    ref_fasta,
    outprefix,
    ref_map_end_allowance=20,
    debug=False,
):
    minimap_out = f"{outprefix}.minimap2.paf"
    consensus_contigs = amplicons_to_consensus_contigs(
        amplicons, ref_fasta, minimap_out, debug=debug
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
    if consensus is None:
        logging.warning("No consensus sequence made!")
    else:
        with open(f"{outprefix}.final_assembly.fa", "w") as f:
            seq = pyfastaq.sequences.Fasta("assembly", consensus)
            print(seq, file=f)
    return consensus
