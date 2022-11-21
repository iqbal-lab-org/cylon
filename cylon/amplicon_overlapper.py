import logging
import os
import re

import pyfastaq
from cylon import utils


def _make_split_contigs_fasta(contigs, outfile):
    with open(outfile, "w") as f:
        for d in contigs:
            length = min(250, int(len(d["seq"]) / 2))
            print(f">{d['name']}.left", file=f)
            print(d["seq"][:length], file=f)
            print(f">{d['name']}.right", file=f)
            print(d["seq"][-length:], file=f)


def _map_split_contigs(to_map_fasta, ref_fasta, end_allowance=20):
    command = f"minimap2 -t 1 -c {ref_fasta} {to_map_fasta}"
    minimap2_out = utils.syscall(command)
    paf_lines = minimap2_out.stdout.strip().split("\n")
    mappings = {}
    for line in paf_lines:
        logging.debug(f"mapping split contigs. paf_line: {line.rstrip()}")
        # if there were no mappings we can get an empty line
        if line == "":
            continue
        fields = line.rstrip().split()
        if fields[4] != "+":
            continue

        name, left_or_right = fields[0].split(".")
        if (left_or_right == "left" and int(fields[2]) > end_allowance) or (
            left_or_right == "right" and int(fields[3]) + end_allowance < int(fields[3])
        ):
            continue

        key = (name, left_or_right)

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
    return bad_contig_indexes


def minimap2_hit_to_nm(hit):
    nm = [x for x in hit if x.startswith("NM:i:")]
    if len(nm) != 1:
        return None
    try:
        nm, i, nm_count = nm[0].split(":")
        return int(nm_count)
    except:
        return None


def split_contigs_on_gaps(contigs_in, gap_length=10):
    contigs_out = []
    regex = re.compile("N{" + str(gap_length) + ",}")
    for contig in contigs_in:
        contigs_out.extend([x for x in re.split(regex, contig) if len(x)>0])
    return contigs_out


def make_trimmed_contigs(contigs_in, window=30):
    contigs_out = []
    contigs_after_gaps = split_contigs_on_gaps(contigs_in)
    for i, seq in enumerate(contigs_after_gaps):
        start = 0
        while start < len(seq) - window and seq[start:start +window].count("N") > 0:
            start += 1
        end = len(seq)
        while end > window and seq[end-window:end].count("N") > 0:
            end -= 1
        if start < end:
            logging.debug(f"contig trim {i}, {len(seq)}, {start}, {end}, {seq}")
            contigs_out.append({"name": str(i), "seq": seq[start:end]})
            logging.debug(f"contig trimmed {i}, {len(seq)}, {start}, {end}, {seq[start:end]}")
    return contigs_out


def consensus_contigs_to_consensus(
    contigs, ref_fasta, outprefix, map_end_allowance=20, debug=False, trim_end_window=30
):
    if contigs is None or len(contigs) == 0:
        logging.warning("No contigs were made. Aborting assembly")
        return None

    fa_to_map = f"{outprefix}.to_map.fa"
    contigs = make_trimmed_contigs(contigs, window=trim_end_window)
    _make_split_contigs_fasta(contigs, fa_to_map)
    mappings = _map_split_contigs(fa_to_map, ref_fasta, end_allowance=map_end_allowance)
    if not debug:
        os.unlink(fa_to_map)
    bad_contig_indexes = _check_mappings(contigs, mappings)
    if len(bad_contig_indexes) > 0:
        contigs = [
            contigs[i] for i in range(len(contigs)) if i not in bad_contig_indexes
        ]
    if len(contigs) == 0:
        logging.warning(
            "Errors aligning contigs to reference to make final sequence. Aborting assembly"
        )
        return None

    consensus = []
    trim_next = 0
    bad_overlap_ns = 10

    for i, contig in enumerate(contigs):
        right_fields = mappings[(contig["name"], "right")]
        consensus.append(contig["seq"][trim_next:])

        if i < len(contigs) - 1:
            next_left_fields = mappings[(contigs[i + 1]["name"], "left")]
            contig_start_in_ref = int(right_fields[7])
            contig_end_in_ref = int(right_fields[8])
            next_contig_start_in_ref = int(next_left_fields[7])
            # Usually we don't expect the contigs to overlap. In which case put
            # Ns between them. But they can sometime overlap. When they do,
            # look for an overlap between them based on the minimap2 mapping.
            # Completely stop if the contigs don't map in the correct
            # order compared to the reference
            if next_contig_start_in_ref < contig_start_in_ref:
                consensus = ""
                break

            # If the contigs do not overlap each other based on ref coords
            if not (contig_start_in_ref < next_contig_start_in_ref < contig_end_in_ref):
                consensus.append("N" * (next_contig_start_in_ref - contig_end_in_ref))
                trim_next = 0
            else:
                overlap_len = contig_end_in_ref - next_contig_start_in_ref
                end_in_contig = len(contig["seq"]) - (
                    int(right_fields[1]) - int(right_fields[3])
                )
                contig_ol_seq = contig["seq"][
                    end_in_contig - overlap_len : end_in_contig
                ]
                next_start = int(next_left_fields[2])
                next_ol_seq = contigs[i + 1]["seq"][
                    next_start : next_start + overlap_len
                ]
                if contig_ol_seq == next_ol_seq:  # good overlap
                    consensus[-1] = contig["seq"][
                        trim_next : end_in_contig - overlap_len
                    ]
                    trim_next = 0
                else:  # non-perfect overlap, pick seq with fewest mismatches
                    this_nm = minimap2_hit_to_nm(right_fields)
                    next_nm = minimap2_hit_to_nm(next_left_fields)
                    if this_nm is None or next_nm is None:
                        consensus = ""
                        break
                    best = "this" if this_nm < next_nm else "next"
                    if this_nm < next_nm:
                        consensus[-1] = contig["seq"][ trim_next : end_in_contig] # FIXME
                        trim_next = overlap_len - 1
                    else:
                        consensus[-1] = contig["seq"][ trim_next : end_in_contig - overlap_len] # FIXME
                        trim_next = 0


    consensus = "".join([x for x in consensus if len(x) > 0])
    if any([x != "N" for x in consensus]):
        return consensus
    else:
        logging.warning(
            "Errors aligning contigs to reference to make final sequence. No consensus made"
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
