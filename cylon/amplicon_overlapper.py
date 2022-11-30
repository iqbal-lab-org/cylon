import logging
from operator import attrgetter, itemgetter
import os
import re

import pyfastaq
import pymummer
from Bio import pairwise2

from cylon import utils


def global_align(seq1, seq2):
    """Returns global alignment strings from N-W alignment of the
    two sequences. Dashes for gaps"""
    alignments = pairwise2.align.globalms(
        seq1,
        seq2,
        1,  # match
        -1,  # mismatch
        -5,  # gap_open
        -3,  # gap_extend
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
            contigs_out.append(seq[start:end])
            logging.debug(
                f"contig trimmed {i}, {len(seq)}, {start}, {end}, {seq[start:end]}"
            )

    return contigs_out


def contigs_list_to_fasta(contigs, outfile):
    with open(outfile, "w") as f:
        for i, seq in enumerate(contigs):
            print(f">{i}", seq, sep="\n", file=f)


def remove_contained_and_bad_order_hits_from_ref_hits(hits):
    if len(hits) <= 1:
        return hits

    # Sort hits in order of longest to shortest. Keep each new hit if
    # it doesn't overlap an existing hit
    hits.sort(reverse=True, key=attrgetter("hit_length_qry"))
    new_hits = [hits[0]]
    for new_hit in hits[1:]:
        can_add = True
        for other_hit in new_hits:
            # if this new hit overlaps other_hit
            if (
                other_hit.qry_start <= new_hit.qry_end
                and new_hit.qry_start <= other_hit.qry_end
            ):
                can_add = False
                break

        if can_add:
            new_hits.append(new_hit)

    new_hits.sort(key=attrgetter("qry_start"))

    # Now have a list of non-overlapping longest hits. Remove hits that have
    # a different reference order or inconsistent gap between contig coords
    # vs ref coords
    i = 0
    while i < len(new_hits) - 1:
        need_to_pop = False

        if new_hits[i + 1].ref_start <= new_hits[i].ref_end:
            need_to_pop = True
        else:
            qry_gap = new_hits[i + 1].qry_start - new_hits[i].qry_end
            assert qry_gap > 0
            ref_gap = new_hits[i + 1].ref_start - new_hits[i].ref_end
            need_to_pop = abs(qry_gap - ref_gap) > 100

        if need_to_pop:
            if new_hits[i].hit_length_qry > new_hits[i + 1].hit_length_qry:
                new_hits.pop(i + 1)
            else:
                new_hits.pop(i)
        else:
            i += 1

    return new_hits


def map_contigs_to_ref(ref_fasta, contigs_fa, outfile):
    runner = pymummer.nucmer.Runner(
        ref_fasta,
        contigs_fa,
        outfile,
        mincluster=5,
        breaklen=500,
        maxmatch=True,
    )
    runner.run()
    mappings = {}
    for hit in pymummer.coords_file.reader(outfile):
        logging.debug(f"nucmer contigs vs ref: {hit}")
        if not hit.ref_start >= hit.ref_end and hit.qry_start >= hit.qry_end:
            continue

        contig_index = int(hit.qry_name)
        if contig_index not in mappings:
            mappings[contig_index] = []
        mappings[contig_index].append(hit)

    mapping_list = []
    for contig_index, hits in mappings.items():
        hits = remove_contained_and_bad_order_hits_from_ref_hits(hits)
        first = sorted(hits, key=attrgetter("qry_start"))[0]
        last = sorted(hits, key=attrgetter("qry_end"))[-1]
        mapping_list.append(
            {
                "contig_index": contig_index,
                "start": first.ref_start,
                "end": last.ref_end,
                "hits": [first, last],
            }
        )
        logging.debug(f"First nucmer match for contig {contig_index}: {first}")
        logging.debug(f"Last nucmer match for contig {contig_index}: {last}")

    mapping_list.sort(key=itemgetter("contig_index"))
    return mapping_list


def remove_out_of_order_or_contained_mappings(mappings):
    i = 0
    while i < len(mappings) - 1:
        this_start = mappings[i]["start"]
        this_end = mappings[i]["end"]
        next_start = mappings[i + 1]["start"]
        next_end = mappings[i + 1]["end"]

        if this_start <= next_start <= next_end <= this_end:
            mappings.pop(i + 1)
            continue
        elif next_start <= this_start <= this_end <= next_end:
            mappings.pop(i)
            continue

        if next_start < this_start or this_end > next_end:
            if this_end - this_start > next_end - next_start:
                mappings.pop(i + 1)
            else:
                mappings.pop(i)
            continue

        i += 1

    return mappings


def self_map_contigs(contigs_fa, outfile, end_allow=20):
    runner = pymummer.nucmer.Runner(
        contigs_fa,
        contigs_fa,
        outfile,
        mincluster=5,
        maxmatch=True,
    )
    runner.run()
    mappings = {}
    for hit in pymummer.coords_file.reader(outfile):
        if (
            hit.ref_start >= hit.ref_end
            or hit.qry_start >= hit.qry_end
            or not int(hit.ref_name) < int(hit.qry_name)
            or hit.ref_end + end_allow < hit.ref_length
            or hit.qry_start > end_allow
        ):
            continue

        key = (int(hit.ref_name), int(hit.qry_name))
        if key not in mappings or hit.hit_length_qry > mappings[key].hit_length_qry:
            mappings[key] = hit

        logging.debug(f"nucmer mapping contigs to each other {hit}")

    return mappings


def consensus_contigs_to_consensus(
    contigs, ref_fasta, outprefix, map_end_allowance=20, debug=False, trim_end_window=20
):
    if contigs is None or len(contigs) == 0:
        logging.warning("No contigs were made. Aborting assembly")
        return None

    contigs = make_trimmed_contigs(contigs, window=trim_end_window)
    contigs_fa = f"{outprefix}.tmp.trimmed_contigs.fa"
    nucmer_contigs_v_ref = f"{outprefix}.tmp.trimmed_contigs_v_ref.coords"
    nucmer_contigs_v_self = f"{outprefix}.tmp.trimmed_contigs_v_self.coords"
    contigs_list_to_fasta(contigs, contigs_fa)
    ref_mappings = map_contigs_to_ref(ref_fasta, contigs_fa, nucmer_contigs_v_ref)
    ref_mappings = remove_out_of_order_or_contained_mappings(ref_mappings)
    self_mappings = self_map_contigs(contigs_fa, nucmer_contigs_v_self)
    if not debug:
        os.unlink(nucmer_contigs_v_ref)
        os.unlink(nucmer_contigs_v_self)
        os.unlink(contigs_fa)
    consensus = []
    to_trim = 0

    for map_index, this_map in enumerate(ref_mappings):
        this_ctg_index = this_map["contig_index"]
        this_contig = contigs[this_ctg_index]
        if map_index >= len(ref_mappings) - 1:  # add last piece of final contig
            consensus.append(this_contig[to_trim:])
            break

        next_map = ref_mappings[map_index + 1]
        next_ctg_index = next_map["contig_index"]
        next_contig = contigs[next_ctg_index]
        self_map_key = (this_ctg_index, next_ctg_index)
        dist_between_contigs = next_map["start"] - this_map["end"]

        # If nucmer found a self-match between this contig and the next
        # one then it gets priority. Also sanity check that we don't expect
        # them to be too far away based on reference mapping
        if dist_between_contigs <= 50 and self_map_key in self_mappings:
            self_hit = self_mappings[self_map_key]
            consensus.append(this_contig[to_trim : self_hit.ref_start])
            this_ns = this_contig[self_hit.ref_start : self_hit.ref_end + 1].count("N")
            next_ns = next_contig[self_hit.qry_start : self_hit.qry_end + 1].count("N")

            if this_ns < next_ns:
                consensus.append(this_contig[self_hit.ref_start : self_hit.ref_end])
                to_trim = self_hit.qry_end
            else:
                to_trim = self_hit.qry_start
        elif this_map["end"] < next_map["start"]:  # no overlap from ref coords
            ref_gap_len = next_map["start"] - this_map["end"] - 1
            this_left_qry_end = this_map["hits"][-1].qry_end
            left_extra = len(this_contig) - this_left_qry_end
            right_extra = next_map["hits"][0].qry_start
            # If the total bases we're adding on from the unmapped contig ends
            # is <= the gap legnth in the reference, then this is ok and we can
            # add all of the contig ends
            if left_extra + right_extra <= ref_gap_len:
                consensus.append(this_contig[to_trim:])
                consensus.append("N" * ref_gap_len)
                to_trim = 0
            else:  # trying to add too many bases from the contig ends. Trim them
                consensus.append(this_contig[to_trim : this_left_qry_end + 1])
                consensus.append("N" * ref_gap_len)
                to_trim = right_extra
        else:  # mapping coords to ref say the contigs overlap
            # Get the sequence from this and next contig that overlap
            overlap_len = this_map["end"] - next_map["start"] + 1
            this_ol_end = this_map["hits"][-1].qry_end
            this_ol_start = this_ol_end - overlap_len + 1
            this_ol_seq = this_contig[this_ol_start : this_ol_end + 1]
            next_ol_start = next_map["hits"][0].qry_start
            next_ol_end = next_ol_start + overlap_len - 1
            next_ol_seq = next_contig[next_ol_start : next_ol_end + 1]

            # Nice case: the two sequences are the same
            if this_ol_seq == next_ol_seq:
                consensus.append(this_contig[to_trim : this_ol_end + 1])
            else:  # different overlapping sequences
                match = global_align(this_ol_seq, next_ol_seq)
                if match is None:
                    consensus.append(this_contig[to_trim:this_ol_start])
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
                        consensus.append(this_contig[to_trim:this_ol_start])
                        consensus.append("N" * overlap_len)
                    else:
                        new_seq = (match[0][:i] + match[1][i:]).replace("-", "")
                        consensus.append(this_contig[to_trim:this_ol_start] + new_seq)

            to_trim = next_ol_end + 1

    consensus = "".join([x for x in consensus if len(x) > 0])
    if any([x != "N" for x in consensus]):
        return consensus
    else:
        logging.warning("Error making consensus, did not get any non-N bases")
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
