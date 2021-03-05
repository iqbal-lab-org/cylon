import collections
import os
import pytest

from viridian import amplicon_overlapper, utils
from viridian import amplicons as amps

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "amplicon_overlapper")

# This is effectively what is made by difflib.SequenceMatcher.find_longest_match()
Match = collections.namedtuple("Match", ("a", "b", "size"))


def test_get_amplicon_overlaps():
    amplicons = [
        amps.Amplicon("a1", 10, 100),
        amps.Amplicon("a2", 10, 100),
        amps.Amplicon("a3", 110, 142),
    ]
    amplicons[1].masked_seq = "AAAAAAAAAAAAAAAAAATGCTGAACAGTCCCCCCC"
    amplicons[2].masked_seq = "CCTGCTGAACGGTTGATGCATCTCATGCTGACNNAGGTGTGGCCAAAAA"

    expect_overlaps = [None, Match(18, 2, 8)]
    got_overlaps = amplicon_overlapper.get_amplicon_overlaps(amplicons, 8)
    assert got_overlaps == expect_overlaps
    assert amplicon_overlapper.get_amplicon_overlaps(amplicons, 9) == [None, None]


def test_amplicons_to_consensus_contigs():
    # ref is 100bp of random sequence
    #                10        20        30        40        50        60        70        80        90
    #      0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    # ref = "GGCAACAAGCCCCGTAACCCAGCTCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGGCTCCTATGCACAGCGCGGACCAACA"
    amplicons = [
        amps.Amplicon("amp1", 9, 39),
        amps.Amplicon("amp2", 24, 75),
        amps.Amplicon("amp3", 63, 99),
    ]
    got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(amplicons, 20)
    assert got_contigs == None

    amplicons[0].masked_seq = "CCCCGTAACCGAGCTCACCAGCGAATCACAA"
    amplicons[0].assemble_success = True
    got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(amplicons, 10)
    assert got_contigs == [amplicons[0].masked_seq]

    amplicons[1].masked_seq = "TCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGG"
    amplicons[1].assemble_success = True
    got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(amplicons, 10)
    expect = ["CCCCGTAACCGAGCTCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGG"]
    assert got_contigs == expect

    amplicons[2].masked_seq = "NNCAGAACACTTTGGCTCCTATGCACAGCGCGGACCAANN"
    amplicons[2].assemble_success = True
    got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(amplicons, 10)
    expect = [
        "CCCCGTAACCGAGCTCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGGCTCCTATGCACAGCGCGGACCAA"
    ]
    assert got_contigs == expect

    amplicons[1].masked_seq = None
    amplicons[1].assemble_success = False
    got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(amplicons, 10)
    expect = [amplicons[0].masked_seq, amplicons[2].masked_seq.strip("N")]
    assert got_contigs == expect

    amplicons[0].masked_seq = None
    amplicons[0].assemble_success = False
    got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(amplicons, 10)
    expect = [amplicons[2].masked_seq.strip("N")]
    assert got_contigs == expect

    amplicons[0].masked_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    amplicons[1].masked_seq = "TCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGG"
    amplicons[2].masked_seq = "NNCAGAACACTTTGGCTCCTATGCACAGCGCGGACCAANN"
    amplicons[0].assemble_success = False
    amplicons[1].assemble_success = False
    amplicons[2].assemble_success = False
    expect = [amplicons[2].masked_seq.strip("N")]
    assert got_contigs == expect

def test_amplicons_to_consensus_contigs_2():
    # This hits case not seen in previous test. Need a combination of amplicons
    # that pass fail pass fail pass. Was a bug where new contig was not being
    # started when we had two amplicons that didn't overlap and got removed
    # ref is 130bp of random sequence
    #                10        20        30        40        50        60        70        80        90        100       110       120
    #      0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    ref = "GGCAACAAGCCCCGTAACCCAGCTCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGGCTCCTATGCACAGCGCGGACCAACAGTCGATCGTTATCCTAACTCTACATACTAA"
    amplicons = [
        amps.Amplicon("amp1", 0, 30),
        amps.Amplicon("amp2", 20, 50),
        amps.Amplicon("amp3", 40, 70),
        amps.Amplicon("amp4", 60, 90),
        amps.Amplicon("amp5", 80, 110),
        amps.Amplicon("amp6", 100, 130),
    ]
    got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(amplicons, 7)
    assert got_contigs == None
    amplicons[0].masked_seq = ref[0:30]
    amplicons[0].assemble_success = True
    amplicons[1].masked_seq = ref[20:50]
    amplicons[1].assemble_success = True
    amplicons[2].masked_seq = ref[40:55] + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    amplicons[2].assemble_success = True
    amplicons[3].masked_seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" + ref[75:90]
    amplicons[3].assemble_success = True
    amplicons[4].masked_seq = ref[80:110]
    amplicons[4].assemble_success = True
    amplicons[5].masked_seq = ref[100:130]
    amplicons[5].assemble_success = True
    got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(amplicons, 7)
    assert got_contigs == [ref[0:50], ref[80:130]]

def test_consensus_contigs_to_consensus():
    ref_fasta = os.path.join(data_dir, "consensus_contigs_to_consensus.fa")
    outprefix = "tmp.consensus_contigs_to_consensus"
    utils.rm_rf(f"{outprefix}.*")
    assert (
        amplicon_overlapper.consensus_contigs_to_consensus(None, ref_fasta, outprefix)
        is None
    )
    assert (
        amplicon_overlapper.consensus_contigs_to_consensus([], ref_fasta, outprefix)
        is None
    )
    # contig is in ref from 1-120
    contig1 = "GGGTCCTCGGCCTACGACTATATCGCATGGCACGGTGCGGCTGTAGGGACACAAGATAATGTTCCGAGCAATTACGCACTTATTTGGTTCAGGAATCAGACTTCCGGTTTCGAACTTTCG"
    contigs = [contig1]
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        contigs, ref_fasta, outprefix
    )
    assert got == contig1
    utils.rm_rf(f"{outprefix}.*")

    # contig2 is in ref from 181-300
    contig2 = "CTATTTGCACCGTTGTAAATGCGCAGTTTGAGCTGTTGTTTCGCGGCACCGTAAGAAAAAAGATGTACTGCCGAACTCGGGGCGTAGTGAGGGGTTCATAGCGAGAAACGTCTTGTACGC"
    contigs = [contig1, contig2]
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        contigs, ref_fasta, outprefix
    )
    assert got == contig1 + "N" * 60 + contig2
    utils.rm_rf(f"{outprefix}.*")

    # contigs in wrong order, should result in aborted assembly
    contigs = [contig2, contig1]
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        contigs, ref_fasta, outprefix
    )
    assert got is None
    utils.rm_rf(f"{outprefix}.*")


def test_assemble_amplicons():
    ref_fasta = os.path.join(data_dir, "assemble_amplicons.ref.fa")
    ref_seq = utils.load_single_seq_fasta(ref_fasta)
    amplicons = [
        amps.Amplicon("a1", 20, 300),
        amps.Amplicon("a2", 240, 550),
        amps.Amplicon("a3", 500, 850),
        amps.Amplicon("a4", 790, 970),
    ]
    outprefix = "tmp.assemble_amplicons"
    utils.rm_rf(f"{outprefix}.*")
    got = amplicon_overlapper.assemble_amplicons(
        amplicons, ref_fasta, outprefix, debug=True
    )
    assert got is None
    utils.rm_rf(f"{outprefix}.*")

    amplicons[0].masked_seq = ref_seq[20:301]
    amplicons[0].assemble_success = True
    got = amplicon_overlapper.assemble_amplicons(
        amplicons, ref_fasta, outprefix, debug=True
    )
    assert got == amplicons[0].masked_seq
    utils.rm_rf(f"{outprefix}.*")

    amplicons[1].masked_seq = ref_seq[250:545]
    amplicons[1].assemble_success = True
    got = amplicon_overlapper.assemble_amplicons(
        amplicons, ref_fasta, outprefix, debug=True
    )
    assert got == ref_seq[20:545]
    utils.rm_rf(f"{outprefix}.*")

    amplicons[3].masked_seq = ref_seq[790:960]
    amplicons[3].assemble_success = True
    got = amplicon_overlapper.assemble_amplicons(
        amplicons, ref_fasta, outprefix, debug=True
    )
    assert got == ref_seq[20:545] + "N" * 245 + ref_seq[790:960]
    utils.rm_rf(f"{outprefix}.*")

    # putting in junk for amplicon 2 means it won't overlap amplicons 1 or 3,
    # and we should only get amplicon 0 back
    amplicons[
        2
    ].masked_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTT"
    amplicons[2].assemble_success = True
    got = amplicon_overlapper.assemble_amplicons(
        amplicons, ref_fasta, outprefix, debug=True
    )
    assert got == ref_seq[20:301]
    utils.rm_rf(f"{outprefix}.*")
