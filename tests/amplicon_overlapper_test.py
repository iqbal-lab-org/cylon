import collections
from unittest import mock
import os
import pytest

from cylon import amplicon_overlapper, utils
from cylon import amplicons as amps

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "amplicon_overlapper")

# This is effectively what is made by difflib.SequenceMatcher.find_longest_match()
Match = collections.namedtuple("Match", ("a", "b", "size"))


def test_remove_contained_and_bad_order_hits_from_ref_hits():
    f = amplicon_overlapper.remove_contained_and_bad_order_hits_from_ref_hits
    assert f([]) == []
    hit1 = mock.Mock()
    hit1.ref_start = 10000
    hit1.ref_end = 11000
    hit1.qry_start = 1000
    hit1.qry_end = 2000
    hit1.hit_length_qry = 1000
    assert f([hit1]) == [hit1]

    hit2 = mock.Mock()
    hit2.ref_start = 10
    hit2.ref_end = 20
    hit2.qry_start = 50
    hit2.qry_end = 60
    hit2.hit_length_qry = 10

    assert f([hit1, hit2]) == [hit1]
    assert f([hit2, hit1]) == [hit1]

    hit3 = mock.Mock()
    hit3.ref_start = 11010
    hit3.ref_end = 11020
    hit3.qry_start = 1998
    hit3.qry_end = 2008
    hit3.hit_length_qry = 10

    assert f([hit3, hit1, hit2]) == [hit1]

    hit3.qry_start = 2010
    hit3.qry_end = 2020
    assert f([hit3, hit1, hit2]) == [hit1, hit3]

    hit3.ref_start = 10995
    hit3.ref_end = 11005
    assert f([hit3, hit1, hit2]) == [hit1]


def test_amplicons_to_consensus_contigs():
    f = amplicon_overlapper.amplicons_to_consensus_contigs
    # ref is 100bp of random sequence
    #                10        20        30        40        50        60        70        80        90
    #      0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    ref = "GGCAACAAGCCCCGTAACCCAGCTCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGGCTCCTATGCACAGCGCGGACCAACA"
    ref_fa = "tmp.amplicons_to_consensus_contigs.ref.fa"
    minimap_out = "tmp.amplicons_to_consensus_contigs.paf"
    utils.rm_rf(minimap_out)

    with open(ref_fa, "w") as f_out:
        print(">ref", ref, sep="\n", file=f_out)

    amplicons = [
        amps.Amplicon("amp1", 9, 39, 1, 2),
        amps.Amplicon("amp2", 24, 75, 3, 4),
        amps.Amplicon("amp3", 63, 99, 5, 6),
    ]
    got_contigs = f(amplicons, ref_fa, minimap_out)
    assert got_contigs == None

    amplicons[0].final_seq = "CCCCGTAACCCAGCTCACCAGCGAATCACAAGT"
    amplicons[0].assemble_success = True
    got_contigs = f(amplicons, ref_fa, minimap_out)
    assert got_contigs == [amplicons[0].final_seq]

    amplicons[1].final_seq = "TCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGG"
    amplicons[1].assemble_success = True
    got_contigs = f(amplicons, ref_fa, minimap_out)
    expect = ["CCCCGTAACCCAGCTCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGG"]
    assert got_contigs == expect

    amplicons[2].final_seq = "CAGAACACTTTGGCTCCTATGCACAGCGCGGACCAA"
    amplicons[2].assemble_success = True
    got_contigs = f(amplicons, ref_fa, minimap_out)
    expect = [
        "CCCCGTAACCCAGCTCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGGCTCCTATGCACAGCGCGGACCAA"
    ]
    assert got_contigs == expect

    amplicons[1].final_seq = None
    amplicons[1].assemble_success = False
    got_contigs = f(amplicons, ref_fa, minimap_out)
    expect = [amplicons[0].final_seq, amplicons[2].final_seq]
    assert got_contigs == expect

    amplicons[0].final_seq = None
    amplicons[0].assemble_success = False
    got_contigs = f(amplicons, ref_fa, minimap_out)
    expect = [amplicons[2].final_seq]
    assert got_contigs == expect

    amplicons[0].final_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    amplicons[1].final_seq = "TCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGG"
    amplicons[2].final_seq = "CAGAACACTTTGGCTCCTATGCACAGCGCGGACCA"
    amplicons[0].assemble_success = False
    amplicons[1].assemble_success = False
    got_contigs = f(amplicons, ref_fa, minimap_out)
    expect = [amplicons[2].final_seq]
    assert got_contigs == expect
    os.unlink(ref_fa)


def test_amplicons_to_consensus_contigs_2():
    # This hits case not seen in previous test. Need a combination of amplicons
    # that pass fail pass fail pass. Was a bug where new contig was not being
    # started when we had two amplicons that didn't overlap and got removed
    # ref is 130bp of random sequence
    #                10        20        30        40        50        60        70        80        90        100       110       120
    #      0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    ref = "GGCAACAAGCCCCGTAACCCAGCTCACCAGCGAATCACAAGTGTTAAGAGACAAAGAAGCGGCAGAACACTTTGGCTCCTATGCACAGCGCGGACCAACAGTCGATCGTTATCCTAACTCTACATACTAAGTGAC"
    ref_fa = "tmp.amplicons_to_consensus_contigs_2.ref.fa"
    minimap_out = "tmp.amplicons_to_consensus_contigs_2.paf"
    utils.rm_rf(minimap_out)
    with open(ref_fa, "w") as f_out:
        print(">ref", ref, sep="\n", file=f_out)
    amplicons = [
        amps.Amplicon("amp1", 0, 34, 1, 1),
        amps.Amplicon("amp2", 20, 57, 1, 2),
        amps.Amplicon("amp3", 40, 70, 1, 1),
        amps.Amplicon("amp4", 60, 90, 1, 1),
        amps.Amplicon("amp5", 78, 119, 2, 1),
        amps.Amplicon("amp6", 100, 135, 1, 3),
    ]
    # got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(amplicons, ref_fa, minimap_out)
    # assert got_contigs == None
    amplicons[0].final_seq = ref[0:34]
    amplicons[0].assemble_success = True
    amplicons[1].final_seq = ref[20:57]
    amplicons[1].assemble_success = True
    amplicons[2].final_seq = ref[40:55] + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    amplicons[2].assemble_success = True
    amplicons[3].final_seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" + ref[75:90]
    amplicons[3].assemble_success = True
    amplicons[4].final_seq = ref[78:119]
    amplicons[4].assemble_success = True
    amplicons[5].final_seq = ref[100:135]
    amplicons[5].assemble_success = True
    got_contigs = amplicon_overlapper.amplicons_to_consensus_contigs(
        amplicons, ref_fa, minimap_out
    )
    assert got_contigs == [ref[0:57], ref[78:135]]
    os.unlink(ref_fa)


def test_split_contigs_on_gaps():
    contigs = ["A", "NNNNANANNCNNNGNNNNTNNNN"]
    got = amplicon_overlapper.split_contigs_on_gaps(contigs, gap_length=3)
    assert got == ["A", "ANANNC", "G", "T"]


def test_make_trimmed_contigs():
    contigs_in = [
        "ACGTG",
        "NNNAGTNNN",
        "NANANAN",
        "ANAATNA",
    ]
    got = amplicon_overlapper.make_trimmed_contigs(contigs_in, window=3)
    assert got == ["ACGTG", "AGT", "AAT"]


def test_consensus_contigs_to_consensus_1():
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

    # contigs in wrong order, should remove smallest wrong contig
    contigs = [contig2, contig1[2:]]  # contig 1 shorter than 2 now
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        contigs, ref_fasta, outprefix
    )
    assert got == contig2
    contigs = [contig2[2:], contig1]  # contig 2 shorter than 1 now
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        contigs, ref_fasta, outprefix
    )
    assert got == contig1
    utils.rm_rf(f"{outprefix}.*")

    # Add one short contig that should get removed because won't map well
    # enough to the ref
    contigs = [contig1, contig2[:10]]
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        contigs, ref_fasta, outprefix
    )
    assert got == contig1
    utils.rm_rf(f"{outprefix}.*")


def test_consensus_contigs_to_consensus_2():
    # Seen in real covid data. Stitching together two contigs that overlap
    # was deleting a chunk of sequence. The overlap is not perfect, so the
    # end of contig 1 and start of contig 2 sequences at the positions in the
    # minimap2 hits are not exactly the same.
    outprefix = "tmp.consensus_contigs_to_consensus_2"
    utils.rm_rf(f"{outprefix}.*")
    worse_overlap_seq = "GCGNGAGACNTGTCACTACAGTTTAAAGACCAATCCAAATNCTACTGACCAGTCTTCNTACATCGTTGATAGTGTTACAGTGANGAATGGTTCCATCC"
    better_overlap_seq = "GCGAGAGACTTGTCACTACAGTTTAAAAGACCAATAAATCCTACTGACCAGTCTTCTTACATCGTTGATAGTGTTACAGTGAAGAATGGTTCCATCC"
    contig1 = "GTCTTAGTGGTTTAGATTCTTTAGACACCTATCCTTCTTTAGAAACTATACAAATTACCATTTCATCTTTTAAATGGGATTTAACTGCTTTTGGCTTAGTTGCAGAGTGGTTTTTGGCATATATTCTTTTCACTAGGTTTTTCTATGTACTTGGATTGGCTGCAATCATGCAATTGTTTTTCAGCTATTTTGCAGTACATTTTATTAGTAATTCTTGGCTTATGTGGTTAATAATTAATCTTGTACAAATGGCCCCGATTTTAGCTATGGTTAGAATGTACATCTTCTTTGCATCATTTTATTATGTATGGAAAAGTTATGTGCATGTTGTAGACGGTTGTAATTCATNNACTTGTATGATGTGTTACAAACGTNATAGAGCAACAAGAGTCGAATGTACAACTATTGTTAATGGTGTTAGAAGGTCCTTTTATGTCTATGCTAATGGAGGTAAAGGCTTTTGCAAACTACACAATTGGAATTGTGTTAATTGTGATACATTCTGTGCTGGTAGTACATTTATTAGTGATGAAGT"
    contig2 = "ATCTTTACTTTGATAAAGCTGGTCAAAAGATTTATGAAAGACATTCTCTCTCTCATTTTGTTAACTTAGACAACCTGAGAGCTAATAACACTAAAGGTTCATTGCCTATTAATGTTATAGTTTTTGATGGTAAATCAAAATGTGAAGAATCATCTGCAAAATCAGCGTCTGTTTACTACAGTCAGCTTATGTGTCAACCTATACTGTTACTAGATCAGGCATTAGTGTCTGATGTTGGTGATAGTGCGGAAGTTGCAGTTAAAATGTTTGATGCTTACGTTAATACGTTTTCATCAACTTTTAACGTACCAATGGAAAAACTCAAAACACTAGTTTCAACTGCAGAAGCTGAATTTGCAAAGAATGTGTCCTTAGACAATGTCTTATCTACTTTTATTTCAGCAGCTCGGCAAGGGTTTGTTGATTCAGATGTAGAAACTAAAGATGTTGTT"
    expect = contig1 + better_overlap_seq + contig2
    ref_fasta = os.path.join(data_dir, "covid.ref.MN908947.3.fa")
    contigs = [contig1 + worse_overlap_seq, better_overlap_seq + contig2]
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        contigs, ref_fasta, outprefix, trim_end_window=25, debug=True
    )
    assert got == expect

    # Swap around which contig has the better/worse overlap sequence, should
    # get the same result. Worth testing because contigs are stitched together
    # as we move along the genome
    contigs = [contig1 + better_overlap_seq, worse_overlap_seq + contig2]
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        contigs, ref_fasta, outprefix
    )
    assert got == expect


def test_consensus_contigs_to_consensus_3():
    # Same as previous test, but the overlap is perfect, just for paranoia
    # checking
    outprefix = "tmp.consensus_contigs_to_consensus_3"
    utils.rm_rf(f"{outprefix}.*")
    overlap = "GCGAGAGACTTGTCACTACAGTTTAAAAGACCAATAAATCCTACTGACCAGTCTTCTTACATCGTTGATAGTGTTACAGTGAAGAATGGTTCCATCC"
    contig1 = "GTCTTAGTGGTTTAGATTCTTTAGACACCTATCCTTCTTTAGAAACTATACAAATTACCATTTCATCTTTTAAATGGGATTTAACTGCTTTTGGCTTAGTTGCAGAGTGGTTTTTGGCATATATTCTTTTCACTAGGTTTTTCTATGTACTTGGATTGGCTGCAATCATGCAATTGTTTTTCAGCTATTTTGCAGTACATTTTATTAGTAATTCTTGGCTTATGTGGTTAATAATTAATCTTGTACAAATGGCCCCGATTTTAGCTATGGTTAGAATGTACATCTTCTTTGCATCATTTTATTATGTATGGAAAAGTTATGTGCATGTTGTAGACGGTTGTAATTCATNNACTTGTATGATGTGTTACAAACGTNATAGAGCAACAAGAGTCGAATGTACAACTATTGTTAATGGTGTTAGAAGGTCCTTTTATGTCTATGCTAATGGAGGTAAAGGCTTTTGCAAACTACACAATTGGAATTGTGTTAATTGTGATACATTCTGTGCTGGTAGTACATTTATTAGTGATGAAGTT"
    contig2 = "ATCTTTACTTTGATAAAGCTGGTCAAAAGATTTATGAAAGACATTCTCTCTCTCATTTTGTTAACTTAGACAACCTGAGAGCTAATAACACTAAAGGTTCATTGCCTATTAATGTTATAGTTTTTGATGGTAAATCAAAATGTGAAGAATCATCTGCAAAATCAGCGTCTGTTTACTACAGTCAGCTTATGTGTCAACCTATACTGTTACTAGATCAGGCATTAGTGTCTGATGTTGGTGATAGTGCGGAAGTTGCAGTTAAAATGTTTGATGCTTACGTTAATACGTTTTCATCAACTTTTAACGTACCAATGGAAAAACTCAAAACACTAGTTTCAACTGCAGAAGCTGAATTTGCAAAGAATGTGTCCTTAGACAATGTCTTATCTACTTTTATTTCAGCAGCTCGGCAAGGGTTTGTTGATTCAGATGTAGAAACTAAAGATGTTGTT"
    contigs = [contig1 + overlap, overlap + contig2]
    ref_fasta = os.path.join(data_dir, "covid.ref.MN908947.3.fa")
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        contigs, ref_fasta, outprefix
    )
    expect = contig1 + overlap + contig2
    assert got == expect


def test_consensus_contigs_to_consensus_4():
    # Somtimes in real data we get a big overlap between two adjacent
    # contigs (esepcially if one of them is small). This is testing we
    # overlap correctly when we have a contig where its right hand end
    # mapping starts after the left end mapping of the next amplicon.
    # Since we stitch from left to right, adding one contig at a time,
    # making sure here that there's no indel error because coords/overlaps
    # positions could get messed up. Data copied+pasted from a real sample,
    # where the danger overlap is between right end of contig 1 and left end
    # of contig 2. Nucmer matches of the ends are:
    #  right_hit 22585	22650	2	67	66	66	100.00	29903	67	1	MN908947.3	1.right
    #   left_hit 22553	22689	1	137	137	137	89.78	29903	175	1	MN908947.3	2.left
    # which means start of 1.right hit is 22585, and start of 2.left hit is 22553.
    # Contig0 is unrelated and expect a gap between it and contig1. Left it in
    # this test anyway just for more paranoia checking
    outprefix = "tmp.consensus_contigs_to_consensus_4"
    utils.rm_rf(f"{outprefix}.*")
    contig0 = "TTGCCAATAGGTATTAACATCACTAGGTTTCAAACTTTACATAGAAGTTATTTGACTCCTGGTGATTCTTCTTCAGGTTGGACAGCTGGT"
    # ----contig2 starts here->
    # -------- 01234567890123456
    contig1 = "AGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATATTACAAACTTGTGCCCTTTTGGTGAAGNTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGAATCAGCAACTGTGTTGC"

    # first 118 bases of contig1 overlap contig 1
    # ------------------------------------------------------------------------------------------------------------contig2 ends here->
    contig2 = "AATCTANTGTTAGATTNNCNAATATTACAAACTTGTGCCCNTTTGGTGAAGTTTTTAANNCCACCAGATTTNNATCTGTTTATGCTTNGAACANGAAGAGAATCAGCAACTGTGNNNNTGATTATTCTGTNCNATATAANTCCGCATCATTTTCCACNTTNNANNGTTANGGNNNNNNNNCTACTAAATTNAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGG"
    ref_fasta = os.path.join(data_dir, "covid.ref.MN908947.3.fa")
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        [contig0, contig1, contig2], ref_fasta, outprefix
    )
    expect = contig0 + "N" * 171 + contig1[:-20] + contig2[98:]
    assert got == expect


def test_consensus_contigs_to_consensus_5():
    # Test when there is a short (nucmer won't find) overlap bewteen the contigs
    outprefix = "tmp.consensus_contigs_to_consensus_5"
    utils.rm_rf(f"{outprefix}.*")
    part_of_ref = "CGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACT"
    contig1 = part_of_ref[:100]
    contig2 = part_of_ref[96:200]
    ref_fasta = os.path.join(data_dir, "covid.ref.MN908947.3.fa")
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        [contig1, contig2], ref_fasta, outprefix
    )
    expect = part_of_ref[:200]
    assert got == expect

    # ref at 100-130 is GTGCTCGTACGTGGCTTTGGAGACTCCGTG
    contig1 = part_of_ref[:100] + "ATGCTCGAACGTGGTTGGACACTCCGTA"
    contig2 = "GTGCTCGTACGTGGCTTTGGAGACTCCGTG" + part_of_ref[131:]
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        [contig1, contig2], ref_fasta, outprefix
    )

    # when it tries to align the overlapping sequences, it gets this:
    # ACATGCTCGAACGTGG--TTGGACACTCCGT
    # --GTGCTCGTACGTGGCTTTGGAGACTCCGT
    # (the first AC is just extra, is ref 99-100, so don't include twice!
    #  and alo need the extra G at the end of middle bit from contig2 )
    expect = part_of_ref[:100] + "ATGCTCGTACGTGGCTTTGGAGACTCCGTG" + part_of_ref[131:]
    assert got == expect


def test_consensus_contigs_to_consensus_6():
    # Test when no nucmer overlap between contigs, but the unmapped parts of
    # their ends are longer than you would expect based on the nucmer coords.
    # Testing here that the contig ends get trimmed
    outprefix = "tmp.consensus_contigs_to_consensus_6"
    utils.rm_rf(f"{outprefix}.*")
    part_of_ref = "CGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACT"
    contig1 = part_of_ref[:100] + "A" * 10
    contig2 = "C" * 10 + part_of_ref[105:200]
    ref_fasta = os.path.join(data_dir, "covid.ref.MN908947.3.fa")
    got = amplicon_overlapper.consensus_contigs_to_consensus(
        [contig1, contig2], ref_fasta, outprefix
    )
    expect = part_of_ref[:100] + "N" * 5 + part_of_ref[105:200]
    assert got == expect


def test_assemble_amplicons():
    ref_fasta = os.path.join(data_dir, "assemble_amplicons.ref.fa")
    ref_seq = utils.load_single_seq_fasta(ref_fasta)
    amplicons = [
        amps.Amplicon("a1", 20, 300, 1, 2),
        amps.Amplicon("a2", 240, 550, 3, 4),
        amps.Amplicon("a3", 500, 850, 5, 6),
        amps.Amplicon("a4", 790, 970, 7, 8),
    ]
    outprefix = "tmp.assemble_amplicons"
    utils.rm_rf(f"{outprefix}.*")
    got = amplicon_overlapper.assemble_amplicons(
        amplicons, ref_fasta, outprefix, debug=True
    )
    assert got is None
    utils.rm_rf(f"{outprefix}.*")

    amplicons[0].final_seq = ref_seq[20:301]
    amplicons[0].assemble_success = True
    got = amplicon_overlapper.assemble_amplicons(
        amplicons, ref_fasta, outprefix, debug=True
    )
    assert got == amplicons[0].final_seq
    utils.rm_rf(f"{outprefix}.*")

    amplicons[1].final_seq = ref_seq[250:545]
    amplicons[1].assemble_success = True
    got = amplicon_overlapper.assemble_amplicons(
        amplicons, ref_fasta, outprefix, debug=True
    )
    assert got == ref_seq[20:545]
    utils.rm_rf(f"{outprefix}.*")

    amplicons[3].final_seq = ref_seq[790:952]
    amplicons[3].assemble_success = True
    got = amplicon_overlapper.assemble_amplicons(
        amplicons, ref_fasta, outprefix, debug=True
    )
    assert got == ref_seq[20:545] + "N" * 245 + ref_seq[790:952]
    utils.rm_rf(f"{outprefix}.*")
