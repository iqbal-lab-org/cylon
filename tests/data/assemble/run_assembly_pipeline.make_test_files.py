#!/usr/bin/env python3

from collections import namedtuple
import random
import pyfastaq


Amplicon = namedtuple("Amplicon", ("name", "start", "end"))
random.seed(42)

ref = random.choices(["A", "C", "G", "T"], k=1000)
ref[200] = "A"
ref_for_amplicons = ref[:900] + ref[901:]
ref_for_amplicons[200] = "T"

amplicons = [
    Amplicon("a1", 10, 400),
    Amplicon("a2", 350, 800),
    Amplicon("a3", 740, 990),
]

amp_seqs = ["".join(ref_for_amplicons[x.start:x.end]) for x in amplicons]

with open("run_assembly_pipeline.reads.fa", "w") as f:
    for i, seq in enumerate(amp_seqs):
        for read_count in range(20):
            to_print = pyfastaq.sequences.Fasta(f"{i}.{read_count}", seq)
            if read_count % 2 == 0:
                to_print.revcomp()
            #print(f">{i}.{read_count}", file=f)
            #print(seq, file=f)
            print(to_print, file=f)

with open("run_assembly_pipeline.ref.fa", "w") as f:
    print(">ref", file=f)
    print("".join(ref), file=f)

with open("run_assembly_pipeline.expect.fa", "w") as f:
    print(">expect", file=f)
    print("".join(ref_for_amplicons), file=f)


with open("run_assembly_pipeline.amplicons.bed", "w") as f:
    for a in amplicons:
        print(a.name, a.start, a.end, sep="\t", file=f)
