#!/usr/bin/env python3

from collections import namedtuple
import json
import random
import pyfastaq

from viridian import amplicons as amps

Amplicon = namedtuple("Amplicon", ("name", "start", "end"))
random.seed(42)

ref = random.choices(["A", "C", "G", "T"], k=1000)
ref[200] = "A"
ref_for_amplicons = ref[:900] + ref[901:]
ref_for_amplicons[200] = "T"

amplicons_json_data = {
  "a1": {
    "start": 10,
    "end": 399,
    "left_primer_end_max": 10,
    "right_primer_start_min": 390
  },
  "a2": {
    "start": 350,
    "end": 799,
    "left_primer_end_max": 355,
    "right_primer_start_min": 790
  },
  "a3": {
    "start": 740,
    "end": 989,
    "left_primer_end_max": 745,
    "right_primer_start_min": 980
  }
}


amplicons = [
    amps.Amplicon("a1", 10, 399, 0, 10),
    amps.Amplicon("a2", 350, 799, 6, 10),
    amps.Amplicon("a3", 740, 989, 6, 10),
]

amp_seqs = ["".join(ref_for_amplicons[x.start:x.end + 1]) for x in amplicons]

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
    print(">ref foo", file=f)
    print("".join(ref), file=f)

with open("run_assembly_pipeline.expect.fa", "w") as f:
    print(">expect", file=f)
    print("".join(ref_for_amplicons), file=f)

with open("run_assembly_pipeline.amplicons.json", "w") as f:
    json.dump(amplicons_json_data, f, indent=2)
