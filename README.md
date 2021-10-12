[![Build Status](https://www.travis-ci.com/iqbal-lab-org/viridian.svg?branch=main)](https://www.travis-ci.com/iqbal-lab-org/viridian)

# viridian
Virus assembler from amplicon sequencing reads

## Install

### From source
These must be installed and in your `$PATH`:
* `racon` (https://github.com/lbcb-sci/racon)
* `minimap2` (https://github.com/lh3/minimap2/).


Clone this repository and run:
```
python3 -m pip install .
```

### Singularity container
Clone this repository and run:
```
sudo singularity build viridian.img Singularity.def
```
to build the container `viridian.img`.


## Example usage

Required input:
1. Reference FASTA file
2. JSON file of amplicons. Described below.
   end end positions of the amplicons
3. Reads, either in a sorted mapped indexed BAM file, or in a FASTA/FASTQ file
(or two FASTA/FASTQ files for paired reads).

The amplicons JSON file must look like this:
```
{
  "amplicon1": {
    "start": 10,
    "end": 399,
    "left_primer_end": 10,
    "right_primer_start": 390
  },
  "amplicon2": { ... etc ...}
}
```
The keys are the amplicon names, and the values are the details for each
amplicon.
All coordinates are 0-based inclusive.
The `start` and `end` entries are the positions of the start and end of the
amplicon.
`left_primer_end` is the rightmost position of the end of the left primer,
and `right_primer_start` is the leftmost position
of the right primer. This means for each amplicon we should have:
`start` < `left_primer_end` < `right_primer_start` < `end`.
(Other key/values can be inside the dictionary
for each amplicon, but will simply be ignored).


Run using a mapped BAM file of ONT reads:
```
viridian assemble --bam reads.bam ont ref.fasta amplicons.json outdir
```

Run using a FASTQ file of ONT reads:
```
viridian assemble --reads_to_map reads.fastq ont ref.fasta amplicons.json outdir
```

Run using two FASTQ files of paired Illumina reads:
```
viridian assemble \
  --reads_to_map reads1.fastq --mates_to_map reads2.fastq \
  illumina ref.fasta amplicons.json outdir
```


## Output

The important files are:
* `consensus.final_assembly.fa`: this contains the consensus sequence.
* `amplicon_data.json`: JSON file containing details of what happened when
  trying to make a consensus sequence of each amplicon.

