![Build Status](https://github.com/iqbal-lab-org/cylon/actions/workflows/build.yaml/badge.svg)

# cylon
Virus assembly module used by viridian

# Important
We recommend that you use
[Viridian](https://github.com/iqbal-lab-org/viridian) instead of
this repository. This repository is intended to be run by
Viridian, not to be used as a stand-alone tool.
Viridian provides a complete end-to-end pipeline for
generating a consensus sequence from reads.

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
sudo singularity build cylon.img Singularity.def
```
to build the container `cylon.img`.


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
    "left_primer_end": 30,
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
cylon assemble --bam reads.bam ont ref.fasta amplicons.json outdir
```

Run using a FASTQ file of ONT reads:
```
cylon assemble --reads_to_map reads.fastq ont ref.fasta amplicons.json outdir
```

Run using two FASTQ files of paired Illumina reads:
```
cylon assemble \
  --reads_to_map reads1.fastq --mates_to_map reads2.fastq \
  illumina ref.fasta amplicons.json outdir
```


## Output

The important files are:
* `consensus.final_assembly.fa`: this contains the consensus sequence.
* `amplicon_data.json`: JSON file containing details of what happened when
  trying to make a consensus sequence of each amplicon.

