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
2. BED file of amplicons. Column 1 = name of amplicon, columns 2 and 3 are start
   end end positions of the amplicons
3. Reads, either in a sorted mapped indexed BAM file, or in a FASTA/FASTQ file
(or two FASTA/FASTQ files for paired reads).


Run using a mapped BAM file of ONT reads:
```
viridian assemble --bam reads.bam ref.fasta amplicons.bed outdir
```

Run using a FASTQ file of ONT reads:
```
viridian assemble --reads_to_map reads.fastq ref.fasta amplicons.bed outdir
```

Run using two FASTQ files of paired Illumina reads:
```
viridian assemble \
  --reads_to_map reads1.fastq --mates_to_map reads2.fastq \
  --minimap_opts "-t 1 -x sr" \
  --min_read_length 50 \
  ref.fasta amplicons.bed outdir
```


## Output

The important files are:
* `consensus.final_assembly.fa`: this contains the consensus sequence.
* `amplicon_data.json`: JSON file containing details of what happened when
  trying to make a consensus sequence of each amplicon.

