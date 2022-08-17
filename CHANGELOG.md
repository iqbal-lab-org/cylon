# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Backwards incompatible change: now uses primer coordinates. The amplicon
  (and primer) information must be given in a JSON file, instead of the previous
  BED file.

- Added command line option `--reads_per_amp_dir`

- When building consensus from contigs, if a contig does not map then remove
  the contig and carry on, instead of aborting the entire assembly

- Fail amplicon if its sequence after stripping Ns is less than 3/4 of
  (length of the amplicon - primer lengths)

### Removed

- Removed the command line function `viridian amplicon_overlap`.

### Fixed

- Handle whitespace in sequence name line of input FASTA file

- Bug fix when using the option `--minimap_opts`. It was used for initial read
  mapping, but not for each racon iteration. Now used for racon iterations as
  well.

- Fixed case viridian would crash when no amplicon contigs map, and so no
  consensus sequence was made. Now it finishes, with `"made_consensus": false`
  and `"consensus": null` in the output JSON file.

- Fixed where Racon could sometimes silently fail due to how reads were sampled
  for an amplicon, resulting in no sequence for that amplicon.


[Unreleased]: https://github.com/iqbal-lab-org/viridian/compare/v0.1.0...HEAD
