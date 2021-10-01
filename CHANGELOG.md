# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- Handle whitespace in sequence name line of input FASTA file

- Bug fix when using the option `--minimap_opts`. It was used for initial read
  mapping, but not for each racon iteration. Now used for racon iterations as
  well.


[Unreleased]: https://github.com/iqbal-lab-org/viridian/compare/v0.10.0...HEAD
