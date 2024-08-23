[![Rust](https://github.com/stela2502/bam_aligner/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/stela2502/bam_aligner/actions/workflows/rust.yml)

# bam_aligner - a tool to construct detected transcripts from 10x bam files

VDJ recombination creates unique Antiboty and T cell receptor genes for each T and B cell.
These transcripts are therefore highly complicated and a valuable research target.

This program tries to collect data from the detected single cells and collides the reads into likely transcripts.

This is work in progress!




# Installation

You need the Rust compiler: https://www.rust-lang.org/tools/install

```
git clone https://github.com/stela2502/bam_aligner
cd bam_aligner
cargo build --release
```

Or:

```
cargo install --git https://github.com/stela2502/bam_aligner
```

# Usage

```
bam_aligner  --help
Simple program to assemble a contig sequence from a BAM file

Usage: bam_aligner [OPTIONS] --bam <BAM> --outfile <OUTFILE>

Options:
  -b, --bam <BAM>                  BAM file path
  -o, --outfile <OUTFILE>          Output FASTA file name
  -n, --num-threads <NUM_THREADS>  number of threads to use
  -h, --help                       Print help
  -V, --version                    Print version
```

Hence you can use it like that:

```
samtools view -hb <source 10x bam file> <chromoisomal area with the genes you are interested in> > bam_subset.bam
bam_aligner --bam bam_subset.bam --outfile ChrRegionTranscripts.fa
```