#!/usr/bin/env python3

import argparse
from pathlib import Path
import gzip
import pysam

def main():
    options = get_options()

    bam_files = parse_bam_config(options.config)
    transcript_sequences = read_transcripts(options.fasta)

    all_results = {}

    for bam_file in bam_files:
        results = parse_bam_file(bam_file, transcript_sequences)

        all_results[bam_file["filename"]] = results


    write_results(all_results,options.outfile)


def parse_bam_config (file):

    bam_files = {}

    with open(file,"rt",encoding="utf8") as infh:
        for line in infh:
            if line.startswith("#"):
                continue

            sections = line.split("\t")

            bamfile = sections[0]

            if not Path(bamfile).exists():
                raise Exception(f"BAM file {bamfile} doesn't exist")
            
            if not bamfile in bam_files:
                bam_files[bamfile] = {}

            seqlen = int(sections[1])
            if seqlen in bam_files[bamfile]:
                raise Exception(f"Duplicate seq length {seqlen} for {bamfile}")
            
            else:
                bam_files[bamfile][seqlen] = int(sections[2])

    if not bam_files:
        raise Exception("Found no valid BAM files in config file")
    
    return bam_files


def read_transcripts(file, keepers):

    transcript_seqs = {}

    id = None
    sequence = ""

    if file.lower().endswith(".gz"):
        infh = gzip.open(file,"rt",encoding="utf8")
    else:
        infh = open(file,"rt",encoding="utf8")

    for line in infh:
        line = line.strip()

        if line.startswith(">"):
            if id is not None and id in keepers:
                transcript_seqs[id] = sequence
            sequence = ""
            id = line.split()[0][1:]
            if "." in id:
                id = id[:id.index(".")]
        else:
            sequence += line


    if id is not None and id in keepers:
        transcript_seqs[id] = sequence

    return transcript_seqs

def get_options():
    parser = argparse.ArgumentParser("Analyse codon usage in RiboSeq BAM files")

    parser.add_argument("codons", type=str, help="File of codon positions prepared by prepare_genome.py script")
    parser.add_argument("fasta", type=str, help="Multi-Fasta File of transcript sequences")
    parser.add_argument("config", type=str,help="Config file for BAM files to analyse")
    parser.add_argument("--outfile",type=str, help="Output file to save to", default="quantitated_codons.txt")

    options = parser.parse_args()

    if not Path(options.codons).exists():
        raise Exception("Codon file path"+options.codons+"doesn't exist")

    if not Path(options.fasta).exists():
        raise Exception("Fasta path"+options.fasta+"doesn't exist")

    if not Path(options.config).exists():
        raise Exception("Config path"+options.config+"doesn't exist")


    return options


if __name__ == "__main__":
    main()