#!/usr/bin/env python3

import argparse
from pathlib import Path
import gzip
import pysam
import sys

def main():
    options = get_options()

    print("Parsing config file",options.config, file=sys.stderr, flush=True)
    bam_files = parse_bam_config(options.config)

    print("Parsing transcripts from",options.fasta, file=sys.stderr, flush=True)
    transcript_sequences = read_transcripts(options.fasta)

    print("Parsing codon information from",options.codons, file=sys.stderr, flush=True)
    codon_data = read_codons(options.codons)

    all_results = {}

    for bam_file in bam_files:
        print("Parsing BAM file",bam_file["filename"], file=sys.stderr, flush=True)

        stats,codon_counts = parse_bam_file(bam_file, codon_data, transcript_sequences)

        write_results(stats,codon_counts,bam_file)

def write_results(stats,codon_counts,bam_file):

    outfile = bam_file["filename"]+"_codonusage.txt"

    with open(outfile,"wt",encoding="utf8") as out:

        print("# Codonusage results", file=out)
        print("# Filename:",bam_file["filename"], file=out)

        for metric in stats:
            print(f"# {metric}: {stats[metric]}", file=out)

        print(f"Codon\tCount_A\tCount_P\tCount_E",file=out)

        for codon in codon_counts["A"]:
            print(codon,codon_counts["A"][codon],codon_counts["P"][codon],codon_counts["E"][codon], sep="\t", file=out)


def parse_bam_file(fileinfo,codons, transcripts):

    file = fileinfo["filename"]
    offsets = fileinfo["offsets"]

    bamfile = pysam.AlignmentFile(file, "rb")

    stats = {
        "total_reads": 0,
        "aligned_reads": 0,
        "good_enough": 0,
        "recognised_sequence": 0,
        "valid_length" : 0,
        "mapped_to_codon": 0,
        "incomplete_codon" : 0,
        "not_at_codon": 0,
        "base_not_aligned": 0
    }

    codon_counts = {"P":{},"A":{},"E":{}}

    for p1 in ["G","A","T","C"]:
        for p2 in ["G","A","T","C"]:
                for p3 in ["G","A","T","C"]:
                        codon_counts["P"][p1+p2+p3] = 0
                        codon_counts["A"][p1+p2+p3] = 0
                        codon_counts["E"][p1+p2+p3] = 0

    # Because the CDS sequences don't have the stop codon
    # included we can't separate those so we'll do a generic
    # STP feature to measure those
    codon_counts["P"]["STP"] = 0
    codon_counts["A"]["STP"] = 0
    codon_counts["E"]["STP"] = 0

    # For the P site we may have the ribosome at the start codon
    # so we need a way to indicate that this was at the promoter
    # we'll only get this for the A site but we'll add it to all
    # dictionaries to make writing results easier.
    codon_counts["P"]["TSS"] = 0
    codon_counts["A"]["TSS"] = 0
    codon_counts["E"]["TSS"] = 0

    for read in bamfile.fetch(until_eof=True):
        
        stats["total_reads"] += 1
        if stats["total_reads"] % 1000000 == 0:
            print(f"Parsed {stats['total_reads']/1000000}M reads from {file}", file=sys.stderr, flush=True)

        # # Just for testing...
        # if stats["total_reads"] == 1000000:
        #     break


        # We don't want reads which aren't aligned at all
        if not read.is_mapped:
            continue
        stats["aligned_reads"] += 1

        # We don't want reads which have MAPQ < 20 as they're not uniquely mapped
        if not read.mapq >= 20:
            continue

        stats["good_enough"] += 1

        # We might be aligning to a reference sequence for which there are no transcripts
        # so we can give up if that's the case - this won't be main chromosomes of course
        if not read.reference_name in codons:
            continue

        stats["recognised_sequence"] += 1


        # Only some read lengths are being used.  If we don't have a recorded offset for
        # this read length then we can move on
        if not read.query_length in offsets:
            continue

        stats["valid_length"] += 1

        # We need to work out the position within the read which we want. This will be the
        # 5' end of the read extended by the offset for this read length.  For forward mapping
        # reads this will just be the offset 
        wanted_read_position = offsets[read.query_length]


        # For reverse reads we need the position subtracted from the end.
        #
        # The positions in aligne
        if read.is_reverse:
            wanted_read_position = (read.query_length-1) - wanted_read_position

        for offset,refpos in read.get_aligned_pairs(matches_only=True):
            if offset == wanted_read_position:
                central_pos = refpos+1

        # If we don't find a position here it's because that base in the read
        # wasn't aligned to a reference position - so would be an insertion
        if not central_pos:
            stats["base_not_aligned"] += 1
            continue

        # Check to see if this is the first base of a codon triplet
        if not central_pos in codons[read.reference_name]:
            # if read.reference_name == "1" and central_pos < 4854381 and central_pos > 4854169:
            #     print(f"No hit at {refpos}")
            stats["not_at_codon"] += 1
            continue

        transcript,position = codons[read.reference_name][central_pos]
        # if read.reference_name == "1" and central_pos < 4854381 and central_pos > 4854169:
        #     print(f"Found hit at {refpos} to {transcript}")


        # We need to check for a stop but we need to keep going so we can still profile the
        # A site.  We shouldn't have stop codons in the CDS, but there are cases where it
        # happens, specifically some of the V-Genes have encoded stops at the end but the CDS
        # includes these for some reason.
        if position == len(transcripts[transcript]):
            # It's the stop codon
            stats["mapped_to_codon"] += 1
            codon_counts["P"]["STP"] += 1

        elif position+2 > len(transcripts[transcript])-1:
            # In some cases we can get incomplete transcripts
            # the sequence just randomly stops in the middle
            # of a codon.
            stats["incomplete_codon"] += 1

        else:
            # The offset is the offset to the P site of the ribosome
            # so we'll pull that out first. The positions start at 1
            # so we can go back 1 for the zero indexes of the string            
            p_codon = transcripts[transcript][position-1:(position+2)]

            if not p_codon in codon_counts["P"]:
                stats["incomplete_codon"] += 1
                continue

            stats["mapped_to_codon"] += 1
            codon_counts["P"][p_codon] += 1

        # Now that we've got that we can go after the A and E sites
        # We may not have them if this is the start or end of the 
        # transcript

        a_position = position-3
        e_position = position+3

        if a_position < 0:
            codon_counts["A"]["TSS"] += 1
        else:
            a_codon = transcripts[transcript][a_position-1:(a_position+2)]
            if a_codon in codon_counts["A"]:
                codon_counts["A"][a_codon] += 1

        if e_position == len(transcripts[transcript]):
            # It's the stop codon
            codon_counts["E"]["STP"] += 1
        elif e_position+2 > len(transcripts[transcript])-1:
            continue
        else:
            e_codon = transcripts[transcript][e_position-1:(e_position+2)]
            if e_codon in codon_counts["E"]:
                codon_counts["E"][e_codon] += 1
            


    return stats,codon_counts


def read_codons(file):
    codon_data = {}

    with open(file,"rt",encoding="utf8") as infh:
        for line in infh:
            chrom,pos,transcript,transcript_pos = line.strip().split("\t")

            if not chrom in codon_data:
                codon_data[chrom] = {}

            codon_data[chrom][int(pos)] = (transcript,int(transcript_pos))


    return codon_data


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
                bam_files[bamfile] = {"filename":bamfile, "offsets": {}}

            for seqlen in range(int(sections[1]),int(sections[2])+1):
                if seqlen in bam_files[bamfile]:
                    raise Exception(f"Duplicate seq length {seqlen} for {bamfile}")
                
                else:
                    bam_files[bamfile]["offsets"][seqlen] = int(sections[3])

    if not bam_files:
        raise Exception("Found no valid BAM files in config file")
    
    return list(bam_files.values())


def read_transcripts(file):

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
            if id is not None:
                transcript_seqs[id] = sequence
            sequence = ""
            id = line.split()[0][1:]
            if "." in id:
                id = id[:id.index(".")]
        else:
            sequence += line


    if id is not None:
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