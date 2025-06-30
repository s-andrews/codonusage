#!/usr/bin/env python3

import argparse
from pathlib import Path
import gzip

def main():
    options = get_options()

    cds_regions = parse_gtf(options.gtf)

    codon_positions = extract_codon_positions(cds_regions)

    save_codons(codon_positions, options.outfile)


def save_codons(codon_positions, outfile):

    with open(outfile,"wt",encoding="utf8") as out:
        for chromosome in codon_positions:
            for position in codon_positions[chromosome]:
                transcript,transcript_position = codon_positions[chromosome][position]

                print("\t".join([str(x) for x in [chromosome,position,transcript,transcript_position]]), file=out)
    



def extract_codon_positions(cds_regions):
    "Extract all of the codon positions from parsed CDS sequences"

    genome_positions = {}

    for transcript in cds_regions:
        forward = cds_regions[transcript][1]
        if not transcript[0] in genome_positions:
            genome_positions[cds_regions[transcript][0]] = {}

        transcript_position = 1
        codon_offset = 0

        for codon in cds_regions[transcript][2]:
            if forward:
                codon_position = codon_offset+codon[0]
                while codon_position <= codon[1]:
                    genome_positions[cds_regions[transcript][0]][codon_position] = (transcript,transcript_position)
                    codon_position += 3
                    transcript_position += 3

                codon_offset = codon_position - (codon[1]+1)


            else:
                codon_position = codon[1]-codon_offset
                while codon_position >= codon[0]:
                    genome_positions[cds_regions[transcript][0]][codon_position] = (transcript,transcript_position)
                    codon_position -= 3
                    transcript_position += 3

                codon_offset = 0 - (codon_position - (codon[0]-1))

    return genome_positions
    



def parse_gtf(file):
    "Parse a GTF file and extract all of the codon positions across the genome"

    cds_regions = {}

    if file.lower().endswith(".gz"):
        infh = gzip.open(file,"rt",encoding="utf8")
    else:
        infh = open(file,"rt",encoding="utf8")

    for line in infh:
        if line.startswith("#"):
            continue

        line = line.strip()
        sections = line.split("\t")

        if not sections[2] == "CDS":
            continue

        chrom = sections[0]
        start = int(sections[3])
        end = int(sections[4])
        forward = sections[6] == "+"

        # We need to find the transcript id
        for tag in sections[8].split(";"):
            if tag.strip().startswith("transcript_id"):
                transcript = tag.strip().replace("transcript_id ","").replace('"',"").strip()
                if not transcript in cds_regions:
                    cds_regions[transcript] = [chrom,forward,[]]

                if not cds_regions[transcript][1] == forward:
                    raise Exception(f"Transcript {transcript} with exons {cds_regions[transcript][2]} had a mix of forward and reverse exons")

                cds_regions[transcript][2].append((start,end))

    infh.close()

    for transcript in cds_regions:
        if cds_regions[transcript][1]:
            # It's forward
            cds_regions[transcript][2].sort(key=lambda x: x[0])
        else:
            # It's revese so we sort backwards
            cds_regions[transcript][2].sort(key=lambda x: x[0], reverse=True)

    return cds_regions


def get_options():
    parser = argparse.ArgumentParser("Prepare a genome for codon analysis")

    parser.add_argument("gtf", type=str, help="GTF file to parse")
    parser.add_argument("fasta", type=str, help="Multi-Fasta File of transcript sequences")
    parser.add_argument("--outfile",type=str, help="Output file to save to", default="codon_data.dat")

    options = parser.parse_args()

    if not Path(options.gtf).exists():
        raise Exception("GTF path"+options.gtf+"doesn't exist")

    if not Path(options.fasta).exists():
        raise Exception("Fasta path"+options.fasta+"doesn't exist")


    return options


if __name__ == "__main__":
    main()