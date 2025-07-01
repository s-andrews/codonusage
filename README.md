# Codonusage
A program to analyse RiboSeq data for the usage of alternate codons.

# Introduction
This program is designed to extract codon usage data from RiboSeq experiments which may contain complex relationships between the length of reads and the position of the P site in the ribosome.

It produces a count table of the number of observations for each codon in the A, P and E sites of the ribosome.

# Installation
To install the program you can either clone the current version of the git repository, or download the latest release.

```
git clone https://github.com/s-andrews/codonusage.git
cd codonusage
```

We reccommend using a virtual python environment to install the program's dependencies.

```
python3 -m venv venv
source venv\bin\activate
pip install -r requirement.txt
```

That should then allow you to run the programs.

# Usage

## Data required 
To run the program you will need four different inputs

1. A GTF file of annotations for the genome you're using.  The table at https://www.ensembl.org/info/data/ftp/index.html has links to GTF files for all of the common model organisms.  It's important that the genome assembly used for the GTF file is the same as the assembly used to align your data

2. A FastA file of CDS sequences for all of the protein coding transcripts described in the GTF file above.  Again, Ensembl provides a CDS file for all of its genomes

3. One or more BAM files aligned to the same reference genome as the GTF and CDS files. The BAM file does not need to be sorted, but the reads should be trimmed (not soft clipped) and any UMI deduplication should have been performed before starting the quantitation.

4. Information on the P site offsets for each file. This should come in the form of the offset from the 5' end of the read to the first base of the P site codon.  Typically this value of 12, but in some library types you might wish to use a program such as [RiboWaltz](https://github.com/LabTranslationalArchitectomics/riboWaltz) to calculate this value.  You can have different offsets for different lengths of read.

## Creating a codon file
The first part of the analysis involves parsing the GTF file to pull out all of the codon positions in the genome.  This operation only requires the GTF file and can be performed once for each genome to produce a file you can re-use.

Running this step is as simple as:

```prepare_genome.py human.gtf.gz```

The parsed genome will be contained in a file called ```codon_data.dat``` by default although you can change this with the ```--outfile``` option.  Your GTF can either be plain text or can be gzip compressed.

## Preparing a config file
For the main analysis you will need to prepare a tab delimited config file.  This lists the BAM files you want to quantitate along with the range of lengths of reads you want to use and the P site offset values you want to use for them.  An example is shown below.

```
#Filename    Min_Len   Max_Len   Offset
file1.bam    28        30        12
file1.bam    31        35        13
file2.bam    28        35        13
```

Only reads whose lengths have a specified offset in the config file will be used in the analysis.

## Running the count
To run the count you use the ```analyse_bam_files.py``` program as follows:

```analyse_bam_files.py codon_data.dat cds_sequences.fa my_files_config.txt```

The codon information is the output of the first step.  Your CDS files can either be plain test fasta, or can be gzip compressed.

# Output
For each BAM file an output file will be created with the same name and location, but with ```_codonusage.txt``` appended to the end.

## Header
At the top of the file will be some statistics about the processing

```
# Codonusage results
# Filename: file1.bam
# total_reads: 35251488
# aligned_reads: 35251488
# good_enough: 20559568
# recognised_sequence: 20557458
# valid_length: 4816676
# mapped_to_codon: 1646830
# incomplete_codon: 1387
# not_at_codon: 3168459
# base_not_aligned: 0
```

The fields here are as follows:

* **Total_Reads**: The total number of reads in the BAM file
* **Aligned_Reads**: The number of reads with an alignment to the reference
* **Good_Enough**: Alignments with a MAPQ value of >=20
* **Recognised_Sequence**: The alignment was to a reference sequence which had at least one coding transcript on it
* **Valid_Length**: Reads whose length matched one of the read lengths specified in the config file
* **Mapped_to_Codon**: Reads whose P site matched the start of a valid codon
* **Incomplete_Codon**: Reads which mapped to an annotated codon but where the transcript truncated before the end of the codon
* **Not_at_a_Codon**: Reads which didn't match an annotated codon
* **Base_not_Aligned**: Reads where there was an insertion such that no direct reference position was available

## Body
The main body of each output file is a list of all codons and their summed counts at the A, P and E sites of the ribosome.

```
Codon   Count_A Count_P Count_E
GGG     31580   34614   32219
GGA     46613   47470   58213
GGT     10774   12205   14845
GGC     14843   16439   18218
```

In addition to standard codons there are two artificial codons:

* ```STP``` Is a generic stop codon, since the stop codon is not generally listed in the CDS files
* ```TSS``` Specifically for the A site this indicates the translation start site - the codon before the initiating methionine