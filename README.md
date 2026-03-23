# 🧬 Sequence Analysis & Phylogeny Pipeline

Bioinformatics pipeline for sequence alignment, consensus generation, and phylogenetic analysis.

## 🎯 Project Goal
Developed an end-to-end bioinformatics pipeline to extract target genomic regions, generate consensus sequences from sequencing reads, and construct phylogenetic relationships between samples.


## Overview

This project implements a complete bioinformatics workflow integrating:

* in silico PCR (iPCR)
* read mapping & consensus generation
* sequence alignment
* phylogenetic tree construction

The pipeline simulates real-world genomic analysis used in microbial and clinical sequencing studies.


## Key Features

* Extract target regions using **in silico PCR**
* Map sequencing reads to reference genomes
* Generate **consensus sequences** from aligned reads
* Perform **global sequence alignment (Needleman–Wunsch)**
* Compute pairwise sequence distances
* Construct phylogenetic trees using **Neighbor-Joining algorithm**


## Project Structure

```
scripts/
│── magop.py            # Main pipeline (entry point)
│── ispcr.py            # In silico PCR (BLAST + seqtk)
│── mapping.py          # Read mapping using minimap2
│── sam.py              # SAM parsing + consensus generation
│── nw.py               # Needleman–Wunsch alignment
│── run_external.py     # Wrapper for external tools
│── map_consensus.py    # CLI for consensus extraction
│── __init__.py         # Package initializer
```


## Workflow

### 1️⃣ In Silico PCR

* Identifies primer binding sites using BLAST
* Pairs forward and reverse primers
* Extracts amplicons from sequences

### 2️⃣ Read Mapping & Consensus

* Maps reads using **minimap2 + samtools**
* Builds consensus sequences from mapped reads
* Performs pileup-based base calling

### 3️⃣ Sequence Alignment

* Implements **Needleman–Wunsch global alignment algorithm**
* Calculates mismatches and alignment scores

### 4️⃣ Phylogenetic Analysis

* Computes pairwise sequence distances
* Builds phylogenetic tree using **Neighbor-Joining algorithm**


## Example Usage

### Run full pipeline

```bash
python magop.py \
  -a assemblies.fasta \
  -p primers.fasta \
  -r sample_R1.fastq sample_R2.fastq \
  -s reference.fasta
```

### Generate consensus from reads

```bash
python map_consensus.py \
  -1 reads_1.fastq \
  -2 reads_2.fastq \
  -r reference.fasta
```



## 📈 Example Output

* Consensus sequences generated from raw sequencing reads
* Extracted amplicon regions from genomes
* Pairwise sequence distance matrix
* Phylogenetic tree in Newick format

Example tree:

```
((Sample1,Sample2),Sample3);
```


##  Requirements

* Python 3.x
* minimap2
* samtools
* BLAST (blastn)
* seqtk
* NumPy
* Biotite


##  Skills Demonstrated

* Bioinformatics pipeline development
* Sequence alignment algorithms
* Genomic data processing
* Command-line tool integration
* Phylogenetic analysis


##  Impact

This pipeline mimics real-world workflows used in:

* microbial genomics
* 16S rRNA analysis
* evolutionary biology
* clinical sequencing pipelines


## 👩‍💻 Author

Divya Reddy
MS Bioinformatics | Georgia Tech
