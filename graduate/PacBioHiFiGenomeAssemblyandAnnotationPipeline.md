---
layout: single
title: "PacBio HiFi Genome Assembly and Annotation Pipeline"
permalink: /code/PacBioHiFiGenomeAssemblyandAnnotationPipeline/
author_profile: true
---

# PacBio HiFi Genome Assembly and Annotation Pipeline

This pipeline describes a complete workflow from raw PacBio HiFi reads through contaminant removal, genome assembly, repeat masking, transcript alignment, and orthology-based annotation. It is designed to be modular, reproducible, and independent of machine-specific paths.

The workflow assumes HiFi reads as input and produces a polished genome assembly with annotated gene models.

A recommended project structure:

```
project/
├── data/
│   ├── raw/
│   ├── reference/
│   └── rna/
├── results/
│   ├── cleaned_reads/
│   ├── alignments/
│   ├── assembly/
│   ├── repeats/
│   ├── annotation/
│   └── expression/
├── scripts/
└── logs/
```

HiFi reads are generated from subreads (if starting from BAM):

```
pbccs data/raw/input.subreads.bam results/cleaned_reads/hifi.bam \
  --report-file results/logs/ccs_report.txt
```

Convert to FASTQ:

```
samtools fastq results/cleaned_reads/hifi.bam > results/cleaned_reads/hifi.fastq
```

Remove contaminant reads by aligning to a contaminant reference:

```
minimap2 -d data/reference/contaminant.mmi data/reference/contaminant.fa

minimap2 -ax map-hifi data/reference/contaminant.mmi results/cleaned_reads/hifi.fastq | \
  samtools sort -o results/alignments/contaminant.bam

samtools index results/alignments/contaminant.bam
```

Extract reads that do not map to contaminants:

```
samtools view -b -f 4 results/alignments/contaminant.bam > results/cleaned_reads/filtered.bam
samtools fastq results/cleaned_reads/filtered.bam > results/cleaned_reads/filtered.fastq
```

Optional length filtering:

```
seqkit seq -m 1000 results/cleaned_reads/filtered.fastq > results/cleaned_reads/filtered_long.fastq
```

Genome assembly using HiFi reads:

```
hifiasm -o results/assembly/genome -t 32 results/cleaned_reads/filtered_long.fastq
```

Optional Hi-C scaffolding:

```
hifiasm -o results/assembly/genome_hic -t 32 \
  --h1 data/reference/hic_R1.fastq.gz \
  --h2 data/reference/hic_R2.fastq.gz \
  results/cleaned_reads/filtered_long.fastq
```

Repeat annotation and masking:

```
BuildDatabase -name genome_db -engine ncbi results/assembly/genome.fa

RepeatModeler -engine ncbi -pa 8 -database genome_db

RepeatMasker -pa 8 -gff -lib genome_db-families.fa \
  -dir results/repeats/ results/assembly/genome.fa
```

Align long-read RNA-seq data:

```
minimap2 -ax splice -uf --secondary=no \
  results/assembly/genome.fa data/rna/long_reads.fastq.gz > results/alignments/rna.sam

samtools view -b results/alignments/rna.sam > results/alignments/rna.bam
samtools sort results/alignments/rna.bam -o results/alignments/rna.sorted.bam
samtools index results/alignments/rna.sorted.bam
```

Transcript assembly and quantification:

```
stringtie results/alignments/rna.sorted.bam \
  -G results/annotation/reference.gtf \
  -o results/expression/transcripts.gtf \
  -A results/expression/gene_abundance.tsv
```

Orthology-based gene annotation via reciprocal BLAST:

```
blastp -query results/annotation/proteins.fa \
  -db data/reference/reference_proteins \
  -outfmt "6 std slen qlen" \
  -max_target_seqs 1 \
  -evalue 1e-5 \
  -num_threads 8 \
  > results/annotation/forward_blast.tsv
```

Extract best hits:

```
cut -f2 results/annotation/forward_blast.tsv | sort | uniq > results/annotation/best_hits.list
```

Retrieve sequences:

```
bioawk -c fastx '
BEGIN { while ((getline k < "results/annotation/best_hits.list") > 0) i[k] = 1 }
{ if (i[$name]) print ">" $name "\n" $seq }
' data/reference/reference_proteins.fa > results/annotation/best_hits.fa
```

Create local BLAST database:

```
makeblastdb -in results/annotation/proteins.fa \
  -dbtype prot \
  -out results/annotation/local_db
```

Reciprocal BLAST:

```
blastp -query results/annotation/best_hits.fa \
  -db results/annotation/local_db \
  -outfmt 6 \
  -max_target_seqs 1 \
  -evalue 1e-5 \
  -num_threads 8 \
  > results/annotation/reciprocal.tsv
```

Identify reciprocal best hits:

```
awk '{print $1 "\t" $2}' results/annotation/forward_blast.tsv > results/annotation/forward_pairs.txt
awk '{print $2 "\t" $1}' results/annotation/reciprocal.tsv > results/annotation/reciprocal_pairs.txt

sort results/annotation/forward_pairs.txt results/annotation/reciprocal_pairs.txt | uniq > results/annotation/rbh.txt
```

Apply gene renaming:

```
python scripts/rename_genes.py
```

This pipeline produces a repeat-masked genome assembly, transcript-supported gene models, ortholog-informed annotations, and expression estimates. It can be extended with additional polishing, scaffolding, or comparative genomics analyses depending on project goals.
