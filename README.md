# Analysis of Genomic Variation in High-Throughput Sequencing Data

This repository contains an integrated analytical framework designed to maximize the utility of affordable whole-genome sequencing (WGS) data. By leveraging replicate sequencing, Mendelian inheritance patterns, and pangenome references, this pipeline resolves low-coverage or repetitive regions that are often problematic in standard 30X WGS data.

## Overview
Standard 30X sequencing often suffers from gaps in repetitive regions or ambiguous calls at low-coverage loci. This project provides a robust methodology to bridge these gaps using a three-generation family-based approach.

Key Features:

- Sample Quality Assessment: Compares blood-derived vs. saliva-derived WGS quality.

- Ambiguous Site Resolution: A custom imputation strategy utilizing replicate support and Mendelian consistency.

- Pangenome Integration: Uses pangenome panels to enhance mapping and variant calling completeness.

- Functional Annotation & Phasing: Clinical variant prioritization (ClinVar) and pedigree-based haplotype resolution.

## Pipeline Architecture
The framework is divided into three primary modules:

1. Quality Control & Alignment
Assesses sequencing metrics, including depth of coverage and mapping quality across different sample sources (blood vs. saliva).

2. Variant Refinement (The Core Engine)
Resolves "ambiguous" genotypes by integrating three layers of evidence:

- Replicate Support: Cross-referencing calls from multiple runs of the same individual.

- Mendelian Consistency: Using family pedigree data to validate or correct inheritance patterns.

- Pangenome Imputation: Utilizing a pangenome reference to resolve complex or repetitive loci.

3. Functional Interpretation
- Clinical Annotation: Filters and identifies "Pathogenic" and "Likely Pathogenic" variants via ClinVar.

- Haplotype Phasing: Gene-wise phasing within coding regions to resolve local haplotypes.

## Data Privacy Notice
Please Note: The genomic data utilized in this thesis is strictly private and is not included in this repository.

This repository provides the source code, scripts, and analytical framework only. 

