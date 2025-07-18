# ORFeome

A pipeline for analyzing ORFeome screening data under different stress like drug pressure.

---

Publication: 

---

## Overview

**ORFeome** is a set of scripts and tools designed to analyze gain-of-function (GoF) screens in _T. brucei_. It identifies open reading frames (ORFs) that are significantly overrepresented in treated samples compared to untreated controls, helping to discover genes associated different stresses like exposure to drug.

This pipeline was initially developed for analyzing screens in kinetoplastid parasites but can be adapted to other organisms. This pipeline revamps, add options, and automatises the data analysis from Carter M, et al. [A Trypanosoma brucei ORFeome-Based Gain-of-Function Library Identifies Genes That Promote Survival during Melarsoprol Treatment](https://journals.asm.org/doi/full/10.1128/msphere.00769-20?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org). mSphere. 2020 Oct 7;5(5):e00769-20. doi: 10.1128/mSphere.00769-20. PMID: 33028684; PMCID: PMC7568655.

---

## Key Features

- Handles multiple replicates per condition  
- Uses DESeq2 for statistical testing
- Computes differential representation of ORFs, set your own Fold Change threshold
- Test the diversity of the control culture
- Highly customizable for different datasets and experimental designs

---


## Conda Environment Setup
```
# Clone the git repository
git clone https://github.com/Franck-Dumetz/ORFeome.git
```
```
# Enter into the ORFeome folder
cd ORFeome
```
# Set up the conda environment
```
# Linux 
conda create -n ORFeome --file conda-linux-64.lock
```
```
# MacOS
conda create -n ORFeome --file conda-osx-64.lock
```
```
# Activate the conda environment
conda activate ORFeome
```

## Running the Pipeline

### Usage

```
./ORF-enrich.sh -A <annotation.gff> -G <genome.fasta> -T <treatments.csv> -C <fold_change> [-R <sra_list.txt> | -F <fastq_directory>] [-u] [-m]
```

### Argument Descriptions

| Flag | Required? | Description |
|------|-----------|-------------|
| `-A` | ✅ | **Genome annotation file** in GFF format. |
| `-G` | ✅ | **Genome file** in FASTA format. |
| `-T` | ✅ | **Treatments CSV file** for grouping samples. <br><br>Must be a **two-column, headerless CSV** with:<br>• Sample names in column 1 (matching SRA IDs or FASTQ names)<br>• Condition/library names in column 2, starting with `treated_` or `untreated_` (e.g., `untreated_NewMP`)<br><br>An example is provided in the `test-data/` folder. |
| `-C` | ✅ | **Fold change** as an integer. |
| `-R` | One of `-R` or `-F` is required | Path to a text file containing SRA accession numbers, one per line:<br>```SRR10846669\nSRR10846670\n...``` |
| `-F` | One of `-R` or `-F` is required | Path to a directory containing `.fastq` files. |
| `-u` | Optional | Process **uniquely aligned** reads only. |
| `-m` | Optional | Process **multi-mapped** reads only.<br>If neither `-u` nor `-m` is used, both types will be processed. |

> **Note:** The order of arguments does not matter, but each flag must come before its corresponding input.

### Running with Test Data
Use the following command to run with the provided test data set:
```
./ORF-enrich.sh -A test-data/test-annot.gff -R test-data/test-sra.txt -G test-data/test-genome.fasta -u -T test-data/test-treatments.csv -C 4
```

