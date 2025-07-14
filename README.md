## Conda Environment Setup
```
# Clone the git repository
git clone https://github.com/Franck-Dumetz/ORFeome.git
# Enter into the ORFeome folder
cd ORFeome
# Set up the conda environment
conda env create -f environment.yml
# Activate the conda environment
conda activate ORFeome
```

## Running the Pipeline
```
./ORF.sh -A <your_annotation_file.gff> -G <your_genome_file.fasta> -T <your_treatments_file.csv> [-R <your_sras_list.txt> | -F <your_fastq_directory_path>] [-u] [-m]
```
Explanation of arguments:
- The genome annotation file. This file must be in GFF format.
- The genome file. This file must be in FASTA format.
- The treatments file. This file is used to help the Deseq script know which groups each sample belongs so that it can do the appropriate comparison. The CSV file must be a two-column file with sample names on the left and condition or library names on the right. The sample names must match either the SRA accession numbers, if you provide an SRA ID list, or the names of the fastq files, if you provide a fastq directory. The names in the right column must start with "treated" or "untreated". So if you have an untreated sample that is from a library names NewMP, you should enter untreated_NewMP in the CSV. An example of a treatments file is in the test-data folder. The file must not contain a header, as the code assumes the first line corresponds to the first sample.
- The SRA ID list file or the FASTQ directory path. This pipeline will either:
    - Take in a txt file containing a list of SRA accession IDs and fetch them from NCBI:
      ```
      SRR10846669
      SRR10846670
      SRR10846671
      ...
      ```
    - Or a FASTQ directory path containing all of the fastq samples.
- The flags (optional). You can choose to only look at unique alignments (-u), only look at multiple alignments (-m), or do both by including both flags. If no flags are included, both types of alignments will be done. 
  
The order of arguments does not matter, as long as the appropriate flag is placed before it. 

### Usage

```
./ORF.sh -A <annotation.gff> -G <genome.fasta> -T <treatments.csv> [-R <sra_list.txt> | -F <fastq_directory>] [-u] [-m]
```

### Argument Descriptions

| Flag | Required? | Description |
|------|-----------|-------------|
| `-A` | ✅ | **Genome annotation file** in GFF format. |
| `-G` | ✅ | **Genome file** in FASTA format. |
| `-T` | ✅ | **Treatments CSV file** for grouping samples. <br><br>Must be a **two-column, headerless CSV** with:<br>• Sample names in column 1 (matching SRA IDs or FASTQ names)<br>• Condition/library names in column 2, starting with `treated_` or `untreated_` (e.g., `untreated_NewMP`)<br><br>An example is provided in the `test-data/` folder. |
| `-R` | One of `-R` or `-F` is required | Path to a text file containing SRA accession numbers, one per line:<br>```\nSRR10846669\nSRR10846670\n...``` |
| `-F` | One of `-R` or `-F` is required | Path to a directory containing `.fastq` files. |
| `-u` | Optional | Process **uniquely aligned** reads only. |
| `-m` | Optional | Process **multi-mapped** reads only.<br>If neither `-u` nor `-m` is used, both types will be processed. |

> **Note:** The order of arguments does not matter, but each flag must come before its corresponding input.

