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

### Running with Test Data
Use the following command to run with the provided test data set:
```
./ORF.sh -A test-data/test-annot.gff -R test-data/test-sra.txt -G test-data/test-genome.fasta -u -T test-data/test-treatments.csv
```

