# Conda Environment Setup
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

# Running the Pipeline
```
./ORF.sh -A <your_annotation_file.gff> -G <your_genome_file.fasta> -T <your_treatments_file.csv> [-R <your_sras_list.txt> | -F <your_fastq_directory_path>] [-u] [-m]
```
Explanation of arguments:
- The genome annotation file. This file must be in GFF format.
- The genome file. This file must be in FASTA format.
- The treatments file. This file is used to help the Deseq script know which groups each sample belongs so that it can do the appropriate comparison. The CSV file must be a two-column file with sample names on the left and condition or library names on the right. The sample names must match either the SRA accession numbers, if you provide an SRA ID list, or the names of the fastq files, if you provide a fastq directory. The names in the right column must start with "treated" or "untreated". So if you have an untreated sample that is from a library names NewMP, you should enter untreated_NewMP in the CSV. An example of a treatments file is in the test-data folder. The file must not contain a header, as the code assumes the first line corresponds to the first sample.
- The SRA ID list file or the FASTQ directory path. This pipeline will either:
-   Take in a txt file containing a list of SRA accession IDs and fetch them from NCBI:
    ```
    SRR10846669
    SRR10846670
    SRR10846671
    ...
    ```
    Or, if no list is given,

The order of arguments does not matter, as long as the appropriate flag is placed before it. 
