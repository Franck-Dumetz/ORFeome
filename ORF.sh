#!/bin/bash

gff=""
sras=""
fasta=""
unique=0
multiple=0
treatments=""

usage="" #change later

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -A)
      gff="$2"
      shift 2
      ;;
    -R)
      sras="$2"
      shift 2
      ;;
    -G)
      fasta="$2"
      shift 2
      ;;
    -u)
      unique=1
      shift 1
      ;;
    -m)
      multiple=1
      shift 1
      ;;
    -T)
      treatments="$2"
      shift 2
      ;;
    *)
      echo "<<Invalid flag: $1>>"
      echo "$usage"
      exit 1
      ;;
  esac
done

#Argument checking
if [[ -z "$sras" ]]; then
  if [[ ! -d fastqs ]]; then
    echo "<<No SRA ID file was provided and no fastqs directory exists>>"
    echo "$usage"
    exit 1
  fi
fi

if [[ -z "$gff" ]]; then
  echo "<<No gff file provided>>"
  echo "$usage"
  exit 1
fi

if [[ -z "$fasta" ]]; then
  echo "<<No fasta file provided>>"
  echo "$usage"
  exit 1
fi

if [[ "$unique" -eq 0 && "$multiple" -eq 0 ]]; then
  echo "<<No alignment flags provided. Defaulting to unique and multiple alignments>>"
  unique=1
  multiple=1
fi

if [[ -z "$treatments" ]]; then
  echo "<<No treatment file provided>>"
  echo "$usage"
  exit 1
fi

#SRA -> FASTQ (populates & creates the fastq directory)
if [[ ! -d fastqs ]]; then
  ./fastq.sh $sras > output.log
  echo "SRAs converted to FASTQs"
fi

#FASTQ -> Trimmed FASTQ (populates & creates the trimmed directory)
mkdir trimmed
trim_galore --output_dir trimmed fastqs/* >> output.log
echo "FASTQs trimmed"

#Trimmed FASTQ -> SAM (populates & creates the sam directory)
mkdir sam

bowtie-build $fasta index

for file in trimmed/*.fq; do
  basename="${file##*/}"
  basename="${basename%_*}"
  if [ "$unique" -eq 1 ]; then
    bowtie --best --strata -t -v 2 -a -m 1 -S index $file > sam/"${basename}_m1.sam"
  fi
  if [ "$multiple" -eq 1 ]; then
    bowtie --best --strata -t -v 2 -a -m 10 -S index $file > sam/"${basename}_m10.sam"
  fi
done
bowtie --best --strata -t -v 2 -a -m 10 -S index trimmed/
echo "Bowtie complete"

# SAM -> BAM
mkdir bam

for file in sam/*; do
  basename="${file##*/}"
  basename=${basename%.*}
  echo "$file"
  samtools view -S -b $file > bam/"${basename}_unsorted.bam"
  samtools sort -o bam/"${basename}.bam" bam/"${basename}_unsorted.bam"
  samtools index bam/"${basename}.bam"
done

rm bam/*_unsorted.bam
echo "BAM files made"

#Trim the gff
python trim-gff.py $gff
echo "GFF trimmed"

#Run seq-monk (generates counts.csv)
python seq-monk.py
echo "SeqMonk complete"

#Create the treatments.csv
if [[ "$unique" -eq 1 && "$multiple" -eq 1 ]]; then
  python make_csv.py -um $treatments
elif [ "$unique" -eq 1 ]; then
  python make_csv.py -u $treatments
else
  python make_csv.py -m $treatments
fi

#echo "Treatments file created"

#Run Deseq2
Rscript Deseq2.R
echo "Deseq analysis complete. Results in overrep_genes.csv"

#Cleaning all created directories
if [[ ! -z "$sras" ]]; then
  rm -r fastqs
fi

rm -r trimmed
rm -r sam
rm -r bam

