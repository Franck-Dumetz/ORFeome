mkdir -p fastqs

while read srr; do
    echo "Processing $srr..."

    prefetch "$srr"

    # Convert to FASTQ format
    fastq-dump --split-3 --outdir fastqs "$srr/$srr.sra"

    rm -r "$srr"

done < $1
