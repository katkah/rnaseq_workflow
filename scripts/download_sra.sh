#!/bin/bash
#PBS -l select=1:ncpus=8:mem=35gb:scratch_local=150gb
#PBS -l walltime=12:00:00
#PBS -N fasterq_dump

# Set umask for proper group permissions
umask 007

# Set all possible temp directories to scratch (before any tool loading)
export TMPDIR=$SCRATCH
export TMP=$SCRATCH
export TEMP=$SCRATCH
export NCBI_HOME=$SCRATCH
export VDB_CONFIG=$SCRATCH

fasterq_dump_path="/path/to/sratoolkit/bin/fasterq-dump"

output_dir="/path/to/output/raw"
srr_list="/path/to/srr_list.txt"

# Create output directory and work in scratch first
mkdir -p $output_dir
cd $SCRATCH
cp $srr_list $SCRATCH/
srr_list="$SCRATCH/$(basename $srr_list)"

echo "Starting SRA downloads in scratch directory: $SCRATCH"
echo "Final output will be synced to: $output_dir"

while read srr; do
    echo "Processing $srr..."

    # Download to scratch for performance
    $fasterq_dump_path $srr --split-files --threads 8

    # Compress the files
    gzip ${srr}_*.fastq

    # Sync files to project directory with proper group permissions
    for file in ${srr}_*.fastq.gz; do
        if [[ -f "$file" ]]; then
            echo "Syncing $file to $output_dir"
            sync_with_group your_group "$file" "$output_dir/$file"
        fi
    done

    # Clean up scratch files
    rm -f ${srr}_*.fastq.gz

    echo "Finished processing $srr"
done < $srr_list

echo "All downloads completed and synced with proper permissions!"

clean_scratch
