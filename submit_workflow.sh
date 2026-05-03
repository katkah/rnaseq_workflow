#!/bin/bash
#PBS -N rnaseq_snakemake
#PBS -l select=1:ncpus=32:mem=135gb:scratch_local=250gb
#PBS -l walltime=24:00:00
#PBS -M your@email.com
#PBS -m abe

# PBS script to run Snakemake workflow on PBS cluster
# This script submits the main Snakemake process, which will then submit individual jobs


# Change to the workflow directory
cd $SCRATCH

# Activate Snakemake environment (Metacentrum way)
module add mambaforge
mamba activate /full/path/to/conda/envs/snakemake

# Create necessary directories
mkdir -p logs results raw
mkdir -p logs/{fastqc,sortmerna,fastqc_trimmed,fastp,star,rsem,multiqc}
mkdir -p results/{fastqc,sortmerna,fastqc_trimmed,fastp,star,rsem,multiqc}

# Set scratch directory for all temporary files
export TMPDIR=$SCRATCH

# Copy necessary files to scratch
cp $PBS_O_WORKDIR/config.yaml .
cp $PBS_O_WORKDIR/Snakefile .
cp $PBS_O_WORKDIR/samples.tsv .



# Run Snakemake with proper resource limits
snakemake \
    --jobs 99 \
    --max-threads 32 \
    --resources mem_gb=135 \
    --shadow-prefix $SCRATCH/ \
    --printshellcmds \
    --keep-going \
    --rerun-incomplete \
    --configfile config.yaml

echo "Snakemake workflow completed!"

cd $PBS_O_WORKDIR
# Create necessary directories
mkdir -p logs results
mkdir -p logs/{fastqc,sortmerna,fastqc_trimmed,fastp,star,rsem,multiqc}
mkdir -p results/{fastqc,sortmerna,fastqc_trimmed,fastp,star,rsem,multiqc}


cd $SCRATCH
# copy results back to original working directory
cp -r results/* $PBS_O_WORKDIR/results/
cp -r logs/* $PBS_O_WORKDIR/logs/

#fix permissions for output files
chmod -R 770 $PBS_O_WORKDIR/results
chmod -R 770 $PBS_O_WORKDIR/logs




# Clean up scratch if needed
if [ -d "$SCRATCHDIR" ]; then
    clean_scratch
fi