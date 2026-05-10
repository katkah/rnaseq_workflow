#!/bin/bash
#PBS -l select=1:ncpus=6:mem=25gb:scratch_local=30gb
#PBS -l walltime=4:00:00
#PBS -N prepare_index_rsem
umask 007

############################################################################################
### Variables
GENOME=/path/to/genome/genome.fa.gz
GTF=/path/to/genome/annotation.gtf
INDEX_OUTPUT_DIR=/path/to/genome/RSEM_index
INDEX_NAME="species_name"   # Base name for RSEM index files (without extension)
THREADS=6

############################################################################################
### Copy inputs to scratch
cp $GTF $SCRATCH/
cp $GENOME $SCRATCH/
cd $SCRATCH/
GENOME=$(basename $GENOME)
gunzip $GENOME
GENOME=${GENOME%.gz}
GTF=$(basename $GTF)

############################################################################################
### Indexing - RUN ONLY ONCE PER GENOME AND RSEM VERSION
RSEM="/path/to/RSEM"
mkdir $SCRATCH/RSEM_index
$RSEM/rsem-prepare-reference --gtf $GTF --num-threads $THREADS $GENOME $SCRATCH/RSEM_index/$INDEX_NAME

############################################################################################
### Copy results
mkdir -p $INDEX_OUTPUT_DIR

cp -r $SCRATCH/RSEM_index/* $INDEX_OUTPUT_DIR/

chmod -R 770 $INDEX_OUTPUT_DIR
