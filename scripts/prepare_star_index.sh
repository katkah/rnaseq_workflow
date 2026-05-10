#!/bin/bash
#PBS -l select=1:ncpus=6:mem=25gb:scratch_local=30gb
#PBS -l walltime=4:00:00
#PBS -N prepare_index_star

############################################################################################
### Variables
GENOME=/path/to/genome/genome.fa.gz
GTF=/path/to/genome/annotation.gtf
INDEX_OUTPUT_DIR=/path/to/genome/STAR_index
# Read length - 1 (e.g. 150bp reads -> 149)
RD_LENGTH=140
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
### STAR create index
source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module add mambaforge
mamba activate /path/to/envs/star_2_7_10b

mkdir -p $SCRATCH/STAR_index

STAR --runMode genomeGenerate \
    --runThreadN $THREADS \
    --genomeDir $SCRATCH/STAR_index \
    --genomeFastaFiles $GENOME \
    --sjdbGTFfile $GTF \
    --sjdbOverhang $RD_LENGTH \
    --genomeSAindexNbases 11
# NOTE: --genomeSAindexNbases 11 is recommended for small genomes (~60 Mb).
# For larger genomes (e.g. human) the default of 14 is appropriate.
# STAR will warn if the value is too large for the genome size.

############################################################################################
### Copy results
mkdir -p $INDEX_OUTPUT_DIR

cp -r $SCRATCH/STAR_index/* $INDEX_OUTPUT_DIR/
chmod -R 770 $INDEX_OUTPUT_DIR
