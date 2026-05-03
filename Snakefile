# RNA-seq Snakemake Workflow
# Author: Learning Snakemake
# Description: FastQC -> MultiQC -> STAR alignment -> RSEM quantification

import pandas as pd

# Load configuration
configfile: "config.yaml"          # Snakemake loads YAML → creates `config` dict
samples_df = pd.read_csv(config["samples"], sep="\t")  # Access `config` dict
SAMPLES = samples_df["sample"].tolist()

# Step 1: Auto-detect single-end vs paired-end data
def detect_data_type(df):
    """Detect if data is single-end or paired-end (all samples must be same type)"""
    has_fq2_column = "fq2" in df.columns
    has_fq2_values = has_fq2_column and not df["fq2"].isna().all()
    return has_fq2_values

def validate_consistent_data_type(df):
    """Ensure all samples are the same type (no mixing single-end and paired-end)"""
    if "fq2" not in df.columns:
        return True  # All single-end (no fq2 column)

    # Check if ALL samples have fq2 values or ALL samples have empty fq2
    all_paired = df["fq2"].notna().all()
    all_single = df["fq2"].isna().all()

    if not (all_paired or all_single):
        # Mixed dataset detected - not allowed
        paired_samples = df[df["fq2"].notna()]["sample"].tolist()
        single_samples = df[df["fq2"].isna()]["sample"].tolist()
        raise ValueError(
            f"Mixed datasets not supported! All samples must be the same type (paired-end or single-end).\n"
            f"Paired-end samples: {paired_samples}\n"
            f"Single-end samples: {single_samples}\n"
            f"Please create separate sample files for different data types."
        )

    return True

# Set global variable for data type
IS_PAIRED = detect_data_type(samples_df)
validate_consistent_data_type(samples_df)  # Enforce consistency
print(f"Detected data type: {'Paired-end' if IS_PAIRED else 'Single-end'} (all {len(SAMPLES)} samples)")

# Step 2: Helper functions for dynamic outputs based on data type
def get_fastqc_outputs():
    """Generate appropriate FastQC output file names based on data type"""
    if IS_PAIRED:
        # Paired-end: separate reports for R1 and R2
        return expand(f"{config['output_dir']}/fastqc/{{sample}}_{{read}}_fastqc.html",
                     sample=SAMPLES, read=[1, 2])
    else:
        # Single-end: one report per sample
        return expand(f"{config['output_dir']}/fastqc/{{sample}}_fastqc.html",
                     sample=SAMPLES)

def get_fastqc_trimmed_outputs():
    """Generate appropriate trimmed FastQC output file names based on data type"""
    if IS_PAIRED:
        # Paired-end: separate reports for R1 and R2 (matching FastQC naming convention)
        return expand(f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R{{read}}_fastqc.html",
                     sample=SAMPLES, read=[1, 2])
    else:
        # Single-end: one report per sample
        return expand(f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_fastqc.html",
                     sample=SAMPLES)

def get_fastp_outputs():
    """Generate appropriate fastp output file names based on data type"""
    if IS_PAIRED:
        # Paired-end: separate trimmed files for R1 and R2
        return expand(f"{config['output_dir']}/fastp/{{sample}}_trimmed_R{{read}}.fastq.gz",
                     sample=SAMPLES, read=[1, 2])
    else:
        # Single-end: one trimmed file per sample
        return expand(f"{config['output_dir']}/fastp/{{sample}}_trimmed.fastq.gz", sample=SAMPLES)

def get_fastqc_zip_outputs():
    """Generate appropriate FastQC zip file names based on data type"""
    if IS_PAIRED:
        # Paired-end: separate zip files for R1 and R2
        return expand(f"{config['output_dir']}/fastqc/{{sample}}_{{read}}_fastqc.zip",
                     sample=SAMPLES, read=[1, 2])
    else:
        # Single-end: one zip file per sample
        return expand(f"{config['output_dir']}/fastqc/{{sample}}_fastqc.zip",
                     sample=SAMPLES)

def get_fastqc_trimmed_zip_outputs():
    """Generate appropriate trimmed FastQC zip file names based on data type"""
    if IS_PAIRED:
        # Paired-end: separate zip files for R1 and R2 (matching FastQC naming convention)
        return expand(f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R{{read}}_fastqc.zip",
                     sample=SAMPLES, read=[1, 2])
    else:
        # Single-end: one zip file per sample
        return expand(f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_fastqc.zip",
                     sample=SAMPLES)

def get_sortmerna_outputs():
    """Generate appropriate SortMeRNA output file names based on data type"""
    if IS_PAIRED:
        # Paired-end: separate non-rRNA files for R1 and R2
        return expand(f"{config['output_dir']}/sortmerna/{{sample}}_non_rRNA_fwd.fq.gz", sample=SAMPLES) + \
               expand(f"{config['output_dir']}/sortmerna/{{sample}}_non_rRNA_rev.fq.gz", sample=SAMPLES)
    else:
        # Single-end: one non-rRNA file per sample
        return expand(f"{config['output_dir']}/sortmerna/{{sample}}_non_rRNA.fq.gz", sample=SAMPLES)

# Define final outputs - this tells Snakemake what we want to produce
rule all:
    input:
        # FastQC reports - raw reads (dynamic based on data type)
        get_fastqc_outputs(),
        # SortMeRNA rRNA-depleted reads (dynamic based on data type)
        get_sortmerna_outputs(),
        # FastQC reports - trimmed reads (dynamic based on data type)
        get_fastqc_trimmed_outputs(),
        # Fastp trimmed files (dynamic based on data type)
        get_fastp_outputs(),
        # MultiQC report (same for both data types)
        f"{config['output_dir']}/multiqc/multiqc_report.html",
        # STAR alignments (same for both data types)
        expand(f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        # STAR transcriptome alignments (same for both data types)
        expand(f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.toTranscriptome.out.bam", sample=SAMPLES),
        # STAR gene counts (same for both data types)
        expand(f"{config['output_dir']}/star/{{sample}}/{{sample}}ReadsPerGene.out.tab", sample=SAMPLES),
        # RSEM quantification - genes (same for both data types)
        expand(f"{config['output_dir']}/rsem/{{sample}}.genes.results", sample=SAMPLES),
        # RSEM quantification - isoforms (same for both data types)
        expand(f"{config['output_dir']}/rsem/{{sample}}.isoforms.results", sample=SAMPLES)



def get_trimmed_fastq_inputs(wildcards):
    """Get appropriate trimmed FASTQ files for each sample based on data type (from fastp output)"""
    if IS_PAIRED:
        # Paired-end: expect R1 and R2 trimmed files from fastp
        return [f"{config['output_dir']}/fastp/{wildcards.sample}_trimmed_R1.fastq.gz",
                f"{config['output_dir']}/fastp/{wildcards.sample}_trimmed_R2.fastq.gz"]
    else:
        # Single-end: expect single trimmed file from fastp
        return [f"{config['output_dir']}/fastp/{wildcards.sample}_trimmed.fastq.gz"]

def get_sortmerna_inputs(wildcards):
    """Get appropriate SortMeRNA output files for each sample based on data type (for fastp input)"""
    if IS_PAIRED:
        # Paired-end: expect non-rRNA files from SortMeRNA
        return [f"{config['output_dir']}/sortmerna/{wildcards.sample}_non_rRNA_fwd.fq.gz",
                f"{config['output_dir']}/sortmerna/{wildcards.sample}_non_rRNA_rev.fq.gz"]
    else:
        # Single-end: expect single non-rRNA file from SortMeRNA
        return [f"{config['output_dir']}/sortmerna/{wildcards.sample}_non_rRNA.fq.gz"]

def get_copy_inputs(wildcards):
    """Get appropriate input raw files for each sample based on data type"""
    sample_data = samples_df[samples_df["sample"] == wildcards.sample].iloc[0]
    inputs = [sample_data["fq1"]]

    # Add fq2 if this is paired-end data
    if IS_PAIRED:
        inputs.append(sample_data["fq2"])

    return inputs

def get_copy_outputs(wildcards):
    """Get appropriate output file names for copied raw files in scratch based on data type"""
    if IS_PAIRED:
        # Paired-end: separate files for R1 and R2
        return [f"raw/{wildcards.sample}_1.fastq.gz",
                f"raw/{wildcards.sample}_2.fastq.gz"]
    else:
        # Single-end: one file per sample
        return [f"raw/{wildcards.sample}.fastq.gz"]

def get_sortmerna_inputs(wildcards):
    if IS_PAIRED:
        return [
            f"{config['output_dir']}/sortmerna/{wildcards.sample}_non_rRNA_fwd.fq.gz",
            f"{config['output_dir']}/sortmerna/{wildcards.sample}_non_rRNA_rev.fq.gz"
        ]
    else:
        return [f"{config['output_dir']}/sortmerna/{wildcards.sample}_non_rRNA.fq.gz"]

# Rule 1: copy_data_to_scratch - Copy raw FASTQ files to scratch for performance (dynamic based on data type)
rule copy_data_to_scratch:
    input:
        get_copy_inputs
    output:
        "raw/{sample}_1.fastq.gz" if IS_PAIRED else "raw/{sample}.fastq.gz",
        "raw/{sample}_2.fastq.gz" if IS_PAIRED else []
    shell:
        """
        # Create scratch directory for raw data if it doesn't exist
        mkdir -p ./raw

        # Copy input files to scratch for performance boost
        cp {input[0]} ./raw/
        """ +
        ("cp {input[1]} ./raw/" if IS_PAIRED else "") +
        """
        """

#Rule 2: FastQC on raw reads - Quality control before any processing (dynamic based on data type)
rule fastqc:
    input:
        get_copy_outputs
    output:
        html=f"{config['output_dir']}/fastqc/{{sample}}_1_fastqc.html" if IS_PAIRED else f"{config['output_dir']}/fastqc/{{sample}}_fastqc.html",
        zip=f"{config['output_dir']}/fastqc/{{sample}}_1_fastqc.zip" if IS_PAIRED else f"{config['output_dir']}/fastqc/{{sample}}_fastqc.zip",
        html2=f"{config['output_dir']}/fastqc/{{sample}}_2_fastqc.html" if IS_PAIRED else [],
        zip2=f"{config['output_dir']}/fastqc/{{sample}}_2_fastqc.zip" if IS_PAIRED else []
    params:
        outdir=f"{config['output_dir']}/fastqc"
    threads: config["fastqc"]["threads"]
    log:
        f"{config['logs_dir']}/fastqc/{{sample}}.log"
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][fastqc]}

        # Create output directory
        mkdir -p {params.outdir}

        # Run FastQC - temp files automatically go to $SCRATCH
        fastqc -t {threads} -o {params.outdir} {input[0]}""" +
        (" {input[1]}" if IS_PAIRED else "") +
        """ 2> {log}
        """

# Rule 2: MultiQC - Aggregate all FastQC reports
rule multiqc:
    input:
        # Include both raw and trimmed FastQC reports (dynamic based on data type)
        raw_fastqc=get_fastqc_zip_outputs(),
        trimmed_fastqc=get_fastqc_trimmed_zip_outputs()
    output:
        f"{config['output_dir']}/multiqc/multiqc_report.html"
    params:
        raw_indir=f"{config['output_dir']}/fastqc",
        trimmed_indir=f"{config['output_dir']}/fastqc_trimmed",
        outdir=f"{config['output_dir']}/multiqc"
    resources:
        mem_gb=config["multiqc"]["mem_gb"]
    log:
        f"{config['logs_dir']}/multiqc/multiqc.log"
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][multiqc_env]}

        # Run MultiQC on both directories - temp files automatically go to $SCRATCH
        multiqc {params.raw_indir} {params.trimmed_indir} -o {params.outdir} 2> {log}

        """

# Rule 2.5: SortMeRNA - rRNA removal
rule sortmerna:
    input:
        get_copy_outputs
    output:
        non_rrna=f"{config['output_dir']}/sortmerna/{{sample}}_non_rRNA_fwd.fq.gz" if IS_PAIRED else f"{config['output_dir']}/sortmerna/{{sample}}_non_rRNA.fq.gz",
        non_rrna2=f"{config['output_dir']}/sortmerna/{{sample}}_non_rRNA_rev.fq.gz" if IS_PAIRED else []
    params:
        outdir=f"{config['output_dir']}/sortmerna",
        database=config["sortmerna"]["database"],
        workdir=f"sortmerna_work_{{sample}}",
        prefix="{sample}_non_rRNA"
    threads: config["sortmerna"]["threads"]
    resources:
        mem_gb=config["sortmerna"]["mem_gb"]
    log:
        f"{config['logs_dir']}/sortmerna/{{sample}}.log"
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][sortmerna_env]}


        # Create unique work directory for this sample to avoid job interference

        mkdir -p {params.workdir}
        mkdir -p {params.workdir}_temp
        cd {params.workdir}

        # Copy database to work directory
        cp {params.database} ./
        database_name=$(basename {params.database})

        # Run SortMeRNA from within unique work directory""" +
        ("""
        sortmerna --ref ./$database_name \\
                  --reads ../{input[0]} \\
                  --reads ../{input[1]} \\
                  --workdir ../{params.workdir}_temp \\
                  --fastx --paired_out --out2 \\
                  --aligned rRNA-reads \\
                  --other {params.prefix} \\
                  --threads {threads} 2> ../{log}

        # Copy all output files explicitly to avoid wildcard issues
        mv {params.prefix}_fwd.fq.gz ../{params.outdir}/
        mv {params.prefix}_rev.fq.gz ../{params.outdir}/
        # Clean up work directories


        """ if IS_PAIRED else """
        sortmerna --ref ./$database_name \\
                  --reads ../{input[0]} \\
                  --workdir ../{params.workdir}_temp \\
                  --fastx \\
                  --aligned rRNA-reads \\
                  --other {params.prefix} \\
                  --threads {threads} 2> ../{log}

        # Copy all output files explicitly to avoid wildcard issues
        mv {params.prefix}.fq.gz ../{params.outdir}/
        # Clean up work directories


        """) +
        """

        """

# Rule 3: FASTP - Quality control and adapter trimming
rule fastp:
    input:
        get_sortmerna_inputs
    output:
        trimmed=f"{config['output_dir']}/fastp/{{sample}}_trimmed_R1.fastq.gz" if IS_PAIRED else f"{config['output_dir']}/fastp/{{sample}}_trimmed.fastq.gz",
        trimmed2=f"{config['output_dir']}/fastp/{{sample}}_trimmed_R2.fastq.gz" if IS_PAIRED else [],
        html=f"{config['output_dir']}/fastp/{{sample}}_fastp.html",
        json=f"{config['output_dir']}/fastp/{{sample}}_fastp.json"
    params:
        outdir=f"{config['output_dir']}/fastp",
        extra=f"{config['fastp']['base_params']} {config['fastp']['paired_params'] if IS_PAIRED else config['fastp']['single_params']}"
    threads: config["fastp"]["threads"]
    resources:
        mem_gb=config["fastp"]["mem_gb"]
    log:
        f"{config['logs_dir']}/fastp/{{sample}}.log"
    shell:
        """
        # Run fastp - temp files automatically go to $SCRATCH""" +
        ("""
        {config[tools][fastp]} \\
            -i {input[0]} \\
            -I {input[1]} \\
            -o {output.trimmed} \\
            -O {output.trimmed2} \\
            -h {output.html} \\
            -j {output.json} \\
            -w {threads} \\
            {params.extra} 2> {log}""" if IS_PAIRED else """
        {config[tools][fastp]} \\
            -i {input[0]} \\
            -o {output.trimmed} \\
            -h {output.html} \\
            -j {output.json} \\
            -w {threads} \\
            {params.extra} 2> {log}""") +
        """
        """

# Rule 4: FastQC on trimmed reads - Quality control after fastp
rule fastqc_trimmed:
    input:
        get_trimmed_fastq_inputs
    output:
        html=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R1_fastqc.html" if IS_PAIRED else f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_fastqc.html",
        zip=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R1_fastqc.zip" if IS_PAIRED else f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_fastqc.zip",
        html2=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R2_fastqc.html" if IS_PAIRED else [],
        zip2=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R2_fastqc.zip" if IS_PAIRED else []
    params:
        outdir=f"{config['output_dir']}/fastqc_trimmed"
    threads: config["fastqc"]["threads"]
    resources:
        mem_gb=config["fastqc"]["mem_gb"]
    log:
        f"{config['logs_dir']}/fastqc_trimmed/{{sample}}.log"
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][fastqc]}

        # Create output directory
        mkdir -p {params.outdir}

        # Run FastQC on trimmed reads - temp files automatically go to $SCRATCH
        fastqc -t {threads} -o {params.outdir} {input[0]}""" +
        (" {input[1]}" if IS_PAIRED else "") +
        """ 2> {log}
        """

# Rule 5: STAR - Align reads to reference transcriptome
rule star_align:
    input:
        fastq=get_trimmed_fastq_inputs,
        genome_dir=config["reference"]["genome_dir"],
        gtf=config["reference"]["gtf"]
    output:
        bam=f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.sortedByCoord.out.bam",
        transcriptome_bam=f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.toTranscriptome.out.bam",
        counts=f"{config['output_dir']}/star/{{sample}}/{{sample}}ReadsPerGene.out.tab"
    params:
        outdir=f"{config['output_dir']}/star/{{sample}}",
        prefix=f"{config['output_dir']}/star/{{sample}}/{{sample}}",
        extra=config["star"]["extra_params"]
    threads: config["star"]["threads"]
    resources:
        mem_gb=config["star"]["mem_gb"]
    log:
        f"{config['logs_dir']}/star/{{sample}}.log"
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][star]}

        # Create unique work directory for this sample to avoid job interference
        mkdir -p ./star_work_{wildcards.sample}
        cd ./star_work_{wildcards.sample}

        # Copy genome index to scratch work dir for performance boost
        echo "Copying genome index to scratch..."
        mkdir -p genome_index
        cp -r {input.genome_dir}/* genome_index/

        # Copy GTF file to scratch for performance boost
        cp {input.gtf} ./


        # Run STAR - outputs created in SCRATCH work directory for performance
        STAR --runThreadN {threads} \\
            --genomeDir ./genome_index/ \\
            --readFilesIn ../{input.fastq[0]}""" +
        (" ../{input.fastq[1]}" if IS_PAIRED else "") +
        """ \\
            --sjdbGTFfile ./$(basename {input.gtf}) \\
            --outFileNamePrefix {wildcards.sample} \\
            {params.extra} 2> ../{log}

        # Copy outputs from SCRATCH work directory to final destination
        mkdir -p ../{params.outdir}
        cp {wildcards.sample}Aligned.sortedByCoord.out.bam ../{params.outdir}/
        cp {wildcards.sample}Aligned.toTranscriptome.out.bam ../{params.outdir}/
        cp {wildcards.sample}ReadsPerGene.out.tab ../{params.outdir}/
        cp {wildcards.sample}SJ.out.tab ../{params.outdir}/
        cp {wildcards.sample}Log.* ../{params.outdir}/

        """

# Rule 6: RSEM - Prepare reference index from genome and GTF
rule rsem_prepare_reference:
    input:
        genome_fasta=config["reference"]["genome_fasta"],
        gtf=config["reference"]["gtf"]
    resources:
        mem_gb=config["rsem"]["mem_gb"]
    output:
        # RSEM creates multiple index files with different extensions
        seq=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.seq",
        grp=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.grp",
        ti=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.ti",
        idx_fa=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.idx.fa",
        transcripts_fa=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.transcripts.fa",
        chrlist=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.chrlist",
        n2g_idx_fa=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.n2g.idx.fa"
    params:
        index_dir=config["reference"]["rsem_index"],
        index_prefix=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}"
    threads: config["rsem"]["threads"]
    log:
        f"{config['logs_dir']}/rsem/prepare_reference.log"
    shell:
        """
        # Load modules and activate environment
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][rsem_env]}

        # Verify RSEM is available
        echo "Checking RSEM availability..."
        ls -la {config[tools][rsem_dir]}/rsem-prepare-reference

        # Create index directory
        mkdir -p {params.index_dir}

        # Copy inputs to scratch for performance
        cp {input.genome_fasta} ./
        cp {input.gtf} ./

        # Extract genome file if compressed
        if [[ {input.genome_fasta} == *.gz ]]; then
            gunzip $(basename {input.genome_fasta})
            GENOME_FILE=$(basename {input.genome_fasta} .gz)
        else
            GENOME_FILE=$(basename {input.genome_fasta})
        fi


        # Prepare RSEM reference
        {config[tools][rsem_dir]}/rsem-prepare-reference \\
            --gtf $(basename {input.gtf}) \\
            --num-threads {threads} \\
            $GENOME_FILE \\
            {params.index_prefix} 2> {log}

        """

# Rule 7: RSEM - Quantify gene expression from transcriptome alignments
rule rsem_quantify:
    input:
        bam=f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.toTranscriptome.out.bam",
        # Index file dependencies (triggers rsem_prepare_reference if needed)
        index_seq=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.seq",
        index_grp=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.grp"
    output:
        genes=f"{config['output_dir']}/rsem/{{sample}}.genes.results",
        isoforms=f"{config['output_dir']}/rsem/{{sample}}.isoforms.results"
    params:
        outdir=f"{config['output_dir']}/rsem",
        prefix=f"{config['output_dir']}/rsem/{{sample}}",
        index_dir=config["reference"]["rsem_index"],
        strandedness=config["rsem"]["strandedness"],
        extra=config["rsem"]["extra_params"],
        seed="12345"  # Fixed seed for reproducibility
    threads: config["rsem"]["threads"]
    resources:
        mem_gb=config["rsem"]["mem_gb"]
    log:
        f"{config['logs_dir']}/rsem/{{sample}}.log"
    shell:
        """
        # Load modules and activate environment
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][rsem_env]}

        # Create unique work directory for this sample to avoid job interference
        mkdir -p ./rsem_work_{wildcards.sample}
        cd ./rsem_work_{wildcards.sample}
        mkdir -p rsem_index
        cp {params.index_dir}/* rsem_index/


        # Run RSEM from unique work directory
        {config[tools][rsem_dir]}/rsem-calculate-expression -p {threads} \\
            """ + ("--paired-end \\" if IS_PAIRED else "") + """
            --seed {params.seed} \\
            --strandedness {params.strandedness} \\
            {params.extra} \\
            ../{input.bam} rsem_index/{config[reference][rsem_index_name]} {wildcards.sample} 2> ../{log}

        # Copy outputs
        mkdir -p ../{params.outdir}
        cp {wildcards.sample}.genes.results ../{params.outdir}/
        cp {wildcards.sample}.isoforms.results ../{params.outdir}/

        """