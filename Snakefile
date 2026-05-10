# RNA-seq Snakemake Workflow
# Author: katkah
# Description: FastQC -> MultiQC -> STAR alignment -> RSEM quantification

import pandas as pd


# Load configuration
configfile: "config.yaml"  # Snakemake loads YAML → creates `config` dict


samples_df = pd.read_csv(config["samples"], sep="\t")  # Access `config` dict
SAMPLES = samples_df["sample"].tolist()
CHUNKS = config["chunks"]


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
print(
    f"Detected data type: {'Paired-end' if IS_PAIRED else 'Single-end'} (all {len(SAMPLES)} samples)"
)


# Constrain {sample} wildcard globally to only match known sample names.
# Without this, the SE split rule pattern (split/{sample}.part_{chunk}.fastq.gz)
# matches PE split files (split/SAMPLE1_1.part_001.fastq.gz) with
# sample=SAMPLE1_1, which doesn't exist in samples.tsv and causes an IndexError.
# Anchoring {sample} to the actual sample list prevents any false wildcard matches.
wildcard_constraints:
    sample="|".join(SAMPLES),


# star_align_pe and star_align_se produce identical output paths so ruleorder
# is still needed to tell Snakemake which rule to use based on data type.
if IS_PAIRED:

    ruleorder: star_align_pe > star_align_se

else:

    ruleorder: star_align_se > star_align_pe


def get_fastqc_outputs(wildcards):
    if IS_PAIRED:
        return expand(
            f"{config['output_dir']}/fastqc/{{sample}}_{{read}}_fastqc.html",
            sample=SAMPLES,
            read=[1, 2],
        )
    else:
        return expand(
            f"{config['output_dir']}/fastqc/{{sample}}_fastqc.html", sample=SAMPLES
        )


def get_fastqc_trimmed_outputs(wildcards):
    if IS_PAIRED:
        return expand(
            f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R{{read}}_fastqc.html",
            sample=SAMPLES,
            read=[1, 2],
        )
    else:
        return expand(
            f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_fastqc.html",
            sample=SAMPLES,
        )


def get_fastqc_zip_outputs(wildcards):
    if IS_PAIRED:
        return expand(
            f"{config['output_dir']}/fastqc/{{sample}}_{{read}}_fastqc.zip",
            sample=SAMPLES,
            read=[1, 2],
        )
    else:
        return expand(
            f"{config['output_dir']}/fastqc/{{sample}}_fastqc.zip", sample=SAMPLES
        )


def get_fastqc_trimmed_zip_outputs(wildcards):
    if IS_PAIRED:
        return expand(
            f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R{{read}}_fastqc.zip",
            sample=SAMPLES,
            read=[1, 2],
        )
    else:
        return expand(
            f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_fastqc.zip",
            sample=SAMPLES,
        )


def get_copy_inputs(wildcards):
    sample_data = samples_df[samples_df["sample"] == wildcards.sample].iloc[0]
    inputs = [sample_data["fq1"]]
    if IS_PAIRED:
        inputs.append(sample_data["fq2"])
    return inputs


# Define final outputs - this tells Snakemake what we want to produce
rule all:
    input:
        # FastQC reports - raw reads (dynamic based on data type)
        get_fastqc_outputs,
        # SortMeRNA rRNA-depleted reads (dynamic based on data type)
        get_fastqc_trimmed_outputs,
        # Fastp trimmed files (dynamic based on data type)
        f"{config['output_dir']}/multiqc/multiqc_report.html",
        # STAR alignments (same for both data types)
        expand(
            f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.sortedByCoord.out.bam",
            sample=SAMPLES,
        ),
        # STAR transcriptome alignments (same for both data types)
        expand(
            f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.toTranscriptome.out.bam",
            sample=SAMPLES,
        ),
        # STAR gene counts (same for both data types)
        expand(
            f"{config['output_dir']}/star/{{sample}}/{{sample}}ReadsPerGene.out.tab",
            sample=SAMPLES,
        ),
        # RSEM quantification - genes (same for both data types)
        expand(f"{config['output_dir']}/rsem/{{sample}}.genes.results", sample=SAMPLES),
        # RSEM quantification - isoforms (same for both data types)
        expand(
            f"{config['output_dir']}/rsem/{{sample}}.isoforms.results", sample=SAMPLES
        ),


rule split_reads_pe:
    input:
        get_copy_inputs,
    output:
        expand("split/{{sample}}_1.part_{c}.fastq.gz", c=CHUNKS),
        expand("split/{{sample}}_2.part_{c}.fastq.gz", c=CHUNKS),
    threads: 1
    params:
        n_chunks=len(CHUNKS),
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][seqkit_env]}

        mkdir -p split
        seqkit split2 -p {params.n_chunks} -O split -1 {input[0]} -2 {input[1]} -f
        """


rule fastqc_pe:
    input:
        get_copy_inputs,
    output:
        html=f"{config['output_dir']}/fastqc/{{sample}}_1_fastqc.html",
        zip=f"{config['output_dir']}/fastqc/{{sample}}_1_fastqc.zip",
        html2=f"{config['output_dir']}/fastqc/{{sample}}_2_fastqc.html",
        zip2=f"{config['output_dir']}/fastqc/{{sample}}_2_fastqc.zip",
    log:
        f"{config['logs_dir']}/fastqc/{{sample}}.log",
    threads: config["fastqc"]["threads"]
    resources:
        mem_gb=config["fastqc"]["mem_gb"],
    params:
        outdir=f"{config['output_dir']}/fastqc",
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][fastqc]}

        mkdir -p {params.outdir}

        fastqc -t {threads} -o {params.outdir} {input} 2> {log}
        """


rule fastqc_se:
    input:
        get_copy_inputs,
    output:
        html=f"{config['output_dir']}/fastqc/{{sample}}_fastqc.html",
        zip=f"{config['output_dir']}/fastqc/{{sample}}_fastqc.zip",
    log:
        f"{config['logs_dir']}/fastqc/{{sample}}.log",
    threads: config["fastqc"]["threads"]
    resources:
        mem_gb=config["fastqc"]["mem_gb"],
    params:
        outdir=f"{config['output_dir']}/fastqc",
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][fastqc]}

        mkdir -p {params.outdir}

        fastqc -t {threads} -o {params.outdir} {input} 2> {log}
        """


# Rule 2: MultiQC - Aggregate all FastQC reports
rule multiqc:
    input:
        # Include both raw and trimmed FastQC reports (dynamic based on data type)
        raw_fastqc=get_fastqc_zip_outputs,
        trimmed_fastqc=get_fastqc_trimmed_zip_outputs,
    output:
        f"{config['output_dir']}/multiqc/multiqc_report.html",
    log:
        f"{config['logs_dir']}/multiqc/multiqc.log",
    resources:
        mem_gb=config["multiqc"]["mem_gb"],
    params:
        raw_indir=f"{config['output_dir']}/fastqc",
        trimmed_indir=f"{config['output_dir']}/fastqc_trimmed",
        outdir=f"{config['output_dir']}/multiqc",
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][multiqc_env]}

        # Run MultiQC on both directories - temp files automatically go to $SCRATCH
        multiqc --force {params.raw_indir} {params.trimmed_indir} -o {params.outdir} 2> {log}

        """


rule sortmerna_pe:
    input:
        r1="split/{sample}_1.part_{chunk}.fastq.gz",
        r2="split/{sample}_2.part_{chunk}.fastq.gz",
    output:
        f"{config['output_dir']}/sortmerna/{{sample}}_{{chunk}}_non_rRNA_fwd.fq.gz",
        f"{config['output_dir']}/sortmerna/{{sample}}_{{chunk}}_non_rRNA_rev.fq.gz",
    log:
        f"{config['logs_dir']}/sortmerna/{{sample}}_{{chunk}}.log",
    shadow:
        "minimal"
    threads: config["sortmerna"]["threads"]
    resources:
        mem_gb=config["sortmerna"]["mem_gb"],
    params:
        database=config["sortmerna"]["database"],
        prefix="{sample}_{chunk}",
        extra=config["sortmerna"]["extra_params"],
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][sortmerna_env]}

        mkdir -p {config[output_dir]}/sortmerna
        cp {params.database} .
        database_name=$(basename {params.database})

        sortmerna --ref $database_name \
                  --reads {input.r1} \
                  --reads {input.r2} \
                  --workdir . \
                  --fastx --paired_out --out2 \
                  --aligned rRNA-reads \
                  --other {params.prefix}_non_rRNA \
                  --threads {threads} \
                  {params.extra} \
                  2> {log}

        mv {params.prefix}_non_rRNA_fwd.fq.gz {config[output_dir]}/sortmerna/
        mv {params.prefix}_non_rRNA_rev.fq.gz {config[output_dir]}/sortmerna/
        """


rule fastp_pe:
    input:
        non_rrna=f"{config['output_dir']}/sortmerna/{{sample}}_{{chunk}}_non_rRNA_fwd.fq.gz",
        non_rrna2=f"{config['output_dir']}/sortmerna/{{sample}}_{{chunk}}_non_rRNA_rev.fq.gz",
    output:
        trimmed=f"{config['output_dir']}/fastp/{{sample}}_{{chunk}}_trimmed_R1.fastq.gz",
        trimmed2=f"{config['output_dir']}/fastp/{{sample}}_{{chunk}}_trimmed_R2.fastq.gz",
        html=f"{config['output_dir']}/fastp/{{sample}}_{{chunk}}_fastp.html",
        json=f"{config['output_dir']}/fastp/{{sample}}_{{chunk}}_fastp.json",
    log:
        f"{config['logs_dir']}/fastp/{{sample}}_{{chunk}}.log",
    shadow:
        "minimal"
    threads: config["fastp"]["threads"]
    resources:
        mem_gb=config["fastp"]["mem_gb"],
    params:
        outdir=f"{config['output_dir']}/fastp",
        extra=f"{config['fastp']['base_params']} {config['fastp']['paired_params']}",
    shell:
        """
        mkdir -p {params.outdir}
        {config[tools][fastp]}\\
            -i {input[0]} \\
            -I {input[1]} \\
            -o {output.trimmed} \\
            -O {output.trimmed2} \\
            -h {output.html} \\
            -j {output.json} \\
            -w {threads} \\
            {params.extra} 2> {log}
        """


rule join_reads_pe:
    input:
        r1=expand(
            f"{config['output_dir']}/fastp/{{{{sample}}}}_{{chunk}}_trimmed_R1.fastq.gz",
            chunk=CHUNKS,
        ),
        r2=expand(
            f"{config['output_dir']}/fastp/{{{{sample}}}}_{{chunk}}_trimmed_R2.fastq.gz",
            chunk=CHUNKS,
        ),
    output:
        r1=f"{config['output_dir']}/joined/{{sample}}_trimmed_R1.fastq.gz",
        r2=f"{config['output_dir']}/joined/{{sample}}_trimmed_R2.fastq.gz",
    log:
        f"{config['logs_dir']}/joined/{{sample}}.log",
    threads: 1
    params:
        outdir=f"{config['output_dir']}/joined",
        logdir=f"{config['logs_dir']}/joined",
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.logdir}
        echo "joining R1 chunks for sample {wildcards.sample}" > {log}
        cat {input.r1} > {output.r1} 2>> {log}
        echo "joining R2 chunks for sample {wildcards.sample}" >> {log}
        cat {input.r2} > {output.r2} 2>> {log}
        echo "Finished joining reads for sample {wildcards.sample}" >> {log}
        """


rule fastqc_trimmed_pe:
    input:
        r1=f"{config['output_dir']}/joined/{{sample}}_trimmed_R1.fastq.gz",
        r2=f"{config['output_dir']}/joined/{{sample}}_trimmed_R2.fastq.gz",
    output:
        html=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R1_fastqc.html",
        zip=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R1_fastqc.zip",
        html2=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R2_fastqc.html",
        zip2=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_R2_fastqc.zip",
    log:
        f"{config['logs_dir']}/fastqc_trimmed/{{sample}}.log",
    threads: config["fastqc"]["threads"]
    resources:
        mem_gb=config["fastqc"]["mem_gb"],
    params:
        outdir=f"{config['output_dir']}/fastqc_trimmed",
        logdir=f"{config['logs_dir']}/fastqc_trimmed",
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][fastqc]}

        # Create output directory
        mkdir -p {params.outdir}
        mkdir -p {params.logdir}

        # Run FastQC on trimmed reads - temp files automatically go to $SCRATCH
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} 2> {log}
        """


rule star_align_pe:
    input:
        fastq1=f"{config['output_dir']}/joined/{{sample}}_trimmed_R1.fastq.gz",
        fastq2=f"{config['output_dir']}/joined/{{sample}}_trimmed_R2.fastq.gz",
        genome_dir=config["reference"]["genome_dir"],
        gtf=config["reference"]["gtf"],
    output:
        bam=f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.sortedByCoord.out.bam",
        transcriptome_bam=f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.toTranscriptome.out.bam",
        counts=f"{config['output_dir']}/star/{{sample}}/{{sample}}ReadsPerGene.out.tab",
        sj=f"{config['output_dir']}/star/{{sample}}/{{sample}}SJ.out.tab",
        log=f"{config['output_dir']}/star/{{sample}}/{{sample}}Log.final.out",
    log:
        f"{config['logs_dir']}/star/{{sample}}.log",
    shadow:
        "minimal"
    threads: config["star"]["threads"]
    resources:
        mem_gb=config["star"]["mem_gb"],
    params:
        outdir=f"{config['output_dir']}/star/{{sample}}",
        extra=config["star"]["extra_params"],
    shell:
        r"""
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][star]}

        # Copy genome index to scratch (shadow dir = scratch via --shadow-prefix $SCRATCH)
        mkdir -p genome_index
        cp -r {input.genome_dir}/* genome_index/
        cp {input.gtf} .

        # Create output directory within shadow dir
        mkdir -p {params.outdir}

        # Run STAR - write outputs directly to declared output paths
        STAR --runThreadN {threads} \
            --genomeDir genome_index \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --sjdbGTFfile $(basename {input.gtf}) \
            --outFileNamePrefix {params.outdir}/{wildcards.sample} \
            {params.extra} \
            2> {log}
        """


rule split_reads_se:
    input:
        get_copy_inputs,
    output:
        expand("split/{{sample}}.part_{c}.fastq.gz", c=CHUNKS),
    threads: 1
    params:
        n_chunks=len(CHUNKS),
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][seqkit_env]}

        mkdir -p split
        seqkit split2 -p {params.n_chunks} -O split {input[0]} -f
        """


rule sortmerna_se:
    input:
        "split/{sample}.part_{chunk}.fastq.gz",
    output:
        f"{config['output_dir']}/sortmerna/{{sample}}_{{chunk}}_non_rRNA.fq.gz",
    log:
        f"{config['logs_dir']}/sortmerna/{{sample}}_{{chunk}}.log",
    shadow:
        "minimal"
    threads: config["sortmerna"]["threads"]
    resources:
        mem_gb=config["sortmerna"]["mem_gb"],
    params:
        database=config["sortmerna"]["database"],
        prefix="{sample}_{chunk}",
        extra=config["sortmerna"]["extra_params"],
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][sortmerna_env]}

        mkdir -p {config[output_dir]}/sortmerna
        cp {params.database} .
        database_name=$(basename {params.database})

        sortmerna --ref $database_name \
                  --reads {input} \
                  --workdir . \
                  --fastx \
                  --aligned rRNA-reads \
                  --other {params.prefix}_non_rRNA \
                  --threads {threads} \
                  {params.extra} \
                  2> {log}

        mv {params.prefix}_non_rRNA.fq.gz {config[output_dir]}/sortmerna/
        """


rule fastp_se:
    input:
        f"{config['output_dir']}/sortmerna/{{sample}}_{{chunk}}_non_rRNA.fq.gz",
    output:
        trimmed=f"{config['output_dir']}/fastp/{{sample}}_{{chunk}}_trimmed.fastq.gz",
        html=f"{config['output_dir']}/fastp/{{sample}}_{{chunk}}_fastp.html",
        json=f"{config['output_dir']}/fastp/{{sample}}_{{chunk}}_fastp.json",
    log:
        f"{config['logs_dir']}/fastp/{{sample}}_{{chunk}}.log",
    shadow:
        "minimal"
    threads: config["fastp"]["threads"]
    resources:
        mem_gb=config["fastp"]["mem_gb"],
    params:
        outdir=f"{config['output_dir']}/fastp",
        extra=f"{config['fastp']['base_params']} {config['fastp']['single_params']}",
    shell:
        """
        mkdir -p {params.outdir}
        {config[tools][fastp]} \
            -i {input} \
            -o {output.trimmed} \
            -h {output.html} \
            -j {output.json} \
            -w {threads} \
            {params.extra} 2> {log}
        """


rule join_reads_se:
    input:
        expand(
            f"{config['output_dir']}/fastp/{{{{sample}}}}_{{chunk}}_trimmed.fastq.gz",
            chunk=CHUNKS,
        ),
    output:
        f"{config['output_dir']}/joined/{{sample}}_trimmed.fastq.gz",
    log:
        f"{config['logs_dir']}/joined/{{sample}}.log",
    threads: 1
    params:
        outdir=f"{config['output_dir']}/joined",
        logdir=f"{config['logs_dir']}/joined",
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.logdir}
        echo "joining chunks for sample {wildcards.sample}" > {log}
        cat {input} > {output} 2>> {log}
        echo "Finished joining reads for sample {wildcards.sample}" >> {log}
        """


rule fastqc_trimmed_se:
    input:
        f"{config['output_dir']}/joined/{{sample}}_trimmed.fastq.gz",
    output:
        html=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_fastqc.html",
        zip=f"{config['output_dir']}/fastqc_trimmed/{{sample}}_trimmed_fastqc.zip",
    log:
        f"{config['logs_dir']}/fastqc_trimmed/{{sample}}.log",
    threads: config["fastqc"]["threads"]
    resources:
        mem_gb=config["fastqc"]["mem_gb"],
    params:
        outdir=f"{config['output_dir']}/fastqc_trimmed",
        logdir=f"{config['logs_dir']}/fastqc_trimmed",
    shell:
        """
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][fastqc]}

        mkdir -p {params.outdir}
        mkdir -p {params.logdir}

        fastqc -t {threads} -o {params.outdir} {input} 2> {log}
        """


rule star_align_se:
    input:
        fastq=f"{config['output_dir']}/joined/{{sample}}_trimmed.fastq.gz",
        genome_dir=config["reference"]["genome_dir"],
        gtf=config["reference"]["gtf"],
    output:
        bam=f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.sortedByCoord.out.bam",
        transcriptome_bam=f"{config['output_dir']}/star/{{sample}}/{{sample}}Aligned.toTranscriptome.out.bam",
        counts=f"{config['output_dir']}/star/{{sample}}/{{sample}}ReadsPerGene.out.tab",
        sj=f"{config['output_dir']}/star/{{sample}}/{{sample}}SJ.out.tab",
        log=f"{config['output_dir']}/star/{{sample}}/{{sample}}Log.final.out",
    log:
        f"{config['logs_dir']}/star/{{sample}}.log",
    shadow:
        "minimal"
    threads: config["star"]["threads"]
    resources:
        mem_gb=config["star"]["mem_gb"],
    params:
        outdir=f"{config['output_dir']}/star/{{sample}}",
        # --quantTranscriptomeBan is intentionally excluded for single-end data.
        # IndelSoftclipSingleend (used in pe_extra_params) bans unpaired mate alignments,
        # which would discard ALL reads from a single-end transcriptome BAM and break RSEM.
        extra=" ".join(
            p
            for p in config["star"]["extra_params"].split()
            if p != "--quantTranscriptomeBan" and p != "IndelSoftclipSingleend"
        ),
    shell:
        r"""
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][star]}

        mkdir -p genome_index
        cp -r {input.genome_dir}/* genome_index/
        cp {input.gtf} .

        mkdir -p {params.outdir}

        STAR --runThreadN {threads} \
            --genomeDir genome_index \
            --readFilesIn {input.fastq} \
            --sjdbGTFfile $(basename {input.gtf}) \
            --outFileNamePrefix {params.outdir}/{wildcards.sample} \
            {params.extra} \
            2> {log}
        """


# Rule 6: RSEM - Prepare reference index from genome and GTF
rule rsem_prepare_reference:
    input:
        genome_fasta=config["reference"]["genome_fasta"],
        gtf=config["reference"]["gtf"],
    output:
        # RSEM creates multiple index files with different extensions
        seq=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.seq",
        grp=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.grp",
        ti=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.ti",
        idx_fa=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.idx.fa",
        transcripts_fa=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.transcripts.fa",
        chrlist=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.chrlist",
        n2g_idx_fa=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.n2g.idx.fa",
    log:
        f"{config['logs_dir']}/rsem/prepare_reference.log",
    shadow:
        "minimal"
    threads: config["rsem"]["threads"]
    resources:
        mem_gb=config["rsem"]["mem_gb"],
    params:
        index_dir=config["reference"]["rsem_index"],
        index_prefix=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}",
    shell:
        """
        # Load modules and activate environment
        {config[modules][loadmodules]}
        {config[modules][mambaforge]}
        {config[modules][rsem_env]}

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
        index_grp=f"{config['reference']['rsem_index']}/{config['reference']['rsem_index_name']}.grp",
    output:
        genes=f"{config['output_dir']}/rsem/{{sample}}.genes.results",
        isoforms=f"{config['output_dir']}/rsem/{{sample}}.isoforms.results",
    log:
        f"{config['logs_dir']}/rsem/{{sample}}.log",
    shadow:
        "minimal"
    threads: config["rsem"]["threads"]
    resources:
        mem_gb=config["rsem"]["mem_gb"],
    params:
        outdir=f"{config['output_dir']}/rsem",
        prefix=f"{config['output_dir']}/rsem/{{sample}}",
        index_dir=config["reference"]["rsem_index"],
        index_name=config["reference"]["rsem_index_name"],
        strandedness=config["rsem"]["strandedness"],
        extra=config["rsem"]["extra_params"],
        seed="12345",  # Fixed seed for reproducibility
        is_paired="--paired-end" if IS_PAIRED else "--single-end",
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
            {params.is_paired} \\
            --seed {params.seed} \\
            --strandedness {params.strandedness} \\
            {params.extra} \\
            ../{input.bam} rsem_index/{params.index_name} {wildcards.sample} 2> ../{log}

        # Copy outputs
        mkdir -p ../{params.outdir}
        cp {wildcards.sample}.genes.results ../{params.outdir}/
        cp {wildcards.sample}.isoforms.results ../{params.outdir}/
        """
