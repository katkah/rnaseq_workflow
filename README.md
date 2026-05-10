# RNA-seq Workflow

Snakemake pipeline for RNA-seq analysis on Metacentrum (Czech HPC). Supports both
paired-end and single-end data. Runs FastQC → SortMeRNA → fastp → STAR → RSEM,
with MultiQC aggregating QC reports.

---

## Pipeline steps

```
raw reads
   │
   ├── FastQC (raw QC)
   │
   ├── seqkit split2 (split into 3 chunks for parallel processing)
   │      │
   │      ├── SortMeRNA (rRNA depletion, per chunk)
   │      │
   │      └── fastp (adapter and quality trimming, per chunk)
   │
   ├── cat (join chunks back)
   │
   ├── FastQC (trimmed QC)
   │
   ├── STAR (alignment)
   │
   ├── RSEM (quantification)
   │
   └── MultiQC (aggregate all QC reports)
```

Paired-end and single-end data are handled by separate rules throughout.
The pipeline detects data type automatically from `samples.tsv` at startup.

### Why reads are split into 3 chunks

SortMeRNA and fastp are the most CPU-intensive steps and also the most
memory-hungry. Running them on a full sample file sequentially would leave most
CPUs idle while one job runs. Splitting each sample into 3 equal chunks with
`seqkit split2` allows SortMeRNA and fastp to run on all 3 chunks in parallel,
making full use of the allocated CPUs.

After trimming, the chunks are concatenated back into a single file per sample
with `cat` before STAR alignment (which handles the full file itself).

The number of chunks is configurable via `chunks` in `config.yaml`. The seqkit
`-p` parameter is derived automatically from the length of that list, so changing
the chunk count requires editing only one line.

---

## Setup

See [INSTALL.md](INSTALL.md) for tool installation instructions.

### 1. Prepare samples.tsv

Tab-separated file with one row per sample:

**Paired-end:**
```
sample  fq1     fq2
SRR001  /path/to/SRR001_R1.fastq.gz     /path/to/SRR001_R2.fastq.gz
SRR002  /path/to/SRR002_R1.fastq.gz     /path/to/SRR002_R2.fastq.gz
```

**Single-end:**
```
sample  fq1
SRR001  /path/to/SRR001.fastq.gz
SRR002  /path/to/SRR002.fastq.gz
```

All samples in one file must be the same type. Mixed PE/SE datasets are not supported.

### 2. Update config.yaml

Set your paths for reference files, conda environments, and tool binaries:

```yaml
reference:
  genome_dir: "/path/to/STAR_index"
  genome_fasta: "/path/to/genome.fa.gz"
  gtf: "/path/to/annotation.gtf"
  rsem_index: "/path/to/RSEM_index"
  rsem_index_name: "species_name"

modules:
  rsem_env: "mamba activate /path/to/envs/rsem_env"
  star: "mamba activate /path/to/envs/star_2_7_10b"
  fastqc: "mamba activate /path/to/envs/fastqc_v0_12_1"
  multiqc_env: "mamba activate /path/to/envs/multiqc_v1_34"
  sortmerna_env: "mamba activate /path/to/envs/sortmerna_v4"
  seqkit_env: "mamba activate /path/to/envs/seqkit"

tools:
  fastp: "/path/to/fastp/fastp"
  rsem_dir: "/path/to/RSEM"
```

Also update `strandedness` under `rsem:` based on your library prep (`forward`,
`reverse`, or `none`). If unknown, use RSeQC `infer_experiment.py` on a test
alignment to determine it.

### 3. Update submit_workflow.sh

Set the path to your Snakemake conda environment:

```bash
mamba activate /path/to/your/envs/snakemake
```

### 4. Dry run

Always do a dry run before submitting to check the DAG and catch any config errors:

```bash
snakemake -n --configfile config.yaml
```

### 5. Submit

```bash
qsub submit_workflow.sh
```

The script requests 32 CPUs, 135 GB RAM, and 250 GB scratch. Walltime is set to 48h.

---

## How scratch and shadow work

Metacentrum provides a fast local scratch disk (`$SCRATCH`) on each compute node.
Using scratch for intermediate files avoids saturating the network storage and speeds
up I/O-heavy steps.

### shadow: "minimal"

Demanding rules (SortMeRNA, fastp, STAR, RSEM) use Snakemake's `shadow: "minimal"`.
When a rule runs with shadow, Snakemake:

1. Creates a temporary working directory inside `$SCRATCH`
2. Sets that as the working directory for the rule
3. Symlinks declared input files into the shadow directory
4. After the rule finishes, copies declared outputs back to their real paths

The `--shadow-prefix $SCRATCH/` in `submit_workflow.sh` tells Snakemake where to
create shadow directories.

This means each demanding rule runs entirely in scratch — temporary files (genome
index copies, sortmerna kvdb/idx, fastp temp files) never touch network storage.

### Rules without shadow

Light rules (split, join, FastQC, MultiQC) run directly in `$PBS_O_WORKDIR` (the
directory from which the job was submitted). These produce small outputs that go
directly to `results/` on network storage.

### Summary

| Rule | Shadow | Working directory |
|------|--------|-------------------|
| split_reads | no | PBS_O_WORKDIR |
| fastqc | no | PBS_O_WORKDIR |
| sortmerna | yes | $SCRATCH |
| fastp | yes | $SCRATCH |
| join_reads | no | PBS_O_WORKDIR |
| fastqc_trimmed | no | PBS_O_WORKDIR |
| star_align | yes | $SCRATCH |
| rsem_prepare_reference | yes | $SCRATCH |
| rsem_quantify | yes | $SCRATCH |
| multiqc | no | PBS_O_WORKDIR |

---

## Resource management

### PBS allocation

`submit_workflow.sh` requests:
- 32 CPUs
- 135 GB RAM
- 250 GB scratch_local
- 48h walltime

These resources are shared across all jobs running concurrently within the PBS job.

### Snakemake resource limits

```bash
snakemake \
    --jobs 2 \
    --cores 32 \
    --max-threads 8 \
    --resources mem_gb=135 \
    ...
```

- `--jobs 2` — run at most 2 rules simultaneously
- `--cores 32` — total CPU budget across all concurrent rules
- `--max-threads 8` — cap threads per rule at 8
- `--resources mem_gb=135` — total memory budget; Snakemake will not start a new
  rule if doing so would exceed this limit

### Per-rule memory

Each rule declares `mem_gb` in its `resources:` block. Snakemake sums these across
concurrent jobs and will not schedule a new job if the total would exceed
`--resources mem_gb=135`.

Current settings:

| Rule | mem_gb | threads |
|------|--------|---------|
| sortmerna | 50 | 4 |
| fastp | 32 | 8 |
| star_align | 32 | 8 |
| rsem_quantify | 48 | 8 |
| rsem_prepare_reference | 48 | 8 |
| fastqc | 8 | 8 |
| multiqc | 4 | 8 |

With `--jobs 2`, two sortmerna jobs can run simultaneously (2 × 50 GB = 100 GB),
staying within the 135 GB budget. Two STAR jobs would use 64 GB — fine. Two RSEM
jobs would use 96 GB — also within budget. Snakemake enforces this automatically.

### Adjusting resources

If you get OOM kills (PBS sends SIGTERM to the whole job), either:
- Reduce concurrent jobs by setting `--jobs 1`
- Lower `mem_gb` for the offending rule in `config.yaml`
- Increase the PBS `mem=` request and `--resources mem_gb=` accordingly

If jobs are too slow, increase `threads` per rule (up to `--max-threads`) and
adjust `--cores` and PBS `ncpus` to match.
