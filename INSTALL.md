# Installation Guide

This pipeline runs on Metacentrum (Czech HPC) using PBS. Each tool runs in its own
conda environment. Some tools required source compilation due to issues with conda
packages — see notes per tool.

## Why not `--use-conda`

Snakemake supports automatic conda environment creation via `--use-conda`, where each
rule declares a `conda:` directive pointing to an environment YAML. This is not used
here because Metacentrum requires PBS modules to be loaded before conda is available
(`source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules`, `module add mambaforge`).
These module commands run inside the shell block of each rule. If `--use-conda` were
used, Snakemake would try to manage conda itself before those modules are loaded,
which fails. Instead, environments are activated manually inside each rule's shell block
after the modules are loaded.

## Why fastp and RSEM are not installed via conda

**fastp:** the conda package for the version needed (v1.3.0) had runtime issues on
Metacentrum. Compiling from the official GitHub source produces a static binary that
works reliably. The binary path is set directly in `config.yaml` under `tools.fastp`.

**RSEM:** similarly, the conda RSEM package had dependency conflicts on Metacentrum.
The solution was to create a conda environment that provides the required runtime
libraries and build tools (gcc, samtools, R, boost, etc.), compile RSEM from source
inside that environment, and then call the resulting binaries directly. The conda
environment is still activated before each RSEM rule runs so the shared libraries
are available at runtime — but the binaries themselves come from the source build,
not from conda. The binary directory is set in `config.yaml` under `tools.rsem_dir`.

## Prerequisites

Log in to a Metacentrum node with internet access and request an interactive session:

```bash
qsub -I -l select=1:ncpus=4:mem=16gb -l walltime=3:00:00
```

Then load mambaforge and set conda package directory:

```bash
module add mambaforge
export CONDA_PKGS_DIRS="/storage/praha1/home/$USER/tools/.conda/pkgs"
export TMPDIR=$SCRATCHDIR
mkdir -p "/storage/praha1/home/$USER/tools/.conda/pkgs"
```

---

## Snakemake environment

The main workflow runs inside a Snakemake conda environment. Create it and note the path —
you will need to update `submit_workflow.sh` to activate it.

```bash
mamba create --prefix /path/to/your/envs/snakemake -c conda-forge -c bioconda \
    snakemake=9.17.1 \
    python=3.13 \
    pandas=2.3.3 \
    -y
```

Update `submit_workflow.sh`:
```bash
mamba activate /path/to/your/envs/snakemake
```

---

## Tool environments

### seqkit v2.13.0

```bash
mamba create --prefix /path/to/your/envs/seqkit -c bioconda seqkit -y
```

### FastQC v0.12.1

```bash
mamba create --prefix /path/to/your/envs/fastqc_v0_12_1 \
    -c bioconda -c conda-forge fastqc -y
```

### MultiQC v1.34

```bash
mamba create --prefix /path/to/your/envs/multiqc_v1_34 \
    -c bioconda -c conda-forge multiqc -y
```

### SortMeRNA v4.3.7

```bash
mamba create --prefix /path/to/your/envs/sortmerna_v4 \
    -c conda-forge sortmerna=4.3.7 -y
```

### STAR v2.7.10b

```bash
mamba create --prefix /path/to/your/envs/star_2_7_10b \
    -c bioconda -c conda-forge bioconda::star==2.7.10b -y
```

### fastp v1.3.0 — compiled from source

The conda package had issues so fastp was compiled from the official GitHub source.

```bash
# Install build dependencies
mamba create --prefix /path/to/your/envs/fastp -c conda-forge \
    gcc gxx cmake zlib -y
mamba activate /path/to/your/envs/fastp

# Clone and build
git clone https://github.com/OpenGene/fastp.git
cd fastp
git checkout v1.3.0
make -j 4
```

The resulting binary is at `fastp/fastp`. Update `config.yaml`:
```yaml
tools:
  fastp: "/path/to/fastp/fastp"
```

### RSEM v1.3.1 — compiled from source inside conda environment

The conda package had issues so RSEM was compiled from source. A conda environment
provides the required build dependencies and runtime libraries.

**Step 1 — create the dependency environment:**

```bash
mamba create --prefix /path/to/your/envs/rsem_env \
    -c conda-forge -c bioconda \
    gcc_linux-64 \
    gxx_linux-64 \
    r-base=4.0 \
    python=3.8 \
    make \
    cmake \
    zlib \
    boost-cpp \
    eigen \
    samtools \
    -y
```

**Step 2 — compile RSEM inside that environment:**

```bash
mamba activate /path/to/your/envs/rsem_env

git clone https://github.com/deweylab/RSEM.git
cd RSEM
git checkout v1.3.1
make -j 4
```

The binaries are in the `RSEM/` directory. Update `config.yaml`:
```yaml
tools:
  rsem_dir: "/path/to/RSEM"
```

---

## config.yaml

After setting up all environments, update `config.yaml` with your paths:

```yaml
modules:
  rsem_env: "mamba activate /path/to/your/envs/rsem_env"
  star: "mamba activate /path/to/your/envs/star_2_7_10b"
  fastqc: "mamba activate /path/to/your/envs/fastqc_v0_12_1"
  multiqc_env: "mamba activate /path/to/your/envs/multiqc_v1_34"
  sortmerna_env: "mamba activate /path/to/your/envs/sortmerna_v4"
  seqkit_env: "mamba activate /path/to/your/envs/seqkit"

tools:
  fastp: "/path/to/fastp/fastp"
  rsem_dir: "/path/to/RSEM"
```

Also update `submit_workflow.sh` with the path to your Snakemake environment.
