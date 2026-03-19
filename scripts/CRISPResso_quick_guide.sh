#!/usr/bin/env bash

###############################################################################
# CRISPResso2 HPC DEBUGGING LOG
#
# Purpose:
#   This script records the installation, testing, and attempted analysis
#   of CRISPResso2 on the HPC, including the troubleshooting steps taken when
#   the real dataset failed to align.
#

###############################################################################

set -u

###############################################################################
# PART 1 — INSTALLATION OF CRISPResso2 ENVIRONMENT
###############################################################################

# Load Conda module available on the HPC
module load python/anaconda/2024.10/3.12.7

# Create a clean Conda environment containing:
#   - python 3.10
#   - cutadapt
#   - fastp
#   - crispresso2
#
# Comment this out if the environment already exists.
conda create -n crispresso_trim \
  -c conda-forge -c bioconda \
  python=3.10 cutadapt fastp crispresso2 -y

# Activate the environment
source activate crispresso_trim

# Check that CRISPResso was installed correctly
CRISPResso --version
CRISPResso --help

###############################################################################
# PART 2 — SYNTHETIC TEST TO CONFIRM CRISPResso2 WORKS
###############################################################################

# Create a dedicated test directory
mkdir -p /gpfs/home/jah14xyu/Postdoc/CRISPResso/crispresso2_test
cd /gpfs/home/jah14xyu/Postdoc/CRISPResso/crispresso2_test || exit 1

# Define a synthetic amplicon and guide
AMPLICONSEQ_TEST="AAACCCAAACCCAAACCCAAACCCAAACCCAAACCCAAACCCAAACCCGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTT"
GUIDESEQ_TEST="CCCAAACCCAAACCCAAACC"

# Construct a minimal FASTQ where the read exactly matches the amplicon
cat > test.fastq <<EOF
@read1
${AMPLICONSEQ_TEST}
+
$(python - <<PY
seq = "${AMPLICONSEQ_TEST}"
print("F" * len(seq))
PY
)
EOF

gzip -f test.fastq

# Run CRISPResso on the synthetic FASTQ
CRISPResso \
  --fastq_r1 test.fastq.gz \
  --amplicon_seq "${AMPLICONSEQ_TEST}" \
  --guide_seq "${GUIDESEQ_TEST}" \
  --output_folder test_output \
  --name synthetic_test \
  --n_processes 1

# Check that output files were generated
ls -lh test_output

###############################################################################
# PART 3 — PREPARATION OF REAL DATA
###############################################################################

# Create a directory to hold the real data
mkdir -p /gpfs/home/jah14xyu/Postdoc/CRISPResso/raw_data

# Clone the GitHub repository containing the example FASTQ files
cd /gpfs/home/jah14xyu/Postdoc/CRISPResso/raw_data || exit 1
git clone https://github.com/asdfjohnny1/CRISPResso2_usage.git

# Reload the Conda module and reactivate the environment
module load python/anaconda/2024.10/3.12.7
source activate crispresso_trim

###############################################################################
# PART 4 — ATTEMPTED ANALYSIS OF 7SK DATASET
###############################################################################

FASTQ_R1="/gpfs/home/jah14xyu/Postdoc/CRISPResso/raw_data/CRISPResso2_usage/data/raw_data/7SK-gRNA_R1_001.fastq.gz"
FASTQ_R2="/gpfs/home/jah14xyu/Postdoc/CRISPResso/raw_data/CRISPResso2_usage/data/raw_data/7SK-gRNA_R2_001.fastq.gz"

AMPLICONSEQ="TTCTGCGATAGCTTTTTCAACACCTGTACGACGCTGTGTGCCATGAACGAGTCCAAACCCGAGGTAGGTTCATCGCAAATTAATAGCGGCGGATCGGTTAACGCCTCCGAGGCGAATGCCAGTCGCTTACGCTCGCCACCCGATAAACCCTTCACCCGACCCGGCACGCCAATCAACGTATTCTGACATTTACCCAGCGAGAGGTCCTGTATAACTTGATCGACACGTTGTACTTTCTGTTTTTGTGTCATATGCCGT"
GUIDESEQ="CCCGATAAACCCTTCACCCG"

FORWARD_PRIMER="TTCTGCGATAGCTTTTTCAACA"
REVERSE_PRIMER="GACAAAAACACAGTATACGGCA"

OUTDIR="/gpfs/home/jah14xyu/Postdoc/CRISPResso/results/7sk"
mkdir -p "$OUTDIR"

# Initial CRISPResso run on the raw paired-end reads
CRISPResso \
  --fastq_r1 "$FASTQ_R1" \
  --fastq_r2 "$FASTQ_R2" \
  --amplicon_seq "$AMPLICONSEQ" \
  --guide_seq "$GUIDESEQ" \
  --output_folder "$OUTDIR" \
  --name "7sk" \
  --n_processes 8

###############################################################################
# PART 5 — TRIMMING ATTEMPTS
###############################################################################

mkdir -p /gpfs/home/jah14xyu/Postdoc/CRISPResso/data/trimmed
cd /gpfs/home/jah14xyu/Postdoc/CRISPResso/data/trimmed || exit 1

# Attempt 1: fastp adapter/quality trimming
fastp \
  -i "$FASTQ_R1" \
  -I "$FASTQ_R2" \
  -o fastp_trimmed_R1.fastq.gz \
  -O fastp_trimmed_R2.fastq.gz \
  --detect_adapter_for_pe \
  --thread 8

# Attempt 2: primer trimming with Cutadapt
cutadapt \
  -g "$FORWARD_PRIMER" \
  -G "$REVERSE_PRIMER" \
  -o primer_trimmed_R1.fastq.gz \
  -p primer_trimmed_R2.fastq.gz \
  "$FASTQ_R1" "$FASTQ_R2"

# Attempt 3: explicit Illumina adapter trimming with Cutadapt
cutadapt \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -o adapter_trimmed_R1.fastq.gz \
  -p adapter_trimmed_R2.fastq.gz \
  "$FASTQ_R1" "$FASTQ_R2"

# Rerun CRISPResso using the adapter-trimmed reads
CRISPResso \
  --fastq_r1 adapter_trimmed_R1.fastq.gz \
  --fastq_r2 adapter_trimmed_R2.fastq.gz \
  --amplicon_seq "$AMPLICONSEQ" \
  --guide_seq "$GUIDESEQ" \
  --output_folder "$OUTDIR" \
  --name "7sk_trimmed" \
  --n_processes 8

###############################################################################
# PART 6 — FORCED SIMPLE PAIRED-END TEST
###############################################################################

CRISPResso \
  --fastq_r1 "$FASTQ_R1" \
  --fastq_r2 "$FASTQ_R2" \
  --amplicon_seq "$AMPLICONSEQ" \
  --output_folder "$OUTDIR/simple_merge_test" \
  --name "7sk_simple_merge" \
  --n_processes 1 \
  --crispresso_merge \
  --min_paired_end_reads_overlap 10

###############################################################################
# PART 7 — DIRECT INSPECTION OF READ CONTENT
###############################################################################

echo "===== FIRST 5 R1 READS ====="
zcat "$FASTQ_R1" | awk 'NR%4==2{print; c++; if(c==5) exit}'

echo "===== FIRST 5 R2 READS ====="
zcat "$FASTQ_R2" | awk 'NR%4==2{print; c++; if(c==5) exit}'

###############################################################################
# PART 8 — TARGETED MATCH CHECKS AGAINST THE SUPPLIED AMPLICON
###############################################################################

echo "===== MATCH CHECKS ====="
echo "R1 reads starting with amplicon start:"
zcat "$FASTQ_R1" | awk 'NR%4==2' | grep -c '^TTCTGCGATAGCTTTTTCAACA'

AMP1="TTCTGCGATAGCTTTTTCAACACCTGTA"
AMP2="GGTTCATCGCAAATTAATAGCGGCGGAT"
AMP3="TGATCGACACGTTGTACTTTCTGTTTTT"

echo "AMP1 matches in R1:"
zcat "$FASTQ_R1" | awk 'NR%4==2' | grep -c "$AMP1"

echo "AMP2 matches in R1:"
zcat "$FASTQ_R1" | awk 'NR%4==2' | grep -c "$AMP2"

echo "AMP3 matches in R2:"
zcat "$FASTQ_R2" | awk 'NR%4==2' | grep -c "$AMP3"

###############################################################################
# PART 9 — READ LENGTH CHECK
###############################################################################

echo "===== READ LENGTHS ====="
echo -n "R1 length: "
zcat "$FASTQ_R1" | awk 'NR%4==2{print length($0); exit}'

echo -n "R2 length: "
zcat "$FASTQ_R2" | awk 'NR%4==2{print length($0); exit}'

###############################################################################
# PART 10 — SECOND DATASET TEMPLATE (U6)
###############################################################################

# Note:
#   This section is a template only and should not be considered validated
#   unless the same sequence-level checks above are also performed.

FASTQ_R1_U6="/gpfs/home/jah14xyu/CRISPResso2/raw_data/raw_fastq/U6-gRNA_R1_001.fastq.gz"
FASTQ_R2_U6="/gpfs/home/jah14xyu/CRISPResso2/raw_data/raw_fastq/U6-gRNA_R2_001.fastq.gz"
AMPLICONSEQ_U6="TTCTGCGATAGCTTTTTCAACACCTGTACGACGCTGTGTGCCATGAACGAGTCCAAACCCGAGGTAGGTTCATCGCAAATTAATAGCGGCGGATCGGTTAACGCCTCCGAGGCGAATGCCAGTCGCTTACGCTCGCCACCCGATAAACCCTTCACCCGACCCGGCACGCCAATCAACGTATTCTGACATTTACCCAGCGAGAGGTCCTGTATAACTTGATCGACACGTTGTACTTTCTGTTTTTGTGTCATATGCCGT"
GUIDESEQ_U6="CCCGATAAACCCTTCACCCG"
OUTDIR_U6="/gpfs/home/jah14xyu/CRISPResso2/results/U6"

mkdir -p "$OUTDIR_U6"

CRISPResso \
  --fastq_r1 "$FASTQ_R1_U6" \
  --fastq_r2 "$FASTQ_R2_U6" \
  --amplicon_seq "$AMPLICONSEQ_U6" \
  --guide_seq "$GUIDESEQ_U6" \
  --output_folder "$OUTDIR_U6" \
  --name "U6" \
  --n_processes 8

