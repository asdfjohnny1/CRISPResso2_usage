###############################################################################
# CRISPResso2 HPC GUIDE
#
# Sections:
#   1) Installation
#   2) Testing
#   3) Running on real data
###############################################################################

###############################################################################
# PART 1 — INSTALLING CRISPResso2
###############################################################################

# Load conda module (try whichever exists on your HPC)
module load  python/anaconda/2024.10/3.12.7

# Create a clean environment (only needed once)
conda create -n crispresso_trim \
 -c conda-forge -c bioconda \
 python=3.10 cutadapt fastp crispresso2

# Activate the environment
source activate crispresso_trim

# Install CRISPResso2 
pip install "git+https://github.com/pinellolab/CRISPResso2@v2.3.3"

# Confirm installation worked
CRISPResso --version
CRISPResso --help

# Response is a bit slow on mine .. not sure why 


###############################################################################
# PART 2 — TESTING THE INSTALLATION
# This uses a tiny synthetic FASTQ to confirm everything runs correctly
###############################################################################

# Make a test directory
mkdir -p /gpfs/home/jah14xyu/Postdoc/CRISPResso/crispresso2_test
cd /gpfs/home/jah14xyu/Postdoc/CRISPResso/crispresso2_test

# create an amplicon (120 bp) and an embedded 20 nt guide
AMPLICONSEQ="AAACCCAAACCCAAACCCAAACCCAAACCCAAACCCAAACCCAAACCCGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTT"
GUIDESEQ="CCCAAACCCAAACCCAAACC"  # 20 nt, present in the amplicon above

# Create a minimal FASTQ containing exactly the amplicon as the read
cat > test.fastq <<EOF
@read1
${AMPLICONSEQ}
+
$(python - <<'PY'
seq = """'"${AMPLICONSEQ}"'"""
print("F"*len(seq))
PY
)
EOF

gzip -f test.fastq

# Run CRISPResso2
CRISPResso \
  --fastq_r1 test.fastq.gz \
  --amplicon_seq "${AMPLICONSEQ}" \
  --guide_seq "${GUIDESEQ}" \
  --output_folder test_output \
  --name synthetic_test \
  --n_processes 1

# Inspect output
ls -lh test_output


###############################################################################
# PART 3 — RUNNING ON REAL DATA
###############################################################################

# uploading data to HPC

# first make directory for data 
mkdir /gpfs/home/nym23vru/CRISPResso2/raw_data

# uplaod data from onedrive to HPC - must use a local terminal - not connected to hpc 
Scp -r /Users/nym23vru/Library/CloudStorage/OneDrive-UniversityofEastAnglia/Chapter\ 2\ -\ 7SK\ Polymerase\ III/Illumina\ Sequencing\ Results/7SK\ and\ U6\ Illumina\ Sequences/raw_fastq  nym23vru@hali.uea.ac.uk:/gpfs/home/nym23vru/CRISPResso2/raw_data 


# Activate environment (if not already active)
module load  python/anaconda/2024.10/3.12.7
source activate crispresso_trim






###############################################################################
# B) PAIRED-END EXAMPLE
###############################################################################


# for 7sK

# Edit these variables
FASTQ_R1="/gpfs/home/nym23vru/CRISPResso2/raw_data/raw_fastq/7SK-gRNA_R1_001.fastq.gz"
FASTQ_R2="/gpfs/home/nym23vru/CRISPResso2/raw_data/raw_fastq/7SK-gRNA_R2_001.fastq.gz"
AMPLICONSEQ="TTCTGCGATAGCTTTTTCAACACCTGTACGACGCTGTGTGCCATGAACGAGTCCAAACCCGAGGTAGGTTCATCGCAAATTAATAGCGGCGGATCGGTTAACGCCTCCGAGGCGAATGCCAGTCGCTTACGCTCGCCACCCGATAAACCCTTCACCCGACCCGGCACGCCAATCAACGTATTCTGACATTTACCCAGCGAGAGGTCCTGTATAACTTGATCGACACGTTGTACTTTCTGTTTTTGTGTCATATGCCGT"
OUTDIR="/gpfs/home/nym23vru/CRISPResso2/results/7sk"
FORWARD_PRIMER="TTCTGCGATAGCTTTTTCAACA"
REVERSE_PRIMER="GACAAAAACACAGTATACGGCA"

# make sure results folder is made
mkdir -p $OUTDIR


# might need to trim the adapters  
mkdir /gpfs/home/nym23vru/CRISPResso2/data/trimmed

cd /gpfs/home/nym23vru/CRISPResso2/data/trimmed

fastp \
 -i $FASTQ_R1 \
 -I $FASTQ_R2 \
 -o trimmed_R1.fastq.gz \
 -O trimmed_R2.fastq.gz \
 --detect_adapter_for_pe \
 --thread 8

# remove primers from sequence 
#


cutadapt \
 -g $FORWARD_PRIMER \
 -G $REVERSE_PRIMER \
 -o trimmed_R1.fastq.gz \
 -p trimmed_R2.fastq.gz \
 $FASTQ_R1 $FASTQ_R2

# Run CRISPResso2
CRISPResso \
  --fastq_r1 "$FASTQ_R1" \
  --fastq_r2 "$FASTQ_R2" \
  --amplicon_seq "$AMPLICONSEQ" \
  --guide_seq "$GUIDESEQ" \
  --output_folder "$OUTDIR" \
  --name "7sk" \
  --n_processes 8

CRISPResso \
 --fastq_r1 trimmed_R1.fastq.gz \
 --fastq_r2 trimmed_R2.fastq.gz \
 --amplicon_seq "$AMPLICONSEQ" \
 --guide_seq "$GUIDESEQ" \
 --output_folder "$OUTDIR" \
 --name "7sk" \
 --n_processes 8

# for 3u6 

# Edit these variables
FASTQ_R1="/gpfs/home/nym23vru/CRISPResso2/raw_data/raw_fastq/U6-gRNA_R1_001.fastq.gz"
FASTQ_R2="/gpfs/home/nym23vru/CRISPResso2/raw_data/raw_fastq/U6-gRNA_R2_001.fastq.gz"
AMPLICONSEQ="TTCTGCGATAGCTTTTTCAACACCTGTACGACGCTGTGTGCCATGAACGAGTCCAAACCCGAGGTAGGTTCATCGCAAATTAATAGCGGCGGATCGGTTAACGCCTCCGAGGCGAATGCCAGTCGCTTACGCTCGCCACCCGATAAACCCTTCACCCGACCCGGCACGCCAATCAACGTATTCTGACATTTACCCAGCGAGAGGTCCTGTATAACTTGATCGACACGTTGTACTTTCTGTTTTTGTGTCATATGCCGT"
GUIDESEQ="CCCGATAAACCCTTCACCCG"
OUTDIR="/gpfs/home/nym23vru/CRISPResso2/results/U6"

# make sure results folder is made
mkdir -p $OUTDIR

CRISPResso \
  --fastq_r1 "$FASTQ_R1" \
  --fastq_r2 "$FASTQ_R2" \
  --amplicon_seq "$AMPLICONSEQ" \
  --guide_seq "$GUIDESEQ" \
  --output_folder "$OUTDIR" \
  --name "7sk" \
  --n_processes 8