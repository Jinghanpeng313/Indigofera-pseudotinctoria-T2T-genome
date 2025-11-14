#!/usr/bin/env bash
# -------------------------------
# Config
# -------------------------------
REF="Hap1.fa"                  # Reference genome
READS_FASTQ="output.fastq"       # Basecalled reads (FASTQ) used by nanopolish
UL_FASTQ="pass.ul.fq"            # (Optional) ultra-long reads; if present, will align too
FAST5_DIR="fast5_pass"           # Directory with FAST5 files for nanopolish indexing
SEQ_SUMMARY=""                   # (Optional) sequencing_summary.txt; leave empty if unknown
THREADS=64
TMP_PREFIX="tmp"

# Paths to tools/scripts
CALC_FREQ="scripts/calculate_methylation_frequency.py"

# Outputs
ALIGN_BAM="output.sorted.bam"
ALIGN_BAM_UL="output.ul.sorted.bam"      # if UL_FASTQ is provided
METH_CALLS="methylation_calls.tsv"
METH_FREQ="methylation_frequency.tsv"

# Hap1 merging (per-site methylation counts)
HAP1_FREQ_LIST=("hap1_methylation_frequency.tsv")  # add more files here if you have multiple cells/runs
HAP1_MERGED="hap1_merged.methy"                    # bed-like: chr  pos  pos  methylated_count  total_count
HAP1_MERGED_SORTED="hap1_merged.methy.sorted"

# Window-level summaries
GENOME_FAI="../01data/hap1.fa.fai"         # fai with contig lengths (chrom \t length ...)
WIN_SIZE=500
HAP1_WINDOWS="hap1.windows.${WIN_SIZE}bp.bed"
HAP1_MERGED_SUM="hap1_merged.methy.${WIN_SIZE}bp.sum.tsv"
HAP1_MERGED_WIN_FREQ="hap1_merged.methy.${WIN_SIZE}bp.freq.tsv"

# -------------------------------
# 0) Nanopolish index (links reads to FAST5)
#    NOTE: nanopolish expects FAST5 directory via -d; optionally provide -s sequencing_summary
# -------------------------------
echo "[Index] nanopolish index"
if [[ -n "${SEQ_SUMMARY}" ]]; then
  nanopolish index -d "${FAST5_DIR}" -s "${SEQ_SUMMARY}" "${READS_FASTQ}"
else
  nanopolish index -d "${FAST5_DIR}" "${READS_FASTQ}"
fi

# -------------------------------
# 1) Align reads (ONT)
#    Option A: use READS_FASTQ
#    Option B (optional): also align UL reads if provided; primary BAM will be ALIGN_BAM
# -------------------------------
echo "[Align] minimap2 -> SAMtools sort (standard reads)"
minimap2 -a -x map-ont -t "${THREADS}" "${REF}" "${READS_FASTQ}" \
  | samtools sort -@ "${THREADS}" -T "${TMP_PREFIX}" -o "${ALIGN_BAM}"

if [[ -s "${UL_FASTQ}" ]]; then
  echo "[Align] minimap2 -> SAMtools sort (ultra-long reads)"
  minimap2 -a -x map-ont -t "${THREADS}" "${REF}" "${UL_FASTQ}" \
    | samtools sort -@ "${THREADS}" -T "${TMP_PREFIX}" -o "${ALIGN_BAM_UL}"
fi

echo "[Index] samtools index"
samtools index "${ALIGN_BAM}"

# (Optional) if you want to replace ALIGN_BAM with the UL BAM, set:
# ALIGN_BAM="${ALIGN_BAM_UL}" && samtools index "${ALIGN_BAM}"

# -------------------------------
# 2) Call methylation with nanopolish
# -------------------------------
echo "[Call] nanopolish call-methylation"
nanopolish call-methylation \
  -t "${THREADS}" \
  -r "${READS_FASTQ}" \
  -b "${ALIGN_BAM}" \
  -g "${REF}" \
  > "${METH_CALLS}"

# -------------------------------
# 3) Calculate per-site methylation frequency
# -------------------------------
echo "[Freq] calculate_methylation_frequency.py"
# -c 2 => group by per-site context; adjust if needed for your script version
"${CALC_FREQ}" -c 2 "${METH_CALLS}" > "${METH_FREQ}"

# -------------------------------
# 4) Merge methylation frequency from multiple runs/cells (hap1 example)
#    Expect input TSVs with columns like:
#    chrom  pos  ...  methylated_count(col5)  unmethylated_count(col6)
#    The merge sums meth/unmeth counts per (chrom,pos) and outputs 5 columns:
#    chrom  pos  pos  meth_sum  total_sum
# -------------------------------
echo "[Merge] merging hap1 per-site counts across runs"
# If you have multiple files, list them in HAP1_FREQ_LIST array above.
# Here we show a robust merge using awk; it accepts a single concatenated stream.
cat "${HAP1_FREQ_LIST[@]}" \
 | grep -v -E '^chro|^chrom' \
 | awk '{a[$1"__"$2]+=$5; b[$1"__"$2]+=$6} END{for(i in a){print i"\t"a[i]"\t"b[i]}}' \
 | sed 's/__/ /' \
 | awk '{print $1"\t"$2"\t"$2"\t"$3"\t"($3+$4)}' \
 | sort -k1,1 -k2,2n \
 > "${HAP1_MERGED}"

# Sanity: also write a clean sorted version (identical if above used sort)
cp "${HAP1_MERGED}" "${HAP1_MERGED_SORTED}"

# -------------------------------
# 5) Compute methylation in fixed windows
#    a) Make windows from genome sizes (fai)
#    b) Sum meth and total counts into windows
#    c) Compute window methylation fraction = meth_sum / total_sum
# -------------------------------
echo "[Windows] make ${WIN_SIZE} bp windows"
bedtools makewindows -g "${GENOME_FAI}" -w "${WIN_SIZE}" > "${HAP1_WINDOWS}"

echo "[Windows] sum counts per window"
# Assumes HAP1_MERGED has: chrom  start  end(=start)  meth_sum  total_sum
# bedtools map will sum columns 4 and 5 into each window
bedtools map -a "${HAP1_WINDOWS}" -b "${HAP1_MERGED_SORTED}" -c 4,5 -o sum \
  > "${HAP1_MERGED_SUM}"

echo "[Windows] compute methylation fraction"
# Input columns from map: win_chr win_start win_end sum_meth sum_total
awk 'BEGIN{OFS="\t"} {meth=$4; tot=$5; frac=(tot>0? meth/tot : "NA"); print $1,$2,$3,meth,tot,frac}' \
  "${HAP1_MERGED_SUM}" \
  > "${HAP1_MERGED_WIN_FREQ}"


