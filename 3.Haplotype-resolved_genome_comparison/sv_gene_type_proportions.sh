#!/usr/bin/env bash
#
# sv_te_gene_partition.sh
#
# Description:
#   Given SyRI SV output, genome GFF3, and TE GFF3 with Class= annotations,
#   this script:
#     1) Extracts structural variants (SVs) from syri.out into BED.
#     2) Builds TE BED from TE GFF.
#     3) Derives gene body / exon / intron-like / upstream2kb / downstream2kb regions.
#     4) Classifies SVs as genic vs intergenic and TE-overlapping vs non-TE.
#     5) Further classifies genic SVs into EXON / INTRON / UPSTREAM2kb / DOWNSTREAM2kb
#        by maximum overlap and a fixed priority.
#     6) Summarizes counts and percentages into summary.tsv.
#
# Usage:
#   ./sv_te_gene_partition.sh syri.out hap1_genome.gff3 hap1.all.TE_repeat.class.gff [OUTDIR]
#
# Input:
#   syri.out                         # SyRI main output
#   hap1_genome.gff3                # Genome annotation with gene/exon features
#   hap1.all.TE_repeat.class.gff    # TE annotation with "Class=" in attributes
#
# Output directory (default: sv_te_gene_stats_out) will contain:
#   - SV.bed
#   - TE.bed
#   - gene_body.bed / exon.bed / intron_like.bed
#   - upstream2kb.bed / downstream2kb.bed
#   - SV_in_genic.bed / SV_in_intergenic.bed
#   - SV_overlap_TE.bed / SV_no_TE.bed
#   - sv_genic_category.tsv
#   - summary.tsv

set -euo pipefail

############################################
# 0. Argument parsing
############################################

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 syri.out hap1_genome.gff3 hap1.all.TE_repeat.class.gff [OUTDIR]" >&2
  exit 1
fi

SV=$1          # syri.out
GFF=$2         # genome GFF3
TEGFF=$3       # TE GFF3
OUTDIR=${4:-sv_te_gene_stats_out}

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "[INFO] Input:"
echo "  SyRI SV file : $SV"
echo "  Genome GFF3   : $GFF"
echo "  TE GFF3       : $TEGFF"
echo "  Output dir    : $(pwd)"

############################################
# 0. Dependency check
############################################
command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not found in PATH" >&2; exit 1; }
command -v python3  >/dev/null 2>&1 || { echo "[ERROR] python3 not found in PATH" >&2; exit 1; }

############################################
# 1) Build SV.bed from SyRI output
############################################
# Keep only "true" SV types (filter out SNP/SYN etc.)
# Allowed SV type keywords (case-sensitive)
SV_KEEP_REGEX='INS|DEL|INV|DUP|TRA|TRANS|INVDP|DUPAL|INVDPAL'

awk -v OFS="\t" -v RGX="$SV_KEEP_REGEX" '
BEGIN{ n=0 }
$1 ~ /^#/ { next }
{
  chr=$1; s=$2; e=$3; id=$9; typ=$11;

  # Skip non-SV types (e.g. SNP / SYN)
  if (typ ~ /SNP/ || typ ~ /SYN/) next;
  if (typ !~ RGX) next;

  # Basic coordinate sanity checks
  if (s=="" || e=="") next;
  if (s ~ /[^0-9]/ || e ~ /[^0-9]/) next;

  # Ensure start <= end
  if (e < s) { t=s; s=e; e=t; }

  # Convert to 0-based BED: [s-1, e)
  bs = s - 1;
  be = e;
  if (bs < 0) bs = 0;

  # For single-point events (e.g. INS), enforce at least 1 bp length
  if (be <= bs) be = bs + 1;

  n++;
  print chr, bs, be, (id=="" ? ("SV"n) : id), typ;
}
' "$SV" | sort -k1,1 -k2,2n > SV.bed

SVN=$(wc -l < SV.bed || echo 0)
echo "[INFO] Total SV in SV.bed: $SVN"

############################################
# 2) Build TE.bed from TE GFF
############################################
# Extract: chr, start, end, class (from Class=... attribute)
awk -v OFS="\t" '
BEGIN{ FS="\t" }
$0 ~ /^#/ { next }
{
  chr  = $1;
  s    = $4;
  e    = $5;
  attr = $9;

  # Default class
  class = "TE";

  # Extract TE class from attribute, e.g. Class=LTR/Gypsy;
  if (attr ~ /Class=[^;]+/) {
    match(attr, /Class=[^;]+/);
    class = substr(attr, RSTART+6, RLENGTH-6);
  }

  bs = s - 1; if (bs < 0) bs = 0;
  be = e;
  if (be <= bs) be = bs + 1;

  print chr, bs, be, class;
}
' "$TEGFF" | sort -k1,1 -k2,2n > TE.bed

############################################
# 3) Gene body / exon / intron-like / upstream2kb / downstream2kb
############################################

# 3.1 Gene body BED: [0-based, with strand]
awk -v OFS="\t" '
BEGIN{ FS="\t" }
$0 ~ /^#/ { next }
$3=="gene" {
  chr    = $1;
  s      = $4;
  e      = $5;
  strand = $7;

  bs = s - 1; if (bs < 0) bs = 0;
  be = e;
  if (be <= bs) be = bs + 1;

  print chr, bs, be, "GENE", strand;
}
' "$GFF" | sort -k1,1 -k2,2n > gene_body.bed

# 3.2 Exon union BED
awk -v OFS="\t" '
BEGIN{ FS="\t" }
$0 ~ /^#/ { next }
$3=="exon" {
  chr = $1;
  s   = $4;
  e   = $5;
  bs  = s - 1; if (bs < 0) bs = 0;
  be  = e;
  if (be <= bs) be = bs + 1;
  print chr, bs, be;
}
' "$GFF" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - > exon.bed

# 3.3 Intron-like BED: gene_body union - exon union (approximation)
cut -f1-3 gene_body.bed \
  | bedtools merge -i - > gene_merged.bed

bedtools subtract -a gene_merged.bed -b exon.bed > intron_like.bed

# 3.4 Upstream2kb / downstream2kb BED (strand-aware)
# First generate raw intervals with UP/DOWN labels
awk -v OFS="\t" '
{
  chr    = $1;
  bs     = $2;
  be     = $3;
  strand = ".";
  if (NF >= 5) strand = $5;

  len = 2000;

  # "left" and "right" relative to coordinate
  left_s  = bs - len; if (left_s < 0) left_s = 0;
  left_e  = bs;
  right_s = be;
  right_e = be + len;

  # For positive strand:
  #   upstream   = left side
  #   downstream = right side
  # For negative strand:
  #   upstream   = right side
  #   downstream = left side
  if (strand == "+") {
    if (left_e > left_s)  print chr, left_s,  left_e,  "UP";
    if (right_e > right_s)print chr, right_s, right_e, "DOWN";
  } else if (strand == "-") {
    if (left_e > left_s)  print chr, left_s,  left_e,  "DOWN";
    if (right_e > right_s)print chr, right_s, right_e, "UP";
  } else {
    # If no strand info, treat as positive
    if (left_e > left_s)  print chr, left_s,  left_e,  "UP";
    if (right_e > right_s)print chr, right_s, right_e, "DOWN";
  }
}
' gene_body.bed \
| sort -k1,1 -k2,2n > updn_raw.bed

# Merge each category separately
awk -v OFS="\t" '$4=="UP"{print $1,$2,$3}'   updn_raw.bed \
  | bedtools merge -i - > upstream2kb.bed

awk -v OFS="\t" '$4=="DOWN"{print $1,$2,$3}' updn_raw.bed \
  | bedtools merge -i - > downstream2kb.bed

############################################
# 4) Genic vs intergenic SV
############################################

# 4.1 Genic union: union of exon / intron-like / upstream2kb / downstream2kb
cat exon.bed \
    intron_like.bed \
    upstream2kb.bed \
    downstream2kb.bed \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > genic_union.bed

# 4.2 SV in genic/intergenic
bedtools intersect -u -a SV.bed -b genic_union.bed > SV_in_genic.bed
bedtools intersect -v -a SV.bed -b genic_union.bed > SV_in_intergenic.bed

############################################
# 5) SV overlap vs non-overlap with TE
############################################
bedtools intersect -u -a SV.bed -b TE.bed > SV_overlap_TE.bed
bedtools intersect -v -a SV.bed -b TE.bed > SV_no_TE.bed

############################################
# 6) Classify genic SVs into EXON / INTRON / UP2kb / DOWN2kb
############################################

# 6.1 Label each genic feature type
{
  awk -v OFS="\t" '{print $1,$2,$3,"EXON"}'          exon.bed
  awk -v OFS="\t" '{print $1,$2,$3,"INTRON"}'        intron_like.bed
  awk -v OFS="\t" '{print $1,$2,$3,"UPSTREAM2kb"}'   upstream2kb.bed
  awk -v OFS="\t" '{print $1,$2,$3,"DOWNSTREAM2kb"}' downstream2kb.bed
} | sort -k1,1 -k2,2n > genic_labeled.sorted.bed

# 6.2 Compute overlap length between genic SVs and labeled features
# SV_in_genic.bed: chr, start, end, id, type
# genic_labeled.sorted.bed: chr, start, end, category
bedtools intersect -wao -a SV_in_genic.bed -b genic_labeled.sorted.bed \
  > genic_overlap.tsv

# 6.3 For each SV ID, keep category with largest overlap, break ties by priority:
#     EXON > INTRON > UPSTREAM2kb > DOWNSTREAM2kb
awk -v OFS="\t" '
function pri(c) {
  if (c=="EXON")         return 4;
  if (c=="INTRON")       return 3;
  if (c=="UPSTREAM2kb")  return 2;
  if (c=="DOWNSTREAM2kb")return 1;
  return 0;
}
{
  sv_id = $4;
  cat   = $8;
  ol    = $9 + 0;

  # Lines with no overlap (cat == "." or ol == 0) are ignored
  if (cat=="." || ol<=0) next;

  if (!(sv_id in best_ol) ||
      ol > best_ol[sv_id] ||
      (ol == best_ol[sv_id] && pri(cat) > pri(best_cat[sv_id]))) {
    best_ol[sv_id]  = ol;
    best_cat[sv_id] = cat;
  }
}
END{
  for (id in best_cat) {
    print id, best_cat[id];
  }
}
' genic_overlap.tsv > sv_genic_category.tsv

############################################
# 7) Summary statistics
############################################

TOTAL=$(wc -l < SV.bed                  || echo 0)
IN_TE=$(wc -l < SV_overlap_TE.bed       || echo 0)
NO_TE=$(wc -l < SV_no_TE.bed            || echo 0)
IN_GENIC=$(wc -l < SV_in_genic.bed      || echo 0)
IN_INTERGENIC=$(wc -l < SV_in_intergenic.bed || echo 0)

EXON_N=$(awk '$2=="EXON"'         sv_genic_category.tsv | wc -l || echo 0)
INTRON_N=$(awk '$2=="INTRON"'     sv_genic_category.tsv | wc -l || echo 0)
UP_N=$(awk '$2=="UPSTREAM2kb"'    sv_genic_category.tsv | wc -l || echo 0)
DOWN_N=$(awk '$2=="DOWNSTREAM2kb"'sv_genic_category.tsv | wc -l || echo 0)

# helper for percentage (2 decimals)
pct() {
  python3 - <<PY
tot = float("$1")
num = float("$2")
print("0.00" if tot == 0 else f"{num*100.0/tot:.2f}")
PY
}

cat > summary.tsv <<TSV
Category\tCount\tPercent
SV_total\t$TOTAL\t100.00
SV_overlap_TE\t$IN_TE\t$(pct "$TOTAL" "$IN_TE")
SV_no_TE\t$NO_TE\t$(pct "$TOTAL" "$NO_TE")
SV_in_genic\t$IN_GENIC\t$(pct "$TOTAL" "$IN_GENIC")
SV_in_intergenic\t$IN_INTERGENIC\t$(pct "$TOTAL" "$IN_INTERGENIC")
SV_in_genic_EXON\t$EXON_N\t$(pct "$IN_GENIC" "$EXON_N")
SV_in_genic_INTRON\t$INTRON_N\t$(pct "$IN_GENIC" "$INTRON_N")
SV_in_genic_UPSTREAM2kb\t$UP_N\t$(pct "$IN_GENIC" "$UP_N")
SV_in_genic_DOWNSTREAM2kb\t$DOWN_N\t$(pct "$IN_GENIC" "$DOWN_N")
TSV

echo
echo "==== summary.tsv ===="
if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' summary.tsv
else
  cat summary.tsv
fi

echo
echo "[INFO] Key intermediate files:"
echo "  SV.bed                         # filtered SV (BED)"
echo "  TE.bed                         # TE intervals with class"
echo "  gene_body.bed / exon.bed / intron_like.bed / upstream2kb.bed / downstream2kb.bed"
echo "  SV_in_genic.bed / SV_in_intergenic.bed"
echo "  SV_overlap_TE.bed / SV_no_TE.bed"
echo "  sv_genic_category.tsv          # category of genic SVs"
echo "  summary.tsv                    # summary counts and percentages"
